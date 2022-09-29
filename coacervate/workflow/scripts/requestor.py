#!/usr/bin/env python3
"""
This file contains the requestor part of our application. There are three areas here:
1. Splitting the data into multiple tasks, each of which can be executed by a provider.
2. Defining what commands must be run within the provider's VM.
3. Scheduling the tasks via a yagna node running locally.
"""
import argparse
import asyncio
from time import sleep
from datetime import timedelta, datetime
from pathlib import Path
from tempfile import gettempdir
from typing import AsyncIterable, Iterator
from uuid import uuid4

from yapapi import Golem, Task, WorkContext
from yapapi.log import enable_default_logger
from yapapi.payload import vm


# CLI arguments definition
arg_parser = argparse.ArgumentParser()
arg_parser.add_argument(
    "--alignments", type=Path, default=Path("/data/results/alignments")
)
arg_parser.add_argument("--vcfs", type=Path, default=Path("/data/results/hc_out"))
arg_parser.add_argument("--subnet", type=str, default="devnet-beta")
arg_parser.add_argument(
    "--image",
    type=str,
    default="635e41034ced5d0622d0760bf6aac8377fdf225154d3f306f4fca805",
)
arg_parser.add_argument("--debug", type=bool, default=False)

# Container object for parsed arguments
args = argparse.Namespace()

PROV_INPATH = Path("/golem/input")
PROV_OUTPATH = Path("/golem/output")
ENTRYPOINT_PATH = "/run/run.sh"
TASK_TIMEOUT = timedelta(minutes=30)


def data(alignments_dir: Path, vcfs_dir: Path) -> Iterator[Task]:
    """Prepare a task object for every region-specific alignment"""

    for reg_dir in alignments_dir.glob("chr*"):
        for almt in reg_dir.glob("*cram"):
            sample = almt.with_suffix("").name
            inputs = {
                "sample": sample,
                "req_align_path": reg_dir.joinpath(f"{sample}.cram"),
                "req_align_index_path": reg_dir.joinpath(f"{sample}.cram.crai"),
                "prov_align_path": PROV_INPATH.joinpath(f"{sample}.cram"),
                "prov_align_index_path": PROV_INPATH.joinpath(f"{sample}.cram.crai"),
                "region_str": reg_dir.name,
                "req_vcf_path": vcfs_dir.joinpath(reg_dir.name, f"{sample}.g.vcf.gz"),
                "req_vcf_index_path": vcfs_dir.joinpath(
                    reg_dir.name, f"{sample}.g.vcf.gz.tbi"
                ),
                "prov_vcf_path": PROV_OUTPATH.joinpath(f"{sample}.g.vcf.gz"),
                "prov_vcf_index_path": PROV_OUTPATH.joinpath(f"{sample}.g.vcf.gz.tbi"),
            }
            yield Task(data=inputs)


async def steps(context: WorkContext, tasks: AsyncIterable[Task]):
    """Prepare a sequence of steps which need to happen for a task to be computed.

    `Script` is a utility which allows us to define a series of commands to interact
    with a provider. It's created using the provided `WorkContext` instance.
    Tasks are provided from a common, asynchronous queue.
    The signature of this function cannot change, as it's used internally by `Executor`.
    """
    script = context.new_script()

    async for task in tasks:
        # Upload input alignments
        script.upload_file(task.data["req_align_path"], task.data["prov_align_path"])
        script.upload_file(
            task.data["req_align_index_path"], task.data["prov_align_index_path"]
        )

        run_args = [
            str(task.data["prov_align_path"]),
            str(task.data["region_str"]),
            str(task.data["prov_vcf_path"]),
        ]

        # future_result = script.run("/usr/bin/which", "tar")
        # future_result = script.run("/usr/bin/tar", "xf", "/run/reference_HG38.tar.zst", "-C", "/golem/entrypoint", "&&", "/bin/ls", "/golem/entrypoint")
        future_result = script.run("/bin/sh", ENTRYPOINT_PATH, *run_args)

        script.download_file(task.data["prov_vcf_path"], task.data["req_vcf_path"])
        script.download_file(
            task.data["prov_vcf_index_path"], task.data["req_vcf_index_path"]
        )

        # Pass the prepared sequence of steps to Executor
        yield script

        # Mark task as accepted and set its result
        task.accept_result(result=await future_result)


async def main():
    # Set of parameters for the VM run by each of the providers
    package = await vm.repo(
        image_hash=args.image,
        min_mem_gib=4.0,
        min_storage_gib=2.0,
    )

    async with Golem(budget=1, subnet_tag=args.subnet) as golem:
        async for completed in golem.execute_tasks(
            steps,
            data(args.alignments, args.vcfs),
            payload=package,
            timeout=TASK_TIMEOUT,
        ):
            print(completed.result.stdout)

    if args.debug:
        sleep(3600)


if __name__ == "__main__":
    args = arg_parser.parse_args()

    loop = asyncio.get_event_loop()
    task = loop.create_task(main())

    # yapapi debug logging to a file
    enable_default_logger(
        log_file=f"/data/workflow/logs/haplotype_caller_{datetime.now().strftime('%Y%m%d-%H%M')}_{args.image}.log"
    )

    try:
        loop.run_until_complete(task)
    except KeyboardInterrupt:
        # Make sure Executor is closed gracefully before exiting
        task.cancel()
        loop.run_until_complete(task)
