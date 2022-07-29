#!/usr/bin/env python3
"""
This file contains the requestor part of our application. There are three areas here:
1. Splitting the data into multiple tasks, each of which can be executed by a provider.
2. Defining what commands must be run within the provider's VM.
3. Scheduling the tasks via a yagna node running locally.
"""
import os
import argparse
import asyncio
from datetime import timedelta
import json
import math
from pathlib import Path, PurePath
from tempfile import gettempdir
from typing import AsyncIterable, Iterator
from uuid import uuid4

from yapapi import Golem, Task, WorkContext
from yapapi.log import enable_default_logger
from yapapi.payload import vm

prov_outpath = Path("/golem/output")
prov_inpath = Path("/golem/input")

# CLI arguments definition
arg_parser = argparse.ArgumentParser()
arg_parser.add_argument("--alignments", type=Path)
arg_parser.add_argument("--vcfs", type=Path)
arg_parser.add_argument("--script", type=Path, default=Path("run.sh"))
arg_parser.add_argument("--subnet", type=str, default="goth")

# Container object for parsed arguments
args = argparse.Namespace()

ENTRYPOINT_PATH = "/golem/entrypoint/run.sh"
TASK_TIMEOUT = timedelta(minutes=10)


def data(alignments_dir: Path, vcfs_dir: Path) -> Iterator[Task]:
    """Prepare a task object for every region-specific alignment
    """

    for reg_dir in alignments_dir.glob('chr*'):
        for almt in reg_dir.glob('*cram'):
            sample = almt.with_suffix('').name
            inputs = {
                'sample': sample,
                'req_align_path': reg_dir.joinpath(f'{sample}.cram'),
                'req_align_index_path': reg_dir.joinpath(f'{sample}.cram.crai'),
                'prov_align_path': prov_inpath.joinpath(f'{sample}.cram'),
                'prov_align_index_path': prov_inpath.joinpath(f'{sample}.cram.crai'),
                'region_str': reg_dir.name,
                'req_vcf_path': vcfs_dir.joinpath(f'{sample}.g.vcf.gz'),
                'req_vcf_index_path': vcfs_dir.joinpath(f'{sample}.g.vcf.gz.tbi'),
                'prov_vcf_path': prov_outpath.joinpath(f'{sample}.g.vcf.gz'),
                'prov_vcf_index_path': prov_outpath.joinpath(f'{sample}.g.vcf.gz.tbi')
            }
            print(inputs)
            yield Task(data=inputs)

async def steps(context: WorkContext, tasks: AsyncIterable[Task]):
    """Prepare a sequence of steps which need to happen for a task to be computed.

    `Script` is a utility which allows us to define a series of commands to interact
    with a provider. It's created using the provided `WorkContext` instance.
    Tasks are provided from a common, asynchronous queue.
    The signature of this function cannot change, as it's used internally by `Executor`.
    """
    script = context.new_script(timeout=timedelta(minutes=5))
    script.upload_file(str(args.script), ENTRYPOINT_PATH)

    async for task in tasks:
        # Upload input alignments
        script.upload_file(task.data['req_align_path'], task.data['prov_align_path'])
        script.upload_file(task.data['req_align_index_path'], task.data['prov_align_index_path'])

        run_args = f"{task.data['prov_align_path']} {task.data['region_str']} {task.data['prov_vcf_path']}"
        future_result = script.run(ENTRYPOINT_PATH, args=run_args)

        script.download_file(task.data['prov_vcf_path'], task.data['req_vcf_path'])
        script.download_file(task.data['prov_vcf_index_path'], task.data['req_vcf_index_path'])

        # Pass the prepared sequence of steps to Executor
        yield script

        # Mark task as accepted and set its result
        task.accept_result(result=await future_result)

        # Re-initialize the script so that `upload_file` is executed only once per worker
        script = context.new_script(timeout=timedelta(minutes=5))


async def main():
    # Set of parameters for the VM run by each of the providers
    package = await vm.repo(
        image_hash="479be4f6c3308dd77014ea07987e2df80bf5eaf59a6a5dc4f1264f94",
        min_mem_gib=1.0,
        min_storage_gib=2.0,
    )

    async with Golem(budget=1, subnet_tag=args.subnet) as golem:
        async for task in golem.execute_tasks(
            steps, data(args.alignments, args.vcfs), payload=package, timeout=TASK_TIMEOUT
        ):
            if task.result:
                print(f"Calling: {task.data['sample']}")
        print("End of `async with Golem`...")


if __name__ == "__main__":
    args = arg_parser.parse_args()

    loop = asyncio.get_event_loop()
    task = loop.create_task(main())

    # yapapi debug logging to a file
    enable_default_logger(log_file="haplotype_caller.log")

    try:
        loop.run_until_complete(task)
    except KeyboardInterrupt:
        # Make sure Executor is closed gracefully before exiting
        task.cancel()
        loop.run_until_complete(task)
