#!/usr/bin/env python3
"""
This file contains the requestor part of our application. There are three areas here:
1. Splitting the data into multiple tasks, each of which can be executed by a provider.
2. Defining what commands must be run within the provider's VM.
3. Scheduling the tasks via a yagna node running locally.
"""
import os
import json
import asyncio
import argparse
import subprocess
from time import sleep
from datetime import timedelta, datetime
from pathlib import Path
from tempfile import gettempdir
from typing import AsyncIterable, Iterator
from uuid import uuid4
from yapapi import Golem, Task, WorkContext
from yapapi.log import enable_default_logger
from yapapi.payload import vm

PROV_INPATH = Path("/golem/input")
PROV_OUTPATH = Path("/golem/output")
ENTRYPOINT_PATH = "/run/run.sh"
TASK_TIMEOUT = timedelta(hours=2)


def data(alignments_dir: Path, vcfs_dir: Path) -> Iterator[Task]:
    """Prepare a task object for every region-specific alignment"""

    for almt in alignments_dir.glob("*cram"):
        sample = almt.with_suffix("").name
        inputs = {
            "sample": sample,
            "req_align_path": alignments_dir.joinpath(f"{sample}.cram"),
            "req_align_index_path": alignments_dir.joinpath(f"{sample}.cram.crai"),
            "prov_align_path": PROV_INPATH.joinpath(f"{sample}.cram"),
            "prov_align_index_path": PROV_INPATH.joinpath(f"{sample}.cram.crai"),
            "region_str": alignments_dir.name,
            "req_vcf_path": vcfs_dir.joinpath(alignments_dir.name, f"{sample}.g.vcf.gz"),
            "req_vcf_index_path": vcfs_dir.joinpath(
                alignments_dir.name, f"{sample}.g.vcf.gz.tbi"
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
    script = context.new_script(timeout=timedelta(minutes=30))

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


async def main(alpath, vcfpath):
    # Set of parameters for the VM run by each of the providers
    package = await vm.repo(
        image_hash='635e41034ced5d0622d0760bf6aac8377fdf225154d3f306f4fca805',
        min_mem_gib=4.0,
        min_storage_gib=2.0,
    )

    async with Golem(budget=5, subnet_tag='public') as golem:
        async for completed in golem.execute_tasks(
            steps,
            data(alpath, vcfpath),
            payload=package,
            timeout=TASK_TIMEOUT,
        ):
            print(completed.result.stdout)


def run(alpath, vcfpath):

    loop = asyncio.get_event_loop()
    task = loop.create_task(main(alpath, vcfpath))

    # yapapi debug logging to a file
    enable_default_logger(
        log_file=f"/tmp/haplotype_caller_{datetime.now().strftime('%Y%m%d-%H%M')}.log"
    )

    # Set app key
    while os.getenv('YAGNA_APPKEY') is None:
        key_list = subprocess.run(["yagna", "app-key", "list", "--json"], capture_output=True)
        os.environ['YAGNA_APPKEY'] = json.loads(key_list.stdout)[0].get('key')
        sleep(5)

    try:
        loop.run_until_complete(task)
    except KeyboardInterrupt:
        # Make sure Executor is closed gracefully before exiting
        task.cancel()
        loop.run_until_complete(task)
