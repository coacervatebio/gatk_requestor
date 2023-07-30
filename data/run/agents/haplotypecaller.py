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
import logging
import argparse
import subprocess
from time import sleep
from typing import List
from datetime import timedelta, datetime
from pathlib import Path, PurePath
from tempfile import gettempdir
from typing import AsyncIterable, Iterator
from uuid import uuid4
from yapapi import Golem, Task, WorkContext
from yapapi.payload import vm
from yapapi.log import enable_default_logger
from run.utils import Alignment

logger = logging.getLogger(__name__)
console_handler = logging.StreamHandler()
console_handler.setFormatter(logging.Formatter("[%(asctime)s %(levelname)s %(name)s] %(message)s"))
logger.addHandler(console_handler)
logger.setLevel(logging.DEBUG)

PROV_INPATH = PurePath("/golem/input")
PROV_OUTPATH = PurePath("/golem/output")
ENTRYPOINT_PATH = "/run/run.sh"
TASK_TIMEOUT = timedelta(hours=2)


def data(pls: List[dict]) -> Iterator[Task]:
    """Prepare a task object for every region-specific alignment"""
    logger.info(f"Processing payloads..")
    logger.debug(f"Payloads: {pls}")

    for pl in pls:
        logger.info(f'Processing alignment {pl["req_align_path"]}')
        sample = pl['sample']
        pl["prov_align_path"] = PROV_INPATH.joinpath(f"{sample}.cram")
        pl["prov_align_index_path"] = PROV_INPATH.joinpath(f"{sample}.cram.crai")
        pl["prov_vcf_path"] = PROV_OUTPATH.joinpath(f"{sample}.g.vcf.gz")
        pl["prov_vcf_index_path"] = PROV_OUTPATH.joinpath(f"{sample}.g.vcf.gz.tbi")
        logger.info(f"Prepared inputs: {pl}")
        yield Task(data=pl)


async def steps(context: WorkContext, tasks: AsyncIterable[Task]):
    """Prepare a sequence of steps which need to happen for a task to be computed.

    `Script` is a utility which allows us to define a series of commands to interact
    with a provider. It's created using the provided `WorkContext` instance.
    Tasks are provided from a common, asynchronous queue.
    The signature of this function cannot change, as it's used internally by `Executor`.
    """
    logger.debug("Executing steps on providers..")
    script = context.new_script(timeout=timedelta(minutes=30))

    async for task in tasks:
        logger.debug(f"Executing task: {task}")

        # Upload input alignments
        script.upload_file(task.data["req_align_path"], task.data["prov_align_path"])
        script.upload_file(
            task.data["req_align_index_path"], task.data["prov_align_index_path"]
        )
        logger.info(f"Uploaded alignment from {task.data['req_align_path']}")

        run_args = [
            str(task.data["prov_align_path"]),
            str(task.data["region_str"]),
            str(task.data["prov_vcf_path"]),
        ]
    
        future_result = script.run("/bin/sh", ENTRYPOINT_PATH, *run_args)

        script.download_file(task.data["prov_vcf_path"], task.data["req_vcf_path"])
        script.download_file(
            task.data["prov_vcf_index_path"], task.data["req_vcf_index_path"]
        )
        logger.info(f"Downloaded VCF to {task.data['req_vcf_path']}")

        # Pass the prepared sequence of steps to Executor
        yield script

        # Mark task as accepted and set its result
        task.accept_result(result=await future_result)


async def main(pls):
    logger.debug('Executing haplotypecaller main..')
    # Set of parameters for the VM run by each of the providers
    image_hash='635e41034ced5d0622d0760bf6aac8377fdf225154d3f306f4fca805'
    logger.info(f'Using provider image: {image_hash}')
    package = await vm.repo(
        image_hash=image_hash,
        min_mem_gib=4.0,
        min_storage_gib=2.0,
    )

    budget = 5
    subnet_tag = 'public'
    logger.info(f'Executing tasks with budget {budget} on {subnet_tag} subnet')
    async with Golem(budget=budget, subnet_tag=subnet_tag) as golem:
        async for completed in golem.execute_tasks(
            steps,
            data(pls),
            payload=package,
            timeout=TASK_TIMEOUT,
        ):
            print(completed.result.stdout)