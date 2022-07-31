#!/usr/bin/env python3
import asyncio
from typing import AsyncIterable

from yapapi import Golem, Task, WorkContext
from yapapi.log import enable_default_logger
from yapapi.payload import vm


async def worker(context: WorkContext, tasks: AsyncIterable[Task]):
    async for task in tasks:
        script = context.new_script()
        # script.upload_file("/home/vagrant/host_shared/golem/haplotype_caller/container/gatk_release/gatk-package-4.2.6.1-local.jar", "/golem/work/gatk.jar")
        # future_result = script.run("/usr/bin/java", "-jar", "/golem/work/gatk.jar", "HaplotypeCaller", "--version")

        # script.upload_file("/home/vagrant/host_shared/golem/haplotype_caller/run.sh", "/golem/entrypoint/run.sh")
        # script.upload_file("/home/vagrant/host_shared/snakemake/results/2_sample/alignments/chr1/HG03633.alt_bwamem_GRCh38DH.20150826.PJL.exome_chr1.cram", "/golem/input/3633_chr1.cram")
        # script.upload_file("/home/vagrant/host_shared/snakemake/results/2_sample/alignments/chr1/HG03633.alt_bwamem_GRCh38DH.20150826.PJL.exome_chr1.cram.crai", "/golem/input/3633_chr1.cram.crai")
        
        # future_result = script.run("/usr/bin/java", "-jar", "/home/gatk-local.jar", "HaplotypeCaller")
        # future_result = script.run("/bin/sh", "-x", "/golem/entrypoint/run.sh")
        future_result = script.run("/bin/ls", "/golem")
        # future_result = script.run("/opt/java/openjdk/bin/java", "-version")
        # future_result = script.run("java", "-jar", "/home/gatk-local.jar", "HaplotypeCaller", "--version")

        yield script

        task.accept_result(result=await future_result)


async def main():
    package = await vm.repo(
        # image_hash="43e2aeeee094fba1f44f5bf5d0d9cb11e9565fd419f14d74b3e8fe83", #My uploaded hello image
        # image_hash="479be4f6c3308dd77014ea07987e2df80bf5eaf59a6a5dc4f1264f94", #Using HC image
        # image_hash="aba6933ce8bfb35654051f8485a9310d9b2f70c09d044c10a50b9cfe", #Just jre
        # image_hash="94642154678189d09ad2fba485efa9344f1aeffa774c4bb71ba346e8", #HC no ref hard vol
        # image_hash="cc4d7097d7a5aa1cd2630e90c482b79e6a1a706165e7438b563b37a5", #steljevoor java image
        # image_hash="b5b66f21a48ab81367ca7f9e5ecba75278b8a5abdeffb939529abfc2", #steljevoor hello-world
        # image_hash="f371fd0344c4d91a0fba69a2b2f96cd49860124ee421d53ee2db76c4", #alpine hc with ref
        image_hash="5e1013bb93f876d38b0fd6c24435bd118a5700c118a64c5d7f61630f", #alpine vols only
    )

    tasks = [Task(data=None)]

    async with Golem(budget=1.0, subnet_tag="goth") as golem:
        async for completed in golem.execute_tasks(worker, tasks, payload=package):
            print(completed.result.stdout)


if __name__ == "__main__":
    enable_default_logger(log_file="hello.log")

    loop = asyncio.get_event_loop()
    task = loop.create_task(main())
    loop.run_until_complete(task)
