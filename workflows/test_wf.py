import os
from datetime import datetime
from pathlib import Path
from typing import List
from flytekit import kwtypes, workflow, dynamic, task, ContainerTask,current_context
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory
from flytekit.extras.tasks.shell import OutputLocation, ShellTask
from flytekitplugins.pod import Pod
from .pod_templates import yagna_requestor_ps
from .utils import get_dir, run_golem, dir_to_alignments, Alignment, VCF
from .flyte_haplotypecaller import main
from .hello_golem import hello

@task(
    container_image='docker.io/coacervate/requestor:latest',
    task_config=Pod(pod_spec=yagna_requestor_ps)
    )
def test_task(als: List[Alignment]) -> List[VCF]:
    # handle alignment and vcf paths before calling golem so it only deals with
    # pathlike objects
    run_golem(main(als=als), log_file=True)
    # return local_dir
    # hello()

@workflow
def wf():
    fd = get_dir(dirpath='s3://my-s3-bucket/input-data/alignments/chr22')
    als = dir_to_alignments(indir=fd)
    # test_task(indir=fd)
