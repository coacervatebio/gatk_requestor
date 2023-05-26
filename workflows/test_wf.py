import os
from datetime import datetime
from pathlib import PurePath, Path
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
def test_task(als: List[Alignment]):# -> List[VCF]:
    payloads = []
    for al in als:
        working_dir = PurePath(al.almt.path).parent
        print(working_dir)
        payload = {
            "sample": al.sample,
            "region_str": al.reg,
            "working_dir": working_dir,
            "req_align_path": Path(al.almt.path),
            "req_align_index_path": Path(al.idx.path),
            "req_vcf_path": working_dir.joinpath(f"{al.sample}.g.vcf.gz"),
            "req_vcf_index_path": working_dir.joinpath(f"{al.sample}.g.vcf.gz.tbi"),
        }
        print(payload)
        payloads.append(payload)
    run_golem(main(pls=payloads), log_file=True)
        # vcf = VCF(
        #     sample=al.sample,
        #     reg=al.reg,
        #     vcf=FlyteFile(path=f'{basename}.g.vcf.gz'),
        #     idx=FlyteFile(path=f'{basename}.g.vcf.gz.tbi')
        # )
        # vcfs.append(vcf)
    # return vcfs
    # return vcf objs
    # hello()

@workflow
def wf():
    fd = get_dir(dirpath='s3://my-s3-bucket/input-data/alignments/chr22')
    als = dir_to_alignments(indir=fd)
    test_task(als=als)
