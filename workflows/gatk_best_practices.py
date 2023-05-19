import os
from pathlib import Path
from typing import List, Tuple
from flytekit import kwtypes, workflow, dynamic, task, ContainerTask, current_context
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory
from flytekitplugins.pod import Pod
from .container_tasks import ic, sc
from .utils import get_dir, run_golem
from .pod_templates import yagna_requestor_ps
from .flyte_haplotypecaller import main

@task(
    container_image='docker.io/coacervate/requestor:latest',
    task_config=Pod(pod_spec=yagna_requestor_ps)
    )
def test_task(indir: FlyteDirectory) -> FlyteDirectory:
    working_dir = current_context().working_directory
    local_dir = Path(os.path.join(working_dir, "vcf_files"))
    local_dir.mkdir(exist_ok=True)
    local_indir = Path(indir)
    run_golem(main(alpath=local_indir, vcfpath=local_dir), log_file=True)
    return local_dir

@dynamic
def process_samples(indir: FlyteDirectory, regs: List[str]) -> FlyteDirectory:

    working_dir = current_context().working_directory
    local_dir = Path(os.path.join(working_dir, "regions"))
    local_dir.mkdir(exist_ok=True)

    per_regs = []
    for f in os.listdir(indir):
        fi = os.path.join(indir, f)
        idx = ic(al_in=fi)
        for r in regs:
            per_reg = sc(al_in=fi, idx_in=idx, reg=r)
            per_reg_idx = ic(al_in=per_reg)
            # per_regs.append(Alignment(almt=per_reg, idx=per_reg_idx))

    return local_dir

@workflow
def wf() -> FlyteDirectory:
    
    fd = get_dir(dirpath='s3://my-s3-bucket/input-data/alignments/full')
    regs = ['chr21', 'chr22']
    out_ = process_samples(indir=fd, regs=regs)
    return out_
