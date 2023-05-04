import os
from pathlib import Path
from typing import List, Tuple
from flytekit import kwtypes, workflow, dynamic, task, ContainerTask, current_context
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory

@dynamic
def process_samples(indir: FlyteDirectory, regs: List[str]) -> FlyteDirectory:

    working_dir = current_context().working_directory
    local_dir = Path(os.path.join(working_dir, "regions"))
    local_dir.mkdir(exist_ok=True)

    for f in os.listdir(indir):
        fi = os.path.join(indir, f)
        idx = ic(al_in=fi)
        for r in regs:
            per_reg = sc(al_in=fi, idx_in=idx, reg=r)
            per_reg_idx = ic(al_in=per_reg)
            os.rename(per_reg, os.path.join(local_dir, os.path.basename(per_reg)))
            os.rename(per_reg, os.path.join(local_dir, os.path.basename(per_reg_idx)))
    return local_dir

@task
def get_dir(dirpath: str) -> FlyteDirectory:
    fd = FlyteDirectory(path=dirpath)
    return fd

@workflow
def wf(dirpath: str) -> FlyteDirectory:
    
    fd = get_dir(dirpath=dirpath)
    regs = ['chr21', 'chr22']
    out_ = process_samples(indir=fd, regs=regs)
    return out_
