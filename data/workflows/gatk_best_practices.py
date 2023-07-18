import os
from pathlib import Path
from typing import List, Tuple
from flytekit import kwtypes, workflow, dynamic, task, ContainerTask, current_context
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory
from flytekitplugins.pod import Pod
from .tasks import index_cram, split_cram, golem_call_variants
from .utils import get_dir, run_golem, Alignment, VCF, return_alignment
from .pod_templates import yagna_requestor_ps
from .config import current_image

@dynamic(
    container_image=current_image,
    task_config=Pod(pod_spec=yagna_requestor_ps)
    )
def process_samples(indir: FlyteDirectory, regs: List[str]) -> List[VCF]:

    per_reg_als = []
    for f in os.listdir(indir):
        s = f.strip('.cram')
        fi = os.path.join(indir, f)
        idx = index_cram(al_in=fi)
        for r in regs:
            per_reg = split_cram(al_in=fi, idx_in=idx, reg=r)
            per_reg_idx = index_cram(al_in=per_reg)
            per_reg_al = return_alignment(sample=s, reg=r, almt=per_reg, idx=per_reg_idx)
            per_reg_als.append(per_reg_al)
    
    vcfs = golem_call_variants(als=per_reg_als)

    return vcfs

@workflow
def wf() -> List[VCF]:
    
    fd = get_dir(dirpath='s3://my-s3-bucket/input-data/alignments/full')
    regs = ['chr21', 'chr22']
    out_ = process_samples(indir=fd, regs=regs)
    return out_
