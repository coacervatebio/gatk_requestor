import os
from datetime import datetime
from pathlib import PurePath, Path
from typing import List, Tuple
from flytekit import kwtypes, workflow, dynamic, task, ContainerTask, current_context
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory
from flytekit.extras.tasks.shell import OutputLocation, ShellTask
from flytekitplugins.pod import Pod
from .pod_templates import yagna_requestor_ps
from .utils import get_dir, run_golem, dir_to_alignments, Alignment, VCF, dir_to_vcfs, prep_db_import, prep_gather_vcfs
from .flyte_haplotypecaller import main
from .hello_golem import hello
from .tasks import golem_call_variants, combine_region, genotype, gather_vcfs
from .config import current_image, reference_location

@workflow
def wf():
    d1 = get_dir(dirpath='s3://my-s3-bucket/input-data/geno_out/chr21')
    d2 = get_dir(dirpath='s3://my-s3-bucket/input-data/geno_out/chr22')
    fmt, fd = prep_gather_vcfs(combi_dirs=[d1, d2])
    gather_vcfs(vnames_fmt=fmt, vdir=fd)