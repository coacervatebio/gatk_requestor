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
from .utils import get_dir, run_golem, dir_to_alignments, Alignment, VCF, dir_to_vcfs
from .flyte_haplotypecaller import main
from .hello_golem import hello
from .tasks import golem_call_variants

@workflow
def wf():
    fd = get_dir(dirpath='s3://my-s3-bucket/data/rp/ask2j94d6nwv9hz2f5kd-n1-0-dn14-0')
    als = dir_to_vcfs(indir=fd)