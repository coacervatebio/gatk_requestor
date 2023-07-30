import os
from pathlib import Path, PurePath
from typing import List, Tuple
from flytekit import kwtypes, workflow, dynamic, task, ContainerTask
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory
from flytekit.extras.tasks.shell import OutputLocation, ShellTask
from flytekitplugins.pod import Pod
from data.workflows.pod_templates import yagna_requestor
from data.workflows.utils import VCF, Alignment, run_golem
from data.workflows.pod_templates import yagna_requestor_ps
from data.workflows.flyte_haplotypecaller import main
from data.workflows.config import current_image
from data.workflows.tasks import index_cram
from data.workflows.utils import get_dir

@workflow
def wf() -> FlyteDirectory:
    return get_dir(dirpath='s3://my-s3-bucket/input-data/alignments/full')
