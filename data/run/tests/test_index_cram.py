import os
from pathlib import Path, PurePath
from typing import List, Tuple
from flytekit import kwtypes, workflow, dynamic, task, ContainerTask
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory
from flytekit.extras.tasks.shell import OutputLocation, ShellTask
from flytekitplugins.pod import Pod
from run.tasks.mapping import index_cram_sh
from run.tasks.utils import get_file

@workflow
def wf():
    # al = get_file(filepath='s3://my-s3-bucket/input-data/alignments/full/HG03633_sub.cram')
    index_cram_sh()
