import os
from pathlib import Path, PurePath
from typing import List, Tuple
from flytekit import kwtypes, workflow, dynamic, task, ContainerTask
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory
from flytekit.extras.tasks.shell import OutputLocation, ShellTask
from flytekitplugins.pod import Pod
from run.tasks.mapping import index_cram
from run.tasks.utils import get_file, compare_files

@workflow
def wf():
    al = get_file(filepath='s3://my-s3-bucket/test-assets/index-cram/HG03633_sub.cram')
    idx_actual = index_cram(al_in=al)
    idx_expected = get_file(filepath='s3://my-s3-bucket/test-assets/index-cram/HG03633_sub.cram.crai')
    equivalent = compare_files(actual=idx_actual, expected=idx_expected)
