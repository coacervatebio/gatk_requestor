import os
from typing import List
from flytekit import task, workflow
from flytekit.types.directory import FlyteDirectory
from flytekit.types.file import FlyteFile


@task
def t1(infile: FlyteFile) -> str:
    content = ''
    with open(infile, 'r') as in_:
        content = in_.read()
    return content


@workflow
def wf(infile: FlyteFile) -> str:
    t1_out = t1(infile=infile)
    return t1_out