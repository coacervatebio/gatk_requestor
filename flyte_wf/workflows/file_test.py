import os
from typing import List
from flytekit import task, workflow
from flytekit.types.directory import FlyteDirectory
from flytekit.types.file import FlyteFile


@task
def t1(inpath: str) -> str:
    ff = FlyteFile(inpath)
    content = ''
    with open(ff, 'r') as in_:
        content = in_.read()
    return content


@workflow
def wf():
    t1(inpath="s3://my-s3-bucket/arbitrary/hello.txt")