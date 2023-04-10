import os
from typing import List
from flytekit import task, workflow
from flytekit.types.directory import FlyteDirectory
from flytekit.types.file import FlyteFile


@task
def t1() -> List[FlyteFile]:
    files = []
    for f in os.listdir(FlyteDirectory("s3://my-s3-bucket/arbitrary")):
        files.append(f)
    return files


@workflow
def wf():
    t1()