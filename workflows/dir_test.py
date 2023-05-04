import os
from typing import List
from flytekit import task, workflow
from flytekit.types.directory import FlyteDirectory
from flytekit.types.file import FlyteFile


@task
def t1(indir: FlyteDirectory) -> List[FlyteFile]:
    files = []
    for f in os.listdir(indir):
        files.append(os.path.join(indir, f))
    return files


@workflow
def wf():
    t1(indir=FlyteDirectory("s3://my-s3-bucket/input-data"))