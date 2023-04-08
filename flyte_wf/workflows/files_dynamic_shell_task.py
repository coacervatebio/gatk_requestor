import os
from pathlib import Path
from typing import List
from flytekit import kwtypes, workflow, dynamic, task
from flytekit.extras.tasks.shell import OutputLocation, ShellTask
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory

s1 = ShellTask(
    name="shorten",
    debug=True,
    script="""
    set -ex
    echo mooooooore >> {inputs.infile}
    """,
    inputs=kwtypes(infile=FlyteFile),
    output_locs=[OutputLocation(var="i", var_type=FlyteFile, location="{inputs.infile}")],
)

@task
def make_flytedir() -> FlyteDirectory:
    local_dir = Path("s3://my-s3-bucket/user-inputs-made")
    local_dir.mkdir(exist_ok=True)
    return FlyteDirectory(path=str(local_dir))

@task
def my_task(indir: FlyteDirectory) -> str:
    for i in os.listdir(indir):
        s1(infile=i)
    return "DYN_DONE"

@workflow
def wf() -> FlyteDirectory:
    flytedir = make_flytedir()
    return my_task(indir=flytedir)
