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

# @task
# def make_flytedir() -> FlyteDirectory:
#     local_dir = Path("s3://my-s3-bucket/user-inputs-made")
#     local_dir.mkdir(exist_ok=True)
#     return FlyteDirectory(path=str(local_dir))

@task
def read_samps(samples: FlyteFile) -> List[str]:
    lst = []
    with open(samples, 'r') as s_in:
        for i in s_in.readlines():
            ffs = f"s3://my-s3-bucket/input-data/{i}"
            lst.append(ffs.strip())
    return lst

@dynamic
def process_samples(infiles: List[str]) -> str:
    for i in infiles:
        ff = FlyteFile(path=i)
        s1(infile=ff)
    # first = infiles[0]
    # s1(infile=FlyteFile(first))
    return "PROCESSED"

@workflow
def wf(samples: FlyteFile) -> str:
    ffs = read_samps(samples=samples)
    out_ = process_samples(infiles=ffs)
    return out_
