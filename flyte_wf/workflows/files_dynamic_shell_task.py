import os
from pathlib import Path
from typing import List
from flytekit import kwtypes, workflow, dynamic, task, ContainerTask
from flytekit.extras.tasks.shell import OutputLocation
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory

calculate_ellipse_area_python = ContainerTask(
    name="ellipse-area-metadata-python",
    input_data_dir="/var/inputs",
    output_data_dir="/var/outputs",
    inputs=kwtypes(a=float, b=float),
    outputs=kwtypes(area=float, metadata=str),
    image="ghcr.io/flyteorg/rawcontainers-python:v2",
    command=[
        "python",
        "calculate-ellipse-area.py",
        "{{.inputs.a}}",
        "{{.inputs.b}}",
        "/var/outputs",
    ],
)

s1 = ContainerTask(
    name="gatk-image-container-task",
    input_data_dir="/var/inputs",
    output_data_dir="/var/outputs",
    inputs=kwtypes(infile=FlyteFile),
    outputs=kwtypes(version=str, headlines=str),
    image="docker.io/coacervate/provider:latest",
    command=[
        "/run/flyte_test.sh",
        "{{.inputs.infile}}",
        "/var/outputs",
    ],
)

@task
def read_samps(samples: FlyteFile) -> List[FlyteFile]:
    lst = []
    with open(samples, 'r') as s_in:
        for i in s_in.readlines():
            ffs = f"s3://my-s3-bucket/input-data/{i}"
            lst.append(FlyteFile(path=ffs.strip()))
    return lst

@dynamic
def process_samples(infiles: List[FlyteFile]) -> str:

    s1_out = []
    for i in infiles:
        s1_out.append(s1(infile=i))

    # first = infiles[0]
    # s1(infile=FlyteFile(first))
    return "PROCESSED"

@workflow
def wf(samples: FlyteFile) -> str:
    ffs = read_samps(samples=samples)
    out_ = process_samples(infiles=ffs)
    return out_
