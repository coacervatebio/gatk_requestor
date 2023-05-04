import os
from typing import List, Tuple
from flytekit import kwtypes, workflow, dynamic, task, ContainerTask
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory
from .pod_templates import yagna_requestor

gt = ContainerTask(
    name="golem-test",
    input_data_dir="/mnt",
    output_data_dir="/var/outputs",
    inputs=kwtypes(),
    outputs=kwtypes(),
    image="docker.io/coacervate/requestor:latest",
    pod_template = yagna_requestor,
    # pod_template_name = "my-pod-template", # Modify vols / svc here for yagna
    command=[
        # "find",
        # "/",
        # "-name",
        # "*cram*"
        # "ls",
        # "-la",
        # "/mnt/indir",
        "python",
        "/data/agents/hello_golem.py"
        # "sleep",
        # "infinity"
    ],
)

bt = ContainerTask(
    name="basic-test",
    input_data_dir="/var/inputs",
    output_data_dir="/var/outputs",
    inputs=kwtypes(indir=FlyteDirectory),
    outputs=kwtypes(),
    image="ghcr.io/flyteorg/rawcontainers-shell:v2",
    # pod_template = pt,
    # pod_template_name = "my-pod-template", # Modify vols / svc here for yagna
    command=[
        # "ls",
        # "-la",
        # "/var/inputs",
        "sleep",
        "infinity"
    ],
)

@task
def get_dir(dirpath: str) -> FlyteDirectory:
    fd = FlyteDirectory(path=dirpath)
    return fd

@task
def get_files(fd: FlyteDirectory) -> List[FlyteFile]:
    infiles = [FlyteFile(os.path.join(fd, i)) for i in os.listdir(fd)]
    return infiles

@workflow
def wf():
    # fd = get_dir(dirpath='s3://my-s3-bucket/cv-in')
    # infiles = get_files(fd=fd)
    gt()