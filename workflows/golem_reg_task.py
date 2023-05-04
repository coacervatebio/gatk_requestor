import os
from typing import List, Tuple
from flytekit import kwtypes, workflow, dynamic, task, ContainerTask
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory
from .pod_templates import yagna_requestor

@task(container_image='docker.io/coacervate/requestor:latest')
def hello_golem():
    print(os.listdir('/root'))

@workflow
def wf():
    hello_golem()