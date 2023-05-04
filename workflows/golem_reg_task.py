import os
import json
import subprocess
from typing import List, Tuple
from flytekit import kwtypes, workflow, dynamic, task, ContainerTask
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory
from flytekitplugins.pod import Pod
from .pod_templates import yagna_requestor_ps

@task(
    container_image='docker.io/coacervate/requestor:latest',
    task_config=Pod(pod_spec=yagna_requestor_ps)
)
def hello_golem():
    key_list = subprocess.run(["yagna", "app-key", "list", "--json"], capture_output=True)
    print(json.loads(key_list.stdout)[0].get('key'))

@workflow
def wf():
    hello_golem()