import os
import json
import subprocess
from typing import List, Tuple
from flytekit import kwtypes, workflow, dynamic, task, ContainerTask
from flytekit.extras.tasks.shell import OutputLocation, ShellTask
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory
from flytekitplugins.pod import Pod
from .pod_templates import yagna_requestor_ps


@task
def get_dir(dirpath: str) -> FlyteDirectory:
    fd = FlyteDirectory(path=dirpath)
    return fd
    
@task
def get_file(filepath: str) -> FlyteFile:
    ff = FlyteFile(path=filepath)
    return ff

@task
def get_files(fd: FlyteDirectory) -> List[FlyteFile]:
    infiles = [FlyteFile(os.path.join(fd, i)) for i in os.listdir(fd)]
    return infiles

@task
def get_file_contents(infile: FlyteFile) -> str:
    content = ''
    with open(infile, 'r') as in_:
        content = in_.read()
    return content

# @task(
#     container_image='docker.io/coacervate/requestor:latest',
#     task_config=Pod(pod_spec=yagna_requestor_ps)
# )
# def get_golem_appkey() -> str:
#     key_list = subprocess.run(["yagna", "app-key", "list", "--json"], capture_output=True)
#     key = json.loads(key_list.stdout)[0].get('key')
#     print(key)
#     return key


hg = ShellTask(
    name="hello_golem",
    debug=True,
    script="""
    set -ex
    python /root/agents/hello_golem.py
    """,
    inputs=kwtypes(),
    # output_locs=[],
    container_image='docker.io/coacervate/requestor:latest',
    # task_config=Pod(pod_spec=yagna_requestor_ps)
)

t1 = ShellTask(
    name="task_1",
    container_image='docker.io/coacervate/requestor:latest',
    debug=True,
    script="""
    set -ex
    echo "Hey there! Let's run some bash scripts using Flyte's ShellTask."
    echo "Showcasing Flyte's Shell Task." >> {inputs.x}
    if grep "Flyte" {inputs.x}
    then
        echo "Found it!" >> {inputs.x}
    else
        echo "Not found!"
    fi
    """,
    inputs=kwtypes(x=FlyteFile),
    output_locs=[OutputLocation(var="i", var_type=FlyteFile, location="{inputs.x}")],
)