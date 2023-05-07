import os
from flytekit import kwtypes, workflow, dynamic, task, ContainerTask
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory
from flytekit.extras.tasks.shell import OutputLocation, ShellTask
from flytekitplugins.pod import Pod
from .pod_templates import yagna_requestor_ps
from .util_tasks import get_dir, hg, get_file, get_file_contents

@task(
    container_image='docker.io/coacervate/requestor:latest',
    task_config=Pod(pod_spec=yagna_requestor_ps)
    )
def local_ls_task(indir: FlyteDirectory):
    print(os.listdir(indir))

@workflow
def wf():
    fd = get_dir(dirpath='s3://my-s3-bucket/input-data/alignments/chr22')
    local_ls_task(indir=fd)
