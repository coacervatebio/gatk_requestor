import os
from flytekit import kwtypes, workflow, dynamic, task, ContainerTask
# from .util_tasks import hg, get_file, get_file_contents
from flytekit.extras.tasks.shell import OutputLocation, ShellTask
from flytekitplugins.pod import Pod
from .pod_templates import yagna_requestor_ps

# @workflow
# def wf() -> str:
#     ff = get_file(filepath='s3://my-s3-bucket/input-data/regions.txt')
#     conts = get_file_contents(infile=ff)
#     return conts
    # hg()

@task(
    container_image='docker.io/coacervate/requestor:latest',
    task_config=Pod(pod_spec=yagna_requestor_ps)
    )
def local_ls_task():
    for d in os.listdir('/root'):
        print(d)

local_ls_shelltask = ShellTask(
    name="local_ls_shelltask",
    debug=True,
    script="""
    set -ex
    ls -la /root
    """,
    inputs=kwtypes(),
    # output_locs=[],
    container_image='docker.io/coacervate/requestor:latest',
    # task_config=Pod(pod_spec=yagna_requestor_ps)
)

@workflow
def wf():
    local_ls_task()
    local_ls_shelltask()