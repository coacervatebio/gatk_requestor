import os
import json
import asyncio
import subprocess
from typing import List, Tuple
from time import sleep
from datetime import datetime
from dataclasses import dataclass
from dataclasses_json import dataclass_json
from flytekit import kwtypes, workflow, dynamic, task, ContainerTask
from flytekit.extras.tasks.shell import OutputLocation, ShellTask
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory
from flytekitplugins.pod import Pod
from yapapi import NoPaymentAccountError
from yapapi.log import enable_default_logger
from .pod_templates import yagna_requestor_ps
from .hello_golem import hello

@dataclass_json
@dataclass
class Alignment:
    sample: str
    reg: str
    almt: FlyteFile
    idx: FlyteFile
    
@dataclass_json
@dataclass
class VCF:
    sample: str
    reg: str
    vcf: FlyteFile
    idx: FlyteFile

@task(container_image='docker.io/coacervate/requestor:latest')
def dir_to_alignments(indir: FlyteDirectory) -> List[Alignment]:
    samps = set()
    for fn in os.listdir(indir):
        sample = fn.split('.')[0]
        reg = sample.split('_')[1]
        samps.add((sample, reg))

    als = []
    for s, r in samps:
        al = Alignment(
            sample=s,
            reg=r,
            almt = os.path.join(indir, f'{s}.cram'),
            idx = os.path.join(indir, f'{s}.cram.crai')
        )
        als.append(al)
    return als

@task(container_image='docker.io/coacervate/requestor:latest')
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

@task(
    container_image='docker.io/coacervate/requestor:latest',
    task_config=Pod(pod_spec=yagna_requestor_ps)
)
def get_golem_appkey() -> str:
    key_list = subprocess.run(["yagna", "app-key", "list", "--json"], capture_output=True)
    key = json.loads(key_list.stdout)[0].get('key')
    print(key)
    return key

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

def run_golem(main, log_file=False):

    if log_file:
        log_file = f"/tmp/haplotype_caller_{datetime.now().strftime('%Y%m%d-%H%M')}.log"
        enable_default_logger(
            log_file=log_file,
            debug_activity_api=True,
            debug_market_api=True,
            debug_payment_api=True,
            debug_net_api=True,
        )

    # Set app key
    while os.getenv('YAGNA_APPKEY') is None:
        key_list = subprocess.run(["yagna", "app-key", "list", "--json"], capture_output=True)
        os.environ['YAGNA_APPKEY'] = json.loads(key_list.stdout)[0].get('key')
        sleep(5)

    loop = asyncio.get_event_loop()
    task = loop.create_task(main())

    try:
        loop.run_until_complete(task)
    except NoPaymentAccountError as e:
        handbook_url = (
            "https://handbook.golem.network/requestor-tutorials/"
            "flash-tutorial-of-requestor-development"
        )
        print(
            f"No payment account initialized for driver `{e.required_driver}` "
            f"and network `{e.required_network}`.\n\n"
            f"See {handbook_url} on how to initialize payment accounts for a requestor node."
        )
        task.cancel()