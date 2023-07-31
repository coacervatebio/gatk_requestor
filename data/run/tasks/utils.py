import os
import json
import asyncio
import subprocess
from pathlib import Path
from typing import List, Tuple
from time import sleep
from datetime import datetime
from dataclasses import dataclass
from dataclasses_json import dataclass_json
from flytekit import kwtypes, workflow, dynamic, task, ContainerTask, current_context
from flytekit.extras.tasks.shell import OutputLocation, ShellTask
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory
from flytekitplugins.pod import Pod
from yapapi import NoPaymentAccountError
from yapapi.log import enable_default_logger
from run.pod.yagna_template import yagna_requestor_ps, yagna_requestor
from run import config

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

@task(container_image=config.current_image)
def return_alignment(sample: str, reg: str, almt: FlyteFile, idx: FlyteFile) -> Alignment:
    return Alignment(sample=sample, reg=reg, almt=almt, idx=idx)

@task(container_image=config.current_image)
def dir_to_vcfs(indir: FlyteDirectory) -> List[VCF]:
    vcfs = {}

    # Fetch FlyteDirectory from object storage and make
    # list of relevant paths
    indir.download()
    all_paths = list(Path(indir.path).rglob('*g.vcf.gz*'))

    for fp in all_paths:
        
        # Parse paths following 'sample_chr.extension' format
        fn = fp.name
        sample = fn.split('.')[0]
        reg = sample.split('_')[-1]
        
        if not vcfs.get(sample):
            vcfs[sample] = VCF(
                sample=sample,
                reg=reg,
                vcf=FlyteFile(path='noop'),
                idx=FlyteFile(path='noop')
            )

        if 'g.vcf.gz.tbi' in fn:
            setattr(vcfs[sample], 'idx', FlyteFile(path=str(fp)))
        elif 'g.vcf.gz' in fn:
            setattr(vcfs[sample], 'vcf', FlyteFile(path=str(fp)))

    return list(vcfs.values())


@task(container_image=config.current_image)
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

@task(container_image=config.current_image)
def get_dir(dirpath: str) -> FlyteDirectory:
    fd = FlyteDirectory(path=dirpath)
    return fd

@task(container_image=config.current_image)
def prep_gather_vcfs(combi_dirs: List[FlyteDirectory]) -> Tuple[str, FlyteDirectory]:
    working_dir = current_context().working_directory
    out_dir = Path(os.path.join(working_dir, "outdir"))
    out_dir.mkdir(exist_ok=True)
    od = FlyteDirectory(path=str(out_dir))
    fnames = []
    for i in combi_dirs:
        i.download()
        for f in os.listdir(i):
            os.rename(os.path.join(i, f), os.path.join(out_dir, os.path.basename(f)))
            if '.tbi' not in f:
                fnames.append(f)
    fnames_fmt = ' '.join([f'-I {i}' for i in fnames])
    return fnames_fmt, od
    
@task(container_image=config.current_image)
def prep_db_import(vcf_objs: List[VCF], region: str) -> Tuple[str, FlyteDirectory]:
    working_dir = current_context().working_directory
    out_dir = Path(os.path.join(working_dir, "outdir"))
    out_dir.mkdir(exist_ok=True)

    vcf_names = []
    for o in [i for i in vcf_objs if i.reg == region]:
        o.vcf.download()
        o.idx.download()
        os.rename(o.vcf.path, os.path.join(out_dir, os.path.basename(o.vcf.path)))
        os.rename(o.idx.path, os.path.join(out_dir, os.path.basename(o.idx.path)))
        vcf_names.append(Path(o.vcf.path).name)

    # Preformat the VCF names for use in GenomicsDBImport
    vcf_names_fmt = ' '.join(f'-V {n}' for n in vcf_names)
    vcf_dir = FlyteDirectory(path=str(out_dir))
    
    return vcf_names_fmt, vcf_dir
    
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
    container_image=config.current_image,
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
    container_image=config.current_image,
    # task_config=Pod(pod_spec=yagna_requestor_ps)
)

t1 = ShellTask(
    name="task_1",
    container_image=config.current_image,
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

def run_golem(main, payloads, log_file=False):

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
    task = loop.create_task(main(payloads))

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

golem_test = ContainerTask(
    name="golem-test",
    input_data_dir="/mnt",
    output_data_dir="/var/outputs",
    inputs=kwtypes(),
    outputs=kwtypes(),
    image=config.current_image,
    pod_template = yagna_requestor,
    # pod_template_name = "my-pod-template", # Modify vols / svc here for yagna
    command=[
        "python",
        "/root/run/agents/hello_golem.py"
    ],
)

sleep_task = ContainerTask(
    name="basic-test",
    input_data_dir="/var/inputs",
    output_data_dir="/var/outputs",
    inputs=kwtypes(),
    outputs=kwtypes(),
    image="ghcr.io/flyteorg/rawcontainers-shell:v2",
    # pod_template = pt,
    # pod_template_name = "my-pod-template", # Modify vols / svc here for yagna
    command=[
        "sleep",
        "infinity"
    ],
)

basic_shell = ShellTask(
    name="basic-shell",
    debug=True,
    script=
    """
    echo hello
    echo world > {outputs.i}
    """,
    inputs=kwtypes(),
    output_locs=[OutputLocation(var="i", var_type=FlyteFile, location="basic_out.txt")],
    container_image=config.current_image
)