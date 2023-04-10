from typing import List, Tuple
from flytekit import kwtypes, workflow, dynamic, task, ContainerTask
from flytekit.types.file import FlyteFile

ic = ContainerTask(
    name="index_cram",
    input_data_dir="/var/inputs",
    output_data_dir="/var/outputs",
    inputs=kwtypes(al_in=FlyteFile),
    outputs=kwtypes(idx_out=FlyteFile),
    image="docker.io/coacervate/requestor:latest",
    command=[
        "samtools",
        "index",
        "/var/inputs/al_in",
        "-o",
        "/var/outputs/idx_out"
    ],
)

sc = ContainerTask(
    name="split_cram",
    input_data_dir="/var/inputs",
    output_data_dir="/var/outputs",
    inputs=kwtypes(al_in=FlyteFile, idx_in=FlyteFile, reg=str),
    outputs=kwtypes(reg_out=FlyteFile),
    image="docker.io/coacervate/requestor:latest",
    command=[
        "samtools",
        "view",
        "-T",
        "/data/resources/reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta",
        "-O",
        "cram",
        "-o",
        "/var/outputs/reg_out",
        "-X",
        "/var/inputs/al_in",
        "/var/inputs/idx_in",
        "{{.inputs.reg}}",
    ],
)

@task
def read_config(samples: FlyteFile, regions: FlyteFile) -> Tuple[List[FlyteFile], List[str]]:
    infiles = []
    with open(samples, 'r') as s_in:
        for i in s_in.readlines():
            ffs = f"s3://my-s3-bucket/input-data/{i}"
            infiles.append(FlyteFile(path=ffs.strip()))

    regs = []
    with open(regions, 'r') as r_in:
        regs = [r.strip() for r in r_in.readlines()]

    return infiles, regs

@dynamic
def process_samples(infiles: List[FlyteFile], regs: List[str]) -> str:

    for i in infiles:
        idx = ic(al_in=i)
        for r in regs:
            per_reg = sc(al_in=i, idx_in=idx, reg=r)

    return "PROCESSED"

@workflow
def wf(samples: FlyteFile, regions: FlyteFile) -> str:
    
    ffs, regs = read_config(samples=samples, regions=regions)
    out_ = process_samples(infiles=ffs, regs=regs)
    return "YUS"
