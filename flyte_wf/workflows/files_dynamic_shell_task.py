from typing import List
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
def read_samps(samples: FlyteFile) -> List[FlyteFile]:
    lst = []
    with open(samples, 'r') as s_in:
        for i in s_in.readlines():
            ffs = f"s3://my-s3-bucket/input-data/{i}"
            lst.append(FlyteFile(path=ffs.strip()))
    return lst

@dynamic
def process_samples(infiles: List[FlyteFile]) -> str:

    for i in infiles:
        idx = ic(al_in=i)
        per_reg = sc(al_in=i, idx_in=idx, reg="chr21")

    return "PROCESSED"

@workflow
def wf(samples: FlyteFile) -> str:
    ffs = read_samps(samples=samples)
    out_ = process_samples(infiles=ffs)
    return out_
