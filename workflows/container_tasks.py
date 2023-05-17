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


cv = ContainerTask(
    name="call_variants",
    input_data_dir="/var/inputs",
    output_data_dir="/var/outputs",
    inputs=kwtypes(als_in=FlyteDirectory),
    outputs=kwtypes(vcf_out=FlyteDirectory),
    image="docker.io/coacervate/requestor:latest",
    command=[
        "ls",
        "/var/inputs"
        # "python",
        # "/data/workflow/scripts/requestor.py",
    ],
)

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
        "/root/reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta",
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