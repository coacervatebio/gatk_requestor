from typing import List, Tuple
from flytekit import kwtypes, workflow, dynamic, task, ContainerTask
from flytekit.types.file import FlyteFile
from flytekit.core.pod_template import PodTemplate

ic = ContainerTask(
    name="index_cram",
    input_data_dir="/var/inputs",
    output_data_dir="/var/outputs",
    inputs=kwtypes(al_in=FlyteFile),
    outputs=kwtypes(idx_out=FlyteFile),
    image="docker.io/coacervate/requestor:latest",
    pod_template = PodTemplate, # Modify vols / svc here for yagna
    command=[
        "samtools",
        "index",
        "/var/inputs/al_in",
        "-o",
        "/var/outputs/idx_out"
    ],
)

@workflow
def wf(samples: FlyteFile, regions: FlyteFile) -> str:
    
    ffs, regs = read_config(samples=samples, regions=regions)
    out_ = process_samples(infiles=ffs, regs=regs)
    return "YUS"