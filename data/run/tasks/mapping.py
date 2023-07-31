from flytekit import kwtypes, ContainerTask
from flytekit.extras.tasks.shell import OutputLocation, ShellTask
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory
from run import config

index_cram_sh = ShellTask(
    name="index_cram",
    debug=True,
    script=
    """
    samtools index {inputs.al} -o {outputs.idx}
    """,
    inputs=kwtypes(al=FlyteFile),
    output_locs=[OutputLocation(var="idx", var_type=FlyteFile, location="/root/results/output/{inputs.al}.crai")],
    container_image=config.current_image
)

index_cram_ct = ContainerTask(
    name="index_cram",
    input_data_dir="/var/inputs",
    output_data_dir="/var/outputs",
    inputs=kwtypes(al_in=FlyteFile),
    outputs=kwtypes(idx_out=FlyteFile),
    image=config.current_image,
    command=[
        "samtools",
        "index",
        "/var/inputs/al_in",
        "-o",
        "/var/outputs/idx_out"
    ],
)

split_cram = ContainerTask(
    name="split_cram",
    input_data_dir="/var/inputs",
    output_data_dir="/var/outputs",
    inputs=kwtypes(al_in=FlyteFile, idx_in=FlyteFile, reg=str),
    outputs=kwtypes(reg_out=FlyteFile),
    image=config.current_image,
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