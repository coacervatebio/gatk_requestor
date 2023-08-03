from flytekit import kwtypes, ContainerTask, task
from flytekit.extras.tasks.shell import OutputLocation, ShellTask
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory
from run import config

index_cram = ShellTask(
    name="index_cram",
    debug=True,
    script=
    """
    samtools index {inputs.al} -o {outputs.idx}
    """,
    inputs=kwtypes(al=FlyteFile),
    output_locs=[OutputLocation(var="idx", var_type=FlyteFile, location="{inputs.al}.crai")],
    container_image=config['current_image']
)

split_cram = ShellTask(
    name="split_cram",
    debug=True,
    script=
    """
    samtools view -T {inputs.ref_loc} \
    -O cram -o {outputs.reg_al} -X {inputs.al} {inputs.idx} {inputs.reg}
    """,
    inputs=kwtypes(al=FlyteFile, idx=FlyteFile, reg=str, ref_loc=str),
    output_locs=[OutputLocation(var="reg_al", var_type=FlyteFile, location="{inputs.al}_{inputs.reg}.cram")],
    container_image=config['current_image']
)

split_cram_ct = ContainerTask(
    name="split_cram_ct",
    input_data_dir="/var/inputs",
    output_data_dir="/var/outputs",
    inputs=kwtypes(al_in=FlyteFile, idx_in=FlyteFile, reg=str),
    outputs=kwtypes(reg_out=FlyteFile),
    image=config['current_image'],
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