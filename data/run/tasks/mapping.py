from flytekit import kwtypes, ContainerTask, task
from flytekit.extras.tasks.shell import OutputLocation, ShellTask
from flytekit.types.file import FlyteFile
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
    inputs=kwtypes(al=FlyteFile, idx=FlyteFile, sample=str, reg=str, ref_loc=str),
    output_locs=[OutputLocation(var="reg_al", var_type=FlyteFile, location="/root/results/output/{inputs.sample}_{inputs.reg}.cram")],
    container_image=config['current_image']
)
