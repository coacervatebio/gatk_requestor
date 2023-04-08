import os
from typing import Tuple

import flytekit
from flytekit import kwtypes, task, workflow
from flytekit.extras.tasks.shell import OutputLocation, ShellTask
from flytekit.types.directory import FlyteDirectory
from flytekit.types.file import FlyteFile

t1 = ShellTask(
    name="task_1",
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
# Passthrough a file and add to it, use it in output locs as well


t2 = ShellTask(
    name="task_2",
    debug=True,
    script="""
    set -ex
    cp {inputs.x} {inputs.y}
    tar -zcvf {outputs.j} {inputs.y}
    """,
    inputs=kwtypes(x=FlyteFile, y=FlyteDirectory),
    output_locs=[
        OutputLocation(var="j", var_type=FlyteFile, location="{inputs.y}.tar.gz")
    ],
)


t3 = ShellTask(
    name="task_3",
    debug=True,
    script="""
    set -ex
    tar -zxvf {inputs.z}
    cat {inputs.y}/$(basename {inputs.x}) | wc -m > {outputs.k}
    """,
    inputs=kwtypes(x=FlyteFile, y=FlyteDirectory, z=FlyteFile),
    output_locs=[OutputLocation(var="k", var_type=FlyteFile, location="output.txt")],
)
# Use arbitrary name for an output and reference it in the `script` arg. Location 
# of output_loc treated like a variable that can be templated in.

# @task
# def create_entities() -> FlyteDirectory:
#     working_dir = flytekit.current_context().working_directory
#     flytedir = os.path.join(working_dir, "testdata")
#     os.makedirs(flytedir, exist_ok=True)
#     flytedir_file = os.path.join(flytedir, ".gitkeep")
#     os.open(flytedir_file, os.O_CREAT)
#     return flytedir

@workflow
def wf() -> FlyteFile:
    x = FlyteFile("s3://my-s3-bucket/input-data/HG03633_short.sam")
    y = FlyteDirectory("s3://my-s3-bucket/input-data/")
    # y = create_entities()
    t1_out = t1(x=x)
    t2_out = t2(x=t1_out, y=y)
    t3_out = t3(x=x, y=y, z=t2_out)
    return t3_out


if __name__ == "__main__":
    print(f"Running wf() {wf()}")