import os
from typing import Tuple, List

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


# t2 = ShellTask(
#     name="task_2",
#     debug=True,
#     script="""
#     set -ex
#     head {inputs.x} > {outputs.j}
#     """,
#     inputs=kwtypes(x=FlyteFile, y=FlyteDirectory),
#     output_locs=[
#         OutputLocation(var="j", var_type=FlyteFile, location="{inputs.y}/flyte_file_short.txt")
#     ],
# )


# t3 = ShellTask(
#     name="task_3",
#     debug=True,
#     script="""
#     set -ex
#     tar -zxvf {inputs.z}
#     cat {inputs.y}/$(basename {inputs.x}) | wc -m > {outputs.k}
#     """,
#     inputs=kwtypes(x=FlyteFile, y=FlyteDirectory, z=FlyteFile),
#     output_locs=[OutputLocation(var="k", var_type=FlyteFile, location="output.txt")],
# )
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

@task
def get_dir() -> List[FlyteFile]:
    files = []
    y = FlyteDirectory(path="my-s3-bucket/data/df/f20bf4f51aa334382a79-n0-0/89608308de764abbe0a6ee9a99dc9fa4")
    for file in [os.path.join(y, lcv) for lcv in sorted(os.listdir(y))]:
        files.append(file)
    return files

@task
def run_tasks(ins: List[FlyteFile]) -> str:
    for i in ins:
        t1(x=i)
    return "RUN TASKS DONE"

@workflow
def wf() -> str:
    files = get_dir()
    done = run_tasks(ins=files)
    # x = FlyteFile("s3://my-s3-bucket/input-data/HG03633_short.sam")
    # y = create_entities()
    # t1_out = t1(x=x)
    # t2_out = t2(x=x, y=y)
    # t2_out = t2(x=t1_out, y=y)
    # t3_out = t3(x=x, y=y, z=t2_out)
    # return t2_out
    # return t3_out
    return done

if __name__ == "__main__":
    print(f"Running wf() {wf()}")