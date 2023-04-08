import os
from typing import List
from flytekit import kwtypes, workflow, dynamic
from flytekit.extras.tasks.shell import OutputLocation, ShellTask
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory

s1 = ShellTask(
    name="shorten",
    debug=True,
    script="""
    set -ex
    echo mooooooore >> {inputs.infile}
    """,
    inputs=kwtypes(infile=FlyteFile),
    output_locs=[OutputLocation(var="i", var_type=FlyteFile, location="{inputs.infile}")],
)

@dynamic
def dyn_task(loops: int, infile: FlyteFile) -> str:

    for i in range(loops):
        s1(infile=infile)
    
    return "DYN_DONE"

@workflow
def wf(infile: FlyteFile) -> str:

    return dyn_task(loops=8, infile=infile)


if __name__ == "__main__":
    print(f"Running wf() {wf()}")