import os
from datetime import datetime
from pathlib import PurePath, Path
from typing import List, Tuple
from flytekit import kwtypes, workflow, dynamic, task, ContainerTask, current_context
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory
from flytekit.extras.tasks.shell import OutputLocation, ShellTask
from flytekitplugins.pod import Pod
from .pod_templates import yagna_requestor_ps
from .utils import get_dir, run_golem, dir_to_alignments, Alignment, VCF, dir_to_vcfs, prep_db_import
from .flyte_haplotypecaller import main
from .hello_golem import hello
from .tasks import golem_call_variants, combine_region, genotype
from .config import current_image, reference_location

gather_vcfs = ShellTask(
    name="gather_vcfs",
    debug=True,
    script=
    """
    cd {inputs.vdir}
    java -jar /usr/local/share/gatk GatherVcfs {inputs.vnames_fmt} -O {outputs.i}
    """,
    inputs=kwtypes(vnames_fmt=str, vdir=FlyteDirectory),
    output_locs=[OutputLocation(var="i", var_type=FlyteFile, location="/root/results/gather_out.vcf.gz")],
    container_image=current_image
)

@task(container_image=current_image)
def prep_gather_vcfs(combi_dirs: List[FlyteDirectory]) -> Tuple[str, FlyteDirectory]:
    working_dir = current_context().working_directory
    out_dir = Path(os.path.join(working_dir, "outdir"))
    out_dir.mkdir(exist_ok=True)
    od = FlyteDirectory(path=str(out_dir))
    fnames = []
    for i in combi_dirs:
        i.download()
        for f in os.listdir(i):
            os.rename(os.path.join(i, f), os.path.join(out_dir, os.path.basename(f)))
            if '.tbi' not in f:
                fnames.append(f)
    fnames_fmt = ' '.join([f'-I {i}' for i in fnames])
    return fnames_fmt, od

@workflow
def wf():
    d1 = get_dir(dirpath='s3://my-s3-bucket/input-data/geno_out/chr21')
    d2 = get_dir(dirpath='s3://my-s3-bucket/input-data/geno_out/chr22')
    fmt, fd = prep_gather_vcfs(combi_dirs=[d1, d2])
    gather_vcfs(vnames_fmt=fmt, vdir=fd)