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
    fnames = []
    for i in combi_dirs:
        i.download()
        # rename files to new dir
        # get fnames and format
        od = FlyteDirectory(path=out_dir.path)
        fnames_fmt = [f'-I {i}' for i in fnames]
    return fnames_fmt, od

@workflow
def wf():
    fd = get_dir(dirpath='s3://my-s3-bucket/input-data/geno_out')
    # vcf_objs = dir_to_vcfs(indir=fd)
    # vnames, vdir = prep_db_import(vcf_objs=vcf_objs, region='chr22')
    # combi = combine_region(vnames_fmt=vnames, vdir=vdir, reg='chr22')
    # genotype(vdir=combi, reg='chr22', refloc=reference_location)
    prep_gather_vcfs()
    gather_vcfs