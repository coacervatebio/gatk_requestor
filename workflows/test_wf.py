import os
from datetime import datetime
from pathlib import PurePath, Path
from typing import List
from flytekit import kwtypes, workflow, dynamic, task, ContainerTask, current_context
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory
from flytekit.extras.tasks.shell import OutputLocation, ShellTask
from flytekitplugins.pod import Pod
from .pod_templates import yagna_requestor_ps
from .utils import get_dir, run_golem, dir_to_alignments, Alignment, VCF, dir_to_vcfs, prep_db_import
from .flyte_haplotypecaller import main
from .hello_golem import hello
from .tasks import golem_call_variants

combine_region = ShellTask(
    name="combine_region",
    debug=True,
    script=
    """
    mkdir /root/results/genomics_db_dir
    cd {inputs.vdir}
    gatk --java-options '-Xmx4g' GenomicsDBImport {inputs.vnames_fmt} -L {inputs.reg} --genomicsdb-workspace-path {outputs.i}
    """,
    inputs=kwtypes(vnames_fmt=str, vdir=FlyteDirectory, reg=str),
    output_locs=[OutputLocation(var="i", var_type=FlyteDirectory, location="/root/results/genomics_db_dir")],
    container_image='docker.io/coacervate/requestor:latest'
)

@workflow
def wf():
    fd = get_dir(dirpath='s3://my-s3-bucket/data/rp/ask2j94d6nwv9hz2f5kd-n1-0-dn14-0')
    vcf_objs = dir_to_vcfs(indir=fd)
    vnames, vdir = prep_db_import(vcf_objs=vcf_objs, region='chr21')
    combine_region(vnames_fmt=vnames, vdir=vdir, reg='chr21')