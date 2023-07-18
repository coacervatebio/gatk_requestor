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
from .utils import get_dir, run_golem, dir_to_alignments, Alignment, VCF, dir_to_vcfs, vcfs_to_dir
from .flyte_haplotypecaller import main
from .hello_golem import hello
from .tasks import golem_call_variants

combine_region = ShellTask(
    name="combine_region",
    debug=True,
    script="ls -lah {inputs.vdir}",
    # script="mkdir /tmp/genomics_db_dir",
    # """
    # gatk --java-options '-Xmx4g' GenomicsDBImport {params} -L {wildcards.reg} --genomicsdb-workspace-path {output}
    # """,
    inputs=kwtypes(vnames=List[str], vdir=FlyteDirectory),
    output_locs=[OutputLocation(var="i", var_type=FlyteDirectory, location="/tmp/genomics_db_dir")],
    container_image='docker.io/coacervate/requestor:latest'
)

@workflow
def wf():
    fd = get_dir(dirpath='s3://my-s3-bucket/data/rp/ask2j94d6nwv9hz2f5kd-n1-0-dn14-0')
    vcf_objs = dir_to_vcfs(indir=fd)
    vnames, vdir = vcfs_to_dir(vcf_objs=vcf_objs)
    combine_region(vnames=vnames, vdir=vdir)