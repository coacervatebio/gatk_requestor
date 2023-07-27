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
from .tasks import golem_call_variants, combine_region
from .config import current_image

# gather_vcfs = ShellTask(
#     name="gather_vcfs",
#     debug=True,
#     script=
#     """
#     cd {inputs.vdir}
#     java -jar /usr/local/share/gatk GatherVcfs {params} -O {output}
#     """,
#     inputs=kwtypes(vnames_fmt=str, vdir=FlyteDirectory, reg=str),
#     output_locs=[OutputLocation(var="i", var_type=FlyteDirectory, location="/root/results/genomics_db_dir")],
#     container_image=current_image
# )

genotype = ShellTask(
    name="genotype",
    debug=True,
    script=
    """
    java -jar /usr/local/share/gatk GenotypeGVCFs -R /root/reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta -V gendb://{inputs.vdir} -O {outputs.v}
    """,
    inputs=kwtypes(vdir=FlyteDirectory, reg=str),
    output_locs=[
        OutputLocation(var="v", var_type=FlyteFile, location="/root/results/combined_{inputs.reg}.g.vcf.gz"),
        OutputLocation(var="i", var_type=FlyteFile, location="/root/results/combined_{inputs.reg}.g.vcf.gz.tbi")
    ],
    container_image=current_image
)


@workflow
def wf():
    fd = get_dir(dirpath='s3://my-s3-bucket/data/lo/ahv88sfgskx7b6f66mjt-n1-0-dn16-0/d3361e88316e93b137f6e07fb98e66ab')
    genotype(vdir=fd, reg='chr21')