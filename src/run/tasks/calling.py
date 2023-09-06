import os
from pathlib import Path, PurePath
from typing import List, Tuple
from flytekit import kwtypes, workflow, dynamic, task, ContainerTask
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory
from flytekit.extras.tasks.shell import OutputLocation, ShellTask
from flytekitplugins.pod import Pod
from run.pod.yagna_template import yagna_requestor_ps
from run.tasks.utils import VCF, Alignment, run_golem
from run.agents.haplotypecaller import main
from run import config

gather_vcfs = ShellTask(
    name="gather_vcfs",
    debug=True,
    script=
    """
    cd {inputs.vdir}
    java -jar /usr/local/share/gatk GatherVcfs {inputs.vnames_fmt} -O {outputs.o}
    """,
    inputs=kwtypes(vnames_fmt=str, vdir=FlyteDirectory),
    output_locs=[OutputLocation(var="o", var_type=FlyteFile, location=f"{config['output_dir']}/gather_out.vcf.gz")],
    container_image=config['current_image']
)

genotype = ShellTask(
    name="genotype",
    debug=True,
    script=
    """
    java -jar /usr/local/share/gatk GenotypeGVCFs -R {inputs.ref_loc} -V gendb://{inputs.db_dir} -O {outputs.o}/combined_{inputs.reg}.g.vcf.gz
    """,
    inputs=kwtypes(db_dir=FlyteDirectory, reg=str, ref_loc=str),
    output_locs=[OutputLocation(var="o", var_type=FlyteDirectory, location=config['output_dir'])],
    container_image=config['current_image']
)

combine_region = ShellTask(
    name="combine_region",
    debug=True,
    script=
    """
    cd {inputs.vdir}
    java -jar /usr/local/share/gatk GenomicsDBImport {inputs.vnames_fmt} -L {inputs.reg} --genomicsdb-workspace-path {outputs.o}
    """,
    inputs=kwtypes(vnames_fmt=str, vdir=FlyteDirectory, reg=str),
    output_locs=[OutputLocation(var="o", var_type=FlyteDirectory, location=f"{config['output_dir']}/genomics_db_dir")],
    container_image=config['current_image']
)

@task(
    container_image=config['current_image'],
    task_config=Pod(pod_spec=yagna_requestor_ps)
    )
def golem_call_variants(als: List[Alignment]) -> List[VCF]:
    payloads = []
    vcfs = []
    for al in als:
        working_dir = PurePath(al.almt.path).parent

        # Download alignment and index from object store to pod for uploading to Golem
        al.almt.download()
        al.idx.download()

        # Prepare payload for passing to requestor agent
        payload = {
            "sample": al.sample,
            "region_str": al.reg,
            "working_dir": working_dir,
            "req_align_path": Path(al.almt.path),
            "req_align_index_path": Path(al.idx.path),
            "req_vcf_path": working_dir.joinpath(f"{al.sample}_{al.reg}.g.vcf.gz"),
            "req_vcf_index_path": working_dir.joinpath(f"{al.sample}_{al.reg}.g.vcf.gz.tbi"),
        }

        # Prepare VCF based on output paths
        vcf = VCF(
            sample=al.sample,
            reg=al.reg,
            vcf=FlyteFile(path=str(payload['req_vcf_path'])),
            idx=FlyteFile(path=str(payload['req_vcf_index_path'])),
        )
        
        vcfs.append(vcf)
        payloads.append(payload)

    # Call requestor agent entrypoint with payloads
    run_golem(main, payloads)
    return vcfs