import os
from datetime import datetime
from pathlib import Path
from flytekit import kwtypes, workflow, dynamic, task, ContainerTask,current_context
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory
from flytekit.extras.tasks.shell import OutputLocation, ShellTask
from flytekitplugins.pod import Pod
from .pod_templates import yagna_requestor_ps
from .util_tasks import get_dir, hg, get_file, get_file_contents
from .flyte_haplotypecaller import main
from .hello_golem import hello
from .golem_utils import run_golem

@task(
    container_image='docker.io/coacervate/requestor:latest',
    task_config=Pod(pod_spec=yagna_requestor_ps)
    )
def test_task(indir: FlyteDirectory) -> FlyteDirectory:
    working_dir = current_context().working_directory
    local_dir = Path(os.path.join(working_dir, "vcf_files"))
    local_dir.mkdir(exist_ok=True)
    local_indir = Path(indir)
    log_file = f"/tmp/haplotype_caller_{datetime.now().strftime('%Y%m%d-%H%M')}.log"
    run_golem(main(alpath=local_indir, vcfpath=local_dir), log_file=log_file)
    return local_dir
    # hello()

@workflow
def wf():
    fd = get_dir(dirpath='s3://my-s3-bucket/input-data/alignments/chr22')
    test_task(indir=fd)
