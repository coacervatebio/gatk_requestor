import os
import filecmp
from pathlib import Path
from typing import List
from flytekit import task, current_context
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory
from run.tasks.utils import VCF
from run import config

@task(container_image=config['current_image'])
def compare_file(actual: FlyteFile, expected: FlyteFile) -> bool:
    actual.download()
    expected.download()
    return filecmp.cmp(actual.path, expected.path, shallow=False)
    
@task(container_image=config['current_image'])
def compare_files(actual: FlyteDirectory, expected: FlyteDirectory, to_compare: List[str]) -> bool:
    actual.download()
    expected.download()

    results = filecmp.cmpfiles(actual.path, expected.path, to_compare, shallow=False)

    # Return True if the 'different' and 'error' lists contain no files
    return len(results[1]) == 0 and len(results[2]) == 0
    
@task(container_image=config['current_image'])
def compare_dirs(actual: FlyteDirectory, expected: FlyteDirectory) -> bool:
    actual.download()
    expected.download()
    return len(filecmp.dircmp(actual, expected).diff_files) == 0
    
@task(container_image=config['current_image'])
def compare_vcf_objs(actual: List[VCF], expected: List[VCF]) -> bool:
    working_dir = current_context().working_directory

    actual_dir = Path(os.path.join(working_dir, "actual"))
    actual_dir.mkdir(exist_ok=True)
    for o in actual:
        o.vcf.download()
        os.rename(o.vcf.path, os.path.join(actual_dir, os.path.basename(o.vcf.path)))
        o.idx.download()
        os.rename(o.idx.path, os.path.join(actual_dir, os.path.basename(o.idx.path)))


    expected_dir = Path(os.path.join(working_dir, "expected"))
    expected_dir.mkdir(exist_ok=True)
    for o in expected:
        o.vcf.download()
        os.rename(o.vcf.path, os.path.join(expected_dir, os.path.basename(o.vcf.path)))
        o.idx.download()
        os.rename(o.idx.path, os.path.join(expected_dir, os.path.basename(o.idx.path)))

    return len(filecmp.dircmp(actual_dir, expected_dir).diff_files) == 0