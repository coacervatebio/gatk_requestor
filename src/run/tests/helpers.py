import os
import filecmp
import pysam as ps
from pathlib import Path
from typing import List
from flytekit import task, current_context
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory
from run.tasks.utils import VCF
from run import config, logger

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

def compare_vcfs(generated_file, expected_file):
    gv = ps.VariantFile(generated_file, "r")
    ev = ps.VariantFile(expected_file, "r")
    g_recs = [str(r) for r in gv.fetch()]
    e_recs = [str(r) for r in ev.fetch()]
    assert g_recs == e_recs

def compare_directories_local(dir1, dir2):
    """
    Compare two directories and their subdirectories recursively.
    
    Args:
        dir1 (str): Path to the first directory.
        dir2 (str): Path to the second directory.
        
    Returns:
        bool: True if the directories and their contents match exactly, False otherwise.
    """
    dcmp = filecmp.dircmp(dir1, dir2, ignore=None)
    logger.debug(f'Made dircmp object between {dir1} and {dir2}')
    
    # Check if files/dirs in the current directory match
    if dcmp.left_only or dcmp.right_only:
        return False
    logger.debug("No uncommon files/dirs at the top level")

    # Additional cmpfiles comparison because dcmp.diff_files does not allow shallow=False
    match, mismatch, errors = filecmp.cmpfiles(dir1, dir2, dcmp.common_files, shallow=False)
    # assert match, f'No files in {dcmp.common_files} matched between {dir1} and {dir2}'
    assert not mismatch, f'The following files do not match: {mismatch}'
    assert not errors, f'The following files threw errors when comparing: {errors}'

    # Recursively compare subdirectories
    for sub_dir in dcmp.common_dirs:
        logger.debug(f'Checking {sub_dir}')
        sub_dir1 = os.path.join(dir1, sub_dir)
        sub_dir2 = os.path.join(dir2, sub_dir)
        if not compare_directories_local(sub_dir1, sub_dir2):
            return False
    
    return True
    
@task(container_image=config['current_image'])
def compare_dirs(actual: FlyteDirectory, expected: FlyteDirectory) -> bool:
    actual.download()
    expected.download()
    assert compare_directories_local(actual, expected), f'Files between {actual.path} and {expected.path} or their subdirectories differ.'
    return True
    
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