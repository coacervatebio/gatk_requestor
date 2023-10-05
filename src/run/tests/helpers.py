import os
import filecmp
import pysam as ps
from time import sleep
from pathlib import Path
from typing import List
from flytekit import task, current_context
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory
from run.tasks.utils import VCF
from run import config, logger

class ComparisonCrawler:
    """
    A unified class for all comparison operations across dirs/files/dataclasses
    """

    def __init__(self, checker, actual_dir, expected_dir, include, ignore) -> None:

        checkers = {
            'simple': SimpleChecker(),
            'vcf': VcfChecker(),
            'cram': CramChecker()
        }

        self.checker = checkers[checker]
        self.actual_dir = actual_dir
        self.expected_dir = expected_dir
        self.include = include
        self.ignore = ignore

    def compare_directories(self):
        """
        Compare two directories and their subdirectories recursively.
        """
        dcmp = filecmp.dircmp(self.actual_dir, self.expected_dir, ignore=None)
        logger.debug(f'Made dircmp object between {self.actual_dir} and {self.expected_dir}')
        
        # Check if files/dirs in the current directory match
        if dcmp.left_only or dcmp.right_only:
            return False
        
        logger.debug("helloooooo")
        logger.debug("No uncommon files/dirs at the top level")

        # Filter files based on 'include' and 'ignore'
        files = dcmp.common_files
        if self.include:
            files = [f for f in files if f in self.include]
        if self.ignore:
            files = [f for f in files if f not in self.ignore]
        logger.debug(f'Files determined as {files} after filtering for include/ignore')

        logger.debug('Checking files at current level')
        for f in files:
            logger.debug(f'Currently comparing {f}')
            f1 = os.path.join(self.actual_dir, f)
            f2 = os.path.join(self.expected_dir, f)
            self.checker.compare_files(f1, f2)

        # Recursively compare subdirectories
        for sub_dir in dcmp.common_dirs:
            logger.debug(f'Checking {sub_dir}')
            sub_dir1 = os.path.join(self.actual_dir, sub_dir)
            sub_dir2 = os.path.join(self.expected_dir, sub_dir)
            if not self.compare_directories(sub_dir1, sub_dir2):
                return False
        
        return True


class SimpleChecker:
    """Byte for byte file comparison checker using `filecmp.cmp`"""

    @staticmethod
    def compare_files(actual_file, expected_file):
        logger.debug(f'Comparing {actual_file} to {expected_file}')
        # Check that cmp returns True (no difference between files)
        assert filecmp.cmp(actual_file, expected_file, shallow=False), \
            f'Deep comparison of {actual_file} and {expected_file} failed.'


class VcfChecker:
    """Compares VCFs based on records, excluding timestamped header"""

    @staticmethod
    def compare_files(actual_file, expected_file):
        logger.debug(f'Comparing {actual_file} to {expected_file}')
        av = ps.VariantFile(actual_file, "r")
        ev = ps.VariantFile(expected_file, "r")
        a_recs = [str(r) for r in av.fetch()]
        e_recs = [str(r) for r in ev.fetch()]
        assert a_recs == e_recs


class CramChecker:
    """Compares CRAMs based on records, ignoring compression differences"""

    @staticmethod
    def compare_files(actual_file, expected_file):
        logger.debug(f'Comparing {actual_file} to {expected_file}')
        aa = ps.AlignmentFile(actual_file, "rc", reference_filename=config['reference_location'])
        ea = ps.AlignmentFile(expected_file, "rc", reference_filename=config['reference_location'])
        a_recs = [str(r) for r in aa.fetch(until_eof=True)]
        e_recs = [str(r) for r in ea.fetch(until_eof=True)]
        assert a_recs == e_recs


@task(container_image=config['current_image'])
def compare_directories(actual: FlyteDirectory, expected: FlyteDirectory, checker: str, include: List[str], ignore: List[str]) -> bool:

    actual.download()
    expected.download()

    crawler = ComparisonCrawler(checker, actual.path, expected.path, include, ignore)

    return crawler.compare_directories()

@task(container_image=config['current_image'])
def compare_files(actual: FlyteFile, expected: FlyteFile, checker: str) -> bool:
    working_dir = current_context().working_directory

    actual_dir = Path(os.path.join(working_dir, "actual"))
    actual_dir.mkdir(exist_ok=True)
    actual.download()
    a_name = os.path.basename(actual.path)
    os.rename(actual.path, os.path.join(actual_dir, a_name))

    expected_dir = Path(os.path.join(working_dir, "expected"))
    expected_dir.mkdir(exist_ok=True)
    expected.download()
    e_name = os.path.basename(expected.path)
    os.rename(expected.path, os.path.join(expected_dir, e_name))

    crawler = ComparisonCrawler(checker, actual_dir, expected_dir, [a_name, e_name], [])

    return crawler.compare_directories()

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

@task(container_image=config['current_image'])
def sleep_task(duration: float):
    sleep(duration)