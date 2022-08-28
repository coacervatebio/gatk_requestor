from pathlib import PurePath, Path
from common import ContainerTester
from runners import SnakemakeRunner
from checkers import VcfChecker


def test_gather_vcfs():

    data_path = PurePath("assets/gather_vcfs/")

    tester = ContainerTester(SnakemakeRunner(), VcfChecker(), data_path)
    tester.run()