from pathlib import PurePath, Path
from tests.common import ContainerTester
from tests.runners import SnakemakeRunner
from tests.checkers import VcfChecker


def test_gather_vcfs():

    data_path = PurePath("assets/gather_vcfs/")

    tester = ContainerTester(SnakemakeRunner(), VcfChecker(), data_path)
    tester.run()
