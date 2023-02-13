from tests.common import ContainerTester
from tests.runners import SnakemakeRunner
from tests.checkers import VcfChecker
from tests.config import unit_assets_root


def test_gather_vcfs():

    data_path = unit_assets_root.joinpath("gather_vcfs")

    tester = ContainerTester(SnakemakeRunner(), VcfChecker(), data_path)
    tester.run()
