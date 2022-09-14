from tests.common import ContainerTester
from tests.runners import SnakemakeRunner
from tests.checkers import SimpleChecker
from tests.config import unit_assets_root

def test_split_cram():

    data_path = unit_assets_root.joinpath("split_cram")

    tester = ContainerTester(SnakemakeRunner(), SimpleChecker(), data_path)
    tester.run()
