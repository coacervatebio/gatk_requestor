from pathlib import PurePath, Path
from common import ContainerTester
from runners import SnakemakeRunner
from checkers import VcfChecker


def test_call_variants():

    data_path = PurePath("assets/call_variants/")

    tester = ContainerTester(SnakemakeRunner(), VcfChecker(), data_path)
    tester.run()