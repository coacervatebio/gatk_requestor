from pathlib import PurePath, Path
from common import ContainerTester
from runners import SnakemakeRunner
from checkers import SimpleChecker


def test_call_variants():

    data_path = PurePath("assets/call_variants/")

    tester = ContainerTester(SnakemakeRunner(), SimpleChecker(), data_path)
    tester.run()