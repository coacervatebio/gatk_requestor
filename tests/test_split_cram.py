from pathlib import PurePath, Path
from common import ContainerTester
from runners import SnakemakeRunner
from checkers import SimpleChecker


def test_split_cram():

    data_path = PurePath("assets/split_cram/")

    tester = ContainerTester(SnakemakeRunner(), SimpleChecker(), data_path)
    tester.run()
