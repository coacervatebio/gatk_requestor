from pathlib import PurePath
from tests.common import ContainerTester
from tests.runners import SnakemakeRunner
from tests.checkers import SimpleChecker


def test_split_cram():

    data_path = PurePath("assets/split_cram/")

    tester = ContainerTester(SnakemakeRunner(), SimpleChecker(), data_path)
    tester.run()
