from pathlib import PurePath
from common import ContainerTester
from runners import SnakemakeRunner


def test_index_cram():

    data_path = PurePath("assets/index_cram/")

    tester = ContainerTester(SnakemakeRunner(), data_path)
    tester.run()
