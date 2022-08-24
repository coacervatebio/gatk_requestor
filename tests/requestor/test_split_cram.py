from pathlib import PurePath
from common import ContainerTester
from runners import SnakemakeRunner


def test_split_cram():

    data_path = PurePath("assets/split_cram/")

    tester = ContainerTester(SnakemakeRunner, data_path)
    tester.run()
