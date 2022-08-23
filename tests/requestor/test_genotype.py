from pathlib import PurePath
from common import ContainerTester
from runners import SnakemakeRunner

def test_genotype():

    data_path = PurePath("assets/genotype/")
    tester = ContainerTester(SnakemakeRunner, data_path)
    tester.run_defaults()
