from tempfile import TemporaryDirectory
from pathlib import Path, PurePath
from common import ContainerTester
from runners import SnakemakeRunner


def test_index_cram():

    data_path = PurePath("assets/index_cram/")
    workdir = Path("assets/tmp_output/")

    tester = ContainerTester(data_path, workdir, SnakemakeRunner)
    tester.run()
    tester.check()
    tester.cleanup()


