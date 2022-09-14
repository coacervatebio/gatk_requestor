from pathlib import Path
from tests.common import ContainerTester, allowed_pattern
from tests.runners import SnakemakeRunner
from tests.checkers import VcfChecker
from tests.config import unit_assets_root

def test_genotype():

    data_path = unit_assets_root.joinpath("genotype")

    tmpdir = Path("/home/vagrant/tmp_guest")
    tester = ContainerTester(SnakemakeRunner(), VcfChecker(), data_path, tmpdir=tmpdir)
    tester.track_unexpected = False

    tester.target_files = filter(allowed_pattern, tester.target_files)

    tester.run()
