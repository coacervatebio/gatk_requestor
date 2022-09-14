from tests.common import ContainerTester, allowed_pattern
from tests.runners import GothRequestorRunner
from tests.checkers import VcfChecker
from tests.config import unit_assets_root


def test_goth_requestor():

    data_path = unit_assets_root.joinpath("call_variants")

    tester = ContainerTester(GothRequestorRunner(), VcfChecker(), data_path)

    tester.track_unexpected = False
    tester.target_files = filter(allowed_pattern, tester.target_files)

    tester.run()
