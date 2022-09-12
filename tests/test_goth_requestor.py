import docker
from pathlib import PurePath
from tests.common import ContainerTester, allowed_pattern
from tests.runners import GothRequestorRunner, DevnetRequestorRunner
from tests.checkers import VcfChecker


def test_goth_requestor():

    data_path = PurePath("assets/call_variants/")

    tester = ContainerTester(GothRequestorRunner(), VcfChecker(), data_path)

    tester.track_unexpected = False
    tester.target_files = filter(allowed_pattern, tester.target_files)

    tester.run()
