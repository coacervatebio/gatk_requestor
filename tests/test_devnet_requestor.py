from pathlib import PurePath
from tests.common import ContainerTester, allowed_pattern
from tests.runners import DevnetRequestorRunner
from tests.checkers import VcfChecker

def test_devnet_requestor():

    data_path = PurePath("assets/call_variants/")

    tester = ContainerTester(DevnetRequestorRunner(), VcfChecker(), data_path)
    
    tester.track_unexpected = False
    tester.target_files = filter(allowed_pattern, tester.target_files)
    
    tester.run()