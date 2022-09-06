from pathlib import PurePath, Path
from common import ContainerTester, allowed_pattern
from runners import SnakemakeRunner
from checkers import VcfChecker
from config import yagna_datadir


def test_call_variants():

    data_path = PurePath("assets/call_variants/")

    tester = ContainerTester(SnakemakeRunner(), VcfChecker(), data_path)
    tester.runner.arb_com = ['-y', 'on'] # Turns on yagna
    tester.runner.vols.append(f'{str(yagna_datadir)}:/home/requestor/.local/share/yagna') # Adds yagna datadir to avoid re-funding test GLM
    
    tester.track_unexpected = False
    tester.target_files = filter(allowed_pattern, tester.target_files)
    
    tester.run()