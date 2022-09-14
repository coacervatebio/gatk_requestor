from tests.common import ContainerTester, allowed_pattern
from tests.runners import SnakemakeRunner
from tests.checkers import VcfChecker
from tests.config import yagna_datadir, unit_assets_root


def test_call_variants():

    data_path = unit_assets_root.joinpath("call_variants")

    tester = ContainerTester(SnakemakeRunner(), VcfChecker(), data_path)
    tester.runner.arb_com = ["-y", "on"]  # Turns on yagna
    tester.runner.vols.append(
        f"{str(yagna_datadir)}:/home/coacervate/.local/share/yagna"
    )  # Adds yagna datadir to avoid re-funding test GLM

    tester.track_unexpected = False
    tester.target_files = filter(allowed_pattern, tester.target_files)

    tester.run()
