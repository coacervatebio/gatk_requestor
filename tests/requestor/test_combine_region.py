from pathlib import PurePath
from common import ContainerTester, allowed_pattern
from runners import SnakemakeRunner
from checkers import SimpleChecker

def test_combine_region():
    # No snakemake rule explicitly states the creation of all the output files
    # Need to specify output from Snakefile but test for more
    data_path = PurePath("assets/combine_region/")
    tester = ContainerTester(SnakemakeRunner(), SimpleChecker(), data_path)

    tester.track_unexpected = False
    tester.target_str = ['/data/results/combi_out/chr21_database']

    # Filter target files against failure patterns
    tester.target_files = filter(allowed_pattern, tester.target_files)

    tester.run()