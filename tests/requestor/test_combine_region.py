from pathlib import PurePath, Path
from common import ContainerTester
from runners import SnakemakeRunner


def allowed_pattern(tf):
    
    expected_cmp_fail_pattern = [
        '__',
        'vcfheader.vcf',
    ]
    
    s_tf = str(tf)
    for ecfp in expected_cmp_fail_pattern:
        if ecfp in s_tf:
            return False
    return True

def test_combine_region():
    # No snakemake rule explicitly states the creation of all the output files
    # Need to specify output from Snakefile but test for more
    data_path = PurePath("assets/combine_region/")
    tester = ContainerTester(SnakemakeRunner(), data_path)

    tester.track_unexpected = False
    tester.target_str = ['/data/results/combi_out/chr21_database']

    # Filter target files against failure patterns
    tester.target_files = filter(allowed_pattern, tester.target_files)

    tester.run()