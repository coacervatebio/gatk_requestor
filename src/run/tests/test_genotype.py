from flytekit import workflow
from run.tasks.calling import genotype
from run.tasks.utils import get_dir
from run.tests.helpers import compare_directories
from run import config

@workflow
def test_genotype_wf():
    db_dir = get_dir(dirpath='s3://my-s3-bucket/test-assets/combine-region-expected')
    actual = genotype(db_dir=db_dir, reg='chr21', ref_loc=config['reference_location'])
    expected = get_dir(dirpath='s3://my-s3-bucket/test-assets/genotype-expected')
    compare_directories(actual=actual, expected=expected, checker='vcf', include=['combined_chr21.g.vcf.gz'], ignore=[])