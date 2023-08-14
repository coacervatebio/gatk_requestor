from flytekit import workflow
from run.tasks.calling import genotype
from run.tasks.utils import get_dir, dir_to_vcfs
from run.tests.helpers import compare_dirs
from run import config

@workflow
def test_genotype_wf():
    db_dir = get_dir(dirpath='s3://my-s3-bucket/test-assets/combine-region-expected')
    actual = genotype(db_dir=db_dir, reg='chr21', ref_loc=config['reference_location'])
    expected = get_dir(dirpath='s3://my-s3-bucket/test-assets/genotype-expected')
    equivalent = compare_dirs(actual=actual, expected=expected)