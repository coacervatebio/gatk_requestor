from flytekit import workflow
from run.tasks.calling import gather_vcfs
from run.tasks.utils import get_dir, get_file
from run.tests.helpers import compare_file
from run import config

@workflow
def test_gather_vcfs_wf():
    gather_dir = get_dir(dirpath='s3://my-s3-bucket/test-assets/gather-vcfs-input')
    actual = gather_vcfs(vnames_fmt='-I combined_chr21.g.vcf.gz -I combined_chr22.g.vcf.gz', vdir=gather_dir)
    expected = get_file(filepath='s3://my-s3-bucket/test-assets/gather_out.vcf.gz')
    equivalent = compare_file(actual=actual, expected=expected)