from flytekit import workflow
from run.tasks.mapping import index_cram
from run.tasks.utils import get_file, compare_files

@workflow
def wf():
    al = get_file(filepath='s3://my-s3-bucket/test-assets/HG03633_sub.cram')
    idx_actual = index_cram(al=al)
    idx_expected = get_file(filepath='s3://my-s3-bucket/test-assets/HG03633_sub.cram.crai')
    equivalent = compare_files(actual=idx_actual, expected=idx_expected)
