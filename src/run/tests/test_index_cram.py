from flytekit import workflow
from run.tasks.mapping import index_cram
from run.tasks.utils import get_file
from run.tests.helpers import compare_files, sleep_task

@workflow
def test_index_cram_wf():
    al = get_file(filepath='s3://my-s3-bucket/test-assets/HG04149_sub.cram')
    idx_actual = index_cram(al=al)
    idx_expected = get_file(filepath='s3://my-s3-bucket/test-assets/HG04149_sub.cram.crai')
    compare_files(actual=idx_actual, expected=idx_expected, checker='simple')
    # sleep_task(duration=float(3600))