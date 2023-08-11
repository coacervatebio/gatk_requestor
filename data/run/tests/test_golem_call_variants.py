from flytekit import workflow
from run.tasks.calling import golem_call_variants
from run.tasks.utils import dir_to_alignments, compare_files, get_dir

@workflow
def test_golem_call_variants_wf():
    aldir = get_dir(dirpath=)
    als = dir_to_alignments(indir=aldir)
    idx_actual = index_cram(al=al)
    idx_expected = get_file(filepath='s3://my-s3-bucket/test-assets/HG03633_sub.cram.crai')
    equivalent = compare_files(actual=idx_actual, expected=idx_expected)