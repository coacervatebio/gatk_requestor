from flytekit import workflow
from run.tasks.calling import golem_call_variants
from run.tasks.utils import dir_to_alignments, get_dir, dir_to_vcfs
from run.tests.helpers import compare_vcf_objs

@workflow
def test_golem_call_variants_wf():
    aldir = get_dir(dirpath='s3://my-s3-bucket/test-assets/indexed-alignments-per-region')
    als = dir_to_alignments(indir=aldir)
    called_actual = golem_call_variants(als=als)
    called_expected = dir_to_vcfs(indir='s3://my-s3-bucket/test-assets/haplotypecaller-expected')
    equivalent = compare_vcf_objs(actual=called_actual, expected=called_expected)
