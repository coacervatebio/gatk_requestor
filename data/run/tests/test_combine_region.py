from flytekit import workflow
from run.tasks.calling import combine_region
from run.tasks.utils import get_dir, dir_to_vcfs

@workflow
def test_golem_call_variants_wf():
    vnames_fmt = '-V HG03633_sub_chr21.g.vcf.gz -V HG04149_sub_chr21.g.vcf.gz'
    vdir = get_dir(dirpath='s3://my-s3-bucket/data/tw/alsqh7dsctc5s6v9w5hz-n1-0-dn15-0/ba734801b6de4b463a2fe53a421ad406')
    actual = combine_region(vnames_fmt=vnames_fmt, vdir=vdir, reg='chr21')
    expected = get_dir(dirpath='s3://my-s3-bucket/test-assets/combine-region-expected')