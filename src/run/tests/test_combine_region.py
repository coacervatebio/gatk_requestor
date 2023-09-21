# from flytekit import workflow
# from run.tasks.calling import combine_region
# from run.tasks.utils import get_dir, dir_to_vcfs
# from run.tests.helpers import compare_files

# @workflow
# def test_combine_region_wf():
#     vnames_fmt = '-V HG03633_sub_chr21.g.vcf.gz -V HG04149_sub_chr21.g.vcf.gz'
#     vdir = get_dir(dirpath='s3://my-s3-bucket/test-assets/combine-region-input')
#     actual = combine_region(vnames_fmt=vnames_fmt, vdir=vdir, reg='chr21')
#     expected = get_dir(dirpath='s3://my-s3-bucket/test-assets/combine-region-expected')
#     equivalent = compare_files(
#         actual=actual,
#         expected=expected,
#         to_compare=[
#             'callset.json',
#             'vidmap.json',
#             'vcfheader.vcf',
#             'chr21$1$46709983/genomicsdb_meta_dir/genomicsdb_column_bounds.json'
#         ]
#     )