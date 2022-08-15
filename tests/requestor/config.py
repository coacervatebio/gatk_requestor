import os
import rootpath
from pathlib import Path

test_tag = 'coacervate_requestor:test'
rpath = Path(rootpath.detect())
yagna_datadir = Path('/home/vagrant/yagna_datadir/')

# Create mount points for input and output to pass to container
almt_dir = Path.joinpath(rpath, 'tests', 'requestor', 'assets', 'alignments')
tmp_out_dir = Path.joinpath(rpath, 'tests', 'requestor', 'assets', 'tmp_output')


# Get names of input samples dynamically
per_chr_samples = []
for reg_dir in almt_dir.glob('chr*'):
    for almt in reg_dir.glob('*cram'):
        per_chr_samples.append(almt.with_suffix('').name)

# Create expected output sample paths
expected_hc_out = []
for sample in per_chr_samples:
    expected_hc_out.append(f"{tmp_out_dir.joinpath(sample)}.g.vcf.gz")
    expected_hc_out.append(f"{tmp_out_dir.joinpath(sample)}.g.vcf.gz.tbi")