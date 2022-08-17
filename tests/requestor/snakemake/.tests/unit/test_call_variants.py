import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_call_variants():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/call_variants/data")
        expected_path = PurePosixPath(".tests/unit/call_variants/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # dbg
        print("/mnt/results/hc_out/chr21/HG03633_chr21.g.vcf.gz /mnt/results/hc_out/chr21/HG03633_chr21.g.vcf.gz.tbi", file=sys.stderr)

        # Run the test job.
        sp.check_output([
            "python",
            "-m",
            "snakemake", 
            "/mnt/results/hc_out/chr21/HG03633_chr21.g.vcf.gz","/mnt/results/hc_out/chr21/HG03633_chr21.g.vcf.gz.tbi",
            "-f", 
            "-j1",
            "--keep-target-files",
            "-s=/mnt/workflow/Snakefile",
            "--config",
            "golem_subnet=goth",
            "--directory",
            workdir,
        ])

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file), 
        # also see common.py.
        common.OutputChecker(data_path, expected_path, workdir).check()
