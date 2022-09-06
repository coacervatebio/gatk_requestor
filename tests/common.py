"""
Common code for unit testing of rules generated with Snakemake 7.12.1.
"""

from cgi import test
import os
import shutil
import subprocess as sp
from pathlib import Path
from icecream import ic
from docker.errors import DockerException, ContainerError, APIError, NotFound
from config import test_name


class ContainerTester:
    def __init__(self, runner, checker, datadir, tmpdir=Path("assets/tmp_output/")):
        self.datadir = datadir
        self.tmpdir = tmpdir
        self.runner = runner
        self.checker = checker
        self.track_unexpected = True #can be set to False when ignoring inputs that change with every run

        # Setup necessary paths, datadir must have input dir (data) and output dir (expected)
        self.input_data = datadir.joinpath('inputs')
        self.expected_output = datadir.joinpath('expected')
        shutil.copytree(self.input_data, self.tmpdir)

        # Determine target files for container execution and comparison
        self.target_files = set(
            (Path(path) / f).relative_to(self.expected_output)
            for path, _, files in os.walk(self.expected_output)
            for f in files
        )
        # Sometimes target_str needs to be overwritten, e.g when it
        # is a dir producing many unspecified targets
        self.runner.target_string = " ".join([f'/{f}' for f in self.target_files])
        
        # Setup default mounts that can be overridden before calling run()
        self.runner.vols = [f"{str(self.tmpdir.joinpath('data', 'results').resolve())}:/data/results"]


    def check_files(self):
        input_files = set(
            (Path(path) / f).relative_to(self.input_data)
            for path, _, files in os.walk(self.input_data)
            for f in files
        )
        unexpected_files = set()
        all_files = []
        for path, _, files in os.walk(self.tmpdir):
            for f in files:
                f = (Path(path) / f).relative_to(self.tmpdir)
                all_files.append(f)
                # print("In focus: ", f)
                if str(f).startswith(".snakemake"):
                    continue
                if f in input_files:
                    # print("Ignoring: ", f)
                    # ignore input files
                    pass
                elif f in self.target_files:
                    print("Comparing: ", f)
                    self.checker.compare_files(self.tmpdir / f, self.expected_output / f)
                elif self.track_unexpected:
                    print("Unexpected: ", f)
                    unexpected_files.add(f)

        missing_targets = []
        for tf in self.target_files:
            if tf not in all_files:
                missing_targets.append(tf)

        if missing_targets:
            raise FileNotFoundError("Missing targets:\n{}".format(
                    "\n".join(sorted(map(str, missing_targets)))
                )
            )
        
        if unexpected_files:
            raise ValueError(
                "Unexpected files:\n{}".format(
                    "\n".join(sorted(map(str, unexpected_files)))
                )
            )

    def clean_con(self):
        try:
            test_con = self.runner.cons.get(test_name)
            test_con.remove(force=True)
        except NotFound:
            pass # Nothing to clean up

    def clean_tmp(self):
        if self.tmpdir.is_dir():
            shutil.rmtree(self.tmpdir)

    def run(self, run_con=True, check=True, clean_con=True, clean_tmp=True):
        try:
            if run_con: self.runner.run()
            if check: self.check_files()
        except ContainerError as ce:
            raise ce
        except APIError as ae:
        # Catches errors relating to getting logs from a container that never started
            raise ae
        # Catch any docker SDK issues
        except NotFound:
            raise
        finally:
            if clean_con: self.clean_con()
            if clean_tmp: self.clean_tmp()


# Small Utils

def allowed_pattern(tf):
    
    expected_cmp_fail_patterns = [
        '__',
        'vcfheader.vcf',
        'gz.tbi'
    ]
    
    s_tf = str(tf)
    for ecfp in expected_cmp_fail_patterns:
        if ecfp in s_tf:
            return False
    return True
