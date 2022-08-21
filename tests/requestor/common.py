"""
Common code for unit testing of rules generated with Snakemake 7.12.1.
"""

from pathlib import Path
import subprocess as sp
import os
import shutil


class ContainerTester:
    def __init__(self, datadir, tmpdir, runner):
        self.datadir = datadir
        self.tmpdir = tmpdir
        self.runner = runner

        # Setup necessary paths, datadir must have input dir (data) and output dir (expected)
        self.input_data = datadir.joinpath('data')
        self.expected_output = datadir.joinpath('expected')
        shutil.copytree(self.input_data, self.tmpdir)

        # Determine target files for container execution and comparison
        self.target_files = set(
            (Path(path) / f).relative_to(self.expected_output)
            for path, _, files in os.walk(self.expected_output)
            for f in files
        )
        print(f"Targets: {self.target_files}")
        
        # Setup default mounts that can be overridden before calling run()
        print(f"self.tmpdir: {self.tmpdir}")
        self.mounts = [f"{str(self.tmpdir.joinpath('mnt', 'results').resolve())}:/mnt/results"]

    def run(self):
        # Run the container with the correct params
        target_str = " ".join([f'/{f}' for f in self.target_files])
        print(f"Mounting: {self.mounts}")
        logs = self.runner.run(target_str, self.mounts)
        print(logs)
    
    def check(self):
        input_files = set(
            (Path(path) / f).relative_to(self.input_data)
            for path, _, files in os.walk(self.input_data)
            for f in files
        )
        unexpected_files = set()
        for path, _, files in os.walk(self.tmpdir):
            for f in files:
                f = (Path(path) / f).relative_to(self.tmpdir)
                if str(f).startswith(".snakemake"):
                    continue
                if f in self.target_files:
                    self.compare_files(self.tmpdir / f, self.expected_output / f)
                elif f in input_files:
                    # ignore input files
                    pass
                else:
                    unexpected_files.add(f)
        if unexpected_files:
            raise ValueError(
                "Unexpected files:\n{}".format(
                    "\n".join(sorted(map(str, unexpected_files)))
                )
            )

    def compare_files(self, generated_file, expected_file):
        # Check that cmp returns no output (no difference between files)
        assert len(sp.check_output(["cmp", generated_file, expected_file])) == 0


    def cleanup(self):
        if self.tmpdir.is_dir():
            shutil.rmtree(self.tmpdir)

