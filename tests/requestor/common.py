"""
Common code for unit testing of rules generated with Snakemake 7.12.1.
"""

import os
import shutil
import subprocess as sp
from pathlib import Path
from docker.errors import DockerException, ContainerError


class ContainerTester:
    def __init__(self, runner, datadir, tmpdir=Path("assets/tmp_output/")):
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
        # print(self.target_files)
        
        # Setup default mounts that can be overridden before calling run()
        self.mounts = [f"{str(self.tmpdir.joinpath('mnt', 'results').resolve())}:/mnt/results"]

    def run(self):
        # Run the container with the correct params
        # A leading slash is included below because the target files are 
        # determined on host and they need to be relative to container fs
        target_str = " ".join([f'/{f}' for f in self.target_files])
        # print(target_str)
        self.runner.run(target_str, self.mounts)

    def check(self, track_unexpected=True):
        # track_unexpected can be set to False when ignoring inputs that change with every run
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
                print("In focus: ", f)
                if str(f).startswith(".snakemake"):
                    continue
                if f in input_files:
                    print("Ignoring: ", f)
                    # ignore input files
                    pass
                elif f in self.target_files:
                    print("Comparing: ", f)
                    self.compare_files(self.tmpdir / f, self.expected_output / f)
                elif track_unexpected:
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

    def compare_files(self, generated_file, expected_file):
        # Check that cmp returns no output (no difference between files)
        assert len(sp.check_output(["cmp", generated_file, expected_file])) == 0

    def cleanup(self):
        if self.tmpdir.is_dir():
            shutil.rmtree(self.tmpdir)

    def run_defaults(self):
        try:
            self.run()
            self.check()
        except ContainerError as ce:
            print(ce.container.logs().decode('utf-8'))
            raise ce
        # Catch any docker SDK issues
        except DockerException as de:
            print(de.container.logs().decode('utf-8'))
            raise de
        finally:
            self.cleanup()
