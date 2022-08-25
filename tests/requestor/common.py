"""
Common code for unit testing of rules generated with Snakemake 7.12.1.
"""

import os
import shutil
import subprocess as sp
from pathlib import Path
from docker.errors import DockerException, ContainerError, APIError


class ContainerTester:
    def __init__(self, runner, datadir, tmpdir=Path("assets/tmp_output/")):
        self.datadir = datadir
        self.tmpdir = tmpdir
        self.runner = runner
        self.track_unexpected = True #can be set to False when ignoring inputs that change with every run


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
        # Sometimes target_str needs to be overwritten, e.g when it
        # is a dir producing many unspecified targets
        self.target_str = [f'/{f}' for f in self.target_files]
        
        # Setup default mounts that can be overridden before calling run()
        self.mounts = [f"{str(self.tmpdir.joinpath('mnt', 'results').resolve())}:/mnt/results"]


    def check(self):
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

    @staticmethod
    def compare_files(generated_file, expected_file):
        # Check that cmp returns no output (no difference between files)
        assert len(sp.check_output(["cmp", generated_file, expected_file])) == 0

    def cleanup(self):
        if self.tmpdir.is_dir():
            shutil.rmtree(self.tmpdir)

    def run(self, container=True, check=True, cleanup=True):
        try:
            if container:
                self.runner.run(self.target_str, self.mounts)
            if check: self.check()
        except APIError as ae:
            raise ae
        except ContainerError as ce:
            print(ce.container.logs().decode('utf-8'))
            raise ce
        # Catch any docker SDK issues
        except DockerException as de:
            raise de
        except:
            raise
        finally:
            if cleanup: self.cleanup()
