import os
import shutil
import logging
from pathlib import Path
from docker.errors import DockerException, ContainerError, APIError, NotFound
from tests.config import test_name

LOGGER = logging.getLogger(__name__)

class ContainerTester:
    """Main composite class for running tests within containers.

    Handles setting up temporary directory, running, checking and cleanup.
    """
    

    def __init__(self, runner, checker, datadir, tmpdir=Path("assets/TEMP/"), track_unexpected=True):
        """
        :param runner: Component Runner class for running container
        :type runner: _type_
        :param checker: Component Checker class for comparing outputs/expected
        :type checker: _type_
        :param datadir: Directory holding mock assets and expected outputs
        :type datadir: PurePath
        :param tmpdir: Directory for holding temporary outputs, defaults to Path("assets/TEMP/")
        :type tmpdir: Path
        :param track_unexpected: Set to False if there are output files not relevant to testing, defaults to True
        :type track_unexpected: Bool
        """
        LOGGER.info(f"Initializing ContainerTester with {runner}, {checker}, and {datadir}")
        self.datadir = datadir
        self.tmpdir = tmpdir
        self.runner = runner
        self.checker = checker
        self.track_unexpected = track_unexpected #can be set to False when ignoring inputs that change with every run

        # Setup necessary paths, datadir must have input dir (data) and output dir (expected)
        LOGGER.info("Copying mock data to temp directory")
        self.input_data = datadir.joinpath('inputs')
        self.expected_output = datadir.joinpath('expected')
        shutil.copytree(self.input_data, self.tmpdir)

        # Determine target files for container execution and comparison
        self.target_files = set(
            (Path(path) / f).relative_to(self.expected_output)
            for path, _, files in os.walk(self.expected_output)
            for f in files
        )
        LOGGER.debug(f"Determined target files as:\n{self.target_files}")

        # Sometimes target_str needs to be overwritten, e.g when it
        # is a dir producing many unspecified targets
        self.runner.target_string = " ".join([f'/{f}' for f in self.target_files])
        
        # Setup default mounts that can be overridden before calling run()
        self.runner.vols = [f"{str(self.tmpdir.joinpath('data', 'results').resolve())}:/data/results"]


    def check_files(self):
        """Iterate through inputs and outputs, comparing and raising issues where appropriate

        :raises FileNotFoundError: Raised when a target file is missing from the outputs
        :raises ValueError: Raised when an unexpected file is encountered
        """
        LOGGER.info("Checking files")
        input_files = set(
            (Path(path) / f).relative_to(self.input_data)
            for path, _, files in os.walk(self.input_data)
            for f in files
        )
        LOGGER.debug(f"Input files determined as: {input_files}")
        unexpected_files = set()
        all_files = []
        for path, _, files in os.walk(self.tmpdir):
            for f in files:
                f = (Path(path) / f).relative_to(self.tmpdir)
                all_files.append(f)
                LOGGER.debug(f"Focusing on {f}")
                if str(f).startswith(".snakemake"):
                    continue
                if f in input_files:
                    LOGGER.debug(f"Ignoring {f} since it's an input file")
                    # ignore input files
                    pass
                elif f in self.target_files:
                    LOGGER.info(f"Comparing {f} against equivalent in {self.expected_output}")
                    self.checker.compare_files(self.tmpdir / f, self.expected_output / f)
                elif self.track_unexpected:
                    LOGGER.warning(f"Found unexpected file: {f}")
                    unexpected_files.add(f)

        missing_targets = []
        for tf in self.target_files:
            if tf not in all_files:
                LOGGER.error(f"Missing expected target file {tf} from all files: \n{all_files}")
                missing_targets.append(tf)

        # Separate if/raise to log all missing files before the test fails/exits
        if missing_targets:
            LOGGER.critical(f"Exiting due to missing targets..")
            raise FileNotFoundError("Missing targets:\n{}".format(
                    "\n".join(sorted(map(str, missing_targets)))
                )
            )
        
        if unexpected_files:
            LOGGER.critical(f"Exiting due to unexpected files..")
            raise ValueError(
                "Unexpected files:\n{}".format(
                    "\n".join(sorted(map(str, unexpected_files)))
                )
            )

    def clean_con(self):
        """Clean up container created by Runner
        """
        try:
            test_con = self.runner.cons.get(test_name)
            LOGGER.info(f"Cleaning up {test_con}")
            test_con.remove(force=True)
        except NotFound:
            pass # Nothing to clean up

    def clean_tmp(self):
        """Remove temporary directory used for test inputs/outputs
        """
        LOGGER.info(f"Cleaning up tmpdir: {self.tmpdir}")
        if self.tmpdir.is_dir():
            shutil.rmtree(self.tmpdir)

    def run(self, run_con=True, check=True, clean_con=True, clean_tmp=True):
        """Main entrypoint for running defined test case, handling errors and cleaning up.

        :param run_con: Call Runner's run() method, defaults to True
        :type run_con: bool, optional
        :param check: Check output of container run, defaults to True
        :type check: bool, optional
        :param clean_con: Remove container after test, set to False to manually run things against container after test exits, defaults to True
        :type clean_con: bool, optional
        :param clean_tmp: Remove temp dir, set to False to manually inspect inputs/outputs, defaults to True
        :type clean_tmp: bool, optional
        :raises ce: ContainerError, logs issues with container process
        :raises ae: APIError, logs issues with Docker daemon
        :raise de: DockerException, catch-all log for Docker issue
        """
        LOGGER.debug("Entering run() method")
        try:
            if run_con: self.runner.run()
            if check: self.check_files()
        except ContainerError as ce:
            LOGGER.critical("Container exited with non-zero code")
            raise ce
        # Catches errors relating to getting logs from a container that never started
        except APIError as ae:
            LOGGER.critical("Docker daemon returned an error")
            raise ae
        # Catch any docker SDK issues
        except DockerException as de:
            LOGGER.critical("Unhandled Docker exception")
            raise de
        finally:
            if clean_con: self.clean_con()
            if clean_tmp: self.clean_tmp()


# Small Utils

def allowed_pattern(tf):
    """
    Filter function to remove files from comparison matching certain patterns
    identifying files that change from run-to-run. E.g. failing comparison
    because the timestamp is different.

    :param tf: Target file being evaluated
    :returns False: File contains failure pattern
    :returns True: File does not contain failure pattern
    """
    expected_cmp_fail_patterns = [
        '__',
        'vcfheader.vcf',
        'gz.tbi'
    ]
    
    s_tf = str(tf)
    for ecfp in expected_cmp_fail_patterns:
        if ecfp in s_tf:
            LOGGER.warn(f"Excluding {s_tf} for containing {ecfp}")
            return False
    return True
