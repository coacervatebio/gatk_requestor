import os
import shutil
import logging
from pathlib import Path

LOGGER = logging.getLogger(__name__)

class TestRunner:
    """Main composite class for running tests within containers.

    Handles setting up temporary directory, running, checking and cleanup.
    """

    def __init__(
        self,
        checker,
        datadir,
        tmpdir,
        track_unexpected=True,
    ):
        """
        :param checker: Component Checker class for comparing outputs/expected
        :type checker: _type_
        :param datadir: Directory holding mock assets and expected outputs
        :type datadir: PurePath
        :param track_unexpected: Set to False if there are output files not relevant to testing, defaults to True
        :type track_unexpected: Bool
        """
        LOGGER.info(
            f"Initializing tester with {checker}, and {datadir}"
        )
        self.checker = checker
        self.datadir = datadir
        self.tmpdir = tmpdir
        self.track_unexpected = track_unexpected  # can be set to False when ignoring inputs that change with every run

        # Setup necessary paths, datadir must have input dir (data) and output dir (expected)
        self.input_data = datadir.joinpath("inputs")
        self.expected_output = datadir.joinpath("expected")

        # Determine target files for container execution and comparison
        self.target_files = set(
            (Path(path) / f).relative_to(self.expected_output)
            for path, _, files in os.walk(self.expected_output)
            for f in files
        )
        LOGGER.debug(f"Determined target files as: {self.target_files}")

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
                print('\n')
                print(type(f))
                print(self.target_files)
                print('\n')
                if str(f).startswith(".snakemake"):
                    continue
                elif f in input_files:
                    LOGGER.debug(f"Ignoring {f} since it's an input file")
                    pass
                elif f in self.target_files:
                    LOGGER.info(
                        f"Comparing {f} against equivalent in {self.expected_output}"
                    )
                    self.checker.compare_files(
                        self.tmpdir / f, self.expected_output / f
                    )
                elif self.track_unexpected and 'test.done' not in str(f):
                    LOGGER.warning(f"Found unexpected file: {f}")
                    unexpected_files.add(f)

        missing_targets = []
        for tf in self.target_files:
            if tf not in all_files:
                LOGGER.error(
                    f"Missing expected target file {tf} from all files: {all_files}"
                )
                missing_targets.append(tf)

        # Separate if/raise to log all missing files before the test fails/exits
        if missing_targets:
            LOGGER.critical(f"Exiting due to missing targets..")
            raise FileNotFoundError(
                "Missing targets:\n{}".format(
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

