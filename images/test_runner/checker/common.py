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

        # Determine input files so they can be excluded from comparison
        self.input_dir = datadir.joinpath("inputs")
        self.input_files = set(
            (Path(path) / f).relative_to(self.input_dir)
            for path, _, files in os.walk(self.input_dir)
            for f in files
        )
        LOGGER.debug(f"Input files determined as: {self.input_files}")

        # Determine expected files
        self.expected_dir = datadir.joinpath("expected")
        self.expected_files = set(
            (Path(path) / f).relative_to(self.expected_dir)
            for path, _, files in os.walk(self.expected_dir)
            for f in files
        )
        LOGGER.debug(f"Determined expected files as: {self.expected_files}")
        
        # Determine comparison files (excludes expected files that should not be compared)
        self.comparison_files = set(filter(self.allowed_pattern, self.expected_files))
        LOGGER.debug(f"Determined comparison files as: {self.comparison_files}")

    def check_files(self):
        """Iterate through inputs and outputs, comparing and raising issues where appropriate

        :raises FileNotFoundError: Raised when a target file is missing from the outputs
        :raises ValueError: Raised when an unexpected file is encountered
        """
        LOGGER.info("Checking files")
        unexpected_files = set()
        all_files = []
        for path, _, files in os.walk(self.tmpdir):
            for f in files:
                f = (Path(path) / f).relative_to(self.tmpdir)
                all_files.append(f)
                LOGGER.debug(f"Focusing on {f}")

                if f in self.comparison_files:
                    LOGGER.info(
                        f"Comparing {f} against equivalent in {self.expected_dir}"
                    )
                    self.checker.compare_files(
                        self.tmpdir / f, self.expected_dir / f
                    )
                elif f in self.expected_files:
                    LOGGER.debug(f"Ignoring {f} since it's expected but not targeted for comparison")
                    continue
                elif ".snakemake" in str(f) or "test.done" in str(f):
                    LOGGER.debug(f"Ignoring {f} since it's a cache/flag file")
                    continue
                elif f in self.input_files:
                    LOGGER.debug(f"Ignoring {f} since it's an input file")
                    continue
                elif self.track_unexpected:
                    LOGGER.warning(f"Found unexpected file: {f}")
                    unexpected_files.add(f)

        missing_targets = []
        for tf in self.expected_files:
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
    
    @staticmethod
    def allowed_pattern(tf):
        """
        Filter function to remove files which may be expected but should
        not be compared, e.g. index files.

        :param tf: Target file being evaluated
        :returns False: File contains failure pattern
        :returns True: File does not contain failure pattern
        """
        expected_cmp_fail_patterns = ["__", "vcfheader.vcf", "gz.tbi"]

        s_tf = str(tf)
        for ecfp in expected_cmp_fail_patterns:
            if ecfp in s_tf:
                LOGGER.warn(f"Excluding {s_tf} for containing {ecfp}")
                return False
        return True

