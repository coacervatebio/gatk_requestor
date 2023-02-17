import sys
import click
import logging
from pathlib import Path
from common import TestRunner
from checkers import SimpleChecker

logging.basicConfig(level=logging.DEBUG)
# LOGGER = logging.getLogger(__name__)
# c_handler = logging.StreamHandler()
# c_handler.setLevel(logging.DEBUG)
# LOGGER.addHandler(c_handler)
# LOGGER.debug("hullo")

runner = TestRunner(SimpleChecker(), Path('/assets'), Path('/temp_test'), track_unexpected=False)
runner.check_files()
