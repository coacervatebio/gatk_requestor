import sys
import click
import time
import logging
from pathlib import Path
from common import TestRunner
from checkers import SimpleChecker, VcfChecker

logging.basicConfig(level=logging.DEBUG)


@click.command()
@click.option()
def check():

    runner = TestRunner(VcfChecker(), Path('/assets'), Path('/temp_test'), track_unexpected=False)

    while not Path('/temp_test/data/results/test.done').is_file():
        time.sleep(2)

    runner.check_files()
