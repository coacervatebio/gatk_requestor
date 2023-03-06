import sys
import click
import time
import logging
from pathlib import Path
from common import TestRunner
from checkers import SimpleChecker, VcfChecker

logging.basicConfig(level=logging.DEBUG)

checkers = {
    'simple': SimpleChecker,
    'vcf': VcfChecker
}

@click.command()
@click.option('-c', '--checker-type', default='simple', help='Checker class as defined in checkers.py for comparing expected to actual oututs')
@click.option('-a', '--assets-root', default='/assets', help='Asset root dir for this test containing INPUTS and EXPECTED dirs')
@click.option('-t', '--temp-dir', default='/temp_test', help='Temp dir where inputs are copied and tests are performed, cleaned up after every run')
@click.option('-u', '--track_unexpected', default=True, help='Test will fail if output is present that is not in the EXPECTED dir')
def check(checker_type, assets_root, temp_dir, track_unexpected):

    runner = TestRunner(checkers[checker_type], Path(assets_root), Path(temp_dir), track_unexpected)

    while not Path(temp_dir).joinpath('data/results/test.done').is_file():
        time.sleep(2)

    runner.check_files()

if __name__ == '__main__':
    check()
