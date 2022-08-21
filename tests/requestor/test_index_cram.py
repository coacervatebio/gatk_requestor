import shutil
import docker
from time import sleep
from tempfile import TemporaryDirectory
from pathlib import Path, PurePosixPath

import common
from config import test_tag, test_sample


def test_index_cram():

    client = docker.from_env()

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath("assets/index_cram/data")
        expected_path = PurePosixPath("assets/index_cram/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # Run the test job.
        logs = client.containers.run(
            test_tag,
            entrypoint="snakemake",
            command=[
                f"/mnt/results/alignments/full/{test_sample}.cram.crai",
                "-f", 
                "-j1",
                "-s=/mnt/workflow/Snakefile",
            ],
            name="test_index_cram",
            auto_remove=True,
            volumes=[f'{str(workdir)}/mnt/results:/mnt/results']
            )

        print(logs.decode('utf-8'))

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file), 
        # also see common.py.
        common.OutputChecker(data_path, expected_path, workdir).check()
        # sleep(240)


