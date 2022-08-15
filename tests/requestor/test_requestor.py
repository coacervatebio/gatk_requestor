import os
import rootpath
from pathlib import Path
from time import sleep
import docker

from config import test_tag, rpath, yagna_datadir

# Create mount points for input and output to pass to container
in_mnt = Path.joinpath(rpath, 'tests', 'requestor', 'assets', 'alignments')
out_mnt = Path.joinpath(rpath, 'tests', 'requestor', 'assets', 'tmp_output')

client = docker.from_env()

# Get names of input samples dynamically
samples = []
for reg_dir in in_mnt.glob('chr*'):
    for almt in reg_dir.glob('*cram'):
        samples.append(almt.with_suffix('').name)

# Create expected output sample paths
expected_out = []
for sample in samples:
    expected_out.append(f"{out_mnt.joinpath(sample)}.g.vcf.gz")
    expected_out.append(f"{out_mnt.joinpath(sample)}.g.vcf.gz.tbi")

def test_goth_requestor():
    # Run container in blocking mode with detach=False
    # Goth must be running and env vars set correctly
    logs = client.containers.run(
        test_tag,
        entrypoint="python",
        command='/mnt/workflow/scripts/requestor.py --subnet goth',
        name="test_goth_requestor",
        environment=[
            f"YAGNA_APPKEY={os.environ['YAGNA_APPKEY']}",
            f"YAGNA_API_URL={os.environ['YAGNA_API_URL']}",
            f"GSB_URL={os.environ['GSB_URL']}",
            ],
        network_mode="host",
        auto_remove=True,
        volumes=[
            f'{str(in_mnt)}:/mnt/results/alignments',
            f'{str(out_mnt)}:/mnt/results/hc_out'
            ]
        )

    print(logs.decode('utf-8'))

    # Test that all output files were correctly created and clean up
    for out_path in expected_out:
        assert Path(out_path).is_file()
        os.remove(out_path)


def test_devnet_requestor():

    logs = client.containers.run(
        test_tag,
        command='-m req_only',
        name="test_devnet_requestor",
        auto_remove=True,
        volumes=[
            f'{str(yagna_datadir)}:/yagna',
            f'{str(in_mnt)}:/mnt/results/alignments',
            f'{str(out_mnt)}:/mnt/results/hc_out'
            ]
        )

    print(logs.decode('utf-8'))

    # Test that all output files were correctly created and clean up
    for out_path in expected_out:
        assert Path(out_path).is_file()
        os.remove(out_path)