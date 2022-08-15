import os
import docker
from pathlib import Path

from config import (
    test_tag,
    yagna_datadir,
    expected_hc_out,
    almt_dir,
    tmp_out_dir
)

client = docker.from_env()

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
            f'{str(almt_dir)}:/mnt/results/alignments',
            f'{str(tmp_out_dir)}:/mnt/results/hc_out'
            ]
        )

    print(logs.decode('utf-8'))

    # Test that all output files were correctly created and clean up
    for out_path in expected_hc_out:
        assert Path(out_path).is_file()
        os.remove(out_path)


def test_devnet_requestor():

    # Run container via usual entrypoint script with requestor-only flag
    # Mounting yagna datadir to avoid re-funding every time
    logs = client.containers.run(
        test_tag,
        command='-m req_only',
        name="test_devnet_requestor",
        auto_remove=True,
        volumes=[
            f'{str(yagna_datadir)}:/yagna',
            f'{str(almt_dir)}:/mnt/results/alignments',
            f'{str(tmp_out_dir)}:/mnt/results/hc_out'
            ]
        )

    print(logs.decode('utf-8'))

    # Test that all output files were correctly created and clean up
    for out_path in expected_hc_out:
        assert Path(out_path).is_file()
        os.remove(out_path)