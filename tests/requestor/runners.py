import os
import docker
from config import test_tag

client = docker.from_env()

class SnakemakeRunner:
    def run(target_files, vols):
    # Run the test job.
        logs = client.containers.run(
            test_tag,
            entrypoint="snakemake",
            # entrypoint="/bin/ls",
            # command="/mnt/results/alignments/full",
            command=[
                target_files,
                "-f", 
                "-j1",
                "-s=/mnt/workflow/Snakefile",
            ],
            name="test_index_cram",
            # auto_remove=True,
            volumes=vols
            )

        return logs.decode('utf-8')


class RequestorRunner():
    def run(target_file, vols):
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
            volumes=vols
            )

        return logs.decode('utf-8')