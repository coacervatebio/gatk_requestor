import os
import docker
from config import test_tag, test_name

client = docker.from_env()

class SnakemakeRunner:

    def __init__(self):
        self.cons = client.containers

    def run(self, target_files, vols):
    # Run the test job.
        logs = self.cons.run(
            test_tag,
            # entrypoint="/bin/ls",
            # command=["-la", "/data/results/alignments/full"],
            entrypoint="snakemake",
            command=[
                *target_files,
                "-f", 
                "-j1",
                "-s=/data/workflow/Snakefile",
            ],
            name=test_name,
            volumes=vols
            ).decode('utf-8')
        print(logs)
        return logs


class RequestorRunner():
    def run(self, target_file, vols):
        self.cons = client.containers
        logs = self.cons.run(
            test_tag,
            entrypoint="python",
            command='/data/workflow/scripts/requestor.py --subnet goth',
            name=test_name,
            environment=[
                f"YAGNA_APPKEY={os.environ['YAGNA_APPKEY']}",
                f"YAGNA_API_URL={os.environ['YAGNA_API_URL']}",
                f"GSB_URL={os.environ['GSB_URL']}",
                ],
            network_mode="host",
            volumes=vols
            )

        return logs.decode('utf-8')