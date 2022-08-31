import os
import docker
from config import test_tag, test_name

client = docker.from_env()

class SnakemakeRunner:

    def __init__(self):
        self.cons = client.containers
        self.target_strings = []
        self.vols = []

    def run(self):
    # Run the test job.
        logs = self.cons.run(
            test_tag,
            # entrypoint="/bin/ls",
            # command=["-la", "/data/results/alignments/full"],
            entrypoint="snakemake",
            command=[
                *self.target_strings,
                "-f", 
                "-j1",
                "-s=/data/workflow/Snakefile",
            ],
            name=test_name,
            volumes=self.vols
            ).decode('utf-8')
        print(logs)
        return logs


class GothRequestorRunner():

    def __init__(self):
        self.cons = client.containers
        self.vols = []

        assert os.environ.get('YAGNA_APPKEY') is not None, 'No appkey set'

    def run(self):
        logs = self.cons.run(
            test_tag,
            entrypoint="python",
            command=f'/data/workflow/scripts/requestor.py --subnet goth',
            name=test_name,
            environment=[
                f"YAGNA_APPKEY={os.environ['YAGNA_APPKEY']}",
                f"YAGNA_API_URL={os.environ['YAGNA_API_URL']}",
                f"GSB_URL={os.environ['GSB_URL']}",
                ],
            network_mode="host",
            volumes=self.vols
            )

        return logs.decode('utf-8')

# class DevnetRequestorRunner():

#     def __init__(self):
#         ...
    
#     def run(self):
#         logs = client.containers.run(
#             test_tag,
#             command='-m req_only', # Default in /data/config/config.yml is devnet-beta
#             name="test_devnet_requestor",
#             auto_remove=True,
#             volumes=[
#                 f'{str(yagna_datadir)}:/yagna',
#                 ]
#             ).decode('utf-8')
#         return logs
