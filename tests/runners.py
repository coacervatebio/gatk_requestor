import os
import logging
import docker
from tests.config import test_tag, test_name, yagna_datadir

LOGGER = logging.getLogger(__name__)
client = docker.from_env()


def log_container_output(con):
    output = con.attach(stdout=True, stream=True, logs=True)
    for line in output:
        LOGGER.info(line)


class SnakemakeRunner:
    """Runner for containers using default entrypoint"""

    def __init__(self):
        """Setup empty defaults and container handle"""
        self.cons = client.containers
        self.target_string = ""
        self.vols = []
        self.arb_com = []

    def run(self):
        """Run container in specific mode with snakemake targets and any arbitrary commands"""
        com = ["-m", "specific", "-o", self.target_string, *self.arb_com]

        LOGGER.debug(
            f"Running container {test_tag} with: "
            f"Command: {' '.join(com)}, "
            f"Volumes: {' -v '.join(self.vols)}"
        )

        container = self.cons.run(
            test_tag, command=com, name=test_name, volumes=self.vols, detach=True
        )
        log_container_output(container)


class GothRequestorRunner:
    """Execute requestor script only using interactive goth on host"""

    def __init__(self):
        self.cons = client.containers
        self.vols = []

        assert os.environ.get("YAGNA_APPKEY") is not None, "No appkey set"

    def run(self):
        """Run container with python entrypoint, sharing necessary config to communicate with host goth"""
        entry = "python"
        com = f"/data/workflow/scripts/requestor.py --subnet goth"
        env = [
            f"YAGNA_APPKEY={os.environ['YAGNA_APPKEY']}",
            f"YAGNA_API_URL={os.environ['YAGNA_API_URL']}",
            f"GSB_URL={os.environ['GSB_URL']}",
        ]
        net = "host"

        LOGGER.debug(
            f"Running container {test_tag} with: "
            f"Entrypoint: {entry}, "
            f"Command: {com}, "
            f"Env vars: {env}, "
            f"Network: {net}, "
            f"Volumes: {self.vols}"
        )

        container = self.cons.run(
            test_tag,
            entrypoint=entry,
            command=com,
            name=test_name,
            environment=env,
            network_mode=net,
            volumes=self.vols,
            detach=True,
        )
        log_container_output(container)


class TestnetRequestorRunner:
    """Execute requestor script using devnet"""

    def __init__(self):
        self.cons = client.containers
        self.vols = []

    def run(self):
        """Run container with yagna on and persistent volume"""
        com = [
            "-m",
            "req_only",  # Default in /data/config/config.yml is devnet-beta
            "-d",
            "on",
        ]
        vols = [
            *self.vols,
            f"{str(yagna_datadir)}:/home/coacervate/.local/share/yagna",
        ]

        LOGGER.debug(
            f"Running container {test_tag} with: "
            f"Command: {com}, "
            f"Volumes: {vols}"
        )

        container = client.containers.run(
            test_tag, command=com, name=test_name, volumes=vols, detach=True
        )
        log_container_output(container)
