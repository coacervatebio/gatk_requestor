import docker
from pathlib import Path

from tests.config import test_tag, ROOTPATH

docker_context = Path.joinpath(ROOTPATH, "coacervate")
client = docker.from_env()


def test_build():
    build = client.images.build(path=str(docker_context), tag=test_tag)
    built = False
    for log in build[1]:
        if "Successfully built" in log.get("stream", "None"):
            built = True
            break
    assert built is True
