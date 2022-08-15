import os
import docker
import rootpath
from time import sleep
from pathlib import Path

from config import test_tag, rpath

docker_context = Path.joinpath(rpath, 'requestor')
client = docker.from_env()

def test_build():
    build = client.images.build(path=str(docker_context), tag=test_tag)
    built = False
    for log in build[1]:
        if 'Successfully built' in log.get('stream', 'None'):
            built = True
            break
    assert built is True