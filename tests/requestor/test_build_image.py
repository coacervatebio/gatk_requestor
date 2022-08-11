import rootpath
from pathlib import Path
import docker

rpath = Path(rootpath.detect())
docker_context = Path.joinpath(rpath, 'requestor')
containerfile = Path.joinpath(rpath, 'requestor', 'Containerfile')

client = docker.from_env()
build = client.images.build(path=str(docker_context), tag='requestor:test')

def test_build():
    built = False
    for log in build[1]:
        if 'Successfully built' in log.get('stream', 'None'):
            built = True
            break
    assert built is True