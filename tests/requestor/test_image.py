import os
import rootpath
from pathlib import Path
import docker

test_tag = 'requestor:test'
rpath = Path(rootpath.detect())
docker_context = Path.joinpath(rpath, 'requestor')
containerfile = Path.joinpath(rpath, 'requestor', 'Containerfile')

client = docker.from_env()
# build = client.images.build(path=str(docker_context), tag=test_tag)

# def test_build():
#     built = False
#     for log in build[1]:
#         if 'Successfully built' in log.get('stream', 'None'):
#             built = True
#             break
#     assert built is True

def test_goth_requestor():

    container = client.containers.run(
        "requestor:test",
        command='/mnt/workflow/scripts/requestor.py --subnet goth',
        environment=[
            f"YAGNA_APPKEY={os.environ['YAGNA_APPKEY']}",
            f"YAGNA_API_URL={os.environ['YAGNA_API_URL']}",
            f"GSB_URL={os.environ['GSB_URL']}",
            ],
        entrypoint="python",
        network_mode="host",
        # detach=True
        )

    print('\n',container)

# def test_devnet_requestor():
#     ...