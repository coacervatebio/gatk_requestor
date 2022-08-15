import os
import rootpath
from pathlib import Path
from time import sleep
import docker

test_tag = 'coacervate_requestor:test'
rpath = Path(rootpath.detect())
docker_context = Path.joinpath(rpath, 'requestor')
containerfile = Path.joinpath(rpath, 'requestor', 'Dockerfile')

client = docker.from_env()
build = client.images.build(path=str(docker_context), tag=test_tag)

def test_build():
    built = False
    for log in build[1]:
        if 'Successfully built' in log.get('stream', 'None'):
            built = True
            break
    assert built is True

def test_goth_requestor():

    # Create mount points for input and output to pass to container
    in_mnt = Path.joinpath(rpath, 'tests', 'requestor', 'assets', 'alignments')
    out_mnt = Path.joinpath(rpath, 'tests', 'requestor', 'assets', 'tmp_output')

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

    # Run container in blocking mode with detach=False
    # Goth must be running and env vars set correctly
    container = client.containers.run(
        test_tag,
        command='/mnt/workflow/scripts/requestor.py --subnet goth',
        environment=[
            f"YAGNA_APPKEY={os.environ['YAGNA_APPKEY']}",
            f"YAGNA_API_URL={os.environ['YAGNA_API_URL']}",
            f"GSB_URL={os.environ['GSB_URL']}",
            ],
        name="test_goth_requestor",
        entrypoint="python",
        network_mode="host",
        auto_remove=True,
        volumes=[
            f'{str(in_mnt)}:/mnt/results/alignments',
            f'{str(out_mnt)}:/mnt/results/hc_out'
            ]
        )

    # Test that all output files were correctly created and clean up
    for out_path in expected_out:
        assert Path(out_path).is_file()
        os.remove(out_path)

# def test_devnet_requestor():
#     ...