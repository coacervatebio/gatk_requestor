import docker
from config import rpath, test_tag

client = docker.from_env()


def test_valid_snakefile():
    logs = client.containers.run(
        test_tag,
        entrypoint='snakemake',
        command='-np -s=/mnt/workflow/rules/hc_golem.smk',
        name="test_valid_snakefile",
        auto_remove=True
    )
    print(logs.decode('utf-8'))