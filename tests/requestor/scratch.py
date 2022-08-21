import shutil
from pathlib import PurePath, Path
from common import ContainerTester
from runners import SnakemakeRunner

data_path = PurePath("assets/index_cram/")
workdir = Path("assets/tmp_output/")

try:
    tester = ContainerTester(data_path, workdir, SnakemakeRunner)
    tester.run()
except Exception as e:
    print(e)
# finally:
    # tester.cleanup()