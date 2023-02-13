from pathlib import Path

ROOTPATH = Path(__file__).parent.parent  # Project root path

test_tag = "coacervate/requestor:test"
test_name = "test_container"
test_sample = "HG03633"
unit_assets_root = ROOTPATH.joinpath('tests', 'unit', 'assets')  
yagna_datadir = "yagna_datadir" # Named volume to persist wallet state
default_tmp = ROOTPATH.joinpath('tests', 'unit', 'assets', 'TEMP')