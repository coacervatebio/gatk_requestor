from pathlib import Path

ROOTPATH = Path(__file__).parent.parent  # Project root path

test_tag = "coacervate_requestor:test"
test_name = "test_container"
test_sample = "HG03633"
unit_assets_root = ROOTPATH.joinpath('tests', 'unit', 'assets')  
yagna_datadir = Path("/home/vagrant/yagna_datadir/")
default_tmp = ROOTPATH.joinpath('test', 'TEMP')