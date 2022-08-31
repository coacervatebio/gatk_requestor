import os
import docker
from pathlib import PurePath
from common import ContainerTester
from runners import GothRequestorRunner, DevnetRequestorRunner
from checkers import VcfChecker

from config import (
    test_tag,
    yagna_datadir,
    expected_hc_out,
    almt_dir,
    tmp_out_dir
)

client = docker.from_env()

# def test_goth_requestor():

#     data_path = PurePath("assets/call_variants/")

#     tester = ContainerTester(GothRequestorRunner(), SimpleChecker(), data_path)
#     tester.run()


def test_devnet_requestor():

    data_path = PurePath("assets/call_variants/")

    tester = ContainerTester(DevnetRequestorRunner(), VcfChecker(), data_path)
    tester.run()