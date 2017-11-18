import pytest
from tempfile import NamedTemporaryFile
from fluff.parse import parse_cmds
from fluff.commands.profile import profile
from fluff.commands.heatmap import heatmap

@pytest.fixture
def bamfile():
    return "tests/data/H3K4me3.bam"

@pytest.fixture
def bwfile():
    return "tests/data/profile.bw"

@pytest.fixture
def bedfile():
    return "tests/data/profile.bed"

@pytest.fixture
def regionfile():
    return "tests/data/profile_region.bed"


def test_profile(bamfile):
    # Only tests of the command runs successfully,
    # doesnt't check the image
    tmp = NamedTemporaryFile(prefix="fluff.", suffix=".png", delete=False)
    args = ["profile",
            "-i", "scaffold_1:44749422-44750067",
            "-d", bamfile,
            "-o", tmp.name]
    args = parse_cmds().parse_args(args)
    profile(args)

def test_heatmap(bamfile, bedfile, bwfile, regionfile):
    # Only tests of the command runs successfully,
    # doesnt't check the image
    tmp = NamedTemporaryFile(prefix="fluff.", suffix=".png", delete=False)
    args = ["heatmap",
            "-f", regionfile,
            "-d", bamfile, bwfile, bedfile,
            "-o", tmp.name]
    args = parse_cmds().parse_args(args)
    heatmap(args)





