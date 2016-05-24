import pytest
from tempfile import NamedTemporaryFile

from fluff.parse import parse_args
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

@pytest.fixture
def profile_args(bamfile):
    tmp = NamedTemporaryFile(prefix="fluff.", suffix=".png", delete=False)
    args = ["profile", 
            "-i", "scaffold_1:44749422-44750067",
            "-d", bamfile,
            "-o", tmp.name]
    return parse_args(args)

@pytest.fixture
def heatmap_args(bamfile, bedfile, bwfile, regionfile):
    tmp = NamedTemporaryFile(prefix="fluff.", suffix=".png", delete=False)
    args = ["heatmap", 
            "-f", regionfile,
            "-d", bamfile, bwfile, bedfile,
            "-o", tmp.name]
    return parse_args(args)

def test_profile(profile_args):
    # Only tests of the command runs successfully, 
    # doesnt't check the image
    profile(profile_args) 

def test_heatmap(heatmap_args):
    # Only tests of the command runs successfully, 
    # doesnt't check the image
    heatmap(heatmap_args) 
