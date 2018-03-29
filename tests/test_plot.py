import pytest
from fluff.plot import profile_screenshot
from tempfile import NamedTemporaryFile 

@pytest.fixture
def region():
    return "scaffold_1:44749422-44750067"

@pytest.fixture
def bamfile():
    return "tests/data/H3K4me3.bam"

def test_profile(region, bamfile):
    with NamedTemporaryFile(suffix=".png") as tmp:
        profile_screenshot(tmp.name, region, [[bamfile]])
