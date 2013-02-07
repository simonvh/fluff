import pytest


@pytest.fixture
def bamfile():
	return "tests/data/H3K4me3.bam"

@pytest.fixture
def bedfile():
	return "tests/data/H3K4me3.bed"

def test_import_trackwrapper():
	from fluff.fluffio import TrackWrapper

def test_bam_window(bamfile):
	from fluff.fluffio import TrackWrapper
	t = TrackWrapper(bamfile)
	assert 0 == len(t[('scaffold_2', 45921140, 45922270, ".")])
	assert 9 == len(t[('scaffold_1', 45921140, 45922270, ".")])
	assert 1 == len(t[('scaffold_1', 45921140, 45922270, "+")])
	assert 8 == len(t[('scaffold_1', 45921140, 45922270, "-")])

def test_bed_window(bedfile):
	from fluff.fluffio import TrackWrapper
	t = TrackWrapper(bedfile)
	assert 0 == len(t[('scaffold_2', 45921140, 45922270, ".")])
	assert 9 == len(t[('scaffold_1', 45921140, 45922270, ".")])
	assert 1 == len(t[('scaffold_1', 45921140, 45922270, "+")])
	assert 8 == len(t[('scaffold_1', 45921140, 45922270, "-")])

