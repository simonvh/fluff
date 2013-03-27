import pytest


@pytest.fixture
def bamfile():
	return "tests/data/H3K4me3.bam"

@pytest.fixture
def bedfile():
	return "tests/data/H3K4me3.bed"

@pytest.fixture
def region():
    return "tests/data/testregion.bed"

@pytest.fixture
def fragments():
    return "tests/data/testfragment.bed"
    

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

def test_fragmentsize(region, fragments):
    from fluff.fluffio import get_binned_stats
    
    result = get_binned_stats(region, fragments, 4)
    assert 1 == len(result)
    
    counts = [int(x) for x in result[0].split("\t")[-4:]]
    assert [1,0,0,1] == counts
    
    result = get_binned_stats(region, fragments, 4, fragmentsize=50)
    counts = [int(x) for x in result[0].split("\t")[-4:]]
    assert [1,0,0,1] == counts
    
    result = get_binned_stats(region, fragments, 4, fragmentsize=100)
    counts = [int(x) for x in result[0].split("\t")[-4:]]
    assert [1,1,1,1] == counts
    
    result = get_binned_stats(region, fragments, 4, fragmentsize=200)
    counts = [int(x) for x in result[0].split("\t")[-4:]]
    assert [2,2,22,2] == counts
