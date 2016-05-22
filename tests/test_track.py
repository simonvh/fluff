import pytest


@pytest.fixture
def bamfile():
    return "tests/data/profile.bam"

@pytest.fixture
def bedfile():
    return "tests/data/profile.bed"

@pytest.fixture
def wigfile():
    return "tests/data/profile.wig"

@pytest.fixture
def bwfile():
    return "tests/data/profile.bw"

@pytest.fixture
def tabixfile():
    return "tests/data/profile.wig.gz"

@pytest.fixture
def interval():
    return "scaffold_1", 44749422, 44750067

@pytest.fixture
def region():
    return "tests/data/profile_region.bed"

def test_bamtrack(bamfile, bedfile, wigfile, bwfile, interval, region):
    from fluff.track import BamTrack, BedTrack, WigTrack, BigWigTrack
    
    # TODO: check if not [2,0,0,0,1,1,2,3,3,3]
    bins = [2,0,0,0,1,1,2,4,3,3]
    tracks =[ BamTrack(bamfile), BedTrack(bedfile)]
    for t in tracks:
        assert 10 == t.count()
        assert 10 == len([i for i in t.fetch(interval)])
        assert 6 == len([i for i in t.fetch(interval, strand="+")])
        assert 4 == len([i for i in t.fetch(interval, strand="-")])
        assert [list(interval) + bins] == [x for x in t.binned_stats(region, 10, split=True)]
    
    tracks += [WigTrack(wigfile), BigWigTrack(bwfile)]
    for track in tracks:
        profile = t.get_profile(interval)
        assert interval[2] - interval[1] == len(profile)
        assert 103 == sum(profile == 2)
        assert 2 == max(profile)
        assert 1 == min(profile)

    bins = [2,0,0,0,1,1,2,2,2,2]
    tracks = [BigWigTrack(bwfile)]
    for t in tracks:
        assert [list(interval) + bins] == [x for x in t.binned_stats(region, 10, split=True)]
