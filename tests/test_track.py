import pytest
from fluff.track import Track
from pytest import approx

@pytest.fixture
def tracks():
    ftypes = ["bam", "bed", "wig", "bg", "bw", "wig.gz", "bg.gz", "bigWig", "bedGraph"]#, "bed.gz"]
    my_tracks = []
    for ftype in ftypes:
        fname = "tests/data/profile." + ftype
        my_tracks.append(Track.load(fname))
    return my_tracks

@pytest.fixture
def interval():
    return ("scaffold_1", 44749422, 44750067)

@pytest.fixture
def region():
    return "tests/data/profile_region.bed"

@pytest.fixture
def region_fs():
    return "tests/data/testregion.bed"

@pytest.fixture
def bamfile():
    return "tests/data/H3K4me3.bam"

@pytest.fixture
def fragments():
    return "tests/data/testfragment.bed"

@pytest.fixture
def bedfile():
    return "tests/data/H3K4me3.bed"

@pytest.fixture
def frag3():
    return "tests/data/testfragments3.bed"

@pytest.fixture
def frag6():
    return "tests/data/testfragments6.bed"

@pytest.fixture
def region3():
    return "tests/data/threeregions3.bed"

@pytest.fixture
def region6():
    return "tests/data/threeregions6.bed"

def test_count(tracks):
    for track in tracks:
        if track.track_type == "feature":
            assert 10 == track.count()

def test_profile(tracks, interval):
    for track in tracks:
        profile = track.get_profile(interval)
        assert interval[2] - interval[1] == len(profile)
        assert 103 == sum(profile == 2)
        assert 2 == max(profile)
        assert 1 == min(profile)

def test_binned_stats(tracks, interval, region):
    bins = {
        "feature": [
             ['scaffold_1', 44749422, 44750067, 2,0,0,0,1,1,2,4,3,3],
             ['scaffold_2', 44749422, 44750067, 0,0,0,0,0,0,0,0,0,0],
             ],
        "profile": [
             ['scaffold_1', 44749422, 44750067,2,0,0,0,1,1,2,2,2,2],
             ['scaffold_2', 44749422, 44750067, 0,0,0,0,0,0,0,0,0,0],
        ]
    }

    for t in tracks:
        print t
        expect = bins[t.track_type]
        if t.track_type == "profile":
            for expect_row, result_row in zip(expect, t.binned_stats(region, 10, split=True, statistic="max")):
                print(expect_row)
                print(result_row)
                print "***"
                assert expect_row[:3] == result_row[:3]
                assert expect_row[4:] == approx(result_row[4:])
        else:
            print(expect)
            print([x for x in t.binned_stats(region, 10, split=True)])
            assert expect == [x for x in t.binned_stats(region, 10, split=True)]
            
def test_fetch(tracks, interval, region):
    
    for t in tracks:
        if t.track_type == "feature":
            assert 10 == len([i for i in t.fetch(interval)])
            assert 6 == len([i for i in t.fetch(interval, strand="+")])
            assert 4 == len([i for i in t.fetch(interval, strand="-")])

def test_fragmentsize(region_fs, fragments):
    
    track = Track.load(fragments)
    result = track.binned_stats(region_fs, 4)
    assert 1 == len(result)

    counts = [int(x) for x in result[0].split("\t")[-4:]]
    assert [1, 0, 0, 1] == counts

    track.fragmentsize = 50
    result = track.binned_stats(region_fs, 4)
    counts = [int(x) for x in result[0].split("\t")[-4:]]
    assert [1, 0, 0, 1] == counts

    track.fragmentsize = 100
    result = track.binned_stats(region_fs, 4)
    counts = [int(x) for x in result[0].split("\t")[-4:]]
    assert [1, 1, 1, 1] == counts

    track.fragmentsize = 200
    result = track.binned_stats(region_fs, 4)
    counts = [int(x) for x in result[0].split("\t")[-4:]]
    assert [2, 2, 2, 2] == counts

def test_fetch_to_counts(frag3, frag6, region3, region6):
    from pybedtools import BedTool
    for f in [frag3, frag6]:
        t = Track.load(f)
        for r in [region3, region6]:
            overlap = [x for x in t.fetch_to_counts(BedTool(r))]
            assert 3 == len(overlap)
            counts = sorted([len(x[1]) + len(x[2]) for x in overlap])
            assert [0, 1, 3] == counts

def test_get_features_by_feature(frag3, frag6, region3, region6):
    from pybedtools import BedTool
    for f in [frag6]:
        t = Track.load(f)
        for r in [region3, region6]:
            result = [len(x[1]) for x in t._get_features_by_feature(BedTool(r))]
            assert [0, 1, 3] == sorted(result)
