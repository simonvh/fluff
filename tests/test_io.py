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


@pytest.fixture
def bedtool():
    import pybedtools
    return pybedtools.BedTool("tests/data/testregion.bed")


def test_import_trackwrapper():
    pass


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
    assert [1, 0, 0, 1] == counts

    result = get_binned_stats(region, fragments, 4, fragmentsize=50)
    counts = [int(x) for x in result[0].split("\t")[-4:]]
    assert [1, 0, 0, 1] == counts

    result = get_binned_stats(region, fragments, 4, fragmentsize=100)
    counts = [int(x) for x in result[0].split("\t")[-4:]]
    assert [1, 1, 1, 1] == counts

    result = get_binned_stats(region, fragments, 4, fragmentsize=200)
    counts = [int(x) for x in result[0].split("\t")[-4:]]
    assert [2, 2, 2, 2] == counts


def test_fetch_to_counts(frag3, frag6, region3, region6):
    from fluff.fluffio import TrackWrapper
    from pybedtools import BedTool

    for f in [frag3, frag6]:
        for r in [region3, region6]:
            t = TrackWrapper(f)
            overlap = [x for x in t.fetch_to_counts(BedTool(r), rmdup=False, rmrepeats=False)]
            assert 3 == len(overlap)
            counts = sorted([len(x[1]) + len(x[2]) for x in overlap])
            assert [0, 1, 3] == counts


def test_get_features_by_feature(frag3, frag6, region3, region6):
    from fluff.fluffio import get_features_by_feature
    from pybedtools import BedTool

    for f in [frag6]:
        for r in [region3, region6]:
            assert [0, 1, 3] == sorted([len(x[1]) for x in get_features_by_feature(BedTool(r), BedTool(f))])
