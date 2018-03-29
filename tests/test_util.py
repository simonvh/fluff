import pytest
import numpy as np

@pytest.fixture
def track_data():
    d = [
        [4, 8, 65, 80, 58, 88, 97, 13, 21, 81],
        [74, 32, 7, 79, 11, 15, 66, 10, 45, 78],
        [91, 7, 12, 5, 69, 74, 32, 83, 60, 9],
    ]

    return d


@pytest.fixture
def track_data_low():
    d = [
        [0, 0, 0, 100, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 5, 0, 0],
    ]

    return d

@pytest.fixture
def cluster_data():
    return (np.array([
            [1,1,1,1,1,100,1,1,1,1,1],
            [2,2,2,2,2,200,2,2,2,2,2],
            [400,300,200,100,0,0,0,100,200,300,400],
            [500,400, 300, 200, 100, 100, 100, 200, 300, 400, 500],
            [1,1,1,1,1,1,1,1,1,1,1],
            [390, 220, 200, 80, -1, 10, 0, 110, 180, 301, 420],
            ], dtype="float"),
            [0,0,1,1,2,1])

def test_get_absolute_scale(track_data):
    from fluff.util import get_absolute_scale
    assert 10.0 == get_absolute_scale(10, track_data)
    assert 99.0 == get_absolute_scale(99.0, track_data)
    assert 99.0 == get_absolute_scale("99.0", track_data)
    assert 83.5 == get_absolute_scale("90.0%", track_data, False)
    assert 7.0 == get_absolute_scale("10%", track_data)
    assert type([]) == type(get_absolute_scale("80%", track_data, True))
    for x, y in zip([82.4, 74.8, 75.8], get_absolute_scale("80%", track_data, True)):
        assert round(x - y, 5) == 0


def test_get_absolute_scale_low(track_data_low):
    from fluff.util import get_absolute_scale
    assert 5.0 == get_absolute_scale("10%", track_data_low)
    assert 5.0 == get_absolute_scale("50%", track_data_low)

def test_cluster_profile(cluster_data):
    from fluff.util import cluster_profile
    cluster_input, cluster_labels = cluster_data
    
    # Test k-means
    ind, labels = cluster_profile(cluster_input, numclusters=3, dist="pearson")
    cluster_labels = np.array(cluster_labels)
    labels = np.array(labels)
    assert len(np.unique(labels)) == len(np.unique(cluster_labels))
    for label in np.unique(cluster_labels):
        l = len(np.unique(labels[np.where(cluster_labels == label)]))
        assert 1 == l

    # Test hierarchical
    results = [
            ("euclidean",  [1,0,4,3,2,5]),
            ("pearson", [5, 2, 3, 4, 0, 1]),
            ]
    
    for dist, result in results:
        ind, labels = cluster_profile(
            cluster_input, 
            cluster_type="h", 
            numclusters=3, 
            dist=dist)
    
        assert result == ind

