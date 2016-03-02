import pytest


@pytest.fixture
def heatmap_data_single_track():
    from numpy import array
    data = {
        "track1": array([
            [0, 0, 0, 0, 0, 10, 10, 10, 10, 10],
            [0, 0, 0, 0, 0, 10, 10, 10, 10, 10],
            [0, 0, 0, 0, 0, 10, 10, 10, 10, 10],
            [0, 0, 0, 0, 0, 10, 10, 10, 10, 10],
            [10, 10, 10, 10, 10, 0, 0, 0, 0, 0],
            [10, 10, 10, 10, 10, 0, 0, 0, 0, 0],
            [10, 10, 10, 10, 10, 0, 0, 0, 0, 0],
            [10, 10, 10, 10, 10, 0, 0, 0, 0, 0],
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        ])
    }
    labels = array([0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2])
    return data, labels


@pytest.fixture
def heatmap_data_multiple_tracks():
    from numpy import array
    data = {
        "track1": array([
            [0, 0, 0, 0, 0, 10, 10, 10, 10, 10],
            [0, 0, 0, 0, 0, 10, 10, 10, 10, 10],
            [0, 0, 0, 0, 0, 10, 10, 10, 10, 10],
            [10, 10, 10, 10, 10, 0, 0, 0, 0, 0],
            [10, 10, 10, 10, 10, 0, 0, 0, 0, 0],
            [10, 10, 10, 10, 10, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 10, 10, 10, 10, 10],
            [0, 0, 0, 0, 0, 10, 10, 10, 10, 10],
            [0, 0, 0, 0, 0, 10, 10, 10, 10, 10],
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        ]),
        "track2": array([
            [10, 10, 10, 10, 10, 10, 10, 10, 10, 10],
            [10, 10, 10, 10, 10, 10, 10, 10, 10, 10],
            [10, 10, 10, 10, 10, 10, 10, 10, 10, 10],
            [10, 10, 10, 10, 10, 0, 0, 0, 0, 0],
            [10, 10, 10, 10, 10, 0, 0, 0, 0, 0],
            [10, 10, 10, 10, 10, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 10, 10, 10, 10, 10],
            [0, 0, 0, 0, 0, 10, 10, 10, 10, 10],
            [0, 0, 0, 0, 0, 10, 10, 10, 10, 10],
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        ])
    }

    labels = array([0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3])
    return data, labels


def test_mirror_single_track(heatmap_data_single_track):
    from fluff.util import mirror_clusters
    data, labels = heatmap_data_single_track
    assert data["track1"][5][0] == 10
    (i, j) = mirror_clusters(data, labels)
    assert i == 0
    assert j == 1


def test_mirror_multiple_track(heatmap_data_multiple_tracks):
    from fluff.util import mirror_clusters
    data, labels = heatmap_data_multiple_tracks
    (i, j) = mirror_clusters(data, labels, 0.05)
    assert i == 1
    assert j == 2
