import pytest


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
