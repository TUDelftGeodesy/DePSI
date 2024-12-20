from depsi.io import read_metadata


def test_read_metadata_lines_pixels():
    metadata = read_metadata("tests/data/example.res")
    assert metadata["n_lines"] == 14065
    assert metadata["n_pixels"] == 49843
