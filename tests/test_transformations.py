import numpy as np
import pytest
import sarxarray
import xarray as xr

from depsi.transformations import radar_to_latlonh, radar_to_time, seconds_of_day


@pytest.mark.parametrize(
    "time, expected_seconds",
    [
        (np.datetime64("2023-03-15T00:00:00"), 0.0),
        (np.datetime64("2023-03-15T12:00:00"), 43200.0),
        (np.datetime64("2023-03-15T23:59:59.999999999"), 86399.999999999),
        (np.datetime64("2023-03-15T12:34:56.123456"), 45296.123456),
    ],
)
def test_seconds_of_day(time, expected_seconds):
    seconds = seconds_of_day(time)
    assert np.isclose(seconds, expected_seconds)


def test_radar_to_time():
    metadata_path = "tests/data/metadata/doris5/20180318/metadata.res"
    metadata_dict = sarxarray.read_metadata(metadata_path, driver="doris5")

    data = xr.open_zarr("tests/data/stm_sparse.zarr")

    azimuths = data.azimuth.values
    ranges = data.range.values
    az_time, rg_time = radar_to_time(azimuths, ranges, metadata_dict)

    az_time_diff = az_time - (
        metadata_dict["first_azimuth_time"]
        + np.array(azimuths / metadata_dict["pulse_repetition_frequency"] * 1e9, dtype="timedelta64[ns]"),
    )
    assert az_time_diff.astype(float).max() < 1  # diff within 1 nano second
    assert np.allclose(rg_time, ranges / metadata_dict["range_sampling_rate"] + metadata_dict["first_range_time"])


def test_radar_to_latlonh():
    metadata_path = "tests/data/metadata/doris5/20180318/metadata.res"
    metadata_dict = sarxarray.read_metadata(metadata_path, driver="doris5")

    data = xr.open_zarr("tests/data/stm_sparse.zarr")

    azimuths = data.azimuth.values
    ranges = data.range.values
    elevations = np.ones(azimuths.shape) * 42
    latlonh = radar_to_latlonh(azimuths, ranges, elevations, metadata_dict)
    assert np.allclose(
        latlonh.T[:5, :],
        np.array(
            [
                [52.74570806, 6.17862461, 41.99994496],
                [52.76170338, 6.27177266, 41.99994498],
                [52.76660854, 6.26778407, 41.99994498],
                [52.76459789, 6.17667703, 41.99994498],
                [52.78162813, 6.2688723, 41.999945],
            ]
        ),
    )
    assert latlonh.shape == (3, azimuths.shape[0])
