import numpy as np
import sarxarray
import xarray as xr

from depsi.transformations import radar_to_latlonh


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
                [52.74080136, 6.17999631, 41.99991671],
                [52.75679655, 6.27313374, 41.99991673],
                [52.76170171, 6.26914562, 41.99991673],
                [52.75969122, 6.178049, 41.99991673],
                [52.77672134, 6.27023375, 41.99991674],
            ]
        ),
    )
    assert latlonh.shape == (3, azimuths.shape[0])
