"""test_classification.py"""

import dask.array as da
import numpy as np
import pytest
import xarray as xr

from pydepsi.classification import _idx_within_distance, _nad_block, _nmad_block, ps_selection

# Create a random number generator
rng = np.random.default_rng(42)


def test_ps_seletion_nad():
    slcs = xr.Dataset(
        data_vars={"amplitude": (("azimuth", "range", "time"), np.ones((10, 10, 10)))},
        coords={"azimuth": np.arange(10), "range": np.arange(10), "time": np.arange(10)},
    )
    res = ps_selection(slcs, 0.5, method="nad", output_chunks=5)
    assert res.sizes["time"] == 10
    assert res.sizes["space"] == 100
    assert "pnt_nad" in res
    assert "azimuth" in res
    assert "range" in res
    assert "space" in res.dims
    assert "time" in res.dims
    assert isinstance(res["pnt_nad"].data, da.core.Array)


def test_ps_seletion_nmad():
    slcs = xr.Dataset(
        data_vars={"amplitude": (("azimuth", "range", "time"), np.ones((10, 10, 10)))},
        coords={"azimuth": np.arange(10), "range": np.arange(10), "time": np.arange(10)},
    )
    res = ps_selection(slcs, 0.5, method="nmad", output_chunks=5)
    assert res.sizes["time"] == 10
    assert res.sizes["space"] == 100
    assert "pnt_nmad" in res
    assert "azimuth" in res
    assert "range" in res
    assert "space" in res.dims
    assert "time" in res.dims
    assert isinstance(res["pnt_nmad"].data, da.core.Array)


def test_ps_seletion_nad_mempersist():
    """When mem_persist=True, results should be a numpy array."""
    slcs = xr.Dataset(
        data_vars={"amplitude": (("azimuth", "range", "time"), np.ones((10, 10, 10)))},
        coords={"azimuth": np.arange(10), "range": np.arange(10), "time": np.arange(10)},
    )
    res = ps_selection(slcs, 0.5, method="nad", output_chunks=5, mem_persist=True)
    assert isinstance(res["pnt_nad"].data, np.ndarray)


def test_ps_seletion_nmad_mempersist():
    """When mem_persist=True, results should be a numpy array."""
    slcs = xr.Dataset(
        data_vars={"amplitude": (("azimuth", "range", "time"), np.ones((10, 10, 10)))},
        coords={"azimuth": np.arange(10), "range": np.arange(10), "time": np.arange(10)},
    )
    res = ps_selection(slcs, 0.5, method="nmad", output_chunks=5, mem_persist=True)
    assert isinstance(res["pnt_nmad"].data, np.ndarray)


def test_ps_seletion_not_implemented():
    slcs = xr.Dataset(
        data_vars={"amplitude": (("azimuth", "range", "time"), np.ones((10, 10, 10)))},
        coords={"azimuth": np.arange(10), "range": np.arange(10), "time": np.arange(10)},
    )
    # catch not implemented method
    with pytest.raises(NotImplementedError):
        ps_selection(slcs, 0.5, method="not_implemented", output_chunks=5)


def test_nad_block_zero_dispersion():
    """NAD for a constant array should be zero."""
    slcs = xr.DataArray(
        data=np.ones((10, 10, 10)),
        dims=("azimuth", "range", "time"),
        coords={"azimuth": np.arange(10), "range": np.arange(10), "time": np.arange(10)},
    )
    res = _nad_block(slcs)
    assert res.shape == (10, 10)
    assert np.all(res == 0)


def test_nmad_block_zero_dispersion():
    """NMAD for a constant array should be zero."""
    slcs = xr.DataArray(
        data=np.ones((10, 10, 10)),
        dims=("azimuth", "range", "time"),
        coords={"azimuth": np.arange(10), "range": np.arange(10), "time": np.arange(10)},
    )
    res = _nmad_block(slcs)
    assert res.shape == (10, 10)
    assert np.all(res == 0)


def test_nad_block_select_two():
    """Should select two pixels with zero dispersion."""
    amp = rng.random((10, 10, 10))  # Random amplitude data
    amp[0, 0:2, :] = 1.0  # Two pixels with constant amplitude
    slcs = xr.Dataset(
        data_vars={"amplitude": (("azimuth", "range", "time"), amp)},
        coords={"azimuth": np.arange(10), "range": np.arange(10), "time": np.arange(10)},
    )
    res = ps_selection(slcs, 1e-10, method="nad", output_chunks=5)  # Select pixels with dispersion lower than 0.00001
    assert res.sizes["time"] == 10
    assert res.sizes["space"] == 2


def test_nmad_block_select_two():
    """Should select two pixels with zero dispersion."""
    amp = rng.random((10, 10, 10))  # Random amplitude data
    amp[0, 0:2, :] = 1.0  # Two pixels with constant amplitude
    slcs = xr.Dataset(
        data_vars={"amplitude": (("azimuth", "range", "time"), amp)},
        coords={"azimuth": np.arange(10), "range": np.arange(10), "time": np.arange(10)},
    )
    res = ps_selection(slcs, 1e-10, method="nmad", output_chunks=5)  # Select pixels with dispersion lower than 0.00001
    assert res.sizes["time"] == 10
    assert res.sizes["space"] == 2


def test__idx_within_distance():
    coords_include = np.array([[1, 1]])
    coords_remain = np.array([[0, 0], [1, 1], [2, 2], [3, 3]])
    idx_within = _idx_within_distance(coords_include, coords_remain, 1)
    assert np.all(idx_within == np.array([1]))

    coords_include = np.array([[1, 1], [2, 2]])
    coords_remain = np.array([[0, 0], [1.1, 1], [2.2, 2], [3, 3]])
    idx_within = _idx_within_distance(coords_include, coords_remain, 1)
    assert np.all(idx_within == np.array([1, 2]))

    coords_include = np.array([[1, 1]])
    coords_remain = np.array([[0, 0], [1, 1], [2, 2], [3, 3]])
    idx_within = _idx_within_distance(coords_include, coords_remain, 2)
    assert np.all(idx_within == np.array([0, 1, 2]))


def test__idx_within_distance_no_drop():
    coords_include = np.array([[100, 100]])
    coords_remain = np.array([[0, 0], [1, 1], [2, 2], [3, 3]])
    idx_within = _idx_within_distance(coords_include, coords_remain, 1)
    assert idx_within is None
