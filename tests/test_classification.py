"""test_classification.py"""

import dask.array as da
import numpy as np
import pytest
import xarray as xr

from depsi.classification import _idx_within_distance, _nad_block, _nmad_block, network_stm_selection, ps_selection

# Create a random number generator
rng = np.random.default_rng(42)


def test_ps_selection_nad():
    slcs = xr.Dataset(
        data_vars={"amplitude": (("azimuth", "range", "time"), np.ones((10, 7, 11)))},
        coords={
            "azimuth": np.arange(10),
            "range": np.arange(7),
            "time": [np.datetime64(f"2015-01-{i:0>2d}") for i in range(1, 12)],
        },
    )
    res = ps_selection(slcs, 0.5, method="nad", output_chunks=5)
    assert res.sizes["time"] == 11
    assert res.sizes["space"] == 70
    assert "time_selection_nad" in res
    assert "azimuth" in res
    assert "range" in res
    assert "space" in res.dims
    assert "time" in res.dims
    assert isinstance(res["time_selection_nad"].data, da.core.Array)


def test_ps_selection_nmad():
    slcs = xr.Dataset(
        data_vars={"amplitude": (("azimuth", "range", "time"), np.ones((10, 9, 11)))},
        coords={
            "azimuth": np.arange(10),
            "range": np.arange(9),
            "time": [np.datetime64(f"2015-01-{i:0>2d}") for i in range(1, 12)],
        },
    )
    res = ps_selection(slcs, 0.5, method="nmad", output_chunks=5)
    assert res.sizes["time"] == 11
    assert res.sizes["space"] == 90
    assert "time_selection_nmad" in res
    assert "azimuth" in res
    assert "range" in res
    assert "space" in res.dims
    assert "time" in res.dims
    assert isinstance(res["time_selection_nmad"].data, da.core.Array)


def test_ps_selection_nad_mempersist():
    """When mem_persist=True, results should be a numpy array."""
    slcs = xr.Dataset(
        data_vars={"amplitude": (("azimuth", "range", "time"), np.ones((10, 7, 11)))},
        coords={
            "azimuth": np.arange(10),
            "range": np.arange(7),
            "time": [np.datetime64(f"2015-01-{i:0>2d}") for i in range(1, 12)],
        },
    )
    res = ps_selection(slcs, 0.5, method="nad", output_chunks=5, mem_persist=True)
    assert isinstance(res["time_selection_nad"].data, np.ndarray)


def test_ps_selection_nmad_mempersist():
    """When mem_persist=True, results should be a numpy array."""
    slcs = xr.Dataset(
        data_vars={"amplitude": (("azimuth", "range", "time"), np.ones((10, 9, 13)))},
        coords={
            "azimuth": np.arange(10),
            "range": np.arange(9),
            "time": [np.datetime64(f"2015-01-{i:0>2d}") for i in range(1, 14)],
        },
    )
    res = ps_selection(slcs, 0.5, method="nmad", output_chunks=5, mem_persist=True)
    assert isinstance(res["time_selection_nmad"].data, np.ndarray)


def test_ps_selection_not_implemented():
    slcs = xr.Dataset(
        data_vars={"amplitude": (("azimuth", "range", "time"), np.ones((10, 5, 7)))},
        coords={
            "azimuth": np.arange(10),
            "range": np.arange(5),
            "time": [np.datetime64(f"2015-01-{i:0>2d}") for i in range(1, 8)],
        },
    )
    # catch not implemented method
    with pytest.raises(NotImplementedError):
        ps_selection(slcs, 0.5, method="not_implemented", output_chunks=5)


def test_network_stm_selection_results():
    stm = xr.Dataset(
        data_vars={
            "amplitude": (("space", "time"), np.ones((100, 10))),
            "time_selection_nad": (("space"), np.linspace(0, 1, 100)),
            "time_selection_nmad": (("space"), np.linspace(0, 1, 100)),
        },
        coords={
            "azimuth": (("space"), np.arange(100)),
            "range": (("space"), np.arange(100)),
            "time": np.arange(10),
            "space": np.arange(100),
        },
    )
    res_nad = network_stm_selection(
        stm, min_dist=20, sortby_var="time_selection_nad", azimuth_spacing=10, range_spacing=10
    )
    res_nmad = network_stm_selection(
        stm, min_dist=20, sortby_var="time_selection_nmad", azimuth_spacing=10, range_spacing=10
    )
    # Fields should remain the same
    assert "pnt_nad" in res_nad
    assert "azimuth" in res_nad
    assert "range" in res_nad
    assert "space" in res_nad.dims
    assert "time" in res_nad.dims
    # Dimensions should be half
    assert res_nad.sizes["space"] == 50
    assert res_nad.sizes["time"] == 10
    assert res_nmad.sizes["space"] == 50
    assert res_nmad.sizes["time"] == 10


def test_network_stm_selection_quality():
    stm = xr.Dataset(
        data_vars={
            "amplitude": (("space", "time"), np.ones((5, 10))),
            "time_selection_nad": (("space"), np.array([0.9, 0.01, 0.9, 0.9, 0.01])),
            "time_selection_nmad": (("space"), np.array([0.01, 0.9, 0.9, 0.9, 0.01])),
        },
        coords={
            "space": np.array([3, 1, 2, 5, 7]),  # non monotonic space coords
            "time": np.arange(10),
            "azimuth": (("space"), np.arange(5)),
            "range": (("space"), np.arange(5)),
        },
    )
    res_nad = network_stm_selection(
        stm, min_dist=3, sortby_var="time_selection_nad", azimuth_spacing=1, range_spacing=1
    )
    res_nmad = network_stm_selection(
        stm, min_dist=3, sortby_var="time_selection_nmad", azimuth_spacing=1, range_spacing=1
    )

    # The two pixels with the lowest NAD should be selected
    assert np.all(np.isclose(res_nad["time_selection_nad"].values, 0.01, rtol=1e-09, atol=1e-09))
    assert np.all(np.isclose(res_nmad["time_selection_nmad"].values, 0.01, rtol=1e-09, atol=1e-09))
    assert np.all(res_nad["space"].values == np.array([1, 7]))
    assert np.all(res_nmad["space"].values == np.array([3, 7]))


def test_network_stm_selection_include_index():
    stm = xr.Dataset(
        data_vars={
            "amplitude": (("space", "time"), np.ones((5, 10))),
            "time_selection_nad": (("space"), np.array([0.01, 0.01, 0.9, 0.9, 0.01])),
            "time_selection_nmad": (("space"), np.array([0.01, 0.01, 0.9, 0.9, 0.01])),
        },
        coords={
            "space": np.array([1, 2, 5, 6, 7]),  # non monotonic space coords
            "time": np.arange(10),
            "azimuth": (("space"), np.arange(5)),
            "range": (("space"), np.arange(5)),
        },
    )
    res_nad = network_stm_selection(
        stm, min_dist=3, include_index=[1], sortby_var="time_selection_nad", azimuth_spacing=1, range_spacing=1
    )
    res_nmad = network_stm_selection(
        stm, min_dist=3, include_index=[1], sortby_var="time_selection_nmad", azimuth_spacing=1, range_spacing=1
    )

    # The two pixels with the lowest NAD should be selected
    assert np.all(np.isclose(res_nad["time_selection_nad"].values, 0.01, rtol=1e-09, atol=1e-09))
    assert np.all(np.isclose(res_nmad["time_selection_nmad"].values, 0.01, rtol=1e-09, atol=1e-09))
    assert np.all(res_nad["space"].values == np.array([2, 7]))
    assert np.all(res_nmad["space"].values == np.array([2, 7]))


def test_network_stm_selection_wrong_csr():
    stm = xr.Dataset(
        data_vars={
            "amplitude": (("space", "time"), np.ones((100, 10))),
            "time_selection_nad": (("space"), np.linspace(0, 1, 100)),
        },
        coords={"azimuth": (("space"), np.arange(100)), "range": (("space"), np.arange(100)), "time": np.arange(10)},
    )
    # catch not implemented method
    with pytest.raises(NotImplementedError):
        network_stm_selection(
            stm,
            min_dist=20,
            sortby_var="time_selection_nad",
            azimuth_spacing=10,
            range_spacing=10,
            crs="not_implemented",
        )


def test_nad_block_zero_dispersion():
    """NAD for a constant array should be zero."""
    slcs = xr.DataArray(
        data=np.ones((10, 7, 9)),
        dims=("azimuth", "range", "time"),
        coords={"azimuth": np.arange(10), "range": np.arange(7), "time": np.arange(9)},
    )
    res = _nad_block(slcs)
    assert res.shape == (10, 7)
    assert np.all(res == 0)


def test_nmad_block_zero_dispersion():
    """NMAD for a constant array should be zero."""
    slcs = xr.DataArray(
        data=np.ones((10, 5, 11)),
        dims=("azimuth", "range", "time"),
        coords={"azimuth": np.arange(10), "range": np.arange(5), "time": np.arange(11)},
    )
    res = _nmad_block(slcs)
    assert res.shape == (10, 5)
    assert np.all(res == 0)


def test_nad_block_select_two():
    """Should select two pixels with zero dispersion."""
    amp = rng.random((10, 7, 11))  # Random amplitude data
    amp[0, 0:2, :] = 1.0  # Two pixels with constant amplitude
    slcs = xr.Dataset(
        data_vars={"amplitude": (("azimuth", "range", "time"), amp)},
        coords={
            "azimuth": np.arange(10),
            "range": np.arange(7),
            "time": [np.datetime64(f"2015-01-{i:0>2d}") for i in range(1, 12)],
        },
    )
    res = ps_selection(slcs, 1e-10, method="nad", output_chunks=5)  # Select pixels with dispersion lower than 1e-10
    assert res.sizes["time"] == 11
    assert res.sizes["space"] == 2


def test_nmad_block_select_two():
    """Should select two pixels with zero dispersion."""
    amp = rng.random((10, 5, 15))  # Random amplitude data
    amp[0, 0:2, :] = 1.0  # Two pixels with constant amplitude
    slcs = xr.Dataset(
        data_vars={"amplitude": (("azimuth", "range", "time"), amp)},
        coords={
            "azimuth": np.arange(10),
            "range": np.arange(5),
            "time": [np.datetime64(f"2015-01-{i:0>2d}") for i in range(1, 16)],
        },
    )
    res = ps_selection(slcs, 1e-10, method="nmad", output_chunks=5)  # Select pixels with dispersion lower than 0.00001
    assert res.sizes["time"] == 15
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
    coords_include = np.array([[100, 105]])
    coords_remain = np.array([[0, 0], [1, 1], [2, 2], [3, 3]])
    idx_within = _idx_within_distance(coords_include, coords_remain, 1)
    assert idx_within is None
