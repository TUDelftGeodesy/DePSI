import numpy as np
import xarray as xr

from depsi.postprocessing import stm_point_filter


def test_stm_point_filter_two_bounds():
    stm = xr.Dataset(
        data_vars={
            "amplitude": (("space", "time"), np.ones((5, 10))),
            "time_selection_nad": (("space"), np.array([0.6, 0.01, 0.4, 0.2, 0.9])),
        },
        coords={
            "space": np.array([3, 1, 2, 5, 7]),  # non monotonic space coords
            "time": np.arange(10),
            "azimuth": (("space"), np.arange(5)),
            "range": (("space"), np.arange(5)),
        },
    )
    filter_layer = "time_selection_nad"
    vmin = 0.2
    vmax = 0.6

    res, res_rejected = stm_point_filter(stm, filter_layer, vmin=vmin, vmax=vmax, return_removed=True)
    assert res.sizes["time"] == 10
    assert res_rejected.sizes["time"] == 10
    assert res.sizes["space"] == 3
    assert res_rejected.sizes["space"] == 2
    assert np.all(vmin <= res[filter_layer].values <= vmax)


def test_stm_point_filter_two_bounds_no_return_removed():
    stm = xr.Dataset(
        data_vars={
            "amplitude": (("space", "time"), np.ones((5, 10))),
            "time_selection_nad": (("space"), np.array([0.3, 0.01, 0.5, 0.2, 0.9])),
        },
        coords={
            "space": np.array([3, 1, 2, 5, 7]),  # non monotonic space coords
            "time": np.arange(10),
            "azimuth": (("space"), np.arange(5)),
            "range": (("space"), np.arange(5)),
        },
    )
    filter_layer = "time_selection_nad"
    vmin = 0.1
    vmax = 0.7

    res = stm_point_filter(stm, filter_layer, vmin=vmin, vmax=vmax, return_removed=False)
    assert res.sizes["time"] == 10
    assert res.sizes["space"] == 3
    assert np.all(vmin <= res[filter_layer].values <= vmax)


def test_stm_point_filter_upper_bound_no_return_removed():
    stm = xr.Dataset(
        data_vars={
            "amplitude": (("space", "time"), np.ones((5, 10))),
            "time_selection_nad": (("space"), np.array([0.3, 0.01, 0.5, 0.2, 0.9])),
        },
        coords={
            "space": np.array([3, 1, 2, 5, 7]),  # non monotonic space coords
            "time": np.arange(10),
            "azimuth": (("space"), np.arange(5)),
            "range": (("space"), np.arange(5)),
        },
    )
    filter_layer = "time_selection_nad"
    vmin = None
    vmax = 0.7

    res = stm_point_filter(stm, filter_layer, vmin=vmin, vmax=vmax, return_removed=False)
    assert res.sizes["time"] == 10
    assert res.sizes["space"] == 4
    assert np.all(res[filter_layer].values <= vmax)


def test_stm_point_filter_lower_bound_no_return_removed():
    stm = xr.Dataset(
        data_vars={
            "amplitude": (("space", "time"), np.ones((5, 10))),
            "time_selection_nad": (("space"), np.array([0.3, 0.01, 0.5, 0.2, 0.9])),
        },
        coords={
            "space": np.array([3, 1, 2, 5, 7]),  # non monotonic space coords
            "time": np.arange(10),
            "azimuth": (("space"), np.arange(5)),
            "range": (("space"), np.arange(5)),
        },
    )
    filter_layer = "time_selection_nad"
    vmin = 0.4
    vmax = None

    res = stm_point_filter(stm, filter_layer, vmin=vmin, vmax=vmax, return_removed=False)
    assert res.sizes["time"] == 10
    assert res.sizes["space"] == 2
    assert np.all(vmin <= res[filter_layer].values)
