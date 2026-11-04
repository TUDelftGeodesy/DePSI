import numpy as np
import xarray as xr

from depsi.postprocessing import stm_point_filter, estimate_reference_point_noise



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
    assert np.all(vmin <= res[filter_layer].values)
    assert np.all(vmax >= res[filter_layer].values)


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
    assert np.all(vmin <= res[filter_layer].values)
    assert np.all(vmax >= res[filter_layer].values)


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


def test_stm_point_filter_no_filter():
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
    vmax = None

    res = stm_point_filter(stm, filter_layer, vmin=vmin, vmax=vmax, return_removed=False)
    assert res.sizes["time"] == 10
    assert res.sizes["space"] == stm.sizes["space"]

# def test_reference_point_noise_estimation():

#     rng = np.random.default_rng(0)
#     n_points, n_epochs = 500, 5
#     sigma = 0.1
#     offset = np.array([0.3, -0.5, 1.0, 0.0, 0.2])   # the "true" RPN
#     da = xr.Dataset(
#         data_vars={
#             "residual_phase": (("space", "time"), rng.normal(0, sigma, (n_points, n_epochs)) + offset),
#         },
#         coords={
#             "space": np.arange(n_points),
#             "time": np.arange(n_epochs),
#         },
#     )


#     # Estimate RPN using mean
#     rpn_mean = estimate_reference_point_noise(da,key_phase="residual_phase", method="mean") 
#     print("Estimated RPN (mean):", rpn_mean["reference_point_noise"].values, "True RPN:", offset)
#     assert np.abs(rpn_mean.values - offset).max() < 0.01  

#     # Estimate RPN using median
#     rpn_median = estimate_reference_point_noise(da, key_phase="residual_phase", method="median")
#     print("Estimated RPN (median):", rpn_median["reference_point_noise"].values, "True RPN:", offset)
#     assert np.abs(rpn_median.values - offset).max() < 0.01
    

def test_reference_point_noise_estimation():
    rng = np.random.default_rng(0)
    n_points, n_epochs = 500, 5
    sigma = 0.1
    offset = np.array([0.3, -0.5, 1.0, 0.0, 0.2])  # the "true" RPN
    stm = xr.Dataset(
        data_vars={
            "residual_phase": (("space", "time"), rng.normal(0, sigma, (n_points, n_epochs)) + offset),
        },
        coords={
            "space": np.arange(n_points),
            "time": np.arange(n_epochs),
        },
    )

    resultmean = estimate_reference_point_noise(stm, key_phase="residual_phase", method="mean")
    assert np.abs(resultmean["reference_point_noise"].values - offset).max() < 0.02
    resultmedian = estimate_reference_point_noise(stm, key_phase="residual_phase", method="median")
    assert np.abs(resultmedian["reference_point_noise"].values - offset).max() < 0.02

# def test_reference_point_noise_invalid_method():
#     stm = xr.Dataset({"residual_phase": (("space", "time"), np.zeros((2, 2)))})
#     with pytest.raises(ValueError):
#         estimate_reference_point_noise(stm, key_phase="residual_phase", method="mode")
