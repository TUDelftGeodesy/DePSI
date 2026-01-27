import numpy as np
import numpy.random as nr
import xarray as xr

from depsi.point_quality import compute_spatiotemporal_consistency


def test_stc():
    rng = nr.default_rng(42)

    stm = xr.Dataset(
        data_vars={
            "amplitude": (("space", "time"), np.ones((5, 10))),
            "unwrapped_phase": (("space", "time"), rng.random((5, 10))),
            "rd_x": (("space"), np.array([150000, 150010, 150100, 150200, 150300])),
            "rd_y": (("space"), np.array([450000, 450010, 450100, 450200, 450300])),
        },
        coords={
            "space": np.array([3, 1, 2, 5, 7]),  # non monotonic space coords
            "time": np.arange(10),
            "azimuth": (("space"), np.arange(5)),
            "range": (("space"), np.arange(5)),
        },
    )
    stm_stc = compute_spatiotemporal_consistency(
        stm=stm,
        min_dist=50,
        max_dist=200,
        x_crd_layer_name="rd_x",
        y_crd_layer_name="rd_y",
        coordinate_type="euclidean",
    )

    expected_stc = [623.54053157, 572.23892935, 572.23892935, 518.81395402, 518.81395402]
    # very high since TS are random 0-1 in meters

    assert np.allclose(stm_stc.stc.values, expected_stc)
    assert stm_stc.stc.values.shape[0] == stm.sizes["space"]
    assert stm_stc.sizes["space"] == stm.sizes["space"]
    assert stm_stc.sizes["time"] == stm.sizes["time"]


def test_stc_geographic():
    rng = nr.default_rng(42)

    stm = xr.Dataset(
        data_vars={
            "amplitude": (("space", "time"), np.ones((5, 10))),
            "unwrapped_phase": (("space", "time"), rng.random((5, 10))),
            "lon": (("space"), np.array([5.31433, 5.31448, 5.31579, 5.31725, 5.31870])),
            "lat": (("space"), np.array([52.03831, 52.03840, 52.03921, 52.04011, 52.04101])),
        },
        coords={
            "space": np.array([3, 1, 2, 5, 7]),  # non monotonic space coords
            "time": np.arange(10),
            "azimuth": (("space"), np.arange(5)),
            "range": (("space"), np.arange(5)),
        },
    )
    stm_stc = compute_spatiotemporal_consistency(
        stm=stm,
        min_dist=50,
        max_dist=200,
        x_crd_layer_name="lon",
        y_crd_layer_name="lat",
        coordinate_type="geographic",
    )

    expected_stc = [623.54053157, 572.23892935, 572.23892935, 518.81395402, 518.81395402]
    # very high since TS are random 0-1 in meters

    assert np.allclose(stm_stc.stc.values, expected_stc)
    assert stm_stc.stc.values.shape[0] == stm.sizes["space"]
    assert stm_stc.sizes["space"] == stm.sizes["space"]
    assert stm_stc.sizes["time"] == stm.sizes["time"]
