import numpy as np
import pytest
import xarray as xr

from depsi.unwrap import periodogram
from depsi.utils import wrap_phase


def get_test_consts(n_obs, n_arcs, velo_min, velo_max, height_min, height_max):
    """function to get constants for testing"""
    wavelength = 0.055465763  # sentinel， in meters
    m2ph = -4 * np.pi / wavelength

    rng = np.random.default_rng(42)  # reset every time for reproducibility
    h2ph = rng.random((n_obs, n_arcs)) * 1e-3  # fixed
    velo = (velo_min - velo_max) * rng.random(
        (n_arcs,)
    ) + velo_max  # velo [-0.02, 0.005), uniform distriution, in meters per year
    height = (height_min - height_max) * rng.random(
        (n_arcs,)
    ) + height_max  # height [-1, 5), uniform distriution, in meters
    return m2ph, n_obs, n_arcs, velo, height, h2ph


def get_arcs_stm(n_obs, n_arcs, velo_min, velo_max, height_min, height_max):
    """function to get arcs stm for testing."""
    TIME_STEP = 15  # Time step in days
    m2ph, n_obs, n_arcs, velo, height, h2ph = get_test_consts(n_obs, n_arcs, velo_min, velo_max, height_min, height_max)
    years = np.linspace(0, (n_obs - 1) * TIME_STEP / 365.25, n_obs)  # years from 0 to n_obs-1
    phs_velo = (np.expand_dims(m2ph * years, axis=1) * velo).T
    phs_height = (m2ph * h2ph * height).T
    phs_obs = phs_velo + phs_height
    phs_obs_wrapped = wrap_phase(phs_obs)
    stm_arcs = xr.Dataset(
        data_vars={
            "phs_obs": (("space", "time"), phs_obs),
            "phs_obs_wrapped": (("space", "time"), phs_obs_wrapped),
            "h2ph_values": (("time", "space"), h2ph),
            "velo": (("space",), velo),
            "height": (("space",), height),
            "years": (("time",), years),
        }
    )

    return stm_arcs


@pytest.mark.parametrize(
    "n_obs, n_arcs, velo_min, velo_max, height_min, height_max",
    [
        (13, 4, -1e-3, 1e-4, -1, 1),  # stable point
        (6, 6, -2.3e-3, 1.5e-4, -0.03, 0.02),  # stable point, extra short time series, low height
        (39, 3, -5e-3, 2e-4, -5, 10),  # non-stable point, long time series, high height
    ],
)
def test_periodogram(n_obs, n_arcs, velo_min, velo_max, height_min, height_max):
    arcs = get_arcs_stm(n_obs, n_arcs, velo_min, velo_max, height_min, height_max)

    results = periodogram(
        stm=arcs,
        key_phs="phs_obs_wrapped",
        key_yeartime="years",
        key_h2ph="h2ph_values",
        std_height=5,
        std_vel=0.01,
        init_step_height=abs(height_max - height_min) / 10,
        init_step_vel=abs(velo_max - velo_min) / 10,
        init_height=(height_min + height_max) / 2,
        init_vel=(velo_min + velo_max) / 2,
    )

    m2ph, n_obs, n_arcs, velo, height, h2ph = get_test_consts(n_obs, n_arcs, velo_min, velo_max, height_min, height_max)
    assert len(results) == 5  # 5 outputs: unwrapped phs, ambiguities, height, vel, coherence
    assert results[0].shape == (n_arcs, n_obs)  # unwrapped phase
    assert results[1].shape == (n_arcs, n_obs)  # ambiguities
    assert results[2].shape == (n_arcs,)  # height
    assert results[3].shape == (n_arcs,)  # velocity
    assert results[4].shape == (n_arcs,)  # coherence

    # Unwrapped phase almost equal
    assert np.allclose(results[0].compute(), arcs["phs_obs"].values, atol=1e-5)
