import dask.array as da
import numpy as np
import pytest
import xarray as xr

from depsi.arc_estimation import _build_periodogram_search_space, _chunk_for_temp_coh_compute, periodogram
from depsi.utils import wrap_phase

WAVELENGTH_S1 = 0.055465763  # m, sentinel-1 wavelength used for testing


def get_test_consts(n_obs, n_arcs, velo_min, velo_max, height_min, height_max):
    """function to get constants for testing"""
    m2ph = -4 * np.pi / WAVELENGTH_S1

    rng = np.random.default_rng(42)  # reset every time for reproducibility
    h2ph = rng.random((n_obs, n_arcs)) * 1e-3  # fixed
    velo = (velo_min - velo_max) * rng.random((n_arcs,)) + velo_max  # uniform distribution, in meters per year
    height = (height_min - height_max) * rng.random((n_arcs,)) + height_max  # uniform distribution, in meters
    return m2ph, n_obs, n_arcs, velo, height, h2ph


def get_arcs_stm(n_obs, n_arcs, velo_min, velo_max, height_min, height_max):
    """function to get arcs stm for testing."""
    TIME_STEP = 15  # Time step in days
    m2ph, n_obs, n_arcs, velo, height, h2ph = get_test_consts(n_obs, n_arcs, velo_min, velo_max, height_min, height_max)
    years = np.linspace(0, (n_obs - 1) * TIME_STEP / 365.25, n_obs)  # years from 0 to n_obs-1
    phs_velo = (np.expand_dims(m2ph * years, axis=1) * velo).T
    phs_height = (m2ph * h2ph @ np.diag(height)).T
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
        },
        attrs={
            "wavelength": WAVELENGTH_S1,  # wavelength in meters
        },
    )

    return stm_arcs


@pytest.mark.parametrize(
    "n_obs, n_arcs, velo_min, velo_max, height_min, height_max",
    [
        (13, 4, -1e-3, 1e-4, -1, 1),  # stable point
        (6, 6, -2.3e-3, 1.5e-4, -0.03, 0.02),  # stable point, extra short time series, low height
        (77, 3, -5e-3, 2e-4, -5, 10),  # non-stable point, long time series, high height
    ],
)
def test_periodogram(n_obs, n_arcs, velo_min, velo_max, height_min, height_max):
    arcs = get_arcs_stm(n_obs, n_arcs, velo_min, velo_max, height_min, height_max)
    std_height = 5  # standard deviation for height
    std_vel = 0.01  # standard deviation for velocity

    results = periodogram(
        stm=arcs,
        key_dphase="phs_obs_wrapped",
        key_Btemporal="years",
        key_h2ph="h2ph_values",
        std_height=std_height,
        std_vel=std_vel,
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

    # Solved height and velocity should be close to the true values, within the standard deviation
    assert np.allclose(results[2].values, arcs["height"].values, atol=std_height)
    assert np.allclose(results[3].values, arcs["velo"].values, atol=std_vel)


@pytest.mark.parametrize("chunk_space, chunk_time", [(10, -1), (10, 5)])
def test_periodogram_chunk(chunk_space, chunk_time):
    n_obs = 21
    n_arcs = 22
    velo_min = -5e-3
    velo_max = 2e-4
    height_min = -7
    height_max = 5
    arcs = get_arcs_stm(n_obs, n_arcs, velo_min, velo_max, height_min, height_max)
    arcs = arcs.chunk({"space": chunk_space, "time": chunk_time})
    std_height = 5  # standard deviation for height
    std_vel = 0.01  # standard deviation for velocity

    results = periodogram(
        stm=arcs,
        key_dphase="phs_obs_wrapped",
        key_Btemporal="years",
        key_h2ph="h2ph_values",
        std_height=std_height,
        std_vel=std_vel,
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

    # Solved height and velocity should be close to the true values, within the standard deviation
    assert np.allclose(results[2].values, arcs["height"].values, atol=std_height)
    assert np.allclose(results[3].values, arcs["velo"].values, atol=std_vel)


def test_periodogram_no_wavelength():
    arcs = get_arcs_stm(13, 4, -1e-3, 1e-4, -1, 1)
    arcs.attrs.pop("wavelength", None)  # remove wavelength to test without it

    # This should raise an error because wavelength is required
    with pytest.raises(ValueError):
        _ = periodogram(
            stm=arcs,
            key_dphase="phs_obs_wrapped",
            key_Btemporal="years",
            key_h2ph="h2ph_values",
        )


def test_periodogram_with_coh_mask_partial_epochs():
    n_obs = 13
    n_arcs = 4
    velo_min = -1e-3
    velo_max = 1e-4
    height_min = -1
    height_max = 1
    arcs = get_arcs_stm(n_obs, n_arcs, velo_min, velo_max, height_min, height_max)

    # Alternate coherent/incoherent epochs for every arc.
    coh_mask = np.zeros((n_arcs, n_obs), dtype=bool)
    coh_mask[:, ::2] = True
    arcs["dt_coh_mask"] = (("space", "time"), coh_mask)

    results = periodogram(
        stm=arcs,
        key_dphase="phs_obs_wrapped",
        key_Btemporal="years",
        key_h2ph="h2ph_values",
        key_coh_mask="dt_coh_mask",
        std_height=5,
        std_vel=0.01,
        init_step_height=abs(height_max - height_min) / 10,
        init_step_vel=abs(velo_max - velo_min) / 10,
        init_height=(height_min + height_max) / 2,
        init_vel=(velo_min + velo_max) / 2,
    )

    unwrapped = results[0].values
    ambiguities = results[1].values
    est_height = results[2].values
    est_vel = results[3].values
    coherence = results[4].values

    assert unwrapped.shape == (n_arcs, n_obs)
    assert ambiguities.shape == (n_arcs, n_obs)

    # Incoherent epochs must remain NaN in outputs with time dimension.
    assert np.all(np.isnan(unwrapped[:, 1::2]))
    assert np.all(np.isnan(ambiguities[:, 1::2]))

    # Coherent epochs should be solved.
    assert np.all(np.isfinite(unwrapped[:, ::2]))
    assert np.all(np.isfinite(ambiguities[:, ::2]))
    assert np.all(np.isfinite(est_height))
    assert np.all(np.isfinite(est_vel))
    assert np.all(np.isfinite(coherence))


def test_periodogram_with_coh_mask_time_only_broadcast_to_space():
    n_obs = 9
    n_arcs = 3
    arcs = get_arcs_stm(n_obs, n_arcs, -1e-3, 1e-4, -1, 1)

    # Time-only mask should broadcast to all arcs in space.
    mask_time = np.ones(n_obs, dtype=bool)
    mask_time[2] = False
    mask_time[7] = False
    arcs["dt_coh_mask_time"] = (("time",), mask_time)

    results = periodogram(
        stm=arcs,
        key_dphase="phs_obs_wrapped",
        key_Btemporal="years",
        key_h2ph="h2ph_values",
        key_coh_mask="dt_coh_mask_time",
    )

    unwrapped = results[0].values
    ambiguities = results[1].values

    assert np.all(np.isnan(unwrapped[:, 2]))
    assert np.all(np.isnan(unwrapped[:, 7]))
    assert np.all(np.isnan(ambiguities[:, 2]))
    assert np.all(np.isnan(ambiguities[:, 7]))


def test_periodogram_with_coh_mask_missing_key_raises():
    arcs = get_arcs_stm(13, 4, -1e-3, 1e-4, -1, 1)

    with pytest.raises(ValueError, match="Coherence mask variable"):
        _ = periodogram(
            stm=arcs,
            key_dphase="phs_obs_wrapped",
            key_Btemporal="years",
            key_h2ph="h2ph_values",
            key_coh_mask="not_present",
        )


def test_periodogram_with_coh_mask_all_false_returns_nan():
    n_obs = 11
    n_arcs = 2
    arcs = get_arcs_stm(n_obs, n_arcs, -1e-3, 1e-4, -1, 1)
    arcs["dt_coh_mask"] = (("space", "time"), np.zeros((n_arcs, n_obs), dtype=bool))

    results = periodogram(
        stm=arcs,
        key_dphase="phs_obs_wrapped",
        key_Btemporal="years",
        key_h2ph="h2ph_values",
        key_coh_mask="dt_coh_mask",
    )

    assert np.all(np.isnan(results[0].values))
    assert np.all(np.isnan(results[1].values))
    assert np.all(np.isnan(results[2].values))
    assert np.all(np.isnan(results[3].values))
    assert np.all(np.isnan(results[4].values))


def test_build_periodogram_search_space():
    """Test the build_search_space function."""

    # Test with a simple case
    height_center = 5
    vel_center = 0.2
    step_height = 1
    step_vel = 0.1
    num_height_search = 2
    num_vel_search = 3

    expect_vels = np.array([-0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
    expect_heights = np.array([3, 4, 5, 6, 7])
    expect_search_space = np.array(np.meshgrid(expect_heights, expect_vels)).T.reshape(-1, 2)

    search_space = _build_periodogram_search_space(
        height_center, vel_center, step_height, step_vel, num_height_search, num_vel_search
    )

    # Sort the search space for comparison
    expect_search_space = expect_search_space[np.lexsort((expect_search_space[:, 1], expect_search_space[:, 0]))]
    search_space = search_space[np.lexsort((search_space[:, 1], search_space[:, 0]))]

    assert np.allclose(search_space, expect_search_space)


@pytest.mark.parametrize("n_arcs, n_obs, n_search", [(1000, 50, 20), (3001, 111, 50)])
def test_chunk_for_temp_coh_compute_dask(n_arcs, n_obs, n_search):
    phase = da.ones((n_arcs, n_obs))  # dummy phase data
    search_space = np.random.rand(n_search, 2)

    phase_chunked, search_space_chunked = _chunk_for_temp_coh_compute(phase, search_space)

    assert isinstance(phase_chunked, da.Array)
    assert isinstance(search_space_chunked, da.Array)
    assert phase_chunked.chunks == phase.chunks  # phase chunking should not change


@pytest.mark.parametrize("n_arcs, n_obs, n_search", [(11, 7, 21)])
def test_chunk_for_temp_coh_compute_np(n_arcs, n_obs, n_search):
    phase = np.ones((n_arcs, n_obs))  # dummy phase data
    search_space = np.random.rand(n_search, 2)

    phase_chunked, search_space_chunked = _chunk_for_temp_coh_compute(phase, search_space)

    # The np arrays should be converted to dask arrays
    assert isinstance(phase_chunked, da.Array)
    assert isinstance(search_space_chunked, da.Array)
