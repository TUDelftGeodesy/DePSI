import numpy as np
import pytest
import xarray as xr

from depsi.model_definition import estimate_model_params

# Test constants
rng = np.random.default_rng(42)
WAVELENGTH = 0.056  # S1 wavelength in meters
m2ph = -4 * np.pi / WAVELENGTH  # meters to phase conversion factor


def stm_linear_height(n_space, n_time) -> xr.Dataset:
    """Create a simple STM dataset with linear and height models."""
    time = np.sort(rng.uniform(-2, 2, size=n_time))
    # ordering time
    space = np.arange(n_space)

    # True parameters
    true_offset = rng.uniform(-np.pi, np.pi, size=n_space)
    true_velocity = rng.uniform(-0.01, 0.01, size=n_space)
    true_height = rng.uniform(0, 100, size=n_space)

    # Height to phase conversion factor
    h2ph = rng.normal(1e-5, 1e-4, (n_space, n_time))

    # Generate unwrapped phase data
    unw_phase = np.zeros((n_space, n_time))
    for idx in range(n_space):
        unw_phase[idx, :] = true_offset[idx] + true_velocity[idx] * time * m2ph + true_height[idx] * h2ph[idx, :] * m2ph

    stm = xr.Dataset(
        {
            "unw_phase": (("space", "time"), unw_phase),
            "h2ph": (("space", "time"), h2ph),
            "true_offset": (("space"), true_offset),
            "true_velocity": (("space"), true_velocity),
            "true_height": (("space"), true_height),
        },
        coords={
            "space": space,
            "time": time,
        },
        attrs={"wavelength": WAVELENGTH},
    )

    return stm


@pytest.mark.parametrize(["n_space", "n_time"], [(11, 15), (41, 31)])
def test_estimate_model_params_linear_height(n_space, n_time):
    """Test model parameter estimation for linear + height model."""
    stm = stm_linear_height(n_space, n_time)

    # Estimate model parameters
    stm_out = estimate_model_params(
        stm,
        models=["height", "linear"],
        key_unw_phase="unw_phase",
        key_h2ph="h2ph",
        key_time="time",
    ).compute()

    # Check if estimated parameters are close to true values
    np.testing.assert_allclose(
        stm_out["pnt_offset"].values,
        stm["true_offset"].values,
        atol=1e-8,
    )
    np.testing.assert_allclose(
        stm_out["pnt_velocity"].values,
        stm["true_velocity"].values,
        atol=1e-8,
    )
    np.testing.assert_allclose(
        stm_out["pnt_height"].values,
        stm["true_height"].values,
        atol=1e-8,
    )


def test_estimate_model_params_invalid_model():
    """Test that an invalid model name raises a ValueError."""
    with pytest.raises(NotImplementedError):
        estimate_model_params(
            stm_linear_height(3, 5),
            models=["invalid_model"],
        )
