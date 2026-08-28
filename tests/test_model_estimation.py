import numpy as np
import pytest
import xarray as xr

from depsi.model_estimation import MODEL_NAMES_PARAMS, estimate_model_params

# Test constants
rng = np.random.default_rng(42)
WAVELENGTH = 0.056  # S1 wavelength in meters
m2ph = -4 * np.pi / WAVELENGTH  # meters to phase conversion factor


def stm_linear_height(n_space, n_time) -> xr.Dataset:
    """Create a simple STM dataset with offset, velocity and height models."""
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


@pytest.mark.parametrize(["n_space", "n_time", "models"], [(11, 15, None), (41, 31, ["height", "velocity", "offset"])])
def test_estimate_model_params_linear_height(n_space, n_time, models):
    """Test model parameter estimation for linear + height model."""
    stm = stm_linear_height(n_space, n_time)

    # Estimate model parameters
    stm_out, param_names = estimate_model_params(
        stm,
        models=models,
        key_observations="unw_phase",
        key_h2ph="h2ph",
        key_time="time",
    )

    # Check if all expected parameters are present
    expected_param_names = ["pnt_offset", "pnt_velocity", "pnt_height"]
    assert set(param_names) == set(expected_param_names)

    # Check if estimated parameters are close to true values
    stm_out = stm_out.compute()
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


def stm_full_models(models, n_space, n_time) -> xr.Dataset:
    """Create an STM dataset with synthetic phase built from exactly the given `models`."""
    time = np.sort(rng.uniform(-2, 2, size=n_time))
    space = np.arange(n_space)

    # True parameters
    true_offset = rng.uniform(-np.pi, np.pi, size=n_space)
    true_velocity = rng.uniform(-0.01, 0.01, size=n_space)
    true_quadratic = rng.uniform(-0.005, 0.005, size=n_space)
    true_cubic = rng.uniform(-0.001, 0.001, size=n_space)
    true_height = rng.uniform(0, 100, size=n_space)
    true_temperature = rng.uniform(-0.05, 0.05, size=n_space)
    true_cross_range = rng.uniform(-0.05, 0.05, size=n_space)
    true_seasonal_a1 = rng.uniform(-0.5, 0.5, size=n_space)
    true_seasonal_a2 = rng.uniform(-0.5, 0.5, size=n_space)

    # Space-and-time varying inputs
    h2ph = rng.normal(1e-5, 1e-4, (n_space, n_time))
    cr2ph = rng.uniform(-1, 1, (n_space, n_time))
    # Time-only varying input (same for every point, like add_stm_attributes builds it)
    temperature = rng.uniform(-10, 30, size=n_time)

    seasonal_a1_col = np.sin(2 * np.pi * time)
    seasonal_a2_col = -(np.cos(2 * np.pi * time) - 1)

    # Per-model phase contribution for one point idx, only summed in for models actually requested
    contributions = {
        "offset": lambda idx: true_offset[idx] * np.ones(n_time),
        "velocity": lambda idx: true_velocity[idx] * time * m2ph,
        "quadratic": lambda idx: true_quadratic[idx] * time**2 * m2ph,
        "cubic": lambda idx: true_cubic[idx] * time**3 * m2ph,
        "height": lambda idx: true_height[idx] * h2ph[idx, :] * m2ph,
        "temperature": lambda idx: true_temperature[idx] * temperature,
        "cross_range": lambda idx: true_cross_range[idx] * cr2ph[idx, :],
        "seasonal": lambda idx: true_seasonal_a1[idx] * seasonal_a1_col + true_seasonal_a2[idx] * seasonal_a2_col,
    }

    unw_phase = np.zeros((n_space, n_time))
    for idx in range(n_space):
        for model in models:
            unw_phase[idx, :] += contributions[model](idx)

    stm = xr.Dataset(
        {
            "unw_phase": (("space", "time"), unw_phase),
            "h2ph": (("space", "time"), h2ph),
            "cr2ph": (("space", "time"), cr2ph),
            "temperature": (("time",), temperature),
            "true_offset": (("space",), true_offset),
            "true_velocity": (("space",), true_velocity),
            "true_quadratic": (("space",), true_quadratic),
            "true_cubic": (("space",), true_cubic),
            "true_height": (("space",), true_height),
            "true_temperature": (("space",), true_temperature),
            "true_cross_range": (("space",), true_cross_range),
            "true_seasonal_a1": (("space",), true_seasonal_a1),
            "true_seasonal_a2": (("space",), true_seasonal_a2),
        },
        coords={
            "space": space,
            "time": time,
        },
        attrs={"wavelength": WAVELENGTH},
    )

    return stm


# Maps each output parameter name to the "true_*" variable holding its ground truth in stm_full_models
TRUE_VALUE_KEYS = {
    "pnt_offset": "true_offset",
    "pnt_velocity": "true_velocity",
    "pnt_quadratic": "true_quadratic",
    "pnt_cubic": "true_cubic",
    "pnt_height": "true_height",
    "pnt_temperature": "true_temperature",
    "pnt_cross_range": "true_cross_range",
    "pnt_seasonal_a1": "true_seasonal_a1",
    "pnt_seasonal_a2": "true_seasonal_a2",
}


@pytest.mark.parametrize(
    ["n_space", "n_time", "models"],
    [
        (11, 15, ["quadratic", "cubic", "temperature", "cross_range", "seasonal"]),
        (
            41,
            31,
            ["offset", "velocity", "quadratic", "cubic", "height", "temperature", "cross_range", "seasonal"],
        ),
    ],
)
def test_estimate_model_params_full_models(n_space, n_time, models):
    """Test model parameter estimation including every model beyond offset/velocity/height."""
    stm = stm_full_models(models, n_space, n_time)

    stm_out, param_names = estimate_model_params(
        stm,
        models=models,
        key_observations="unw_phase",
        key_h2ph="h2ph",
        key_time="time",
        key_cross_range="cr2ph",
        key_temperature="temperature",
    )

    # Check if all expected parameters are present
    expected_param_names = [name for model in models for name in MODEL_NAMES_PARAMS[model]]
    assert set(param_names) == set(expected_param_names)

    # Check if estimated parameters are close to true values
    stm_out = stm_out.compute()
    for param_name in param_names:
        np.testing.assert_allclose(
            stm_out[param_name].values,
            stm[TRUE_VALUE_KEYS[param_name]].values,
            atol=1e-8,
        )


def test_estimate_model_params_rejects_invalid_model_name():
    """Test that an unsupported model name raises a NotImplementedError."""
    with pytest.raises(NotImplementedError):
        estimate_model_params(
            stm_linear_height(3, 5),
            models=["invalid_model"],
        )


def test_estimate_model_params_rejects_empty_model_list():
    """Test that an empty model list raises a ValueError."""
    with pytest.raises(ValueError):
        estimate_model_params(
            stm_linear_height(3, 5),
            models=[],
        )


def test_estimate_model_params_dataset_behavior():
    """Test the output dataset's contract: expected variables, shapes, and coordinate preservation."""
    n_space, n_time = 2, 5
    stm = stm_linear_height(n_space, n_time)

    stm_out, param_names = estimate_model_params(
        stm,
        models=["offset", "velocity", "height"],
        key_observations="unw_phase",
        key_h2ph="h2ph",
        key_time="time",
    )
    stm_out = stm_out.compute()

    # Fast guards: expected variables are present, with the expected shapes
    assert set(param_names).issubset(set(stm_out.data_vars))
    assert {"expected_phases_yhat", "phase_residuals"}.issubset(set(stm_out.data_vars))
    for param_name in param_names:
        assert stm_out[param_name].shape == (n_space,)
    assert stm_out["expected_phases_yhat"].shape == (n_space, n_time)
    assert stm_out["phase_residuals"].shape == (n_space, n_time)

    # Coordinate preservation
    np.testing.assert_array_equal(stm_out["space"].values, stm["space"].values)
    np.testing.assert_array_equal(stm_out["time"].values, stm["time"].values)

    # Semantic invariant tied to the documented contract: phase_residuals = unw_phase - expected_phases_yhat
    np.testing.assert_allclose(
        stm_out["phase_residuals"].values,
        stm_out["unw_phase"].values - stm_out["expected_phases_yhat"].values,
        atol=1e-8,
    )
