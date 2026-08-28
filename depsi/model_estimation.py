"""Estimation of model parameters as defined by the models in depsi.model_definition.

This module contains functions to estimate model parameters for observations based on specified model components.
"""

import numpy as np
import xarray as xr

from depsi.model_definition import MODEL_NAMES_PARAMS, construct_design_matrix

# MODEL_NAMES_PARAMS now lives in depsi.model_definition (re-exported here for backwards compatibility,
# e.g. scripts/depsi_main.py imports it from this module).


def estimate_model_params(
    stm: xr.Dataset,
    models: list[str] = None,
    key_observations: str = "unw_phase",
    key_h2ph: str = "h2ph",
    key_time: str = "time",
    key_cross_range: str = "cr2ph",
    key_temperature: str = "temperature",
    wavelength: float | None = None,
) -> tuple[xr.Dataset, list]:
    """Estimate model parameters for all points in the Space-Time Matrix.

    The model estimation is performed based on observation values and specified model components.
    The function supports multiple models, which can be combined to form a comprehensive model for each point.
    By default, we use a model which combines an offset, velocity and height model.
    The estimated model parameters are added to the input STM as new variables.

    Currently supported models are:
    - "offset": point offset at reference epoch
    - "velocity": point velocity
    - "height": point height
    - "temperature": point thermal expansion coefficient
    - "cross_range": point cross-range component
    - "seasonal": point seasonal (sin/cos) component

    Parameters
    ----------
    stm : xr.Dataset
        Space-time dataset
    models : list[str], optional
        List of model names to be used for estimation.
        By default this argument is None, which means using the default models: ["velocity", "height"].
    key_observations : str, optional
        Key for the observations in the STM to be modeled, by default "unw_phase"
    key_h2ph : str, optional
        Key for height-to-phase conversion factor in the STM, by default "h2ph"
    key_time : str, optional
        Key for time coordinate in the STM, by default "time"
    key_cross_range : str, optional
        Key for the cross-range component in the STM, by default "cr2ph". Only used if "cross_range"
        is in `models`.
    key_temperature : str, optional
        Key for the temperature in the STM, by default "temperature". Only used if "temperature" is in
        `models`.
    wavelength : float | None, optional
        Wavelength used for phase conversion, by default None

    Returns
    -------
    xr.Dataset
        Dataset containing estimated model parameters for each point.
        All estimated parameters are assumed to only have "space" dimension.
        In case a model has multiple parameters, they will be split into separate variables.
        The layers `expected_phases_yhat` (y_hat=Ax_hat) and `phase_residuals` (e_hat=y-y_hat) are also
        added into the dataset, with dimensions ("space", "time")
    list
        Parameter layer names
    """
    # Get model list with a standard order as in MODEL_NAMES_PARAMS
    if models is None:
        # Default model components
        models = [
            "offset",
            "velocity",
            "height",
        ]
    else:
        # model list should not be empty
        if len(models) == 0:
            raise ValueError("At least one model must be specified.")
        # Check models should be in the keys of MODEL_NAMES_PARAMS
        for model in models:
            if model not in MODEL_NAMES_PARAMS.keys():
                raise NotImplementedError(
                    f"Model '{model}' is not supported. Available models are: {list(MODEL_NAMES_PARAMS.keys())}"
                )

        # Reorder models according to MODEL_NAMES_PARAMS
        ordered_models = []
        for model_name in MODEL_NAMES_PARAMS.keys():
            if model_name in models:
                ordered_models.append(model_name)
        models = ordered_models

    # Get wavelength from metadata if not provided
    if wavelength is None:
        if "wavelength" in stm.attrs:
            wavelength = stm.attrs["wavelength"]
        else:
            raise ValueError("Wavelength must be provided either as an argument or in the dataset attributes.")
    m2ph = -4 * np.pi / wavelength

    # Prepare space-time arguments
    # These are args with (space, time) dimensions, and "time" as the core dimension
    # It will be split into corresponding points of the observations, each with (time,) dimension
    st_args = []
    st_args_keys = []  # keep track of keys for args
    if "height" in models:
        st_args.append(stm[key_h2ph])
        st_args_keys.append("h2ph")
    if "cross_range" in models:
        st_args.append(stm[key_cross_range])
        st_args_keys.append("cross_range")

    # Prepare keyword args
    # These are arguments which will be broadcasted to each point
    kwargs = {
        "models": models,
        "m2ph": m2ph,
        "time": stm[
            key_time
        ].values,  # Time coordinate is broadcasted to each point. It's evaluated since it's 1d and small.
        "st_args_keys": st_args_keys,  # pass the keys for st_args to identify them in the function
    }
    if "temperature" in models:
        # Temperature only varies over time (not per point), so it's broadcast like "time" rather than
        # split per point via st_args.
        kwargs["temperature"] = stm[key_temperature].values

    # Organize output parameter names
    param_names = []
    for model in models:
        param_names.extend(MODEL_NAMES_PARAMS[model])

    # Apply model estimation for each point
    params, yhat, ehat = xr.apply_ufunc(
        _estimate_model_params_one_point,
        stm[key_observations],
        *st_args,
        input_core_dims=[["time"]] * (len(st_args) + 1),  # all args plus key_observations have core dim "time"
        output_core_dims=[["params"], ["time"], ["time"]],
        vectorize=True,
        dask="parallelized",
        output_dtypes=[float, float, float],
        dask_gufunc_kwargs={"output_sizes": {"params": len(param_names), "time": len(stm.time)}},
        kwargs=kwargs,
    )

    # Assign parameter names to the output dataset
    stm_out = stm.copy()
    for i, param_name in enumerate(param_names):
        stm_out = stm_out.assign({param_name: (["space"], params.isel(params=i).data)})
    stm_out = stm_out.assign(
        {"expected_phases_yhat": (["space", "time"], yhat.data), "phase_residuals": (["space", "time"], ehat.data)}
    )

    return stm_out, param_names


def _estimate_model_params_one_point(
    obs: np.ndarray,
    *st_args,
    models: list[str],
    m2ph: float,
    time: np.ndarray,
    st_args_keys: list[str],
    temperature: np.ndarray | None = None,
) -> tuple:
    """Estimate model parameters for a single point based on unwrapped phase values."""
    # st_args holds the space-time varying inputs (h2ph and/or cross_range, when requested), keyed by
    # st_args_keys. n_epochs/time cover the models that only need the (shared) time coordinate.
    # temperature is likewise shared across points (not split via st_args), so it's passed separately.
    model_inputs = dict(zip(st_args_keys, st_args, strict=True))
    model_inputs["n_epochs"] = len(time)
    model_inputs["time"] = time
    if temperature is not None:
        model_inputs["temperature"] = temperature

    A = construct_design_matrix(models, m2ph, **model_inputs)

    # Estimate model parameters using least squares
    # currently without Qyy
    # params, _, _, _ = np.linalg.lstsq(A, obs)  # rcond=None
    params = np.linalg.inv(A.T @ A) @ A.T @ obs
    y_hat = A @ params
    e_hat = obs - y_hat.reshape(obs.shape)

    return params.T, y_hat, e_hat
