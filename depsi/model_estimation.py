"""Definition of models used in deformation analysis.

This module contains functions to estimate model parameters for observations based on specified model components.
"""

import numpy as np
import xarray as xr

from depsi.model_definition import a_height, a_linear

# Look-up dictionary from model names to output model parameter names
# Only models defined here can be used in the model estimation function
# Note that the order of output parameters will follow the order in this list
MODEL_NAMES_PARAMS = {
    "linear": ["pnt_offset", "pnt_velocity"],  # offset + linear velocity
    "height": ["pnt_height"],  # point height
}


def estimate_model_params(
    stm: xr.Dataset,
    models: list[str] = None,
    key_observations: str = "unw_phase",
    key_h2ph: str = "h2ph",
    key_time: str = "time",
    wavelength: float | None = None,
) -> xr.Dataset:
    """Estimate model parameters for all points in the Space-Time Matrix.

    The model estimation is performed based on observation values and specified model components.
    The function supports multiple models, which can be combined to form a comprehensive model for each point.
    By default, we use a model which combines an offset, velocity and height model.
    The estimated model parameters are added to the input STM as new variables.

    Parameters
    ----------
    stm : xr.Dataset
        Space-time dataset
    models : list[str], optional
        List of model names to be used for estimation.
        By default this argument is None, which means using the default models: ["linear", "height"].
    key_observations : str, optional
        Key for the observations in the STM to be modeled, by default "unw_phase"
    key_h2ph : str, optional
        Key for height-to-phase conversion factor in the STM, by default "h2ph"
    key_time : str, optional
        Key for time coordinate in the STM, by default "time"
    wavelength : float | None, optional
        Wavelength used for phase conversion, by default None

    Returns
    -------
    xr.Dataset
        Dataset containing estimated model parameters for each point.
        All estimated parameters are assumed to only have "space" dimension.
        In case a model has multiple parameters, they will be splitted into separate variables.
    """
    # Get model list with a standard order as in MODEL_NAMES_PARAMS
    if models is None:
        # Default model components
        models = [
            "linear",
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

    # Prepare keyword args
    # These are arguments which will be broadcasted to each point
    kwargs = {
        "models": models,
        "m2ph": m2ph,
        "time": stm[key_time],  # Time coordinate is broadcasted to each point
        "st_args_keys": st_args_keys,  # pass the keys for st_args to identify them in the function
    }

    # Organize output parameter names
    param_names = []
    for model in models:
        param_names.extend(MODEL_NAMES_PARAMS[model])

    # Apply model estimation for each point
    params = xr.apply_ufunc(
        _estimate_model_params_one_point,
        stm[key_observations],
        *st_args,
        input_core_dims=[["time"]] * (len(st_args) + 1),  # all args plus key_observations have core dim "time"
        output_core_dims=[[]] * len(param_names),  # per point per parameter, the output is a scalar
        vectorize=True,
        dask="parallelized",
        output_dtypes=[float] * len(param_names),
        kwargs=kwargs,
    )

    # Assign parameter names to the output dataset
    stm_out = stm.copy()
    for i, param_name in enumerate(param_names):
        stm_out[param_name] = params[i]

    return stm_out, param_names


def _estimate_model_params_one_point(
    obs: np.ndarray,
    *st_args,
    models: list[str],
    m2ph: float,
    time: np.ndarray,
    st_args_keys: list[str],
) -> tuple:
    """Estimate model parameters for a single point based on unwrapped phase values."""
    # Loop though models and build A matrix
    for model in models:
        if model == "linear":
            A_model = a_linear(time, m2ph)
        elif model == "height":
            # Get index of "h2ph" in st_args based on st_args_keys
            h2ph_idx = st_args_keys.index("h2ph")
            A_model = a_height(
                st_args[h2ph_idx],
                m2ph,
            )

        # Stack A matrices horizontally
        if model == models[0]:
            A = A_model
        else:
            A = np.hstack((A, A_model))

    # Estimate model parameters using least squares
    params, _, _, _ = np.linalg.lstsq(A, obs, rcond=None)

    return tuple(params)
