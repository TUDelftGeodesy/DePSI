"""Definition of models used in deformation analysis.

This module contains functions to define various models to describe unwrapped phase data.
After unwrapped phase is available in a Space-Time Matrix (STM),
the function `estimate_model_params` can be used to estimate model parameters for each point in the STM.
This function assembles multiple model components, by calling corresponding A matrix construction functions
to build relevant columns in the A matrix, and performs least squares estimation of model parameters.

The parameters estimation is performed per point by applying the point-wise function `_estimate_model_params_one_point`
using `xarray.apply_ufunc`, which allows efficient processing of large datasets with Dask support.
"""

import numpy as np
import xarray as xr

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
    key_unw_phase: str = "unw_phase",
    key_h2ph: str = "h2ph",
    key_time: str = "time",
    key_temperature: str = "temperature",
    key_cr2ph: str = "cr2ph",
    wavelength: float | None = None,
) -> xr.Dataset:
    """Estimate model parameters for all points in the Space-Time Matrix.

    The model estimation is performed based on unwrapped phase values and specified model components.
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
    key_unw_phase : str, optional
        Key for unwrapped phase in the STM, by default "unw_phase"
    key_h2ph : str, optional
        Key for height-to-phase conversion factor in the STM, by default "h2ph"
    key_time : str, optional
        Key for time coordinate in the STM, by default "time"
    key_temperature : str, optional
        Key for temperature in the STM, by default "temperature"
    key_cr2ph : str, optional
        Key for conversion factor from temperature to phase in the STM, by default "cr2ph"
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
        # model should not be empty
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
    # It will be split into corresponding points of the unwrapped phase, each with (time,) dimension
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
        stm[key_unw_phase],
        *st_args,
        input_core_dims=[["time"]] * (len(st_args) + 1),  # all args plus unw_phase have core dim "time"
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
    unw_phase: np.ndarray,
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
    params, _, _, _ = np.linalg.lstsq(A, unw_phase, rcond=None)

    return tuple(params)


def a_height(h2ph, m2ph):
    """Create the A matrix to estimate point height.

    Parameters
    ----------
    h2ph : np.ndarray
        Array of height-related values. Shape (n_epochs, ).
    m2ph : float
        Meter to phase conversion factor.

    Returns
    -------
    np.ndarray
        The A matrix, which is the reshaped height array multiplied by the given multiplier. Shape (n, 1).

    """
    A_h2ph = (np.reshape(h2ph, (len(h2ph), 1))) * m2ph

    return A_h2ph


def a_linear(time, m2ph):
    """Generate the A matrix for a linear model based on time.

    Parameters
    ----------
    time : np.ndarray
        Array of time values (in years). Shape (n, ).
    m2ph : float
        Meter to phase conversion factor.

    Returns
    -------
    np.ndarray
        The A matrix, a 2-column matrix with ones in the first column and the time in the second column. Shape (n, 2).
    """
    A = np.ones((len(time), 2))
    A[:, 1] = time * m2ph

    return A


def a_seasonal(t):
    """Construct the A matrix to estimate a seasonal component, a linear velocity, and an offset.

    This function constructs an A matrix based on a seasonal model. The seasonal component is modeled as
    a sine and cosine function of time (in years).
    The construction is based on Eq. 4.17 from the thesis by Freek.

    The function estimates the seasonal component using the following equation:
    y = a1 * D1 + a2 * D2
    where:
        a1 = sin(2*pi*t)
        D1 = A * cos(2*pi*t_0) (t_0 is the offset of the phase)
        a2 = -(cos(2*pi*t) - 1)
        D2 = A * sin(2*pi*t_0)

    Parameters
    ----------
    t : np.ndarray
        Array of time values (in years). Shape (n, ).

    Returns
    -------
    np.ndarray
        The A matrix, which is a 2-column matrix used in the regression. Shape (n, 2).
    """
    A_1 = np.sin(2 * np.pi * t)
    A_2 = -(np.cos(2 * np.pi * t) - 1)

    A = np.zeros((len(t), 2))
    A[:, 0] = A_1
    A[:, 1] = A_2

    return A


def a_temperature(temp):
    """Create the A matrix for temperature data.

    Parameters
    ----------
    temp : np.ndarray
        Array of temperature values. Shape (n, ).

    Returns
    -------
    np.ndarray
        The A matrix, which is a reshaped column vector of temperature values. Shape (n, 1).

    """
    A_temp = np.reshape(temp, (len(temp), 1))

    return A_temp


def a_cross_range(b_cr):
    """Create the A matrix to estimate the cross-range component (in radians).

    This function constructs the A matrix for the cross-range component by reshaping the input `b_cr` array
    into a column vector.

    Parameters
    ----------
    b_cr : np.ndarray
        Array of cross-range values (in radians). Shape (n, ).

    Returns
    -------
    np.ndarray
        The A matrix, which is the reshaped cross-range array. Shape (n, 1).

    """
    A_cr = np.reshape(b_cr, (len(b_cr), 1))

    return A_cr
