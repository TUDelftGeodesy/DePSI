"""Definition of models used in deformation analysis.

This module contains functions generating columns of the A matrix for different model components, a
look-up dictionary of the models these support, and a function to construct the combined A matrix for a
chosen set of models.
"""

import numpy as np

# Look-up dictionary from model names to output model parameter names.
# Only models defined here can be used by construct_design_matrix (and, by extension, by
# depsi.model_estimation.estimate_model_params and depsi.arc_estimation).
# Note that the order of output parameters follows the order in this dictionary.
MODEL_NAMES_PARAMS = {
    "offset": ["pnt_offset"],  # point offset at reference epoch
    "velocity": ["pnt_velocity"],  # point velocity
    "quadratic": ["pnt_quadratic"],  # point quadratic term
    "cubic": ["pnt_cubic"],  # point cubic term
    "height": ["pnt_height"],  # point height
    "temperature": ["pnt_temperature"],  # point thermal expansion coefficient
    "cross_range": ["pnt_cross_range"],  # point cross-range component
    # a1/a2 coefficients of the sin/cos seasonal terms, matching a_seasonal's own notation below
    "seasonal": ["pnt_seasonal_a1", "pnt_seasonal_a2"],
}

# Required keys in `model_inputs` for each model supported by construct_design_matrix.
MODEL_REQUIRED_INPUTS = {
    "offset": ["n_epochs"],
    "velocity": ["time"],
    "quadratic": ["time"],
    "cubic": ["time"],
    "height": ["h2ph"],
    "temperature": ["temperature"],
    "cross_range": ["cross_range"],
    "seasonal": ["time"],
}


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


def a_offset(n_epochs):
    """Generate the A matrix for an offset model.

    Parameters
    ----------
    n_epochs : int
        Number of epochs (observations).

    Returns
    -------
    np.ndarray
        The A matrix, which is a column vector of ones. Shape (n_epochs, 1).
    """
    return np.ones((n_epochs, 1))


def a_velocity(Btemporal, m2ph):
    """Generate the A matrix for a velocity model based on temporal baseline.

    Parameters
    ----------
    Btemporal : np.ndarray
        Array of temporal baseline values (in years). Shape (n, ).
    m2ph : float
        Meter to phase conversion factor.

    Returns
    -------
    np.ndarray
        The A matrix, which is a reshaped column vector of temporal baseline values multiplied by
        the meter to phase conversion factor. Shape (n, 1).
    """
    return np.reshape(Btemporal * m2ph, (len(Btemporal), 1))


def a_quadratic(Btemporal, m2ph):
    """Generate the A matrix for a quadratic (acceleration) model based on temporal baseline.

    Parameters
    ----------
    Btemporal : np.ndarray
        Array of temporal baseline values (in years). Shape (n, ).
    m2ph : float
        Meter to phase conversion factor.

    Returns
    -------
    np.ndarray
        The A matrix, which is a reshaped column vector of squared temporal baseline values multiplied
        by the meter to phase conversion factor. Shape (n, 1).
    """
    return np.reshape(Btemporal**2 * m2ph, (len(Btemporal), 1))


def a_cubic(Btemporal, m2ph):
    """Generate the A matrix for a cubic model based on temporal baseline.

    Parameters
    ----------
    Btemporal : np.ndarray
        Array of temporal baseline values (in years). Shape (n, ).
    m2ph : float
        Meter to phase conversion factor.

    Returns
    -------
    np.ndarray
        The A matrix, which is a reshaped column vector of cubed temporal baseline values multiplied by
        the meter to phase conversion factor. Shape (n, 1).
    """
    return np.reshape(Btemporal**3 * m2ph, (len(Btemporal), 1))


def a_seasonal(Btemporal):
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
    Btemporal : np.ndarray
        Array of time values (in years). Shape (n, ).

    Returns
    -------
    np.ndarray
        The A matrix, which is a 2-column matrix used in the regression. Shape (n, 2).
    """
    A_1 = np.sin(2 * np.pi * Btemporal)
    A_2 = -(np.cos(2 * np.pi * Btemporal) - 1)

    A = np.zeros((len(Btemporal), 2))
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


def construct_design_matrix(models: list[str], m2ph: float, **model_inputs) -> np.ndarray:
    """Construct the combined design (A) matrix for a chosen set of models.

    This centralizes the per-model A matrix construction and horizontal stacking that was previously
    hardcoded separately in depsi.model_estimation and depsi.arc_estimation. Each requested model
    contributes one block of columns, via its corresponding `a_*` function in this module, in the order
    given in `models`.

    Parameters
    ----------
    models : list[str]
        Model names to include, in the order they should be stacked. Must be a subset of the keys of
        MODEL_NAMES_PARAMS.
    m2ph : float
        Meter to phase conversion factor. Used by the "velocity", "quadratic", "cubic" and "height"
        models.
    **model_inputs
        Named arrays/values required by the requested models:
        - "offset" requires `n_epochs` (int): number of epochs (observations).
        - "velocity" requires `time` (np.ndarray): temporal baseline values (in years).
        - "quadratic" requires `time` (np.ndarray): temporal baseline values (in years). Shares the same
          input as "velocity".
        - "cubic" requires `time` (np.ndarray): temporal baseline values (in years). Shares the same
          input as "velocity".
        - "height" requires `h2ph` (np.ndarray): height-to-phase conversion factor values.
        - "temperature" requires `temperature` (np.ndarray): temperature values.
        - "cross_range" requires `cross_range` (np.ndarray): cross-range values (in radians).
        - "seasonal" requires `time` (np.ndarray): time values (in years). Shares the same input as
          "velocity".
        Only the inputs required by the requested `models` need to be provided.

    Returns
    -------
    np.ndarray
        The stacked design matrix A, with one block of columns per requested model, in the order given
        in `models`. Shape (n_epochs, n_params).

    Raises
    ------
    NotImplementedError
        If a requested model is not in MODEL_NAMES_PARAMS.
    ValueError
        If `model_inputs` is missing a key required by one of the requested `models`.
    """
    # Check that all requested models are supported
    for model in models:
        if model not in MODEL_NAMES_PARAMS:
            raise NotImplementedError(
                f"Model '{model}' is not supported. Available models are: {list(MODEL_NAMES_PARAMS.keys())}"
            )

    # Check for missing required inputs for each model
    missing_inputs = {
        model: [key for key in MODEL_REQUIRED_INPUTS[model] if key not in model_inputs] for model in models
    }
    missing_inputs = {model: keys for model, keys in missing_inputs.items() if keys}
    if missing_inputs:
        details = ", ".join(f"'{model}' requires {keys}" for model, keys in missing_inputs.items())
        raise ValueError(f"Missing required model_inputs: {details}")

    A_blocks = []
    for model in models:
        if model == "offset":
            A_blocks.append(a_offset(model_inputs["n_epochs"]))
        elif model == "velocity":
            A_blocks.append(a_velocity(model_inputs["time"], m2ph))
        elif model == "quadratic":
            A_blocks.append(a_quadratic(model_inputs["time"], m2ph))
        elif model == "cubic":
            A_blocks.append(a_cubic(model_inputs["time"], m2ph))
        elif model == "height":
            A_blocks.append(a_height(model_inputs["h2ph"], m2ph))
        elif model == "temperature":
            A_blocks.append(a_temperature(model_inputs["temperature"]))
        elif model == "cross_range":
            A_blocks.append(a_cross_range(model_inputs["cross_range"]))
        elif model == "seasonal":
            A_blocks.append(a_seasonal(model_inputs["time"]))

    return np.hstack(A_blocks)
