"""Definition of models used in deformation analysis.

This module contains functions generating columns of the A matrix for different model components.
"""

import numpy as np


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


def a_velocity(time, m2ph):
    """Generate the A matrix for a velocity model based on time.

    Parameters
    ----------
    time : np.ndarray
        Array of time values (in years). Shape (n, ).
    m2ph : float
        Meter to phase conversion factor.

    Returns
    -------
    np.ndarray
        The A matrix, which is a reshaped column vector of time values multiplied by the given multiplier. Shape (n, 1).
    """
    return np.reshape(time * m2ph, (len(time), 1))


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
