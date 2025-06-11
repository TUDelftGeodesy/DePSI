"""A module to estimate atmosphere signal from network STMs with unwrapped
phases.

This estimation is done in two steps in general:
1. A temporal filtering applied per point to extract the temporal high frequency
signal.
2. A spatial kriging filtering per epoch to estimate the atmosphere signal per
epoch.
"""


import numpy as np
import xarray as xr
from scipy import signal
from scipy.ndimage import convolve1d


def estimate_non_linear_deformation(
        residual_phase: xr.DataArray,
        baseline_years: xr.DataArray,
        filter_length: int,
        method='block'
    ):
    """Estimation of non-linear deformation by filtering

    Function to apply a low-pass filter to the time series to remove the
    non-linear deformation.

    Parameters
    ----------
    residual_phase: xr.DataArray
        The residual phase time series to apply the filter to.
    baseline_years: xr.DataArray
        The baseline years corresponding to the time series.
    filter_length: int
        Length of the filter (year) to apply a low-pass filter to the time series.
    method: str
        Method to use for the estimation, e.g. 'block', 'triangle', or 'gaussian'.
    Returns
    -------

    """
    # TODO: check why * 1000.0
    # TODO: check why half width is used

    # Build low_pass filter
    if method == 'block':
        low_pass_filter = np.zeros(baseline_years.size)
        no_points = int(round(0.5 * filter_length)) + 1
        low_pass_filter[:no_points] = signal.windows.boxcar(no_points)
    elif method == 'triangle':
        low_pass_filter = np.zeros(baseline_years.size)
        no_points = int(round(0.5 * filter_length)) + 1
        low_pass_filter[:no_points] = signal.windows.triang(no_points)
    elif method == 'gaussian':
        no_points = baseline_years.size
        low_pass_filter = signal.windows.gaussian(no_points, std=0.5 * filter_length / 3)
    else:
        raise NotImplementedError(
            f"Method {method} is not implemented. "
            "Available methods are: 'block', 'triangle', 'gaussian'."
        )

    # TODO check baseline_years is monotonic
    # TODO check baseline_years size is equal to residual_phase size
    # Distances in time
    distances_matrix = np.abs(
        baseline_years.values[:, None] - baseline_years.values[None, :]
    )

    # Create a low-pass filter and compute weights matrix
    weights_matrix = low_pass_filter[distances_matrix.astype(int)]
    weights_matrix /= weights_matrix.sum(axis=1, keepdims=True)  # Normalize weights

    # Apply the low-pass filter and return the non-linear deformation
    return np.dot(weights_matrix, residual_phase.values)


def krige_in_space():
    """Kriging in space to estimate the atmosphere signal
    """


def estimate_atmosphere_phase(stm, atmosphere_master):
    """Estimate the atmosphere phase from the network STMs with unwrapped phases.

    This function applies a temporal filter to extract the high-frequency
    atmospheric signal and then uses spatial kriging to estimate the atmospheric
    phase per epoch.
    """
    # Step 1: Apply temporal filtering to extract high-frequency atmospheric signal
    non_linear = estimate_non_linear_deformation(
        residual_phase=stm["d_phase"],  # Placeholder for actual data
        baseline_years=stm["time"],  # Placeholder for actual data
        filter_length=5,  # Example filter length in years
        method='block'  # Example method
    )
    atmosphere_slave = stm["d_phase"] - non_linear + atmosphere_master

    # Step 2: Apply spatial kriging to estimate atmospheric phase per epoch
    krige_in_space()
