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
    half_width = 0.5 * filter_length
    if method == 'block':
        low_pass_filter = np.zeros_like(baseline_years, dtype=float)
        no_points = int(round(half_width)) + 1
        low_pass_filter[:no_points] = signal.windows.boxcar(no_points)
    elif method == 'triangle':
        low_pass_filter = np.zeros_like(baseline_years, dtype=float)
        no_points = int(round(half_width)) + 1
        low_pass_filter[:no_points] = signal.windows.triang(no_points)
    elif method == 'gaussian':
        low_pass_filter = signal.windows.gaussian(baseline_years.size, std=half_width / 3)
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
    weights_matrix = weights_matrix / weights_matrix.sum(axis=1, keepdims=True)  # Normalize weights

     # Apply the low-pass filter and return the non-linear deformation
    def apply_filter(residual_vector, weights):
        return np.dot(weights, residual_vector)

    # Function to apply filter to residual phase across time
    return xr.apply_ufunc(
        apply_filter,
        residual_phase,
        kwargs={'weights': weights_matrix},
        input_core_dims=[['time']],
        output_core_dims=[['time']],
        vectorize=True,
        dask='parallelized',
        output_dtypes=[residual_phase.dtype]
    )



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
