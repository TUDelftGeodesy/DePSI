"""A module to estimate atmosphere signal from network STMs with unwrapped phases.

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
        psc_phase: xr.DataArray,
        baseline_years: xr.DataArray,
        filter_length: int,
        method='block',
        mode='mirror'
    ):
    """Apply a low-pass filter to the time series to remove the non-linear deformation.

    Parameters
    ----------
    psc_phase: xr.DataArray
        The residual phase time series to apply the filter to.
    baseline_years: xr.DataArray
        The baseline years corresponding to the time series.
    filter_length: int
        Length of the filter (year) to apply a low-pass filter to the time series.
    method: str, optional
        Method to use for building the window , e.g. 'block', 'triangle', or
        'gaussian', default is 'block', see `scipy.signal.windows` for more
    mode: str, optional
        The mode to use for the convolution, default is 'mirror', see
        `scipy.ndimage.convolve1d` for more.

    Returns
    -------
    xr.DataArray
        The non-linear deformation estimated from the time series.
    """
    # Check if baseline_years size is equal to psc_phase size
    if baseline_years.size != psc_phase["time"].size:
        raise ValueError(
            "The size of baseline_years must match the time dimension of psc_phase."
        )
    # Check baseline_years is monotonic
    is_monotonic_increasing = np.all(np.diff(baseline_years.values) >= 0)
    is_monotonic_decreasing = np.all(np.diff(baseline_years.values) <= 0)
    if not (is_monotonic_increasing or is_monotonic_decreasing):
        raise ValueError("baseline_years must be monotonic.")


    # Build the window for the low-pass filter
    std_dev = filter_length / 3
    window_size = int(6 * std_dev) | 1
    if method == 'block':
        window = signal.windows.boxcar(window_size)

    elif method == 'triangle':
        window = signal.windows.triang(window_size)

    elif method == 'gaussian':
        window = signal.windows.gaussian(window_size, std=std_dev)

    else:
        raise NotImplementedError(
            f"Method {method} is not implemented. "
            "Available methods are: 'block', 'triangle', 'gaussian'."
        )

    # normalize the window
    window = window / window.sum()

     # Apply the low-pass filter and return the non-linear deformation
    def apply_filter(data):
        return convolve1d(data, window, mode=mode)

    # Function to apply filter to residual phase across time
    return xr.apply_ufunc(
        apply_filter,
        psc_phase,
        input_core_dims=[['time']],
        output_core_dims=[['time']],
        vectorize=True,
        dask='parallelized',
        output_dtypes=[psc_phase.dtype]
    )



def krige_in_space():
    """Kriging in space to estimate the atmosphere signal."""


def estimate_atmosphere_phase(stm):
    """Estimate the atmosphere phase.

    This function applies a temporal filter to extract the high-frequency
    atmospheric signal and then uses spatial kriging to estimate the atmospheric
    phase per epoch.
    """
    # Step 1: Apply temporal filtering to extract high-frequency atmospheric signal
    non_linear = estimate_non_linear_deformation(
        psc_phase=stm["d_phase"],
        baseline_years=stm["time"],
        filter_length=5,
        method='block'
    )
    stm["atmosphere_estimated"] = stm["d_phase"] - non_linear + stm["atmosphere_base"]

    # Step 2: Apply spatial kriging to estimate atmospheric phase per epoch
