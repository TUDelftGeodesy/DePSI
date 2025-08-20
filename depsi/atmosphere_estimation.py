"""A module to estimate atmosphere signal from network STMs with unwrapped phases.

This estimation is done in two steps:
1. A temporal filtering applied per point to extract the unmodeled deformation (low-frequency signal).
signal.
2. A spatial least-squares prediction based on the residuals per epoch to estimate the atmosphere signal per
epoch.
"""

import numpy as np
import xarray as xr
from scipy import signal


def _get_signal_window_with_zero_padding(
    type,
    timespan,
    filter_length,
    sampling_rate
) -> np.ndarray:

    # Determine window of size to cover the full range of time differences
    window_size = int(timespan) + 1  # Ensure window size is an odd integer

    # Determine the core window size based on the filter length and sampling rate
    core_window_size = int(filter_length * sampling_rate) + 1

    window = np.zeros(window_size)
    start = (window_size - core_window_size) // 2
    end = start + core_window_size

    window[start:end] = signal.windows.get_window(type, core_window_size, fftbins=False)

    return window


def estimate_unmodeled_displacement(
        psc_phase_residuals: xr.DataArray,
        baseline_years: xr.DataArray,
        filter_length: int,
        sampling_rate: int = 1000,
        filter_type='gaussian',
    ) -> xr.DataArray:
    """Apply a low-pass filter to the time series to estimate the unmodeled deformation.

    Parameters
    ----------
    psc_phase_residuals: xr.DataArray
        The PSC time series residuals to apply the filter to.
    baseline_years: xr.DataArray
        The baseline years corresponding to the time series.
    filter_length: int
        Length of the filter (year) to apply a low-pass filter to the time series.
    sampling_rate: int, optional.
        Sampling rate of the time series, default is 1000.
    filter_type: str, optional
        Method to use for building the window , e.g. 'block', 'triangle', or
        'gaussian', default is 'gaussian', see `scipy.signal.windows` for more.

    Returns
    -------
    xr.DataArray
        The unmodeled deformation estimated from the time series.
    """
    # Check if baseline_years size is equal to psc_phase_residuals size
    if baseline_years.size != psc_phase_residuals["time"].size:
        raise ValueError(
            "The size of baseline_years must match the time dimension of psc_phase_residuals."
        )
    # Check baseline_years is monotonic
    is_monotonic_increasing = np.all(np.diff(baseline_years.values) >= 0)
    is_monotonic_decreasing = np.all(np.diff(baseline_years.values) <= 0)
    if not (is_monotonic_increasing or is_monotonic_decreasing):
        raise ValueError("baseline_years must be monotonic.")

    # Build the window for the low-pass filter
    baseline_scaled = baseline_years * sampling_rate
    timespan = baseline_scaled.max() - baseline_scaled.min() + 1

    if (filter_length * sampling_rate) > timespan:
        raise ValueError(
            "Filter length * sampling_rate is too large compared to the temporal span. "
            "Adjust the filter_length or sampling_rate."
        )

    window_size = int(timespan) + 1 # Ensure window size is an odd integer

    if filter_type == 'block':
        window = _get_signal_window_with_zero_padding(
            type='boxcar',
            timespan=timespan,
            filter_length=filter_length,
            sampling_rate=sampling_rate
        )

    elif filter_type == 'triangle':
        window = _get_signal_window_with_zero_padding(
            type='triang',
            timespan=timespan,
            filter_length=filter_length,
            sampling_rate=sampling_rate
        )
    elif filter_type == 'gaussian':
        std_dev = filter_length * sampling_rate / 6  # ±3σ covers the window
        window = signal.windows.gaussian(window_size, std=std_dev)
    else:
        raise NotImplementedError(
            f"Filter type {filter_type} is not implemented. "
            "Available types are: 'block', 'triangle', 'gaussian'."
        )

    # normalize the window
    window = window / window.sum()

    # Extract weights for each time difference
    time_diffs = np.subtract.outer(baseline_scaled.values, baseline_scaled.values)
    center_index = window_size // 2
    weight_indices = np.clip(center_index + np.round(time_diffs).astype(int), 0, window_size - 1)
    weight_matrix = window[weight_indices]

    # Normalize weights along axis=1 (per row)
    weight_matrix /= np.sum(weight_matrix, axis=1, keepdims=True)

    # Apply the low-pass filter and return the non-linear deformation
    def apply_filter(data):
        return np.einsum('ij,j->i', weight_matrix, data)

    return xr.apply_ufunc(
        apply_filter,
        psc_phase_residuals,
        input_core_dims=[['time']],
        output_core_dims=[['time']],
        vectorize=True,
        dask='parallelized',
        output_dtypes=[psc_phase_residuals.dtype]
    )


def krige_in_space():
    """Kriging in space to estimate the atmosphere signal."""


def estimate_atmosphere_phase(stm: xr.Dataset, **unmodeled_displacement_kwargs) -> xr.Dataset:
    """Estimate the atmosphere phase.

    This function applies a temporal filter to extract the high-frequency
    atmospheric signal and then uses spatial kriging to estimate the atmospheric
    phase per epoch.

    Parameters
    ----------
    stm: xr.Dataset
        The network STMs with unwrapped phases.
    unmodeled_displacement_kwargs: dict
        Additional keyword arguments for the `estimate_unmodeled_displacement` function.

    Returns
    -------
    xr.Dataset
        The STM with the estimated atmospheric phase.
    """
    # Step 1: Apply temporal filtering to extract high-frequency atmospheric signal
    unmodeled_disp = estimate_unmodeled_displacement(
        psc_phase_residuals=stm["psc_phase_residuals"],
        baseline_years=stm["time"],
        **unmodeled_displacement_kwargs,
    )
    stm["atmosphere_estimates"] = stm["psc_phase_residuals"] - unmodeled_disp + stm["atmosphere_mother"]

    # Step 2: Apply spatial kriging to estimate atmospheric phase per epoch
