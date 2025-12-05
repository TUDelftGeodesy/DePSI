import numpy as np
import pytest
import xarray as xr
from numpy.testing import assert_almost_equal
from scipy import signal

from depsi.atmosphere_estimation import estimate_unmodeled_displacement


def _calculate_weigths(baseline_years, filter_length, sampling_rate, filter_type):
    baseline_scaled = baseline_years * sampling_rate
    timespan = baseline_scaled.max() - baseline_scaled.min() + 1

    # Create a window of size `timespan` to cover the full range of time differences
    window_size = int(timespan) + 1 # Ensure window size is an odd integer
    window = np.zeros(window_size)

    # Determine the core window size based on the filter length and sampling rate
    core_window_size = int(filter_length * sampling_rate) + 1
    start = (window_size - core_window_size) // 2
    end = start + core_window_size

    if filter_type == 'block':
        window[start:end] = signal.windows.boxcar(core_window_size)
    elif filter_type == 'triangle':
        window[start:end] = signal.windows.triang(core_window_size)
    elif filter_type == 'gaussian':
        std_dev = filter_length * sampling_rate / 6
        window = signal.windows.gaussian(window_size, std=std_dev)
    else:
        raise NotImplementedError(f"Method '{filter_type}' is not implemented.")
    window = window / window.sum()
    time_diffs = (baseline_scaled.values[:, None] - baseline_scaled.values[None, :])
    center_index = window_size // 2
    weight_indices = np.clip(center_index + np.round(time_diffs).astype(int), 0, window_size - 1)
    weight_matrix = window[weight_indices]
    weight_matrix /= np.sum(weight_matrix, axis=1, keepdims=True)
    return weight_matrix

class TestEstimateUnmodeledDisplacement:
    def test_estimate_unmodeled_displacement_block(self):
        """Test the block filter_type for estimating non-linear deformation."""
        psc_phase_residuals = xr.DataArray(np.random.rand(5, 10), dims=('space', 'time'))
        baseline_years = xr.DataArray(np.sort(np.random.rand(10)), dims='time')
        filter_length = 1
        sampling_rate = 2

        # Actual
        result = estimate_unmodeled_displacement(
            psc_phase_residuals=psc_phase_residuals,
            baseline_years=baseline_years,
            filter_length=filter_length,
            sampling_rate=sampling_rate,
            filter_type='block'
        )
        actual = result.isel(space=0).data

        # Expected
        weight_matrix = _calculate_weigths(baseline_years, filter_length, sampling_rate, 'block')
        expected = np.sum(weight_matrix * psc_phase_residuals.isel(space=0).values[None, :], axis=1)

        assert_almost_equal(actual, expected)

    def test_estimate_unmodeled_displacement_triangle(self):
        """Test the triangle filter_type for estimating non-linear deformation."""
        psc_phase_residuals = xr.DataArray(np.random.rand(5, 10), dims=('space', 'time'))
        baseline_years = xr.DataArray(np.sort(np.random.rand(10)), dims='time')
        filter_length = 1
        sampling_rate = 2

        # Actual
        result = estimate_unmodeled_displacement(
            psc_phase_residuals=psc_phase_residuals,
            baseline_years=baseline_years,
            filter_length=filter_length,
            sampling_rate=sampling_rate,
            filter_type='triangle'
        )
        actual = result.isel(space=0).data

        # Expected
        weight_matrix = _calculate_weigths(baseline_years, filter_length, sampling_rate, 'triangle')
        expected = np.sum(weight_matrix * psc_phase_residuals.isel(space=0).values[None, :], axis=1)

        assert_almost_equal(actual, expected)

    def test_estimate_unmodeled_displacement_gaussian(self):
        """Test the gaussian filter_type for estimating non-linear deformation."""
        psc_phase_residuals = xr.DataArray(np.random.rand(5, 10), dims=('space', 'time'))
        baseline_years = xr.DataArray(np.sort(np.random.rand(10)), dims='time')
        filter_length = 1
        sampling_rate = 2

        # Actual
        result = estimate_unmodeled_displacement(
            psc_phase_residuals=psc_phase_residuals,
            baseline_years=baseline_years,
            filter_length=filter_length,
            sampling_rate = sampling_rate,
            filter_type='gaussian'
        )
        actual = result.isel(space=0).data

        # Expected
        weight_matrix = _calculate_weigths(baseline_years, filter_length, sampling_rate, 'gaussian')
        expected = np.sum(weight_matrix * psc_phase_residuals.isel(space=0).values[None, :], axis=1)

        assert_almost_equal(actual, expected)

    def test_estimate_unmodeled_displacement_invalid_method(self):
        """Test that an error is raised for an invalid filter_type."""
        psc_phase_residuals = xr.DataArray(np.random.rand(10, 5), dims=('time', 'space'))
        baseline_years = xr.DataArray(np.sort(np.random.rand(10)), dims='time')

        with pytest.raises(NotImplementedError):
            estimate_unmodeled_displacement(
                psc_phase_residuals=psc_phase_residuals,
                baseline_years=baseline_years,
                filter_length=1,
                sampling_rate=2,
                filter_type='invalid_method'
            )

    def test_estimate_unmodeled_displacement_mismatched_sizes(self):
        """Test that an error is raised for mismatched sizes of psc_phase_residuals and baseline_years."""
        psc_phase_residuals = xr.DataArray(np.random.rand(10, 5), dims=('time', 'space'))
        baseline_years = xr.DataArray(np.arange(5), dims='time')

        with pytest.raises(ValueError):
            estimate_unmodeled_displacement(
                psc_phase_residuals=psc_phase_residuals,
                baseline_years=baseline_years,
                filter_length=1,
                sampling_rate=2,
                filter_type='block'
            )

    def test_estimate_unmodeled_displacement_non_monotonic_years(self):
        """Test that an error is raised for non-monotonic baseline_years."""
        psc_phase_residuals = xr.DataArray(np.random.rand(10, 5), dims=('time', 'space'))
        baseline_years = xr.DataArray(np.array([0, 2, 1, 3]), dims='time')

        with pytest.raises(ValueError):
            estimate_unmodeled_displacement(
                psc_phase_residuals=psc_phase_residuals,
                baseline_years=baseline_years,
                filter_length=1,
                sampling_rate=2,
                filter_type='block'
            )

    def test_estimate_unmodeled_displacement_large_sampling_rate(self):
        """Test that an error is raised for non-monotonic baseline_years."""
        psc_phase_residuals = xr.DataArray(np.random.rand(10, 5), dims=('time', 'space'))
        baseline_years = xr.DataArray(np.sort(np.random.rand(10)), dims='time')

        with pytest.raises(ValueError):
            estimate_unmodeled_displacement(
                psc_phase_residuals=psc_phase_residuals,
                baseline_years=baseline_years,
                filter_length=2,
                sampling_rate = 10,
                filter_type='block'
            )
