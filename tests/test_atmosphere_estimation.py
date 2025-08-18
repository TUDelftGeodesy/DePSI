import numpy as np
import pytest
import xarray as xr
from numpy.testing import assert_almost_equal
from scipy import signal

from depsi.atmosphere import estimate_non_linear_deformation


def _calculate_weigths(baseline_years, filter_length, temporal_scale, method):
    baseline_scaled = baseline_years * temporal_scale
    timespan = baseline_scaled.max() - baseline_scaled.min()
    window_size = int(2 * timespan) + 1
    if method == 'block':
        window = signal.windows.boxcar(window_size)
    elif method == 'triangle':
        window = signal.windows.triang(window_size)
    elif method == 'gaussian':
        std_dev = filter_length * temporal_scale / 6
        window = signal.windows.gaussian(window_size, std=std_dev)
    else:
        raise NotImplementedError(f"Method '{method}' is not implemented.")
    window = window / window.sum()
    time_diffs = (baseline_scaled.values[:, None] - baseline_scaled.values[None, :])
    center_index = window_size // 2
    weight_indices = np.clip(center_index + np.round(time_diffs).astype(int), 0, window_size - 1)
    weight_matrix = window[weight_indices]
    weight_matrix /= np.sum(weight_matrix, axis=1, keepdims=True)
    return weight_matrix

class TestEstimateNonLinearDeformation:
    def test_estimate_non_linear_deformation_block(self):
        """Test the block method for estimating non-linear deformation."""
        psc_phase = xr.DataArray(np.random.rand(5, 10), dims=('space', 'time'))
        baseline_years = xr.DataArray(np.sort(np.random.rand(10)), dims='time')
        filter_length = 1
        temporal_scale = 10

        # Actual
        result = estimate_non_linear_deformation(
            psc_phase=psc_phase,
            baseline_years=baseline_years,
            filter_length=filter_length,
            temporal_scale=temporal_scale,
            method='block'
        )
        actual = result.isel(space=0).data

        # Expected
        weight_matrix = _calculate_weigths(baseline_years, filter_length, temporal_scale, 'block')
        expected = np.sum(weight_matrix * psc_phase.isel(space=0).values[None, :], axis=1)

        assert_almost_equal(actual, expected)

    def test_estimate_non_linear_deformation_triangle(self):
        """Test the triangle method for estimating non-linear deformation."""
        psc_phase = xr.DataArray(np.random.rand(5, 10), dims=('space', 'time'))
        baseline_years = xr.DataArray(np.sort(np.random.rand(10)), dims='time')
        filter_length = 1
        temporal_scale = 10

        # Actual
        result = estimate_non_linear_deformation(
            psc_phase=psc_phase,
            baseline_years=baseline_years,
            filter_length=filter_length,
            temporal_scale=temporal_scale,
            method='triangle'
        )
        actual = result.isel(space=0).data

        # Expected
        weight_matrix = _calculate_weigths(baseline_years, filter_length, temporal_scale, 'triangle')
        expected = np.sum(weight_matrix * psc_phase.isel(space=0).values[None, :], axis=1)

        assert_almost_equal(actual, expected)

    def test_estimate_non_linear_deformation_gaussian(self):
        """Test the gaussian method for estimating non-linear deformation."""
        psc_phase = xr.DataArray(np.random.rand(5, 10), dims=('space', 'time'))
        baseline_years = xr.DataArray(np.sort(np.random.rand(10)), dims='time')
        filter_length = 2
        temporal_scale = 10

        # Actual
        result = estimate_non_linear_deformation(
            psc_phase=psc_phase,
            baseline_years=baseline_years,
            filter_length=filter_length,
            temporal_scale = temporal_scale,
            method='gaussian'
        )
        actual = result.isel(space=0).data

        # Expected
        weight_matrix = _calculate_weigths(baseline_years, filter_length, temporal_scale, 'gaussian')
        expected = np.sum(weight_matrix * psc_phase.isel(space=0).values[None, :], axis=1)

        assert_almost_equal(actual, expected)

    def test_estimate_non_linear_deformation_invalid_method(self):
        """Test that an error is raised for an invalid method."""
        psc_phase = xr.DataArray(np.random.rand(10, 5), dims=('time', 'space'))
        baseline_years = xr.DataArray(np.sort(np.random.rand(10)), dims='time')

        with pytest.raises(NotImplementedError):
            estimate_non_linear_deformation(
                psc_phase=psc_phase,
                baseline_years=baseline_years,
                filter_length=2,
                method='invalid_method'
            )

    def test_estimate_non_linear_deformation_mismatched_sizes(self):
        """Test that an error is raised for mismatched sizes of psc_phase and baseline_years."""
        psc_phase = xr.DataArray(np.random.rand(10, 5), dims=('time', 'space'))
        baseline_years = xr.DataArray(np.arange(5), dims='time')

        with pytest.raises(ValueError):
            estimate_non_linear_deformation(
                psc_phase=psc_phase,
                baseline_years=baseline_years,
                filter_length=2,
                method='block'
            )

    def test_estimate_non_linear_deformation_non_monotonic_years(self):
        """Test that an error is raised for non-monotonic baseline_years."""
        psc_phase = xr.DataArray(np.random.rand(10, 5), dims=('time', 'space'))
        baseline_years = xr.DataArray(np.array([0, 2, 1, 3]), dims='time')

        with pytest.raises(ValueError):
            estimate_non_linear_deformation(
                psc_phase=psc_phase,
                baseline_years=baseline_years,
                filter_length=2,
                method='block'
            )
