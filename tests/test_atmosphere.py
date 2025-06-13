import numpy as np
import pytest
import xarray as xr
from numpy.testing import assert_array_equal
from scipy import signal

from depsi.atmosphere import estimate_non_linear_deformation


class TestEstimateNonLinearDeformation:
    def test_estimate_non_linear_deformation_block(self):
        """Test the block method for estimating non-linear deformation."""
        psc_phase = xr.DataArray(np.random.rand(10, 5), dims=('time', 'space'))
        baseline_years = xr.DataArray(np.arange(10), dims='time')

        # Actual
        result = estimate_non_linear_deformation(
            psc_phase=psc_phase,
            baseline_years=baseline_years,
            filter_length=2,
            method='block'
        )
        actual = result.isel(space=0).data

        # Expected
        low_pass_filter = np.zeros_like(baseline_years, dtype=float)
        low_pass_filter[:2] = signal.windows.boxcar(2)
        distances_matrix = np.abs(
            baseline_years.values[:, None] - baseline_years.values[None, :]
        )
        weights_matrix = low_pass_filter[distances_matrix.astype(int)]
        weights_matrix = weights_matrix / weights_matrix.sum(axis=1, keepdims=True)
        expected = np.dot(weights_matrix, psc_phase.isel(space=0).values)

        assert_array_equal(actual, expected)


    def test_estimate_non_linear_deformation_triangle(self):
        """Test the triangle method for estimating non-linear deformation."""
        psc_phase = xr.DataArray(np.random.rand(10, 5), dims=('time', 'space'))
        baseline_years = xr.DataArray(np.arange(10), dims='time')

        # Actual
        result = estimate_non_linear_deformation(
            psc_phase=psc_phase,
            baseline_years=baseline_years,
            filter_length=2,
            method='triangle'
        )
        actual = result.isel(space=0).data

        # Expected
        low_pass_filter = np.zeros_like(baseline_years, dtype=float)
        low_pass_filter[:2] = signal.windows.triang(2)
        distances_matrix = np.abs(
            baseline_years.values[:, None] - baseline_years.values[None, :]
        )
        weights_matrix = low_pass_filter[distances_matrix.astype(int)]
        weights_matrix = weights_matrix / weights_matrix.sum(axis=1, keepdims=True)
        expected = np.dot(weights_matrix, psc_phase.isel(space=0).values)

        assert_array_equal(actual, expected)

    def test_estimate_non_linear_deformation_gaussian(self):
        """Test the gaussian method for estimating non-linear deformation."""
        psc_phase = xr.DataArray(np.random.rand(10, 5), dims=('time', 'space'))
        baseline_years = xr.DataArray(np.arange(10), dims='time')

        # Actual
        result = estimate_non_linear_deformation(
            psc_phase=psc_phase,
            baseline_years=baseline_years,
            filter_length=2,
            method='gaussian'
        )
        actual = result.isel(space=0).data

        # Expected
        low_pass_filter = signal.windows.gaussian(baseline_years.size, std=1 / 3)
        distances_matrix = np.abs(
            baseline_years.values[:, None] - baseline_years.values[None, :]
        )
        weights_matrix = low_pass_filter[distances_matrix.astype(int)]
        weights_matrix = weights_matrix / weights_matrix.sum(axis=1, keepdims=True)
        expected = np.dot(weights_matrix, psc_phase.isel(space=0).values)

        assert_array_equal(actual, expected)

    def test_estimate_non_linear_deformation_invalid_method(self):
        """Test that an error is raised for an invalid method."""
        psc_phase = xr.DataArray(np.random.rand(10, 5), dims=('time', 'space'))
        baseline_years = xr.DataArray(np.arange(10), dims='time')

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
