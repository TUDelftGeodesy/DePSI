import numpy as np
import pytest
import xarray as xr
from numpy.testing import assert_almost_equal
from scipy import signal
from scipy.ndimage import convolve1d

from depsi.atmosphere import estimate_non_linear_deformation


class TestEstimateNonLinearDeformation:
    def test_estimate_non_linear_deformation_block(self):
        """Test the block method for estimating non-linear deformation."""
        psc_phase = xr.DataArray(np.random.rand(5, 10), dims=('space', 'time'))
        baseline_years = xr.DataArray(np.sort(np.random.rand(10)), dims='time')
        filter_length = 2

        # Actual
        result = estimate_non_linear_deformation(
            psc_phase=psc_phase,
            baseline_years=baseline_years,
            filter_length=filter_length,
            method='block'
        )
        actual = result.isel(space=0).data

        # Expected
        std_dev = filter_length / 3
        window_size = int(6 * std_dev) | 1
        window = signal.windows.boxcar(window_size)
        window = window / window.sum()
        result = convolve1d(psc_phase, window, mode='mirror')
        expected = result[0, :]

        assert_almost_equal(actual, expected)


    def test_estimate_non_linear_deformation_triangle(self):
        """Test the triangle method for estimating non-linear deformation."""
        psc_phase = xr.DataArray(np.random.rand(5, 10), dims=('space', 'time'))
        baseline_years = xr.DataArray(np.sort(np.random.rand(10)), dims='time')
        filter_length = 2

        # Actual
        result = estimate_non_linear_deformation(
            psc_phase=psc_phase,
            baseline_years=baseline_years,
            filter_length=filter_length,
            method='triangle'
        )
        actual = result.isel(space=0).data

        # Expected
        std_dev = filter_length / 3
        window_size = int(6 * std_dev) | 1
        window = signal.windows.triang(window_size)
        window = window / window.sum()
        result = convolve1d(psc_phase, window, mode='mirror')
        expected = result[0, :]

        assert_almost_equal(actual, expected)

    def test_estimate_non_linear_deformation_gaussian(self):
        """Test the gaussian method for estimating non-linear deformation."""
        psc_phase = xr.DataArray(np.random.rand(5, 10), dims=('space', 'time'))
        baseline_years = xr.DataArray(np.sort(np.random.rand(10)), dims='time')
        filter_length = 2

        # Actual
        result = estimate_non_linear_deformation(
            psc_phase=psc_phase,
            baseline_years=baseline_years,
            filter_length=2,
            method='gaussian'
        )
        actual = result.isel(space=0).data

        # Expected
        std_dev = filter_length / 3
        window_size = int(6 * std_dev) | 1
        window = signal.windows.gaussian(window_size, std=std_dev)
        window = window / window.sum()
        result = convolve1d(psc_phase, window, mode='mirror')
        expected = result[0, :]

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
