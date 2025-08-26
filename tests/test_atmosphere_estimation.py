import numpy as np
import pykrige
import pytest
import xarray as xr
from numpy.testing import assert_almost_equal, assert_allclose
from scipy import signal
from scipy.optimize import curve_fit

from depsi.atmosphere_estimation import calculate_empirical_variogram, calculate_variogram_cloud, estimate_unmodeled_displacement, fit_variogram, setup_kriging_system


@pytest.fixture
def get_test_data():
    return xr.DataArray(
        np.array(
            [
                0.51615017,  0.49838693,  1.14898968,  0.73637701,  0.98923508,
                0.29524828, -0.03234519,  0.36147726,  1.42115834,  0.82486138
            ]
        ),
        dims=('space',),
        coords={
            'x': ('space', np.array([ 13.9,  69.5, 166.8, 180.7, 194.6, 194.6, 194.6, 208.5, 208.5, 222.4])),
            'y': ('space', np.array([ 8524, 13224, 13084, 13448,  6356,  6804,  8924,  2288, 11804, 5792])),
        }
    )



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


class TestCalculateVariogramCloud:
    def test_calculate_variogram_cloud(self):
        da = xr.DataArray(
            np.random.rand(5),
            dims=('space',),
            coords={
                'x': ('space', np.arange(5)),
                'y': ('space', np.arange(5)),
            }
        )
        actual_distances, actual_variances = calculate_variogram_cloud(da)

        x, y, z = da.x.values, da.y.values, da.values

        # full matrix
        distances = np.hypot(x - x[:, np.newaxis], y - y[:, np.newaxis])
        variances = (z - z[:, np.newaxis])**2

        # select upper triangle
        mask = np.triu(np.ones(distances.shape), k=1).astype(bool)
        pairwise_distances = distances[mask]
        pairwise_variances = variances[mask]

        # apply cutoff
        mask = pairwise_distances < 10000
        expected_distances = pairwise_distances[mask]
        expected_variances = pairwise_variances[mask]

        assert_almost_equal(actual_variances, expected_variances)
        assert_almost_equal(actual_distances, expected_distances)


class TestCalculateEmpiricalVariogram:
    def test_calculate_empirical_variogram_standard(self):
        da = xr.DataArray(
            np.random.rand(5),
            dims=('space',),
            coords={
                'x': ('space', np.arange(5)),
                'y': ('space', np.arange(5)),
            }
        )

        actual_lags, actual_variances = calculate_empirical_variogram(da, method='standard')

        distances, variances = calculate_variogram_cloud(da)
        nlags = 50
        bins = np.linspace(distances.min(), distances.max() + 1e-3, nlags + 1)
        lags, semivariances = [], []
        for left, right in zip(bins[:-1], bins[1:]):
            mask = (distances >= left) & (distances < right)
            if mask.any():
                lags.append(distances[mask].mean())
                semivariances.append(np.mean(variances[mask]))

        expected_lags, expected_variances = np.array(lags), np.array(semivariances)
        assert_almost_equal(actual_lags, expected_lags)
        assert_almost_equal(actual_variances, expected_variances)

    def test_calculate_empirical_variogram_unbiased(self):
        da = xr.DataArray(
            np.random.rand(5),
            dims=('space',),
            coords={
                'x': ('space', np.arange(5)),
                'y': ('space', np.arange(5)),
            }
        )

        actual_lags, actual_variances = calculate_empirical_variogram(da, method='unbiased')

        distances, variances = calculate_variogram_cloud(da)
        nlags = 50
        bins = np.linspace(distances.min(), distances.max() + 1e-3, nlags + 1)
        lags, semivariances = [], []
        for left, right in zip(bins[:-1], bins[1:]):
            mask = (distances >= left) & (distances < right)
            if mask.any():
                lags.append(distances[mask].mean())
                ch = 0.457 + 0.494 / len(variances[mask]) + 0.045 / len(variances[mask]) ** 2
                semivariances.append(1 / ch * np.mean(variances[mask] ** 0.25) ** 4)

        expected_lags, expected_variances = np.array(lags), np.array(semivariances)
        assert_almost_equal(actual_lags, expected_lags)
        assert_almost_equal(actual_variances, expected_variances)


    def test_calculate_empirical_variogram_unbiased_robust(self):
        da = xr.DataArray(
            np.random.rand(5),
            dims=('space',),
            coords={
                'x': ('space', np.arange(5)),
                'y': ('space', np.arange(5)),
            }
        )

        actual_lags, actual_variances = calculate_empirical_variogram(da, method='unbiased_robust')

        distances, variances = calculate_variogram_cloud(da)
        nlags = 50
        bins = np.linspace(distances.min(), distances.max() + 1e-3, nlags + 1)
        lags, semivariances = [], []
        for left, right in zip(bins[:-1], bins[1:]):
            mask = (distances >= left) & (distances < right)
            if mask.any():
                lags.append(distances[mask].mean())
                semivariances.append(1 / 0.457 * np.median(variances[mask] ** 0.25) ** 4)

        expected_lags, expected_variances = np.array(lags), np.array(semivariances)
        assert_almost_equal(actual_lags, expected_lags)
        assert_almost_equal(actual_variances, expected_variances)

class TestFitVariogram:
    def test_fit_variogram_gaussian(self, get_test_data):
        da = get_test_data

        _, lags_variances = fit_variogram(da, variogram_model='gaussian')
        lags, estimated_semivariances, semivariances = lags_variances
        actual_residual = semivariances - estimated_semivariances

        def gaussian_model(h, psill, range_, nugget):
            return psill * (1.0 - np.exp(-(h**2.0) / (range_ * 4.0 / 7.0) ** 2.0)) + nugget

        lags, semivariances = calculate_empirical_variogram(da)
        initial_guess = [0.4, 2000, 0.15] # psill, range, nugget
        popt, _ = curve_fit(gaussian_model, lags, semivariances, p0=initial_guess, bounds=(0, np.inf))
        estimated_semivariances = gaussian_model(lags, *popt)
        expected_residual = semivariances - estimated_semivariances

        assert_allclose(actual_residual, expected_residual, atol=1e-1)

    def test_fit_variogram_gaussian_params(self, get_test_data):
        da = get_test_data

        params, _ = fit_variogram(da, variogram_model='gaussian')
        assert len(params) == 3
        assert "sill" in params
        assert "range" in params
        assert "nugget" in params

    def test_fit_variogram_without_kwrags(self, get_test_data):
        da = get_test_data

        _, (lags, _, _) = fit_variogram(da, variogram_model='gaussian')

        # test default values
        assert lags.shape[0] < 50
        assert lags.max() < 10000

    def test_fit_variogram_with_kwrags(self, get_test_data):
        da = get_test_data

        kwrgs_empirical_variogram = {
            'method': 'unbiased_robust',
            'nlags': 30,
            'cutoff': 5000,
        }
        _, (lags, _, empirical_var) = fit_variogram(da, variogram_model='gaussian', **kwrgs_empirical_variogram)

        # test default values
        assert lags.shape[0] < 30
        assert lags.max() < 5000
        assert_allclose(empirical_var[0], 0.3839, atol=1e-3)


class TestSetupKrigingSystem:
    def test_setup_kriging_system(self, get_test_data):
        da = get_test_data
        krige_obj = setup_kriging_system(da)

        assert isinstance(krige_obj, pykrige.uk.UniversalKriging)
        assert krige_obj.variogram_model == 'gaussian'
        assert len(krige_obj.variogram_model_parameters) == 3
        assert krige_obj.variogram_model_parameters[2] > 0  # nugget
        assert krige_obj.regional_linear_drift == True

    def test_setup_kriging_system_other_method(self, get_test_data):
        da = get_test_data
        with pytest.raises(NotImplementedError):
            setup_kriging_system(da, method="ordinary")