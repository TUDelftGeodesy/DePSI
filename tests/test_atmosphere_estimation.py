import numpy as np
import pykrige
import pytest
import xarray as xr
from numpy.random import default_rng
from numpy.testing import assert_allclose, assert_almost_equal
from scipy import signal
from scipy.optimize import curve_fit

from depsi.atmosphere_estimation import (
    calculate_empirical_variogram,
    calculate_variogram_cloud,
    estimate_atmosphere_phase,
    estimate_unmodeled_displacement,
    fit_variogram,
    setup_kriging_system,
    solve_kriging,
    solve_kriging_per_single_time,
)

# Create a random number generator for consistent random data
rng = default_rng(seed=42)


@pytest.fixture
def get_test_data():
    return xr.DataArray(
        np.array(
            [
                0.51615017,
                0.49838693,
                1.14898968,
                0.73637701,
                0.98923508,
                0.29524828,
                -0.03234519,
                0.36147726,
                1.42115834,
                0.82486138,
            ]
        ),
        dims=("space",),
        coords={
            "x": ("space", np.array([13.9, 69.5, 166.8, 180.7, 194.6, 194.6, 194.6, 208.5, 208.5, 222.4])),
            "y": ("space", np.array([8524, 13224, 13084, 13448, 6356, 6804, 8924, 2288, 11804, 5792])),
        },
    )


def _calculate_weigths(baseline_years, filter_length, sampling_rate, filter_type):
    baseline_scaled = baseline_years * sampling_rate
    timespan = baseline_scaled.max() - baseline_scaled.min() + 1

    # Create a window of size `timespan` to cover the full range of time differences
    window_size = int(timespan) + 1  # Ensure window size is an odd integer
    window = np.zeros(window_size)

    # Determine the core window size based on the filter length and sampling rate
    core_window_size = int(filter_length * sampling_rate) + 1
    start = (window_size - core_window_size) // 2
    end = start + core_window_size

    if filter_type == "block":
        window[start:end] = signal.windows.boxcar(core_window_size)
    elif filter_type == "triangle":
        window[start:end] = signal.windows.triang(core_window_size)
    elif filter_type == "gaussian":
        std_dev = filter_length * sampling_rate / 6
        window = signal.windows.gaussian(window_size, std=std_dev)
    else:
        raise NotImplementedError(f"Method '{filter_type}' is not implemented.")
    window = window / window.sum()
    time_diffs = baseline_scaled.values[:, None] - baseline_scaled.values[None, :]
    center_index = window_size // 2
    weight_indices = np.clip(center_index + np.round(time_diffs).astype(int), 0, window_size - 1)
    weight_matrix = window[weight_indices]
    weight_matrix /= np.sum(weight_matrix, axis=1, keepdims=True)
    return weight_matrix


class TestEstimateUnmodeledDisplacement:
    def test_estimate_unmodeled_displacement_block(self):
        """Test the block filter_type for estimating non-linear deformation."""
        psc_phase_residuals = xr.DataArray(rng.random((5, 10)), dims=("space", "time"))
        baseline_years = xr.DataArray(np.sort(rng.random(10)), dims="time")
        filter_length = 1
        sampling_rate = 2

        # Actual
        result = estimate_unmodeled_displacement(
            psc_phase_residuals=psc_phase_residuals,
            baseline_years=baseline_years,
            filter_length=filter_length,
            sampling_rate=sampling_rate,
            filter_type="block",
        )
        actual = result.isel(space=0).data

        # Expected
        weight_matrix = _calculate_weigths(baseline_years, filter_length, sampling_rate, "block")
        expected = np.sum(weight_matrix * psc_phase_residuals.isel(space=0).values[None, :], axis=1)

        assert_almost_equal(actual, expected)

    def test_estimate_unmodeled_displacement_triangle(self):
        """Test the triangle filter_type for estimating non-linear deformation."""
        psc_phase_residuals = xr.DataArray(rng.random((5, 10)), dims=("space", "time"))
        baseline_years = xr.DataArray(np.sort(rng.random(10)), dims="time")
        filter_length = 1
        sampling_rate = 2

        # Actual
        result = estimate_unmodeled_displacement(
            psc_phase_residuals=psc_phase_residuals,
            baseline_years=baseline_years,
            filter_length=filter_length,
            sampling_rate=sampling_rate,
            filter_type="triangle",
        )
        actual = result.isel(space=0).data

        # Expected
        weight_matrix = _calculate_weigths(baseline_years, filter_length, sampling_rate, "triangle")
        expected = np.sum(weight_matrix * psc_phase_residuals.isel(space=0).values[None, :], axis=1)

        assert_almost_equal(actual, expected)

    def test_estimate_unmodeled_displacement_gaussian(self):
        """Test the gaussian filter_type for estimating non-linear deformation."""
        psc_phase_residuals = xr.DataArray(rng.random((5, 10)), dims=("space", "time"))
        baseline_years = xr.DataArray(np.sort(rng.random(10)), dims="time")
        filter_length = 1
        sampling_rate = 2

        # Actual
        result = estimate_unmodeled_displacement(
            psc_phase_residuals=psc_phase_residuals,
            baseline_years=baseline_years,
            filter_length=filter_length,
            sampling_rate=sampling_rate,
            filter_type="gaussian",
        )
        actual = result.isel(space=0).data

        # Expected
        weight_matrix = _calculate_weigths(baseline_years, filter_length, sampling_rate, "gaussian")
        expected = np.sum(weight_matrix * psc_phase_residuals.isel(space=0).values[None, :], axis=1)

        assert_almost_equal(actual, expected)

    def test_estimate_unmodeled_displacement_invalid_method(self):
        """Test that an error is raised for an invalid filter_type."""
        psc_phase_residuals = xr.DataArray(rng.random((10, 5)), dims=("time", "space"))
        baseline_years = xr.DataArray(np.sort(rng.random(10)), dims="time")

        with pytest.raises(NotImplementedError):
            estimate_unmodeled_displacement(
                psc_phase_residuals=psc_phase_residuals,
                baseline_years=baseline_years,
                filter_length=1,
                sampling_rate=2,
                filter_type="invalid_method",
            )

    def test_estimate_unmodeled_displacement_mismatched_sizes(self):
        """Test that an error is raised for mismatched sizes of psc_phase_residuals and baseline_years."""
        psc_phase_residuals = xr.DataArray(rng.random((10, 5)), dims=("time", "space"))
        baseline_years = xr.DataArray(np.arange(5), dims="time")

        with pytest.raises(ValueError):
            estimate_unmodeled_displacement(
                psc_phase_residuals=psc_phase_residuals,
                baseline_years=baseline_years,
                filter_length=1,
                sampling_rate=2,
                filter_type="block",
            )

    def test_estimate_unmodeled_displacement_non_monotonic_years(self):
        """Test that an error is raised for non-monotonic baseline_years."""
        psc_phase_residuals = xr.DataArray(rng.random((10, 5)), dims=("time", "space"))
        baseline_years = xr.DataArray(np.array([0, 2, 1, 3]), dims="time")

        with pytest.raises(ValueError):
            estimate_unmodeled_displacement(
                psc_phase_residuals=psc_phase_residuals,
                baseline_years=baseline_years,
                filter_length=1,
                sampling_rate=2,
                filter_type="block",
            )

    def test_estimate_unmodeled_displacement_large_sampling_rate(self):
        """Test that an error is raised for non-monotonic baseline_years."""
        psc_phase_residuals = xr.DataArray(rng.random((10, 5)), dims=("time", "space"))
        baseline_years = xr.DataArray(np.sort(rng.random(10)), dims="time")

        with pytest.raises(ValueError):
            estimate_unmodeled_displacement(
                psc_phase_residuals=psc_phase_residuals,
                baseline_years=baseline_years,
                filter_length=2,
                sampling_rate=10,
                filter_type="block",
            )


class TestCalculateVariogramCloud:
    def test_calculate_variogram_cloud(self, get_test_data):
        da = get_test_data
        actual_distances, actual_variances = calculate_variogram_cloud(da)

        x, y, z = da.x.values, da.y.values, da.values

        # full matrix
        distances = np.hypot(x - x[:, np.newaxis], y - y[:, np.newaxis])
        variances = (z - z[:, np.newaxis]) ** 2

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
            rng.random(5),
            dims=("space",),
            coords={
                "x": ("space", np.arange(5)),
                "y": ("space", np.arange(5)),
            },
        )

        actual_lags, actual_variances = calculate_empirical_variogram(da, method="standard")

        distances, variances = calculate_variogram_cloud(da)
        nlags = 50
        bins = np.linspace(distances.min(), distances.max() + 1e-3, nlags + 1)
        lags, semivariances = [], []
        for left, right in zip(bins[:-1], bins[1:], strict=False):
            mask = (distances >= left) & (distances < right)
            if mask.any():
                lags.append(distances[mask].mean())
                semivariances.append(np.mean(variances[mask]))

        expected_lags, expected_variances = np.array(lags), np.array(semivariances)
        assert_almost_equal(actual_lags, expected_lags)
        assert_almost_equal(actual_variances, expected_variances)

    def test_calculate_empirical_variogram_unbiased(self):
        da = xr.DataArray(
            rng.random(5),
            dims=("space",),
            coords={
                "x": ("space", np.arange(5)),
                "y": ("space", np.arange(5)),
            },
        )

        actual_lags, actual_variances = calculate_empirical_variogram(da, method="unbiased")

        distances, variances = calculate_variogram_cloud(da)
        nlags = 50
        bins = np.linspace(distances.min(), distances.max() + 1e-3, nlags + 1)
        lags, semivariances = [], []
        for left, right in zip(bins[:-1], bins[1:], strict=False):
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
            rng.random(5),
            dims=("space",),
            coords={
                "x": ("space", np.arange(5)),
                "y": ("space", np.arange(5)),
            },
        )

        actual_lags, actual_variances = calculate_empirical_variogram(da, method="unbiased_robust")

        distances, variances = calculate_variogram_cloud(da)
        nlags = 50
        bins = np.linspace(distances.min(), distances.max() + 1e-3, nlags + 1)
        lags, semivariances = [], []
        for left, right in zip(bins[:-1], bins[1:], strict=False):
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

        _, lags_variances = fit_variogram(da, variogram_model="gaussian")
        lags, estimated_semivariances, semivariances = lags_variances
        actual_residual = semivariances - estimated_semivariances

        def gaussian_model(h, psill, range_, nugget):
            return psill * (1.0 - np.exp(-(h**2.0) / (range_ * 4.0 / 7.0) ** 2.0)) + nugget

        lags, semivariances = calculate_empirical_variogram(da)
        initial_guess = [0.4, 2000, 0.15]  # psill, range, nugget
        popt, _ = curve_fit(gaussian_model, lags, semivariances, p0=initial_guess, bounds=(0, np.inf))
        estimated_semivariances = gaussian_model(lags, *popt)
        expected_residual = semivariances - estimated_semivariances

        assert_allclose(actual_residual, expected_residual, atol=1e-1)

    def test_fit_variogram_gaussian_params(self, get_test_data):
        da = get_test_data

        params, _ = fit_variogram(da, variogram_model="gaussian")
        assert len(params) == 3
        assert "sill" in params
        assert "range" in params
        assert "nugget" in params

    def test_fit_variogram_without_kwrags(self, get_test_data):
        da = get_test_data

        _, (lags, _, _) = fit_variogram(da, variogram_model="gaussian")

        # test default values
        assert lags.shape[0] < 50
        assert lags.max() < 10000

    def test_fit_variogram_with_kwrags(self, get_test_data):
        da = get_test_data

        _, (lags, _, empirical_var) = fit_variogram(
            da,
            variogram_model="gaussian",
            empirical_variogram_method="unbiased_robust",
            empirical_variogram_nlags=30,
            empirical_variogram_cutoff=5000
        )

        # test default values
        assert lags.shape[0] < 30
        assert lags.max() < 5000
        assert_allclose(empirical_var[0], 0.3839, atol=1e-3)

class TestSetupKrigingSystem:
    def test_setup_kriging_system_default(self, get_test_data):
        da = get_test_data
        krige_obj = setup_kriging_system(da)

        assert isinstance(krige_obj, pykrige.uk.UniversalKriging)
        assert krige_obj.variogram_model == "gaussian"
        assert len(krige_obj.variogram_model_parameters) == 3
        assert krige_obj.variogram_model_parameters[2] > 0  # nugget
        assert krige_obj.regional_linear_drift

    def test_setup_kriging_system_other_method(self, get_test_data):
        da = get_test_data
        with pytest.raises(NotImplementedError):
            setup_kriging_system(da, method="ordinary")

    def test_setup_kriging_system_kwargs(self, get_test_data):
        da = get_test_data
        variogram_args = {
            "variogram_model": "power",
            "variogram_parameters": {"scale": 0.5, "exponent": 1.5, "nugget": 0.1},
            "drift_terms": None,
        }
        krige_obj = setup_kriging_system(da, variogram_args=variogram_args)

        assert isinstance(krige_obj, pykrige.uk.UniversalKriging)
        assert krige_obj.variogram_model == "power"
        assert len(krige_obj.variogram_model_parameters) == 3
        assert_almost_equal(
            krige_obj.variogram_model_parameters, list(variogram_args["variogram_parameters"].values())
        )
        assert not krige_obj.regional_linear_drift


class TestSolveKrigingPerSingleTime:
    def test_solve_kriging_per_single_time_grid(self, get_test_data):
        da = get_test_data
        x_min, x_max = da.x.min(), da.x.max()
        y_min, y_max = da.y.min(), da.y.max()

        grid_resolution = 500  # in meter
        x_grid = np.arange(x_min, x_max + grid_resolution, grid_resolution)
        y_grid = np.arange(y_min, y_max + grid_resolution, grid_resolution)
        grid = xr.Dataset(coords={"x": x_grid, "y": y_grid})

        zvalues, sigmasq = solve_kriging_per_single_time(da, grid)

        assert zvalues.shape[0] == len(grid.y) and zvalues.shape[1] == len(grid.x)
        assert not np.any(np.isnan(zvalues))
        assert not np.any(np.isnan(sigmasq))

    def test_solve_kriging_per_single_time_points(self, get_test_data):
        da = get_test_data
        points = xr.Dataset(coords=da.coords)

        zvalues, _ = solve_kriging_per_single_time(da, points)

        assert zvalues.shape[0] == len(points.space)
        assert not np.any(np.isnan(zvalues))

    def test_solve_kriging_per_single_time_n_neighbours(self, get_test_data):
        da = get_test_data
        points = xr.Dataset(coords=da.coords)

        zvalues, _ = solve_kriging_per_single_time(da, points, n_nearest_neighbors=5)

        assert zvalues.shape[0] == len(points.space)
        assert not np.any(np.isnan(zvalues))

    def test_solve_kriging_per_single_time_n_neighbours_grid(self, get_test_data):
        da = get_test_data
        x_min, x_max = da.x.min(), da.x.max()
        y_min, y_max = da.y.min(), da.y.max()

        grid_resolution = 500  # in meter
        x_grid = np.arange(x_min, x_max + grid_resolution, grid_resolution)
        y_grid = np.arange(y_min, y_max + grid_resolution, grid_resolution)
        grid = xr.Dataset(coords={"x": x_grid, "y": y_grid})

        with pytest.raises(NotImplementedError):
            solve_kriging_per_single_time(da, grid, n_nearest_neighbors=5)


class TestSolveKriging:
    def test_solve_kriging_no_time(self, get_test_data):
        da = get_test_data

        with pytest.raises(ValueError) as excinfo:
            solve_kriging(da, da)
        assert "ps_atmosphere must have a 'time' dimension." in str(excinfo.value)

    def test_solve_kriging_chunked_space(self, get_test_data):
        da = get_test_data
        da.chunk({"space": 5})

        with pytest.raises(ValueError) as excinfo:
            solve_kriging(da, da)
            assert "ps_atmosphere must not be chunked" in str(excinfo.value)

    def test_solve_kriging_with_time(self, get_test_data):
        da = get_test_data
        points = xr.Dataset(coords=da.coords)

        da = xr.concat([da] * 3, dim="time")
        da = da.assign_coords(time=np.array([1, 2, 3]))

        results = solve_kriging(da, points)
        assert isinstance(results, xr.Dataset)
        assert "predicted" in results and "sigmasq" in results
        assert len(results.time) == len(da.time)

    def test_solve_kriging_with_time_in_grid(self, get_test_data):
        da = get_test_data
        da = xr.concat([da] * 3, dim="time")
        da = da.assign_coords(time=np.array([1, 2, 3]))

        points = xr.Dataset(coords=da.coords)
        with pytest.raises(ValueError) as excinfo:
            solve_kriging(da, points)

        assert "Grid must not have 'time' dimension" in str(excinfo.value)


class TestEstimateAtmospherePhase:
    def test_estimate_atmosphere_phase_defaults(self, get_test_data):
        da = get_test_data
        da = xr.concat([da] * 20, dim="time")
        da = da.assign_coords(time=np.sort(rng.random(20) * 10))
        stm = da.to_dataset(name="psc_phase_residuals")
        stm["atmosphere_mother"] = (("time", "space"), rng.random((len(da.time), len(da.space))))

        results = estimate_atmosphere_phase(stm)
        assert isinstance(results, xr.Dataset)
        assert "psc_phase_residuals" in results
        assert "atmosphere_mother" in results
        assert "unmodeled_disp" in results
        assert "atmosphere_estimates" in results
        assert "atmosphere_predicted" in results
        assert "atmosphere_sigmasq" in results

    def test_estimate_atmosphere_phase_kwargs(self, get_test_data):
        da = get_test_data
        da = xr.concat([da] * 20, dim="time")
        da = da.assign_coords(time=np.sort(rng.random(20) * 10))
        stm = da.to_dataset(name="psc_phase_residuals")
        stm["atmosphere_mother"] = (("time", "space"), rng.random((len(da.time), len(da.space))))

        results = estimate_atmosphere_phase(
            stm,
            unmodeled_displacement_args={
            "filter_length": 2,
            "sampling_rate": 1,
            "filter_type": "triangle",
        },
            kriging_args={
                "variogram_args": {
            "variogram_model": "power",
            "variogram_parameters": {"scale": 0.5, "exponent": 1.5, "nugget": 0.1},
            "drift_terms": None,
            }
        }
        )

        assert isinstance(results, xr.Dataset)
        assert "psc_phase_residuals" in results
        assert "atmosphere_mother" in results
        assert "unmodeled_disp" in results
        assert "atmosphere_estimates" in results
        assert "atmosphere_predicted" in results
        assert "atmosphere_sigmasq" in results
