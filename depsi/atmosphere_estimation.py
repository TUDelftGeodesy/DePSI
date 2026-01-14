"""A module to estimate atmosphere signal from network STMs with unwrapped phases.

This estimation is done in two steps:
1. A temporal filtering applied per point to extract the unmodeled deformation (low-frequency signal).
signal.
2. A spatial least-squares prediction based on the residuals per epoch to estimate the atmosphere signal per
epoch.
Using: https://geostat-framework.readthedocs.io/projects/pykrige/en/stable/generated/pykrige.uk.UniversalKriging.html#pykrige.uk.UniversalKriging
"""

from logging import getLogger

import numpy as np
import pykrige
import xarray as xr
from scipy import signal
from scipy.spatial import KDTree
from scipy.spatial.distance import pdist

logger = getLogger(__name__)


def _get_signal_window_with_zero_padding(type, timespan, filter_length, sampling_rate) -> np.ndarray:
    """Build a signal window with zero padding.

    Parameters
    ----------
    type: str
        Type of window to build, e.g. 'boxcar', 'triang', see `scipy.signal.windows` for more.
    timespan: float
        The total timespan to cover with the window.
    filter_length: int
        Length of the filter (year) to apply a low-pass filter to the time series.
    sampling_rate: int
        Sampling rate of the time series.

    Returns
    -------
    np.ndarray
        The signal window with zero padding.
    """
    # Determine window of size to cover the full range of time differences
    window_size = int(timespan) | 1  # Ensure window size is an odd integer

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
    filter_length: int = 1,
    sampling_rate: int = 1000,
    filter_type="gaussian",
) -> xr.DataArray:
    """Apply a low-pass filter to the time series to estimate the unmodeled deformation.

    Parameters
    ----------
    psc_phase_residuals: xr.DataArray
        The PSC time series residuals to apply the filter to.
    baseline_years: xr.DataArray
        The baseline years corresponding to the time series.
    filter_length: int
        Length of the filter (year) to apply a low-pass filter to the time series. Default is 1.
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
        raise ValueError("The size of baseline_years must match the time dimension of psc_phase_residuals.")
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

    window_size = int(timespan) | 1  # Ensure window size is an odd integer

    if filter_type == "block":
        window = _get_signal_window_with_zero_padding(
            type="boxcar", timespan=timespan, filter_length=filter_length, sampling_rate=sampling_rate
        )

    elif filter_type == "triangle":
        window = _get_signal_window_with_zero_padding(
            type="triang", timespan=timespan, filter_length=filter_length, sampling_rate=sampling_rate
        )
    elif filter_type == "gaussian":
        std_dev = filter_length * sampling_rate / 6  # ±3σ covers the window
        window = signal.windows.gaussian(window_size, std=std_dev)
    else:
        raise NotImplementedError(
            f"Filter type {filter_type} is not implemented. Available types are: 'block', 'triangle', 'gaussian'."
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
    def _apply_filter(data):
        return np.einsum("ij,j->i", weight_matrix, data)

    unmodeled_disp = xr.apply_ufunc(
        _apply_filter,
        psc_phase_residuals,
        input_core_dims=[["time"]],
        output_core_dims=[["time"]],
        vectorize=True,
        dask="parallelized",
        output_dtypes=[psc_phase_residuals.dtype],
    )
    return unmodeled_disp.transpose(*psc_phase_residuals.dims)  # align dims order


def calculate_variogram_cloud(da: xr.DataArray, cutoff: float = 10000.0):
    """Calculate the variogram cloud for a DataArray.

    Parameters
    ----------
    da: xr.DataArray
        The DataArray containing the variable of interest (including 'x' and 'y'
        coordinates) to calculate the variogram cloud.
    cutoff: float
        The maximum distance to consider for the variogram cloud. Default is
        10000.0.

    Returns
    -------
    pairwise_distances: np.ndarray
        The pairwise distances between the points in the DataArray.
    variances: np.ndarray
        The variances corresponding to the pairwise distances.
    """
    # Check if there x, y coords
    if "x" not in da.coords or "y" not in da.coords:
        raise ValueError("DataArray must have coordinates 'x' and 'y'.")

    pairwise_distances = pdist(np.column_stack((da.coords["x"].values, da.coords["y"].values)), metric="euclidean")

    z_values = da.values.flatten()[:, None]
    variances = pdist(z_values, metric="sqeuclidean")

    # apply cutoff
    mask = pairwise_distances < cutoff
    return pairwise_distances[mask], variances[mask]


def _calculate_binned_variances(variances, method: str = "standard"):
    """Calculate the binned variances based on the method specified.

    standard: Returns Experimental variogram.
    unbiased: Returns unbiased variogram mentioned in (Cressie-Hawkins, 1980).
    unbiased_robust: Returns unbiased robust variogram mentioned in (Cressie, 1993).
    """
    if method == "standard":
        return np.mean(variances)
    elif method == "unbiased":
        ch = 0.457 + 0.494 / len(variances) + 0.045 / len(variances) ** 2
        return 1 / ch * np.mean(variances**0.25) ** 4
    elif method == "unbiased_robust":
        return 1 / 0.457 * np.median(variances**0.25) ** 4


def calculate_empirical_variogram(da: xr.DataArray, method: str = "standard", nlags: int = 50, cutoff=10000.0):
    """Calculate the empirical variogram of a DataArray.

    Parameters
    ----------
    da: xr.DataArray
        The DataArray containing the data to calculate the empirical variogram.
    method: str
        The method to use for calculating the empirical variogram. Options are:
        'standard', 'unbiased', 'unbiased_robust'. Default is 'standard'.
    nlags: int
        The number of lags to use for the empirical variogram. Default is 50.
    cutoff: float
        The maximum distance to consider for the empirical variogram. Default is 10000.0.

    Returns
    -------
    lags: np.ndarray
        The lags for the empirical variogram.
    semivariance: np.ndarray
        The semivariance for the empirical variogram.
    """
    distances, variances = calculate_variogram_cloud(da, cutoff=cutoff)

    # Bin edges (ensure last bin includes max distance)
    bins = np.linspace(distances.min(), distances.max() + 1e-3, nlags + 1)

    lags, semivariances = [], []

    for left, right in zip(bins[:-1], bins[1:], strict=False):
        mask = (distances >= left) & (distances < right)
        if mask.any():
            lags.append(distances[mask].mean())
            semivariances.append(_calculate_binned_variances(variances[mask], method=method))

    return np.array(lags), np.array(semivariances)


def _check_variogram_args(kwargs):
    valid_keys = {
        "variogram_model",
        "variogram_parameters",
        "drift_terms",
    }
    for key in kwargs.keys():
        if key not in valid_keys:
            raise ValueError(f"Invalid keyword argument: {key}")


def _check_empirical_variogram_args(kwargs):
    valid_keys = {
        "method",
        "nlags",
        "cutoff",
    }
    for key in kwargs.keys():
        if key not in valid_keys:
            raise ValueError(f"Invalid keyword argument: {key}")


def fit_variogram(
    da: xr.DataArray,
    lags: np.ndarray = None,
    semivariances: np.ndarray = None,
    variogram_model: str = "gaussian",
    empirical_variogram_method: str = "standard",
    empirical_variogram_nlags: int = 50,
    empirical_variogram_cutoff: float = 10000.0,
):
    """Fit a variogram model to the empirical variogram.

    If the arguments `lags` or `semivariances` are None, the empirical variogram
    can be calculated using arguments `empirical_variogram_method`,
    `empirical_variogram_nlags`, and `empirical_variogram_cutoff`.

    Parameters
    ----------
    da: xr.DataArray
        The DataArray containing the variable of interest (including 'x' and 'y'
        coordinates) to fit the variogram model.
    lags: np.ndarray, optional
        The lags for the empirical variogram. If None, the empirical variogram
        will be calculated.
    semivariances: np.ndarray, optional
        The semivariances for the empirical variogram. If None, the empirical
        variogram will be calculated.
    variogram_model: str
        The variogram model to fit. Options are: 'linear', 'power', 'gaussian',
        'spherical', 'exponential', 'hole-effect'. Default is 'gaussian'.
    empirical_variogram_method: str
        The `method` argument can be one of 'standard', 'unbiased', 'unbiased_robust'.
        Default is 'standard'.
    empirical_variogram_nlags: int
        The number of lags to use for the empirical variogram. Default is 50.
    empirical_variogram_cutoff: float
        The maximum distance to consider for the empirical variogram. Default is 10000.0.

    Returns
    -------
    variogram_parameters: dict
        The fitted variogram parameters.
    (lags, estimated_semivariances, semivariances): tuple
        A tuple containing the lags, the estimated semivariances from the fitted
        model, and the empirical semivariances.
    """
    if lags is None or semivariances is None:
        lags, semivariances = calculate_empirical_variogram(
            da, empirical_variogram_method, empirical_variogram_nlags, empirical_variogram_cutoff
        )

    # see equations and reference in
    # https://geostat-framework.readthedocs.io/projects/pykrige/en/stable/variogram_models.html
    # Gaussian model uses "effective range" introduced in  Pebesma, E.J. &
    # Wesseling, C.G. (1998). "Gstat: a program for geostatistical modelling,
    # prediction and simulation." Computers & Geosciences, 24(1), 17–31
    variogram_dict = {
        "linear": pykrige.variogram_models.linear_variogram_model,
        "power": pykrige.variogram_models.power_variogram_model,
        "gaussian": pykrige.variogram_models.gaussian_variogram_model,
        "spherical": pykrige.variogram_models.spherical_variogram_model,
        "exponential": pykrige.variogram_models.exponential_variogram_model,
        "hole-effect": pykrige.variogram_models.hole_effect_variogram_model,
    }

    variogram_function = variogram_dict.get(variogram_model)

    # see reference in
    # https://github.com/GeoStat-Framework/PyKrige/blob/e02baad442ac99b22f038b09b6290e7abacc17ae/src/pykrige/core.py#L582
    estimated_model_parameters = pykrige.core._calculate_variogram_model(
        lags, semivariances, variogram_model, variogram_function, weight=False
    )

    # Prepare the parameters
    variogram_parameters = {}
    if variogram_model == "linear":
        variogram_parameters["slope"] = estimated_model_parameters[0]
        variogram_parameters["nugget"] = estimated_model_parameters[1]
    elif variogram_model == "power":
        variogram_parameters["scale"] = estimated_model_parameters[0]
        variogram_parameters["exponent"] = estimated_model_parameters[1]
        variogram_parameters["nugget"] = estimated_model_parameters[2]
    else:
        variogram_parameters["sill"] = estimated_model_parameters[0] + estimated_model_parameters[2]
        variogram_parameters["range"] = estimated_model_parameters[1]
        variogram_parameters["nugget"] = estimated_model_parameters[2]

    estimated_semivariances = variogram_function(estimated_model_parameters, lags)
    return variogram_parameters, (lags, estimated_semivariances, semivariances)


def setup_kriging_system(
    da: xr.DataArray,
    method="universal",
    empirical_variogram_args: dict = None,
    variogram_args: dict = None,
):
    """Kriging in space for a single time step.

    Make sure that coordinates 'x' and 'y' are present in the DataArray and they
    are in metric units.

    Parameters
    ----------
    da: xr.DataArray
        The DataArray containing the data to be interpolated. It must have
        coordinates 'x' and 'y'.
    method: str
        The kriging method to use, e.g. 'universal'. Default is 'universal'.
        Other methods are not implemented yet.
        See https://geostat-framework.readthedocs.io/projects/pykrige/en/stable/generated/pykrige.uk.UniversalKriging.html#pykrige.uk.UniversalKriging
        for more details.
    empirical_variogram_args: dict
        Additional keyword arguments to pass to the function
        `calculate_empirical_variogram`. Allowed keys are: 'method', 'nlags',
        'cutoff'. If it left empty, default parameters will be used as
        {"method":"standard", "nlags": 50, "cutoff"=10000.0}. The `method`
        argument can be one of 'standard', 'unbiased', 'unbiased_robust'.
        Default is 'standard'. See the documentation of
        `calculate_empirical_variogram` for more details.
    variogram_args: dict
        Additional keyword arguments to pass to the kriging method. Valid keys are:
        'variogram_model', 'variogram_parameters', 'drift_terms'.

        The `variogram_model` can be one of: 'linear', 'power', 'gaussian',
        'spherical', 'exponential', 'hole-effect'. Default is 'gaussian'.

        The variogram parameters can be provided as a dictionary, for example,
        # linear
            {'slope': slope, 'nugget': nugget}
        # power
            {'scale': scale, 'exponent': exponent, 'nugget': nugget}
        # gaussian, spherical, exponential and hole-effect:
            {'sill': s, 'range': r, 'nugget': n}
            # OR
            {'psill': p, 'range': r, 'nugget': n}
        If `variogram_parameters` are not provided, they will be estimated from
        the empirical variogram using the specified values in `empirical_variogram_args`.

        The `drift_terms` argument is only used for universal kriging. Supported drift
        terms are currently 'regional_linear', 'point_log', 'external_Z',
        'specified', and 'functional'. Default is 'regional_linear', which activates
        a first-order drift model.

    Returns
    -------
    kriging_obj: Any
        A kriging object that can be used to perform kriging interpolation.
    """
    # Check that da.data shape is 2d
    if len(da.data.shape) > 2:
        raise ValueError("DataArray must be 2D with coordinates 'x' and 'y'.")

    # Check if there x, y coords
    if "x" not in da.coords or "y" not in da.coords:
        raise ValueError("DataArray must have coordinates 'x' and 'y'.")

    # Calculate empirical variogram
    if not empirical_variogram_args:
        logger.info("Estimating variogram with default parameters.")
        empirical_variogram_args = {}

    _check_empirical_variogram_args(empirical_variogram_args)

    lags, semivariances = calculate_empirical_variogram(da, **empirical_variogram_args)

    # Check if variogram parameters are provided
    # if not, estimate them
    if not variogram_args:
        variogram_args = {}

    _check_variogram_args(variogram_args)

    variogram_model = variogram_args.get("variogram_model", "gaussian")
    variogram_parameters = variogram_args.get("variogram_parameters", None)
    if variogram_parameters is None:
        variogram_parameters, _ = fit_variogram(da, lags, semivariances, variogram_model)

    # Create a kriging instance
    if method == "universal":
        # see input arguments in
        # https://github.com/GeoStat-Framework/PyKrige/blob/e02baad442ac99b22f038b09b6290e7abacc17ae/src/pykrige/uk.py#L220
        kriging_obj = pykrige.uk.UniversalKriging(
            da.coords["x"],
            da.coords["y"],
            da,
            variogram_model=variogram_model,
            variogram_parameters=variogram_parameters,
            exact_values=False,  #  If True, results would be input values at input locations
            drift_terms=variogram_args.get("drift_terms", "regional_linear"),  # this activates drift of order 1
        )

        # Adjust some variables
        kriging_obj.lags = lags
        kriging_obj.semivariance = semivariances

        return kriging_obj
    else:
        raise NotImplementedError(f"{method} is not implemented yet.")


def solve_kriging_per_single_time(
    da: xr.DataArray,
    prediction_coords: xr.Dataset | xr.DataArray,
    method="universal",
    n_nearest_neighbors: int | None = None,
    kriging_backend: str = "vectorized",
    empirical_variogram_args: dict = None,
    variogram_args: dict = None,
):
    """Kriging in space for a single time step.

    Make sure that coordinates 'x' and 'y' are present in the DataArray and they
    are in metric units.

    Parameters
    ----------
    da: xr.DataArray
        The DataArray containing the variable of interest to be interpolated. It
        must have coordinates 'x' and 'y'.
    prediction_coords: xr.Dataset | xr.DataArray | None
        The prediction coordinates on which to interpolate the data. It should
        have coordinates 'x' and 'y'.
    method: str
        The kriging method to use, e.g. 'universal'. Default is 'universal'.
        Other methods are not implemented yet.
    n_nearest_neighbors: int | None
        The number of nearest neighbors to use for interpolation. If None, all
        points in the DataArray are used for interpolation.
    kriging_backend: str
        Specifies which approach to use in kriging. Specifying "vectorized" will solve
        the entire kriging problem at once in a vectorized operation. This approach is
        faster but also can consume a significant amount of memory for large grids
        and/or large datasets. Specifying "loop" will loop through each point at which
        the kriging system is to be solved. This approach is slower but also less
        memory-intensive. Default is "vectorized".
    empirical_variogram_args: dict
        Additional keyword arguments to pass to the function
        `calculate_empirical_variogram`. Allowed keys are: 'method', 'nlags',
        'cutoff'. If it left empty, default parameters will be used as
        {"method": "standard", "nlags": 50, "cutoff": 10000.0}. The `method`
        argument can be one of 'standard', 'unbiased', 'unbiased_robust'.
        Default is 'standard'. See the documentation of
        `calculate_empirical_variogram` for more details.
    variogram_args: dict
        Additional keyword arguments to pass to the kriging method. Valid keys
        are: 'variogram_model', 'variogram_parameters', 'drift_terms'.

        The `variogram_model` can be one of: 'linear', 'power', 'gaussian',
        'spherical', 'exponential', 'hole-effect'. Default is 'gaussian'.

        The variogram parameters can be provided as a dictionary, for example, #
        linear
            {'slope': slope, 'nugget': nugget}
        # power
            {'scale': scale, 'exponent': exponent, 'nugget': nugget}
        # gaussian, spherical, exponential and hole-effect:
            {'sill': s, 'range': r, 'nugget': n} # OR {'psill': p, 'range': r,
            'nugget': n}
        If `variogram_parameters` are not provided, they will be estimated from
        the empirical variogram using the specified values in
        `empirical_variogram_args`.

        The `drift_terms` argument is only used for universal kriging. Supported
        drift terms are currently 'regional_linear', 'point_log', 'external_Z',
        'specified', and 'functional'. Default is 'regional_linear', which
        activates a first-order drift model.

    Returns
    -------
    zvalues: np.ndarray
        The interpolated values at the prediction coordinates.
    sigmasq: np.ndarray
        The associated variance (sigmasq) for the interpolated values.
    """
    # setup kriging system
    kriging_obj = setup_kriging_system(da, method, empirical_variogram_args, variogram_args)

    # Check if prediction_coords has x, y coords
    if "x" not in prediction_coords.coords or "y" not in prediction_coords.coords:
        raise ValueError("Prediction coordinates must have coordinates 'x' and 'y'.")

    # if there is "space" in dimension,
    # style is points, otherwise it is a grid
    if "space" in prediction_coords.dims:
        interpolation_style = "points"
    else:
        interpolation_style = "grid"

    if n_nearest_neighbors is None:
        # Calculates a kriged grid and the associated variance
        # result has shape (M, N): M y coords and N x coords
        # all grid points are used for interpolation, more efficient
        zvalues, sigmasq = kriging_obj.execute(
            interpolation_style,
            prediction_coords.coords["x"],  # shape (N,)
            prediction_coords.coords["y"],  # shape (M,)
            backend=kriging_backend,
        )
        return zvalues.data, sigmasq.data  # numpy.ndarray
    else:
        if "space" not in prediction_coords.dims:
            raise NotImplementedError(
                "Kriging with nearest neighbors is not implemented for grid interpolation. "
                "Because this method can be very slow and memory intensive for a large grid. "
            )

        if n_nearest_neighbors > da.size:
            raise ValueError(
                f"n_nearest_neighbors ({n_nearest_neighbors}) cannot be larger than "
                f"the number of points in da ({da.size})."
            )

        # Find the nearest neighbors
        tree = KDTree(np.stack((da.coords["x"], da.coords["y"]), axis=1))
        _, indices = tree.query(
            np.stack((prediction_coords.coords["x"], prediction_coords.coords["y"]), axis=1), k=n_nearest_neighbors
        )

        neighbor_x = np.take(da.coords["x"].values, indices)
        neighbor_y = np.take(da.coords["y"].values, indices)
        neighbor_z = np.take(da.values, indices)

        def _apply_kriging_one_point(index):
            # This is a workaround to count for nearest neighbors
            # because pykrige does not support nearest neighbors
            kriging_obj.X_ADJUSTED = neighbor_x[index]
            kriging_obj.Y_ADJUSTED = neighbor_y[index]
            kriging_obj.Z = neighbor_z[index]

            zvalues, sigmasq = kriging_obj.execute(
                "points",
                prediction_coords.coords["x"].data[index],
                prediction_coords.coords["y"].data[index],
                backend="loop",  # use 'loop' backend for single point
            )
            return np.concatenate([zvalues, sigmasq])

        # Loop over each point in prediction_coords and apply krige
        zvalues = np.empty(indices.shape[0])
        sigmasq = np.empty(indices.shape[0])

        # TODO: Parallelize this loop if needed
        for index, _ in enumerate(indices):
            zvalues[index], sigmasq[index] = _apply_kriging_one_point(index)

        return zvalues, sigmasq


def solve_kriging(
    ps_atmosphere: xr.DataArray, prediction_coords: xr.Dataset | xr.DataArray, method="universal", **kwargs
):
    """Kriging in space to estimate the atmosphere signal time series.

    Parameters
    ----------
    ps_atmosphere: xr.DataArray
        The DataArray containing the atmosphere signal with coordinates 'x' and 'y'.
        It must have a time dimension.
    prediction_coords: xr.DataArray
        The prediction coordinates on which to interpolate the atmosphere
        signal. It should have coordinates 'x' and 'y'.
    method: str
        The kriging method to use, e.g. 'universal'. Default is 'universal'.
        Other methods are not implemented yet.
        see https://geostat-framework.readthedocs.io/projects/pykrige/en/stable/generated/pykrige.uk.UniversalKriging.html#pykrige.uk.UniversalKriging
        for more details.
    kwargs: dict
        Additional keyword arguments to pass to the function
        `solve_kriging_per_single_time`. Allowed keys are:
        "n_nearest_neighbors", "kriging_backend", "empirical_variogram_args",
        "variogram_args", See the documentation of
        `solve_kriging_per_single_time` for more details.

    Returns
    -------
    xr.Dataset
        A dataset containing the interpolated atmosphere signal and the associated
        variance (sigmasq) for each time step.
    """
    # Check if ps_atmosphere has time dimension
    if "time" not in ps_atmosphere.dims:
        raise ValueError(
            "ps_atmosphere must have a 'time' dimension. Otherwise, use `krige_per_single_time` function directly."
        )

    # Check if ps_atmosphere has 'x' and 'y' coordinates
    if "x" not in ps_atmosphere.coords or "y" not in ps_atmosphere.coords:
        raise ValueError("ps_atmosphere must have coordinates 'x' and 'y'.")

    # Remove "time" because we will apply kriging per time step
    input_core_dims = list(ps_atmosphere.sizes)
    if "time" in input_core_dims:
        input_core_dims.remove("time")

    coords_no_time = {k: v for k, v in ps_atmosphere.coords.items() if "time" not in v.dims}

    # Check if prediction_coords has "x" and "y" coordinates
    if "x" not in prediction_coords.coords or "y" not in prediction_coords.coords:
        raise ValueError("Prediction coordinates must have coordinates 'x' and 'y'.")

    # Check if 'time' in prediction_coords dims
    if "time" in prediction_coords.dims:
        raise ValueError("Prediction coordinates must not have 'time' dimension.")

    # Check if `input_core_dims` are not chunked
    if ps_atmosphere.chunks is not None:
        chunk_sizes = dict(zip(list(ps_atmosphere.sizes), ps_atmosphere.chunks, strict=False))
        if any(len(chunk_sizes[dim]) != 1 for dim in input_core_dims):
            logger.warning(
                "Rechunking ps_atmosphere to have non-chunked core dimensions for kriging. "
                f"Current chunks: {ps_atmosphere.chunks}"
            )
            ps_atmosphere = ps_atmosphere.chunk({dim: -1 for dim in input_core_dims})

    def apply_kriging_per_single_time(data: np.ndarray):
        """Apply kriging for a single time step."""
        da = xr.DataArray(data=data, coords=coords_no_time, dims=input_core_dims)
        predicted, sigmasq = solve_kriging_per_single_time(da, prediction_coords, method=method, **kwargs)
        return predicted, sigmasq

    predicted, sigmasq = xr.apply_ufunc(
        apply_kriging_per_single_time,
        ps_atmosphere,
        input_core_dims=[input_core_dims],
        output_core_dims=[
            list(prediction_coords.sizes)[::-1],
            list(prediction_coords.sizes)[::-1],
        ],  # result has shape (y, x)
        dask="parallelized",
        vectorize=True,
        output_dtypes=[ps_atmosphere.dtype, ps_atmosphere.dtype],
        dask_gufunc_kwargs={
            "output_sizes": dict(prediction_coords.sizes)
        },  # this is needed when prediction_coords has "x" and "y" coordinates
    )

    # Update time values
    predicted = predicted.assign_coords(time=ps_atmosphere["time"].data)
    sigmasq = sigmasq.assign_coords(time=ps_atmosphere["time"].data)

    return xr.Dataset({"predicted": predicted, "sigmasq": sigmasq})


def estimate_atmosphere_phase(
    stm: xr.Dataset,
    prediction_coords: xr.Dataset | xr.DataArray = None,
    psc_phase_residuals="psc_phase_residuals",
    atmosphere_mother: int | str = "atmosphere_mother",
    unmodeled_displacement_args: dict = None,
    kriging_args: dict = None,
) -> xr.Dataset:
    """Estimate the atmosphere phase.

    This function applies a temporal filter to extract the high-frequency
    atmospheric signal and then uses spatial kriging to predict the atmospheric
    phase per epoch.

    Parameters
    ----------
    stm: xr.Dataset
        The STM with unwrapped phases.
    prediction_coords: xr.Dataset | xr.DataArray
        The prediction coordinates on which to interpolate the atmosphere phase.
        It should have coordinates 'x' and 'y'. If None, the coordinates from
        the stm Dataset will be used.
    psc_phase_residuals: str
        The name of the variable in the stm Dataset that contains
        the PSC phase residuals. Default is "psc_phase_residuals".
    atmosphere_mother: int | str
        A string indicating the name of the variable in the stm Dataset
        that contains the atmosphere mother or an integer indicating the time
        index of the atmosphere.
    unmodeled_displacement_args: dict
        Keyword arguments for the `estimate_unmodeled_displacement` function.
        The allowed keys are: `filter_length`, `sampling_rate`, `filter_type`.
        For example: {"filter_length": 1, "sampling_rate": 1000, "filter_type":
        "gaussian"}. See the documentation of `estimate_unmodeled_displacement`
        for more details. The argument `unmodeled_displacement_args` can be left
        empty. Then default parameters i.e. {"filter_length": 1,
        "sampling_rate": 1000, "filter_type": "gaussian"} will be used.
    kriging_args: dict
        Additional keyword arguments to pass to the function
        `solve_kriging_per_single_time`. Allowed keys are:
        "n_nearest_neighbors", "kriging_backend", "empirical_variogram_args",
        "variogram_args", See the documentation of
        `solve_kriging_per_single_time` for more details.

    Returns
    -------
    xr.Dataset
        The predicted atmosphere phase and the associated variance (sigmasq) for
        each time step.
    """
    # Step 1: Apply temporal filtering to extract high-frequency atmospheric signal
    if not unmodeled_displacement_args:
        unmodeled_displacement_args = {}
    unmodeled_disp = estimate_unmodeled_displacement(
        psc_phase_residuals=stm[psc_phase_residuals],
        baseline_years=stm["time"],
        **unmodeled_displacement_args,
    )

    # Add results to stm
    stm["unmodeled_disp"] = unmodeled_disp

    # Get atmosphere mother
    if isinstance(atmosphere_mother, str):
        if atmosphere_mother not in stm:
            raise ValueError(f"atmosphere_mother '{atmosphere_mother}' not found in stm Dataset.")
        atmosphere_mother = stm[atmosphere_mother]
    else:
        atmosphere_mother = stm.isel(time=atmosphere_mother)[psc_phase_residuals]

    # Estimate atmosphere phase
    stm["atmosphere_estimates"] = stm[psc_phase_residuals] - stm["unmodeled_disp"] + atmosphere_mother

    # Step 2: Apply spatial kriging to predict atmospheric phase per epoch
    if prediction_coords is None:
        prediction_coords = xr.Dataset(coords=stm.isel(time=0).coords)

    if not kriging_args:
        kriging_args = {}
    predicted_atmosphere = solve_kriging(
        ps_atmosphere=stm["atmosphere_estimates"],
        prediction_coords=prediction_coords,
        **kriging_args,
    )
    # rename variables "predicted" to "atmosphere_predicted" and "sigmasq" to "atmosphere_sigmasq"
    predicted_atmosphere = predicted_atmosphere.rename(
        {"predicted": "atmosphere_predicted", "sigmasq": "atmosphere_sigmasq"}
    )
    return predicted_atmosphere
