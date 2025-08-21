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

logger = getLogger(__name__)



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


def krige_per_single_time(
        da: xr.DataArray,
        grid: xr.Dataset | xr.DataArray | None = None,
        method='universal',
        n_nearest_neighbors: int | None = None,
        **kwargs
    ):
    """Kriging in space for a single time step.

    Make sure that coordinates 'x' and 'y' are present in the DataArray and they
    are in metric units.

    Parameters
    ----------
    da: xr.DataArray
        The DataArray containing the data to be interpolated. It must have
        coordinates 'x' and 'y'.
    grid: xr.Dataset | xr.DataArray | None
        The grid on which to interpolate the data. It should have coordinates
        'x' and 'y'. If None, the function returns a kriging object.
    method: str
        The kriging method to use, e.g. 'universal'. Default is 'universal'.
    n_nearest_neighbors: int | None
        The number of nearest neighbors to use for interpolation. If None, all
        points in the DataArray are used for interpolation.
    kwargs: dict
        Additional keyword arguments to pass to the kriging method, such as
        'variogram_model', 'variogram_parameters', etc.
    Returns
    -------
    zvalues: np.ndarray
        The interpolated values at the grid points.
    sigmasq: np.ndarray
        The associated variance (sigmasq) for the interpolated values.
    """
    # Check that da.data shape is 2d
    if len(da.data.shape) > 2:
        raise ValueError("DataArray must be 2D with coordinates 'x' and 'y'.")

    # Check if there x, y coords
    if 'x' not in da.coords or 'y' not in da.coords:
        raise ValueError("DataArray must have coordinates 'x' and 'y'.")

    # Define default variogram parameters
    # as implemented in Matlab version
    variogram_model = kwargs.get('variogram_model', 'gaussian')
    if variogram_model == 'gaussian':
        default_variogram_parameters = {
            'sill': 0.8 * da.var(),
            'range': 1000.0,
            'nugget': 0.2 * da.var()
        }
    else:
        default_variogram_parameters = None

    # Create a kriging instance
    if method == 'universal':
        kriging_obj = pykrige.uk.UniversalKriging(
            da.coords['x'],
            da.coords['y'],
            da,
            variogram_model=variogram_model,
            variogram_parameters=kwargs.get('variogram_parameters', default_variogram_parameters),
            nlags=kwargs.get('nlags', 50),
            exact_values=False,  #  If True, results would be input values at input locations
            drift_terms=['regional_linear']  # this activates drift of order 1
        )
    else:
        raise NotImplementedError(f"{method} is not implemented yet.")

    if grid is None:
        return kriging_obj

    # Check if grid has x, y coords
    if 'x' not in grid.coords or 'y' not in grid.coords:
        raise ValueError("Grid must have coordinates 'x' and 'y'.")

    # if there is "space" in dimension,
    # style is points, otherwise it is a grid
    if 'space' in grid.dims:
        interpolation_style = 'points'
    else:
        interpolation_style = 'grid'

    if n_nearest_neighbors is None:
        # Calculates a kriged grid and the associated variance
        # result has shape (M, N): M y coords and N x coords
        # all grid points are used for interpolation, more efficient
        zvalues, sigmasq = kriging_obj.execute(
            interpolation_style,
            grid.coords['x'],  # shape (N,)
            grid.coords['y'],  # shape (M,)
            backend=kwargs.get('backend', 'vectorized'),
        )
        return zvalues.data, sigmasq.data  # numpy.ndarray
    else:
        if 'space' not in da.dims and 'space' not in grid.dims:
            raise NotImplementedError(
                "Kriging with nearest neighbors is not implemented for grid interpolation."
            )

        if n_nearest_neighbors > da.size:
            raise ValueError(
                f"n_nearest_neighbors ({n_nearest_neighbors}) cannot be larger than "
                f"the number of points in da ({da.size})."
            )

        # Find the nearest neighbors
        tree = KDTree(np.stack((da.coords['x'], da.coords['y']), axis=1))
        _, indices = tree.query(
            np.stack((grid.coords['x'], grid.coords['y']), axis=1),
            k=n_nearest_neighbors
        )

        neighbor_x = np.take(da.coords['x'].values, indices)
        neighbor_y = np.take(da.coords['y'].values, indices)
        neighbor_z = np.take(da.values, indices)

        def _apply_kriging_one_point(index):
            # This is a workaround to count for nearest neighbors
            # because pykrige does not support nearest neighbors
            kriging_obj.X_ADJUSTED = neighbor_x[index]
            kriging_obj.Y_ADJUSTED = neighbor_y[index]
            kriging_obj.Z = neighbor_z[index]

            zvalues, sigmasq = kriging_obj.execute(
                'points',
                grid.coords['x'].data[index],
                grid.coords['y'].data[index],
                backend='loop'
                )
            return np.concatenate([zvalues, sigmasq])

        # Loop over each point in grid and krige
        zvalues = np.empty(indices.shape[0])
        sigmasq = np.empty(indices.shape[0])
        for index, _ in enumerate(indices):
            zvalues[index], sigmasq[index] = _apply_kriging_one_point(index)

        return zvalues, sigmasq


def krige_in_space(
    ps_atmosphere: xr.DataArray,
    grid: xr.Dataset | xr.DataArray ,
    method='universal',
    **kwargs
):
    """Kriging in space to estimate the atmosphere signal.

    Parameters
    ----------
    ps_atmosphere: xr.DataArray
        The DataArray containing the atmosphere signal with coordinates 'x' and 'y'.
        It must have a time dimension.
    grid: xr.DataArray
        The grid on which to interpolate the atmosphere signal. It should have
        coordinates 'x' and 'y'.
    method: str
        The kriging method to use, e.g. 'universal'. Default is 'universal'.
    kwargs: dict
        Additional keyword arguments to pass to the kriging method, such as
        'variogram_model', 'variogram_parameters', etc.

    Returns
    -------
    xr.Dataset
        A dataset containing the interpolated atmosphere signal and the associated
        variance (sigmasq) for each time step.
    """
    # Check if ps_atmosphere has time dimension
    if 'time' not in ps_atmosphere.dims:
        raise ValueError(
            "ps_atmosphere must have a 'time' dimension. "
            "Otherwise, use `krige_per_single_time` function directly."
        )

    # Check if ps_atmosphere has 'x' and 'y' coordinates
    if 'x' not in ps_atmosphere.coords or 'y' not in ps_atmosphere.coords:
        raise ValueError("ps_atmosphere must have coordinates 'x' and 'y'.")

    # Remove "time" because we will apply kriging per time step
    input_core_dims = list(ps_atmosphere.sizes)
    if 'time' in input_core_dims:
        input_core_dims.remove('time')

    # Check if grid has 'x' and 'y' coordinates
    if 'x' not in grid.coords or 'y' not in grid.coords:
        raise ValueError("Grid must have coordinates 'x' and 'y'.")

    # Check if `input_core_dims` are not chunked
    if ps_atmosphere.chunks is not None:
        chunk_sizes = dict(zip(list(ps_atmosphere.sizes), ps_atmosphere.chunks, strict=False))
        if any(len(chunk_sizes[dim]) !=1 for dim in input_core_dims):
            raise ValueError(
                "ps_atmosphere must not be chunked in the core dimensions "
                f"{input_core_dims}."
            )

    def apply_krige_per_single_time(data: np.ndarray):
        """Apply kriging for a single time step."""
        da = xr.DataArray(
            data=data,
            coords=ps_atmosphere.coords,
            dims=input_core_dims
        )

        interpolated, sigmasq = krige_per_single_time(da, grid, method=method, **kwargs)
        return interpolated, sigmasq

    interpolated, sigmasq = xr.apply_ufunc(
        apply_krige_per_single_time,
        ps_atmosphere,
        input_core_dims=[input_core_dims],
        output_core_dims=[list(grid.sizes)[::-1], list(grid.sizes)[::-1]], # result has shape (y, x)
        dask="parallelized",
        vectorize=True,
        output_dtypes=[ps_atmosphere.dtype, ps_atmosphere.dtype],
        dask_gufunc_kwargs = {"output_sizes": dict(grid.sizes)},  # this is needed when grid has "x" and "y" coordinates
    )

    # Update time values
    interpolated = interpolated.assign_coords(
        time = ps_atmosphere["time"].data
    )
    sigmasq = sigmasq.assign_coords(
        time = ps_atmosphere["time"].data
    )

    return xr.Dataset({
        'interpolated': interpolated,
        'sigmasq': sigmasq
    })


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
