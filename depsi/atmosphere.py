"""A module to estimate atmosphere signal from network STMs with unwrapped phases.

This estimation is done in two steps in general:
1. A temporal filtering applied per point to extract the temporal high frequency
signal.
2. A spatial kriging filtering per epoch to estimate the atmosphere signal per
epoch.
"""

from logging import getLogger

import numpy as np
import pykrige
import xarray as xr
from scipy import signal
from scipy.ndimage import convolve1d

logger = getLogger(__name__)


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
        The PSC time series to apply the filter to.
    baseline_years: xr.DataArray
        The baseline years corresponding to the time series.
    filter_length: int
        Length of the filter (year) to apply a low-pass filter to the time series.
        This is used to determine the size of the window i.e. 2 * filter_length + 1.
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
    std_dev = filter_length / 3  # Standard deviation defined based on a Gaussian function
    window_size = int(6 * std_dev) | 1  # Window size is ±3 standard deviations

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

    return xr.apply_ufunc(
        apply_filter,
        psc_phase,
        input_core_dims=[['time']],
        output_core_dims=[['time']],
        vectorize=True,
        dask='parallelized',
        output_dtypes=[psc_phase.dtype]
    )

def krige_per_single_time(
        da: xr.DataArray,
        grid: xr.Dataset | xr.DataArray | None = None,
        method='universal',
        **kwargs
    ):
    """Kriging in space for a single time step.

    Make sure that coordinates 'x' and 'y' are present in the DataArray and they
    are in metric units.
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

    # Calculates a kriged grid and the associated variance
    # result has shape (M, N): M y coords and N x coords
    return kriging_obj.execute(
        interpolation_style,
        grid.coords['x'],  # shape (N,)
        grid.coords['y'],  # shape (M,)
        backend=kwargs.get('backend', 'vectorized'),
    )

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
        return interpolated.data, sigmasq.data

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

    # # Update time values
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


def estimate_atmosphere_phase(stm, grid=None):
    """Estimate the atmosphere phase.

    This function applies a temporal filter to extract the high-frequency
    atmospheric signal and then uses spatial kriging to estimate the atmospheric
    phase per epoch.
    """
    # Step 1: Apply temporal filtering to extract high-frequency atmospheric signal
    non_linear = estimate_non_linear_deformation(
        psc_phase=stm["d_phase"],
        baseline_years=stm["time"],
        filter_length=30,
        method='gaussian'
    )
    stm["atmosphere_estimated"] = stm["d_phase"] - non_linear + stm["atmosphere_base"]

    # Step 2: Apply spatial kriging to estimate atmospheric phase per epoch
    interpolated = krige_in_space(
        ps_atmosphere=stm["atmosphere_estimated"],
        grid=grid,
        grid_size=100,
    )

    return interpolated
