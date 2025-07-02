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

    # Calculates a kriged grid and the associated variance
    # result has shape (M, N): M y coords and N x coords
    return kriging_obj.execute(
        'grid',
        grid.coords['x'],  # shape (N,)
        grid.coords['y'],  # shape (M,)
        backend=kwargs.get('backend', 'vectorized'),
    )


def _create_grid(bbox, grid_size: int = 100):
    """Create a grid based on the bounding box and grid size."""
    if isinstance(bbox, list | tuple) or len(bbox) == 4:
        x_min, y_min, x_max, y_max = bbox
    else:
        raise ValueError("Bounding box must be a list or tuple of four elements: [x_min, y_min, x_max, y_max].")

    new_x = np.arange(x_min, x_max + grid_size, grid_size)
    new_y = np.arange(y_min, y_max + grid_size, grid_size)
    return xr.Dataset(coords={'x': new_x, 'y': new_y})


def krige_in_space(
    ps_atmosphere: xr.DataArray,
    grid: xr.Dataset | None = None,
    grid_size: int = 100,
    method='universal',
    **kwargs
):
    """Kriging in space to estimate the atmosphere signal."""
    # Check if ps_atmosphere has time dimension
    if 'time' not in ps_atmosphere.dims:
        raise ValueError(
            "ps_atmosphere must have a 'time' dimension. "
            "Otherwise, use `krige_per_single_time` function directly."
        )

    # Check if ps_atmosphere has 'x' and 'y' coordinates
    if 'x' not in ps_atmosphere.coords or 'y' not in ps_atmosphere.coords:
        raise ValueError("ps_atmosphere must have coordinates 'x' and 'y'.")

    # Check if ps_atmosphere is chunked in time
    if ps_atmosphere.chunks is not None and 'time' not in ps_atmosphere.chunks:
        logger.info(
            "It is better if `ps_atmosphere` is chunked in time. "
            "Use `chunk` method to chunk it in time."
        )

    if grid is None:
        bbox = [
            ps_atmosphere["x"].min(),
            ps_atmosphere["y"].min(),
            ps_atmosphere["x"].max(),
            ps_atmosphere["y"].max()
        ]
        grid = _create_grid(bbox, grid_size)

    def apply_krige_per_single_time(data: np.ndarray):
        """Apply kriging for a single time step."""
        da = xr.DataArray(
            data,
            dims=['space'],
            coords={'x': ('space', ps_atmosphere["x"].data), 'y': ('space', ps_atmosphere["y"].data)}
        )

        interpolated, sigmasq = krige_per_single_time(da, grid, method=method, **kwargs)
        return interpolated.data, sigmasq.data

    # setting output_sizes gives `FutureWarning`, but it is needed when
    # ps_atmosphere is chunked, see
    # https://docs.xarray.dev/en/stable/generated/xarray.apply_ufunc.html
    interpolated, sigmasq = xr.apply_ufunc(
        apply_krige_per_single_time,
        ps_atmosphere,
        input_core_dims=[['space']],
        output_core_dims=[['y', 'x'], ['y', 'x']],
        dask='parallelized',
        vectorize=True,
        output_dtypes=[ps_atmosphere.dtype, ps_atmosphere.dtype],
        output_sizes={'y': len(grid["y"].data), 'x': len(grid["x"].data)}
    )

    # Update coords values
    interpolated = interpolated.assign_coords(
        x = grid["x"].data,
        y = grid["y"].data,
        time = ps_atmosphere["time"].data
    )
    sigmasq = sigmasq.assign_coords(
        x = grid["x"].data,
        y = grid["y"].data,
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
