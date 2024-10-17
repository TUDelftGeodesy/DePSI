"""Functions for scatterer selection related operations."""

from typing import Literal

import numpy as np
import xarray as xr
from scipy.spatial import KDTree


def ps_selection(
    slcs: xr.Dataset,
    threshold: float,
    method: Literal["nad", "nmad"] = "nad",
    output_chunks: int = 10000,
    mem_persist: bool = False,
) -> xr.Dataset:
    """Select Persistent Scatterers (PS) from an SLC stack, and return a Space-Time Matrix.

    The selection method is defined by `method` and `threshold`.
    The selected pixels will be reshaped to (space, time), where `space` is the number of selected pixels.
    The unselected pixels will be discarded.
    The original `azimuth` and `range` coordinates will be persisted.
    The computed NAD or NMAD will be added to the output dataset as a new variable. It can be persisted in
    memory if `mem_persist` is True.

    Parameters
    ----------
    slcs : xr.Dataset
        Input SLC stack. It should have the following dimensions: ("azimuth", "range", "time").
        There should be a `amplitude` variable in the dataset.
    threshold : float
        Threshold value for selection.
    method : Literal["nad", "nmad"], optional
        Method of selection, by default "nad".
        - "nad": Normalized Amplitude Dispersion
        - "nmad": Normalized median absolute deviation
    output_chunks : int, optional
        Chunk size in the `space` dimension, by default 10000
    mem_persist : bool, optional
        If true persist the NAD or NMAD in memory, by default False.


    Returns
    -------
    xr.Dataset
        Selected STM, in form of an xarray.Dataset with two dimensions: (space, time).

    Raises
    ------
    NotImplementedError
        Raised when an unsupported method is provided.
    """
    # Make sure there is no temporal chunk
    # since later a block function assumes all temporal data is available in a spatial block
    slcs = slcs.chunk({"time": -1})

    # Calculate selection mask
    match method:
        case "nad":
            nad = xr.map_blocks(
                _nad_block, slcs["amplitude"], template=slcs["amplitude"].isel(time=0).drop_vars("time")
            )
            nad = nad.compute() if mem_persist else nad
            slcs = slcs.assign(pnt_nad=nad)
            mask = nad < threshold
        case "nmad":
            nmad = xr.map_blocks(
                _nmad_block, slcs["amplitude"], template=slcs["amplitude"].isel(time=0).drop_vars("time")
            )
            nmad = nmad.compute() if mem_persist else nmad
            slcs = slcs.assign(pnt_nmad=nmad)
            mask = nmad < threshold
        case _:
            raise NotImplementedError

    # Get the 1D index on space dimension
    mask_1d = mask.stack(space=("azimuth", "range")).drop_vars(["azimuth", "range", "space"])  # Drop multi-index coords
    index = mask_1d["space"].where(mask_1d.compute(), other=0, drop=True)  # Evaluate the 1D mask to index

    # Reshape from Stack ("azimuth", "range", "time") to Space-Time Matrix  ("space", "time")
    stacked = slcs.stack(space=("azimuth", "range"))

    # Drop multi-index coords for space coordinates
    # This will also azimuth and range coordinates, as they are part of the multi-index coordinates
    stm = stacked.drop_vars(["space", "azimuth", "range"])

    # Assign a continuous index the space dimension
    # Assign azimuth and range back as coordinates
    stm = stm.assign_coords(
        {
            "space": (["space"], range(stm.sizes["space"])),
            "azimuth": (["space"], stacked["azimuth"].values),
            "range": (["space"], stacked["range"].values),
        }
    )  # keep azimuth and range as coordinates

    # Apply selection
    stm_masked = stm.sel(space=index)

    # Re-order the dimensions to community preferred ("space", "time") order
    stm_masked = stm_masked.transpose("space", "time")

    # Rechunk is needed because after apply maksing, the chunksize will be inconsistant
    stm_masked = stm_masked.chunk(
        {
            "space": output_chunks,
            "time": -1,
        }
    )

    # Reset space coordinates
    stm_masked = stm_masked.assign_coords(
        {
            "space": (["space"], range(stm_masked.sizes["space"])),
        }
    )

    # Compute NAD or NMAD if mem_persist is True
    # This only evaluate a very short task graph, since NAD or NMAD is already in memory
    if mem_persist:
        match method:
            case "nad":
                stm_masked["pnt_nad"] = stm_masked["pnt_nad"].compute()
            case "nmad":
                stm_masked["pnt_nmad"] = stm_masked["pnt_nmad"].compute()

    return stm_masked


def network_stm_seletcion(
    stm: xr.Dataset,
    min_dist: int | float,
    include_index: list[int] = None,
    sortby_var: str = "pnt_nmad",
    crs: int | str = "radar",
    x_var: str = "azimuth",
    y_var: str = "range",
    azimuth_spacing: float = None,
    range_spacing: float = None,
):
    """Select a Space-Time Matrix (STM) from a candidate STM for network processing.

    The selection is based on two criteria:
    1. A minimum distance between selected points.
    2. A sorting metric to select better points.

    The candidate STM will be sorted by the sorting metric.
    The selection will be performed iteratively, starting from the best point.
    In each iteration, the best point will be selected, and points within the minimum distance will be removed.
    The process will continue until no points are left in the candidate STM.

    Parameters
    ----------
    stm : xr.Dataset
        candidate Space-Time Matrix (STM).
    min_dist : int | float
        Minimum distance between selected points.
    include_index : list[int], optional
        Index of points in the candidate STM that must be included in the selection, by default None
    sortby_var : str, optional
        Sorting metric for selecting points, by default "pnt_nmad"
    crs : int | str, optional
        EPSG code of Coordinate Reference System of `x_var` and `y_var`, by default "radar".
        If crs is "radar", the distance will be calculated based on radar coordinates, and
        azimuth_spacing and range_spacing must be provided.
    x_var : str, optional
        Data variable name for x coordinate, by default "azimuth"
    y_var : str, optional
        Data variable name for y coordinate, by default "range"
    azimuth_spacing : float, optional
        Azimuth spacing, by default None. Required if crs is "radar".
    range_spacing : float, optional
        Range spacing, by default None. Required if crs is "radar".

    Returns
    -------
    xr.Dataset
        Selected network Space-Time Matrix (STM).

    Raises
    ------
    ValueError
        Raised when `azimuth_spacing` or `range_spacing` is not provided for radar coordinates.
    NotImplementedError
        Raised when an unsupported Coordinate Reference System is provided.
    """
    match crs:
        case "radar":
            if (azimuth_spacing is None) or (range_spacing is None):
                raise ValueError("Azimuth and range spacing must be provided for radar coordinates.")
        case _:
            raise NotImplementedError

    # Get coordinates and sorting metric, load them into memory
    stm_select = None
    stm_remain = stm[[x_var, y_var, sortby_var]].compute()

    # Select the include_index if provided
    if include_index is not None:
        stm_select = stm_remain.isel(space=include_index)

        # Remove points within min_dist of the included points
        coords_include = np.column_stack(
            [stm_select["azimuth"].values * azimuth_spacing, stm_select["range"].values * range_spacing]
        )
        coords_remain = np.column_stack(
            [stm_remain["azimuth"].values * azimuth_spacing, stm_remain["range"].values * range_spacing]
        )
        idx_drop = _idx_within_distance(coords_include, coords_remain, min_dist)
        if idx_drop is not None:
            stm_remain = stm_remain.where(~(stm_remain["space"].isin(idx_drop)), drop=True)

    # Reorder the remaining points by the sorting metric
    stm_remain = stm_remain.sortby(sortby_var)

    while stm_remain.sizes["space"] > 0:
        print(f"Remaining points: {stm_remain.sizes['space']}")

        # Select one point with best sorting metric
        stm_now = stm_remain.isel(space=0)

        # Add the selected point to the selection
        if stm_select is None:
            stm_select = stm_now.copy()
        else:
            stm_select = xr.concat([stm_select, stm_now], dim="space")

        # Remove the selected point from the remaining points
        stm_remain = stm_remain.isel(space=slice(1, None)).copy()

        # Remove points in stm_remain within min_dist of stm_now
        coords_remain = np.column_stack(
            [stm_remain["azimuth"].values * azimuth_spacing, stm_remain["range"].values * range_spacing]
        )
        coords_stmnow = np.column_stack(
            [stm_now["azimuth"].values * azimuth_spacing, stm_now["range"].values * range_spacing]
        )
        idx_drop = _idx_within_distance(coords_stmnow, coords_remain, min_dist)
        if idx_drop is not None:
            stm_drop = stm_remain.isel(space=idx_drop)
            stm_remain = stm_remain.where(~(stm_remain["space"].isin(stm_drop["space"])), drop=True)

    # Get the selected points by space index from the original stm
    stm_out = stm.sel(space=stm_select["space"].values)

    return stm_out


def _nad_block(amp: xr.DataArray) -> xr.DataArray:
    """Compute Normalized Amplitude Dispersion (NAD) for a block of amplitude data.

    Parameters
    ----------
    amp : xr.DataArray
        Amplitude data, with dimensions ("azimuth", "range", "time").
        This can be extracted from an SLC xr.Dataset.

    Returns
    -------
    xr.DataArray
        Normalized Amplitude Dispersion (NAD) data, with dimensions ("azimuth", "range").
    """
    # Compute amplitude dispersion
    # By defalut, the mean and std function from Xarray will skip NaN values
    # However, if there is NaN value in time series, we want to discard the pixel
    # Therefore, we set skipna=False
    # Adding epsilon to avoid zero division
    nad_da = amp.std(dim="time", skipna=False) / (amp.mean(dim="time", skipna=False) + np.finfo(amp.dtype).eps)

    return nad_da


def _nmad_block(amp: xr.DataArray) -> xr.DataArray:
    """Compute Normalized Median Absolute Deviation(NMAD) for a block of amplitude data.

    Parameters
    ----------
    amp : xr.DataArray
        Amplitude data, with dimensions ("azimuth", "range", "time").
        This can be extracted from an SLC xr.Dataset.

    Returns
    -------
    xr.DataArray
        Normalized Median Absolute Dispersion (NMAD) data, with dimensions ("azimuth", "range").
    """
    # Compoute NMAD
    median_amplitude = amp.median(dim="time", skipna=False)
    mad = (np.abs(amp - median_amplitude)).median(dim="time")  # Median Absolute Deviation
    nmad = mad / (median_amplitude + np.finfo(amp.dtype).eps)  # Normalized Median Absolute Deviation

    return nmad


def _idx_within_distance(coords_ref, coords_others, min_dist):
    """Get the index of points in coords_others that are within min_dist of coords_ref.

    Parameters
    ----------
    coords_ref : np.ndarray
        Coordinates of reference points. Shape (n, 2).
    coords_others : np.ndarray
        Coordinates of other points. Shape (m, 2).
    min_dist : int, float
        distance threshold.

    Returns
    -------
    np.ndarray
        Index of points in coords_others that are within `min_dist` of `coords_ref`.
    """
    kd_ref = KDTree(coords_ref)
    kd_others = KDTree(coords_others)
    sdm = kd_ref.sparse_distance_matrix(kd_others, min_dist)
    if len(sdm) > 0:
        idx = np.array(list(sdm.keys()))[:, 1]
        return idx
    else:
        return None
