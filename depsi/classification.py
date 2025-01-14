"""Functions for scatterer selection related operations."""

from datetime import datetime
from typing import Literal

import dask.array as da
import numpy as np
import pyproj
import xarray as xr
from scipy.spatial import KDTree

from depsi.point_quality import _estimate_breakpoints
from depsi.utils import _npdatetime64_to_datetime, crop_slc_spacetime

REQUIRED_BREAKPOINT_KEYS = ["db_segmentation", "search_method", "cost_function", "min_obs_partition"]


def ps_selection(
    slcs: xr.Dataset,
    threshold: float,
    method: Literal["nad", "nmad"] = "nad",
    output_chunks: int = 10000,
    mem_persist: bool = False,
    ps_selection_start_date: datetime | str | None = None,
    ps_selection_end_date: datetime | str | int | None = None,
    recalibration_jump_size: int = 10,
    do_rd_coordinate_conversion: bool = False,
    do_breakpoint_analysis: bool = False,
    breakpoint_kwargs: dict | None = None,
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
        Threshold value for selection for "nad" / "nmad".
    method : Literal["nad", "nmad"], optional
        Method of selection, by default "nad".
        - "nad": Normalized Amplitude Dispersion
        - "nmad": Normalized median absolute deviation
    output_chunks : int, optional
        Chunk size in the `space` dimension, by default 10000
    mem_persist : bool, optional
        If true persist the NAD or NMAD in memory, by default False.
    ps_selection_start_date : datetime | str | None, optional
      the start date of the time window to be used for the ps_selection, in one of three formats:
      - datetime object
      - str object, formatted as YYYYMMDD
      - None, no cropping in time requested for the ps_selection (default)
    ps_selection_end_date : datetime | str | int | None, optional
      the end date of the time window to be used for the ps_selection, in one of four formats:
      - datetime object
      - str object, formatted as YYYYMMDD
      - int object, which is interpreted as the number of images intended in the crop (including the start date). If
        more images are requested than exist since the start date, all images from start_date until the last image
        are provided.
      - None, no cropping in time requested for the ps_selection (default)
    recalibration_jump_size: int, optional
      the number of images that the recalibration NAD / NMAD variables remains constant. Will start after the
      initialization epoch (if ps_selection_start_date and ps_selection_end_date are not None). Defaults to 10.
    do_rd_coordinate_conversion: bool, optional
      boolean to trigger coordinate conversion from latitude/longitude (WGS84) to RD_X/RD_Y (Rijksdriehoek). This only
      makes sense for AoIs located in the Netherlands. Defaults to False.
    do_breakpoint_analysis: bool, optional
      boolean to trigger the breakpoint analysis. Defaults to False.
    breakpoint_kwargs: dict | None, optional
      the keyword arguments required for the breakpoint analysis. Required if do_breakpoint_analysis is set to True.
      Formatted as a dictionary with required keys:
      - db_segmentation: True or False, whether or not to do partitioning in dB. Advised False
      - search_method: 'pelt' or 'binseg'. Advised 'pelt'
      - cost_function: 'l#' with # replaced by 0-3. Advised 'l2'
      - min_obs_partition: integer. Advised min 0.5 years converted to # images, for Sentinel-1 27 (6 day interval)


    Returns
    -------
    xr.Dataset
        Selected STM, in form of an xarray.Dataset with two dimensions: (space, time).

    Raises
    ------
    NotImplementedError
        Raised when an unsupported method is provided.
    """
    if do_breakpoint_analysis:
        assert breakpoint_kwargs is not None, "Breakpoint analysis requested without keyword arguments!"
        assert isinstance(breakpoint_kwargs, dict), f"breakpoint_kwargs should be dict but is {type(breakpoint_kwargs)}"
        assert np.all(
            [key in breakpoint_kwargs.keys() for key in REQUIRED_BREAKPOINT_KEYS]
        ), f"Keys {REQUIRED_BREAKPOINT_KEYS} are required but received {breakpoint_kwargs.keys()}!"

    # Make sure there is no temporal chunk
    # since later a block function assumes all temporal data is available in a spatial block
    slcs = slcs.chunk({"time": -1})

    # Calculate selection mask
    match method:
        case "nad":
            if ps_selection_start_date is not None:
                selection_crop = crop_slc_spacetime(
                    slcs, start_date=ps_selection_start_date, end_date=ps_selection_end_date
                )
                nad = xr.map_blocks(
                    _nad_block,
                    selection_crop["amplitude"],
                    template=selection_crop["amplitude"].isel(time=0).drop_vars("time"),
                )
                ps_selection_times = selection_crop["time"].values
            else:
                nad = xr.map_blocks(
                    _nad_block, slcs["amplitude"], template=slcs["amplitude"].isel(time=0).drop_vars("time")
                )
                ps_selection_times = []
            nad = nad.compute() if mem_persist else nad
            slcs = slcs.assign(pnt_nad=nad)
            mask = nad < threshold
        case "nmad":
            if ps_selection_start_date is not None:
                selection_crop = crop_slc_spacetime(
                    slcs, start_date=ps_selection_start_date, end_date=ps_selection_end_date
                )
                nmad = xr.map_blocks(
                    _nmad_block,
                    selection_crop["amplitude"],
                    template=selection_crop["amplitude"].isel(time=0).drop_vars("time"),
                )
                ps_selection_times = selection_crop["time"].values
            else:
                nmad = xr.map_blocks(
                    _nmad_block, slcs["amplitude"], template=slcs["amplitude"].isel(time=0).drop_vars("time")
                )
                ps_selection_times = []
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

    # add incremental and recalibration NAD / NMAD
    incremental_imgs = []
    recalibration_imgs = []
    recalibration_idx = -1  # start at -1 so that the first addition will trigger a reset of the data layer
    recalibration_data_layer = None
    for date in stm_masked["time"].values:
        if date in ps_selection_times:  # only gets triggered in case there is an initialization epoch for the duration
            # of the initialization epoch
            start_date = ps_selection_start_date
            end_date = ps_selection_end_date
            recalibration_idx = 0
        else:
            start_date = _npdatetime64_to_datetime(stm_masked["time"].values[0])
            end_date = _npdatetime64_to_datetime(date)
            recalibration_idx += 1  # add, and do modulo the jump size, so that it will be 0 every time a new image
            # should be loaded
            recalibration_idx %= recalibration_jump_size
        current_crop = crop_slc_spacetime(stm_masked, start_date=start_date, end_date=end_date)
        match method:
            case "nad":
                nad = xr.map_blocks(
                    _nad_block,
                    current_crop["amplitude"],
                    template=current_crop["amplitude"].isel(time=0).drop_vars("time"),
                )
                incremental_imgs.append(nad)
                if recalibration_idx == 0:
                    recalibration_data_layer = nad.copy()
            case "nmad":
                nmad = xr.map_blocks(
                    _nmad_block,
                    current_crop["amplitude"],
                    template=current_crop["amplitude"].isel(time=0).drop_vars("time"),
                )
                incremental_imgs.append(nmad)
                if recalibration_idx == 0:
                    recalibration_data_layer = nmad.copy()
            case _:
                raise NotImplementedError

        recalibration_imgs.append(recalibration_data_layer.copy())

    # format the images, and add them to the data array
    incremental_nad_nmad = da.vstack(incremental_imgs).T
    recalibration_nad_nmad = da.vstack(recalibration_imgs).T
    match method:
        case "nad":
            stm_masked_inc = stm_masked.assign({"incremental_nad": (["space", "time"], incremental_nad_nmad)})
            stm_masked_inc = stm_masked_inc.assign({"recalibration_nad": (["space", "time"], recalibration_nad_nmad)})
        case "nmad":
            stm_masked_inc = stm_masked.assign({"incremental_nmad": (["space", "time"], incremental_nad_nmad)})
            stm_masked_inc = stm_masked_inc.assign({"recalibration_nmad": (["space", "time"], recalibration_nad_nmad)})
        case _:
            raise NotImplementedError

    # Rechunk is needed because after calculating incremental NAD/NMAD, the chunksize will be inconsistant
    stm_masked_inc = stm_masked_inc.chunk(
        {
            "space": output_chunks,
            "time": -1,
        }
    )

    # Add RD coordinates if requested
    if do_rd_coordinate_conversion:
        wgs84 = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:28992", always_xy=True).transform
        # Convert Lat and Lon to RD-coordinates
        rd_x, rd_y = wgs84(stm_masked_inc["lon"], stm_masked_inc["lat"])

        # Add RD coordinates to the dataset
        stm_masked_inc = stm_masked_inc.assign({"rd_x": (["space"], rd_x)})
        stm_masked_inc = stm_masked_inc.assign({"rd_y": (["space"], rd_y)})

    # Add extra time coordinate variables for time intervals since first image
    days = np.array(
        [
            (_npdatetime64_to_datetime(date) - _npdatetime64_to_datetime(stm_masked["time"].values[0])).days
            for date in stm_masked["time"].values
        ]
    )
    stm_masked_inc = stm_masked_inc.assign({"days_since_first_img": (["time"], days)})
    stm_masked_inc = stm_masked_inc.assign({"years_since_first_img": (["time"], days / 365.2425)})

    if do_breakpoint_analysis:
        breakpoints, partition_identifiers = _estimate_breakpoints(
            stm_masked_inc["amplitude"],
            breakpoint_kwargs["db_segmentation"],
            breakpoint_kwargs["search_method"],
            breakpoint_kwargs["cost_function"],
            breakpoint_kwargs["min_obs_partition"],
        )
        stm_masked_inc = stm_masked_inc.assign({"breakpoints": (["space", "time"], breakpoints)})
        stm_masked_inc = stm_masked_inc.assign({"partition_id": (["space", "time"], partition_identifiers)})

    # Compute NAD or NMAD if mem_persist is True
    # This only evaluate a very short task graph, since NAD or NMAD is already in memory
    if mem_persist:
        match method:
            case "nad":
                for key in ["pnt_nad", "inremental_nad"]:
                    stm_masked_inc[key] = stm_masked[key].compute()
            case "nmad":
                for key in ["pnt_nmad", "inremental_nmad"]:
                    stm_masked_inc[key] = stm_masked[key].compute()

    return stm_masked_inc


def network_stm_selection(
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

    # Build a list of the index of selected points
    if stm_select is None:
        space_idx_sel = []
    else:
        space_idx_sel = stm_select["space"].values.tolist()

    while stm_remain.sizes["space"] > 0:
        # Select one point with best sorting metric
        stm_now = stm_remain.isel(space=0)

        # Append the selected point index
        space_idx_sel.append(stm_now["space"].values.tolist())

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
    stm_out = stm.sel(space=space_idx_sel)

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
