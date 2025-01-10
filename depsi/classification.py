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


def detect_side_lobes(stm: xr.Dataset, max_pixel_dist: float, min_correlation: float) -> tuple[np.ndarray]:
    """Detect and mask side-lobe points based on the phase correlation between points.

    It first finds points on the same range and azimuth and only considers points close by. Then it
    computes the complex DD phase and computes the correlation.

    Parameters
    ----------
    stm : xarray.Dataset
      An input stm must include 'range', 'azimuth', 'pnt_idx', 'sd_complex', and 'nmad_full'.
    max_pixel_dist : float
      The maximum allowed spatial distance (in pixels) between points to be considered potential side-lobes.
    min_correlation : float
      The minimum correlation threshold to classify points as side-lobes. 0 means no correlation, 1 is maximum
      correlation

    Returns
    -------
    side_lobes_array : np.ndarray
      An array containing the indices of the detected side-lobe points.
    mask_side_lobes : np.ndarray
      A boolean mask where 'False' indicates detected side-lobe points.
    """
    # Lazy load variables
    range_vals = stm["range"].data
    azimuth_vals = stm["azimuth"].data
    sd_complex = stm["sd_complex"].data
    amplitude_vals = stm["sd_amplitude"].data
    # nmad_full_vals = stm["nmad_full"].data
    nr_epochs = len(stm.time)

    point_idx = stm["pnt_idx"].values

    # Define an empty set where the sidelobes will be stored
    side_lobes = set()

    for point in point_idx:
        # Skip the point if it is already detected as a sidelobe
        if point in side_lobes:
            continue

        # Get the range and azimuth coordinates of the point
        range_i = range_vals[point]
        azimuth_i = azimuth_vals[point]

        # Search for points close by with same range and azimuth coordinates
        idx_range = np.where(
            np.logical_and(
                range_vals == range_i,  # Same range value
                np.abs(azimuth_vals - azimuth_i) < max_pixel_dist,  # Within pixel distance
            )
        )[0]

        idx_azimuth = np.where(
            np.logical_and(
                azimuth_vals == azimuth_i,  # Same azimuth value
                np.abs(range_vals - range_i) < max_pixel_dist,  # Within pixel distance
            )
        )[0]

        potential_side_lobe_idx = np.union1d(idx_range, idx_azimuth)
        potential_side_lobe_idx = potential_side_lobe_idx.compute() if isinstance(potential_side_lobe_idx, da.Array) else potential_side_lobe_idx

        for point2 in potential_side_lobe_idx:
            if point2 != point and point2 not in side_lobes:  # Skip the current and already detected side-lobe points
                dd_complex = _compute_dd_for_correlation(
                    sd_complex[point, :], sd_complex[point2, :]
                )  # Compute DD between the two points
                corr = _calculate_phase_correlation(
                    dd_complex, nr_epochs
                )  # Check the phase difference between two points and compute correlation

                if (
                    corr >= min_correlation
                ):  # The  point with the lowest mean amplitude will be detected as the side-lobe
                    mean_ampl_p1 = np.mean(amplitude_vals[point, :])
                    mean_ampl_p2 = np.mean(amplitude_vals[point2, :])

                    if mean_ampl_p2 < mean_ampl_p1:
                        side_lobes.add(
                            point2
                        )  # point 2 has the lowest mean amplitude, so it is dected as the side-lobe
                    else:
                        side_lobes.add(point)  # point 1 has the lowest mean amplitude, so it is dected as the side-lobe

    # Make an array of the set
    side_lobes_array = np.array(list(side_lobes))

    mask_side_lobes = np.ones(len(point_idx), dtype=bool)  # Create a mask
    mask_side_lobes[side_lobes_array] = False

    return side_lobes_array, mask_side_lobes


def _calculate_phase_correlation(dd_complex, nr_epochs):
    """Compute correlation between phase time series of two pixels based on their double-difference (DD) phasors.

    This function calculates the phase similarity between two complex-valued time series
    by analyzing the angular differences in their double-difference (DD) phasors.
    The correlation is normalized over the number of epochs to produce a value between 0 and 1,
    where 1 indicates perfect correlation.

    Parameters
    ----------
    dd_complex : array-like
      A complex-valued array representing the double-difference phasors
                             between two pixels over multiple epochs.
    nr_epochs : int
      The number of time epochs (observations) over which the correlation is calculated.

    Returns
    -------
    float: A correlation value between 0 and 1, representing the phase similarity of the two time series.
    """
    corr = np.abs(np.sum(np.exp(1j * (np.angle(dd_complex))))) / nr_epochs
    return corr


def _compute_dd_for_correlation(complex_p1, complex_p2):
    """Compute the complex double-difference (DD) between two complex-valued time series.

    This function calculates the element-wise product of the complex conjugate of the first time series (`complex_p1`)
    and the second time series (`complex_p2`). The result represents the phase difference
    between the two series, which is useful for detecting similarities in phase behavior.

    Parameters
    ----------
    complex_p1 : (np.ndarray)
      A complex-valued array representing the first time series.
    complex_p2 : (np.ndarray)
      A complex-valued array representing the second time series.

    Returns
    -------
    np.ndarray: An array of complex values representing the phase differences (double differences)
                between the two input time series.
    """
    complex_conj_p1 = np.conj(complex_p1)
    dd_complex = complex_conj_p1 * complex_p2

    return dd_complex
