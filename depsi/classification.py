"""Functions for scatter selection related operations."""

from typing import Literal

import numpy as np
import xarray as xr


def ps_selection(
    slcs: xr.Dataset,
    threshold: float,
    method: Literal["nad", "nmad"] = "nad",
    output_chunks: int = 10000,
    mem_persist: bool = False,
) -> xr.Dataset:
    """Select Persistent Scatterers (PS) from a SLC stack, and return a Space-Time Matrix.

    The selection method is defined by `method` and `threshold`.
    The selected pixels will be reshaped to (space, time), where `space` is the number of selected pixels.
    The unselected pixels will be discarded.
    The original `azimuth` and `range` coordinates will be persisted.
    The computed NAD or NMAD will be added to the output dataset as a new variable. It can be persisted in
    memory if `mem_persist` is True.

    Parameters
    ----------
    slcs : xr.Dataset
        nput SLC stack. It should have the following dimensions: ("azimuth", "range", "time").
        There should be a `amplitude` variable in the dataset.
    threshold : float
        Threshold value for selection.
    method : Literal["nad", "nmad"], optional
        Method of selection, by default "nad".
        - "nad": Normalized Amplitude Dispersion
        - "nmad": Normalized Median Amplitude Deviation
    output_chunks : int, optional
        Chunk size in the `space` dimension, by default 10000
    mem_persist : bool, optional
        If true persist the NAD or NMAD in memory, by default False.


    Returns
    -------
    xr.Dataset
        Selected PS, in form of an xarray.Dataset with two dimensions: (space, time).

    Raises
    ------
    NotImplementedError
        Raised when an unsupported method is provided.
    """
    # Make sure there is not temporal chunk
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


def _nad_block(amp: xr.DataArray) -> xr.DataArray:
    """Compute Normalized Amplitude Dispersion (NAD) for a block of amplitude data.

    Parameters
    ----------
    amp : xr.DataArray
        Amplitude data, with dimensions ("azimuth", "range", "time").
        This can be extracted from a SLC xr.Dataset.

    Returns
    -------
    xr.DataArray
        Normalized Amplitude Dispersion (NAD) data, with dimensions ("azimuth", "range").
    """
    # Compoute amplitude dispersion
    # By defalut, the mean and std function from Xarray will skip NaN values
    # However, if there is NaN value in time series, we want to discard the pixel
    # Therefore, we set skipna=False
    # Adding epsilon to avoid zero division
    nad_da = amp.std(dim="time", skipna=False) / (amp.mean(dim="time", skipna=False) + np.finfo(amp.dtype).eps)

    return nad_da


def _nmad_block(amp: xr.DataArray) -> xr.DataArray:
    """Compute Normalized Median Absolute Dispersion (NMAD) for a block of amplitude data.

    Parameters
    ----------
    amp : xr.DataArray
        Amplitude data, with dimensions ("azimuth", "range", "time").
        This can be extracted from a SLC xr.Dataset.

    Returns
    -------
    xr.DataArray
        Normalized Median Absolute Dispersion (NMAD) data, with dimensions ("azimuth", "range").
    """
    # Compoute NMAD
    median_amplitude = amp.median(dim="time", skipna=False)
    mad = (np.abs(amp - median_amplitude)).median(dim="time")  # Median Absolute Dispersion
    nmad = mad / (median_amplitude + np.finfo(amp.dtype).eps)  # Normalized Median Absolute Dispersion

    return nmad
