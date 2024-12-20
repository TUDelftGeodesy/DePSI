"""Functions for scatter selection related operations."""

from typing import Literal

import numpy as np
import xarray as xr


def ps_selection(
    slcs: xr.Dataset, threshold: float, method: Literal["nad", "nmad"] = "nad", output_chunks: int = 10000
) -> xr.Dataset:
    """Select Persistent Scatterers (PS) from a SLC stack, and return a Space-Time Matrix.

    The selection method is defined by `method` and `threshold`.
    The selected pixels will be reshaped to (space, time), where `space` is the number of selected pixels.
    The unselected pixels will be discarded.
    The original `azimuth` and `range` coordinates will be persisted.

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

    Returns
    -------
    xr.Dataset
        Selected PS, in form of an xarray.Dataset with two dimensions: (space, time).

    Raises
    ------
    NotImplementedError
        Raised when an unsupported method is provided.
    """
    match method:
        case "nad":
            nad = xr.map_blocks(
                _nad_block, slcs["amplitude"], template=slcs["amplitude"].isel(time=0).drop_vars("time")
            )
            mask = nad < threshold
        case _:
            raise NotImplementedError

    # Get the 1D index on space dimension
    mask_1d = mask.stack(space=("azimuth", "range")).drop_vars(["azimuth", "range"])
    index = mask_1d["space"].where(mask_1d.compute(), other=0, drop=True)  # Evaluate the 1D mask to index

    # Reshape from Stack ("azimuth", "range", "time") to Space-Time Matrix  ("space", "time")
    stacked = slcs.stack(space=("azimuth", "range"))

    # Drop multi-index coords for space coordinates
    # This will also azimuth and range coordinates, as they are part of the multi-index coordinates
    stm = stacked.drop_vars(["space"])

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

    return stm_masked


def _nad_block(amp: xr.DataArray) -> xr.DataArray:
    # Time dimension order
    t_order = list(amp.dims).index("time")

    # Compoute amplitude dispersion
    # By defalut, the mean and std function from Xarray will skip NaN values
    # However, if there is NaN value in time series, we want to discard the pixel
    # Therefore, we set skipna=False
    # Adding epsilon to avoid zero division
    nad_da = amp.std(axis=t_order, skipna=False) / (amp.mean(axis=t_order, skipna=False) + np.finfo(amp.dtype).eps)

    return nad_da
