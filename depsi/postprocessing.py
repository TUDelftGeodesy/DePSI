import numpy as np
import xarray as xr


def stm_point_filter(
    stm: xr.Dataset,
    layer_to_filter: str,
    vmin: float | int | None = None,
    vmax: float | int | None = None,
    return_removed: bool = False,
) -> xr.Dataset | tuple[xr.Dataset, xr.Dataset]:
    """Filter the points in an STM based on a given layer name and an allowed value range.

    To filter with just an upper or lower limit, leave `vmin` to `None` (just upper limit) or `vmax` to `None` (just
    lower limit). Leaving both `vmin` and `vmax` to `None` will return the original STM.

    Parameters
    ----------
    stm: xr.Dataset
        STM with at least the layer `layer_to_filter` and dimension `space`
    layer_to_filter: str
        Name of the layer on which the filter should be applied
    vmin: float | int | None
        Minimum allowed value. Values equal to `vmin` are included. No lower bound is used if set to `None`.
    vmax: float | int | None
        Maximum allowed value. Values equal to `vmax` are included. No upper bound is used if set to `None`.
    return_removed: bool, default False
        If True, a second STM with the removed points is returned as second output alongside the STM with the retained
        points

    Returns
    -------
    xr.Dataset (if `return_removed` is `False`)
        STM containing all points where `layer_to_filter` is within the given bounds
    tuple[xr.Dataset, xr.Dataset] (if `return_removed` is `True`)
        STM containing all points where `layer_to_filter` is within the given bounds
        STM containing all points where `layer_to_filter` is outside the given bounds
    """
    if vmin is not None:
        if vmax is not None:
            mask = (stm[layer_to_filter] >= vmin) & (stm[layer_to_filter] <= vmax)
        else:
            mask = stm[layer_to_filter] >= vmin
        mask = mask.values
    else:
        if vmax is not None:
            mask = stm[layer_to_filter] <= vmax
            mask = mask.values
        else:
            mask = np.array([True] * stm.sizes["space"])  # no filter applied, so return everything

    retained_stm = stm.sel(space=stm["space"].values[mask])

    if return_removed:
        removed_stm = stm.sel(space=stm["space"].values[~mask])
        return retained_stm, removed_stm

    return retained_stm
