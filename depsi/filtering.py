import dask.array as da
import xarray as xr


def stm_point_filter(
    stm: xr.Dataset, layer_to_filter: str, filter_bounds: tuple, return_removed: bool = False
) -> xr.Dataset | tuple[xr.Dataset, xr.Dataset]:
    """Filter the points in an STM based on a given layer name and an allowed value range.

    To filter with just an upper or lower limit, use `filter_bounds=(None, upper value)` (just upper limit) or
    `filter_bounds=(lower value, None)` (just lower limit).

    Parameters
    ----------
    stm: xr.Dataset
        STM with at least the layer `layer_to_filter` and dimension `space`
    layer_to_filter: str
        Name of the layer on which the filter should be applied
    filter_bounds: tuple
        (minimum allowed value, maximum allowed value). Values equal to the given value are included. If either bound
        is set to `None`, no lower or upper bound is used.
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

    Raises
    ------
    AssertionError
        - if length of `filter_bounds` is not 2
    """
    assert len(filter_bounds) == 2, f"Expected filter_bounds=(lower, upper), got {filter_bounds} instead."

    if filter_bounds[0] is not None:
        if filter_bounds[1] is not None:
            mask = (stm[layer_to_filter] >= filter_bounds[0]) & (stm[layer_to_filter] <= filter_bounds[1])
        else:
            mask = stm[layer_to_filter] >= filter_bounds[0]
    else:
        if filter_bounds[1] is not None:
            mask = stm[layer_to_filter] <= filter_bounds[1]
        else:
            mask = da.array([True] * len(list(stm["space"].values)))  # no filter applied, so return everything

    mask = mask.values  # compute it

    retained_stm = stm.sel(space=stm["space"].values[mask])

    if return_removed:
        removed_stm = stm.sel(space=stm["space"].values[~mask])
        return retained_stm, removed_stm

    return retained_stm
