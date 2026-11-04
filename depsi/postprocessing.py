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

def estimate_reference_point_noise(
    stm: xr.Dataset,
    key_phase: str="phase_residuals",
    method: str = "mean",
    space_dim: str = "space",
) -> xr.DataArray:
    """Estimate the reference point noise (RPN) per epoch.

    All points are relative to a common reference point, so the noise of the
    reference point appears as the same offset in every point at each epoch.
    Following Chang and Hanssen (2016), Eq. 15, this offset is estimated as the
    spatial mean (or median) of the residual phases per epoch.

    Parameters
    ----------
    stm: xr.Dataset
        The input dataset containing the residual phases.
    key_phase: str
        The key of the residual phase data in the dataset. Default is "residual_phase".
    method: str
        Statistic used to estimate the RPN per epoch: "mean" or "median".
        Default is "mean". anything else will be treated as "median".

    Returns
    -------
    xr.DataArray
        The estimated reference point noise, one value per epoch.

    References
    ----------
    Chang, L. and Hanssen, R. F. (2016). A probabilistic approach for InSAR
    time-series postprocessing. IEEE Transactions on Geoscience and Remote
    Sensing, 54(1), 421-430.
    """
    stm_rpn = stm.copy()  # Avoid modifying the original dataset

    residual_phase = stm_rpn[key_phase]

    if method == "mean":
        reference_point_noise = residual_phase.mean(dim=space_dim, skipna=True)
    else:
        # A median needs all points of an epoch in one chunk when using dask
        reference_point_noise = residual_phase.median(dim=space_dim, skipna=True)

    stm_rpn["reference_point_noise"] = reference_point_noise

    return stm_rpn