"""functions for densification."""

import numpy as np
import xarray as xr
from scipy.spatial import KDTree

from depsi.arc_estimation import periodogram
from depsi.utils import compute_phase_difference, concatenate_stms, generate_pnt_uids


def densification(
    stm_dens_pnts: xr.Dataset,
    stm_network_pnts: xr.Dataset,
    idx_refpnt: int | None = None,
    n_connections: int = 1,
    key_xcoord: str = "azimuth",
    key_ycoord: str = "range",
    key_h2ph: str = "h2ph_values",
    key_sdphase: str = "sd_phase",
    key_Btemporal: str = "years",
    **kwargs_arc_estimation,
) -> xr.Dataset:
    """Densification ambiguity estimation w.r.t. network points.

    Parameters
    ----------
    stm_dens_pnts : xarray.Dataset
        Dataset containing densification points.
    stm_network_pnts : xarray.Dataset
        Dataset containing network points.
    idx_refpnt : int or None, optional
        Index of the reference point in the network points (default is None).
        The phase of this point will be used to calculate the unwrapped phases of densification points.
        If None, this function will look for the "idx_refpnt" attribute in stm_network_pnts.
        If not found, an error will be raised.
    n_connections : int, optional
        Number of nearest neighbors to connect each densification point to (default is 1). Should be an odd number.
    key_xcoord : str, optional
        Name of the coordinate representing the x-axis (default is 'azimuth').
    key_ycoord : str, optional
        Name of the coordinate representing the y-axis (default is 'range').
    key_h2ph : str, optional
        Name of the data variable containing height-to-phase values (default is 'h2ph_values').
        Assuming the same key in densification and network points.
    key_sdphase : str, optional
        Name of the data variable containing phase data (default is 'sd_phase').
        Assuming the same key in densification and network points.
    key_Btemporal : str, optional
        Name of the data variable containing temporal baselines for arc estimation (default is 'years').
    **kwargs_arc_estimation :
        Additional keyword arguments to pass to the arc estimation function.

    Returns
    -------
    xarray.Dataset
        Densified points with estimated ambiguities, including both the original network points and
        the densification points.
    """
    # Make pnt_uid if not present
    if "pnt_uid" not in stm_dens_pnts:
        stm_dens_pnts = generate_pnt_uids(stm_dens_pnts)
    if "pnt_uid" not in stm_network_pnts:
        stm_network_pnts = generate_pnt_uids(stm_network_pnts)

    # Get wavelength from stm_network_pnts if not provided
    if "wavelength" not in stm_network_pnts.attrs:
        raise ValueError(
            "Wavelength is not provided and not found in the network STM attributes."
            "Please make sure it is provided as an attribute of stm_network_pnts."
            "For example: stm_network_pnts = stm_network_pnts.assign_attrs({'wavelength': wavelength})"
        )
    wavelength = stm_network_pnts.attrs["wavelength"]

    # Determine reference index
    if idx_refpnt is None:
        if "idx_refpnt" in stm_network_pnts.attrs:
            idx_refpnt = stm_network_pnts.attrs["idx_refpnt"]
        else:
            raise ValueError(
                "The attribute 'idx_refpnt' not found in stm_network_pnts.\n"
                "Please provide this attribute or set idx_refpnt input argument explicitly."
            )

    # Check number of connections
    if n_connections % 2 == 0:
        raise ValueError("n_connections should be an odd number.")
    if n_connections < 1:
        raise ValueError("n_connections should be at least 1.")
    if n_connections > 1:
        raise NotImplementedError("Currently only n_connections=1 is supported.")

    # Make sure stm_dens_pnts and stm_network_pnts have the same time coordinates
    if not np.array_equal(stm_dens_pnts["time"].data, stm_network_pnts["time"].data):
        raise ValueError("stm_dens_pnts and stm_network_pnts must have the same 'time' coordinates.")

    # Remove network points that are also in densification points
    mask = np.isin(stm_dens_pnts["pnt_uid"].data, stm_network_pnts["pnt_uid"].data)
    stm_dens_pnts = stm_dens_pnts.isel(space=np.where(~mask)[0])

    # Make a copy of densification points
    stm_dens_pnts_output = stm_dens_pnts.copy()

    # Query densification connections
    idx_dens_pnts, idx_network_pnts = _determine_dens_connections(
        stm_dens_pnts, stm_network_pnts, n_connections, key_xcoord, key_ycoord
    )

    # form densification arcs
    h2ph = (
        stm_dens_pnts[key_h2ph].isel(space=idx_dens_pnts).data
        + stm_network_pnts.isel(space=idx_network_pnts)[key_h2ph].data
    ) / 2  # take the mean for arc h2ph
    dd_phase = compute_phase_difference(
        stm_network_pnts.isel(space=idx_network_pnts)[key_sdphase].data,
        stm_dens_pnts[key_sdphase].isel(space=idx_dens_pnts).data,
        "subtract",
    )  # double difference phase, method should be subtract otherwise ambiguity check fails
    Btemp = stm_dens_pnts[key_Btemporal].data  # time baselines
    stm_densification_arcs = xr.Dataset(
        coords={
            "idx_dens": (("space",), idx_dens_pnts),
            "idx_network": (("space",), idx_network_pnts),
            "Btemp": (("time",), Btemp),
        },
        data_vars={
            "h2ph": (("space", "time"), h2ph),
            "dd_phase": (("space", "time"), dd_phase),
        },
        attrs={"wavelength": wavelength},
    )

    # unwrap densification arcs phase
    _, ambiguities, _, _, temporal_coh_arc = periodogram(
        stm_densification_arcs,
        key_dphase="dd_phase",
        key_h2ph="h2ph",
        key_Btemporal="Btemp",
        **kwargs_arc_estimation,
    )

    stm_densification_arcs["ambiguities"] = ambiguities

    if n_connections == 1:
        # When single connection, drop indexes of densification points
        # Use idx_network as the index to join with network points
        stm_densification_arcs = stm_densification_arcs.drop_vars("idx_dens")

    # Calculate estimated ambiguities for densification points
    # Because arc ambiguities = densification ambiguities - network ambiguities
    # => densification ambiguities = network ambiguities + arc ambiguities
    estimated_ambiguities = (
        stm_network_pnts.isel(space=stm_densification_arcs["idx_network"])["ambiguities"].data
        + stm_densification_arcs["ambiguities"].data
    )
    stm_dens_pnts_output["ambiguities"] = (("space", "time"), estimated_ambiguities)

    # Compute double-difference unwrapped phase
    ref_phase = stm_network_pnts[key_sdphase].isel(space=idx_refpnt).data
    unwrapped_phase_densification = (
        stm_dens_pnts_output[key_sdphase].data
        + estimated_ambiguities * 2 * np.pi
        - np.tile(ref_phase, (stm_dens_pnts_output.sizes["space"], 1))
    )
    stm_dens_pnts_output["unwrapped_phase"] = (("space", "time"), unwrapped_phase_densification)

    # Local temporal coherence with one arc connection is the same as arc temporal coherence
    stm_dens_pnts_output["local_temp_coh"] = (("space",), temporal_coh_arc.data)

    # Attach network points to the output
    # Concat in space dimension
    # For variables with (space,) or (space, time) dimensions, missing values will be filled with NaNs
    # For time-only variables, the values will be directly copied from stm_network_pnts or stm_dens_pnts_output
    # If a time-only variable exists in both STMs, it is assumed to be identical across those STMs.
    stm_dens_pnts_output = concatenate_stms([stm_network_pnts, stm_dens_pnts_output])

    return stm_dens_pnts_output


def _determine_dens_connections(
    stm_dens_pnts, stm_network_pnts, n_connections, key_xcoord="azimuth", key_ycoord="range"
) -> tuple[np.ndarray, np.ndarray]:
    """Determine densification connections between densification points and network points using KDTree.

    Parameters
    ----------
    stm_dens_pnts : xarray.Dataset
        Dataset containing densification points.
    stm_network_pnts : xarray.Dataset
        Dataset containing network points.
    n_connections : int
        Number of nearest neighbors to connect each densification point to.
    key_xcoord : str, optional
        Name of the coordinate representing the x-axis (default is 'azimuth').
    key_ycoord : str, optional
        Name of the coordinate representing the y-axis (default is 'range').

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        Tuple containing arrays of indices for densification points and network points.
    """
    coords_network = np.stack((stm_network_pnts[key_xcoord].data, stm_network_pnts[key_ycoord].data), axis=-1)
    coords_densification = np.stack((stm_dens_pnts[key_xcoord].data, stm_dens_pnts[key_ycoord].data), axis=-1)
    tree = KDTree(coords_network)

    distances, indices = tree.query(coords_densification, k=1)

    if indices.ndim == 1:
        indices = indices[:, np.newaxis]  # Make it 2D for uniformity

    # Build densification arcs from indices
    idx_dens_pnts = np.repeat(np.arange(stm_dens_pnts.sizes["space"]), indices.shape[1])  # Point to densify
    idx_network_pnts = indices.flatten()  # Network point

    return idx_dens_pnts, idx_network_pnts
