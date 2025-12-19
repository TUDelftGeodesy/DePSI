"""functions for densification."""

import numpy as np
import xarray as xr
from scipy.spatial import KDTree

from depsi.arc_estimation import periodogram
from depsi.utils import compute_phase_difference, generate_pnt_uids


def densification(
    stm_densification: xr.Dataset,
    stm_network_pnts: xr.Dataset,
    wavelength: float,
    n_connections: int = 1,
    key_xcoord: str = "azimuth",
    key_ycoord: str = "range",
    key_h2ph: str = "h2ph_values",
    key_sdphase: str = "sd_phase",
    key_btemp: str = "years",
    phase_diff_method: str = "subtract",
    **kwargs_arc_estimation,
) -> xr.Dataset:
    """Densification ambiguity estimation w.r.t. network points.

    Parameters
    ----------
    stm_densification : xarray.Dataset
        Dataset containing densification points.
    stm_network_pnts : xarray.Dataset
        Dataset containing network points.
    wavelength : float
        Wavelength used for phase unwrapping and ambiguity estimation.
    n_connections : int, optional
        Number of nearest neighbors to connect each densification point to (default is 1).
    key_xcoord : str, optional
        Name of the coordinate representing the x-axis (default is 'azimuth').
    key_ycoord : str, optional
        Name of the coordinate representing the y-axis (default is 'range').
    key_h2ph : str, optional
        Name of the data variable containing height-to-phase values (default is 'h2ph_values').
    key_sdphase : str, optional
        Name of the data variable containing phase data (default is 'sd_phase').
    key_btemp : str, optional
        Name of the data variable containing time baselines for arc estimation (default is 'years').
    phase_diff_method : str, optional
        Method for computing phase differences (default is 'subtract').
    **kwargs_arc_estimation :
        Additional keyword arguments to pass to the arc estimation function.

    Returns
    -------
    xarray.Dataset
        Densified points with estimated ambiguities, including both network points and
        densification points with computed ambiguity estimates.
    """
    # Make pnt_uid if not present
    if "pnt_uid" not in stm_densification:
        stm_densification = generate_pnt_uids(stm_densification)
    if "pnt_uid" not in stm_network_pnts:
        stm_network_pnts = generate_pnt_uids(stm_network_pnts)

    # Check number of connections
    if n_connections > 1:
        raise NotImplementedError("Currently only n_connections=1 is supported.")

    # Make sure stm_densification and stm_network_pnts have the same time coordinates
    if not np.array_equal(stm_densification["time"].values, stm_network_pnts["time"].values):
        raise ValueError("stm_densification and stm_network_pnts must have the same 'time' coordinates.")

    # Remove network points that are also in densification points
    mask = np.isin(stm_densification["pnt_uid"].values, stm_network_pnts["pnt_uid"].values)
    stm_densification = stm_densification.isel(space=np.where(~mask)[0])

    # Make a copy of densification points
    stm_densification_output = stm_densification.copy()

    # Query densification connections
    idx_dens_pnts, idx_network_pnts = _query_dens_connections(
        stm_densification, stm_network_pnts, n_connections, key_xcoord, key_ycoord
    )

    # form densification arcs
    h2ph = (
        stm_densification[key_h2ph].isel(space=idx_dens_pnts).values
        + stm_network_pnts.isel(space=idx_network_pnts)[key_h2ph].values
    ) / 2  # take the mean for arc h2ph
    dd_phase = compute_phase_difference(
        stm_densification[key_sdphase].isel(space=idx_dens_pnts).values,
        stm_network_pnts.isel(space=idx_network_pnts)[key_sdphase].values,
        phase_diff_method,
    )  # double difference phase
    Btemp = stm_densification[key_btemp].values  # time baselines
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
    )

    # unwrap densification arcs phase
    _, ambiguities, _, _, _ = periodogram(
        stm_densification_arcs,
        key_dphase="dd_phase",
        key_h2ph="h2ph",
        key_Btemp="Btemp",
        wavelength=wavelength,
        **kwargs_arc_estimation,
    )

    stm_densification_arcs["ambiguities"] = ambiguities

    if n_connections == 1:
        stm_densification_arcs = stm_densification_arcs.drop_vars("idx_dens").rename_vars(
            {"idx_network": "idx_network"}
        )

    # Calculate estimated ambiguities for densification points
    # Because arc ambiguities = network ambiguities - densification ambiguities
    # => densification ambiguities = network ambiguities - arc ambiguities
    estimated_ambiguities = (
        stm_network_pnts.isel(space=stm_densification_arcs["idx_network"])["ambiguities"].values
        - stm_densification_arcs["ambiguities"].values
    )
    stm_densification_output["ambiguities"] = (("space", "time"), estimated_ambiguities)

    # Attach network points to the output
    # Join in space dimension, keep all data variables
    stm_densification_output = xr.concat([stm_network_pnts, stm_densification_output], dim="space", data_vars="all")

    return stm_densification_output


def _query_dens_connections(
    stm_densification, stm_network_pnts, n_connections, key_xcoord="azimuth", key_ycoord="range"
) -> tuple[np.ndarray, np.ndarray]:
    """Query densification connections between densification points and network points using KDTree.

    Parameters
    ----------
    stm_densification : xarray.Dataset
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
    coords_network = np.stack((stm_network_pnts[key_xcoord].values, stm_network_pnts[key_ycoord].values), axis=-1)
    coords_densification = np.stack(
        (stm_densification[key_xcoord].values, stm_densification[key_ycoord].values), axis=-1
    )
    tree = KDTree(coords_network)

    distances, indices = tree.query(coords_densification, k=1)

    if indices.ndim == 1:
        indices = indices[:, np.newaxis]  # Make it 2D for uniformity

    # Build densification arcs from indices
    idx_dens_pnts = np.repeat(np.arange(stm_densification.sizes["space"]), indices.shape[1])  # Point to densify
    idx_network_pnts = indices.flatten()  # Network point

    return idx_dens_pnts, idx_network_pnts
