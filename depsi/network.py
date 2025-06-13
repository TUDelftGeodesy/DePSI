"""Module for creating networks from STM points."""

import logging
import math
from typing import Literal

import numpy as np
import xarray as xr
from scipy.spatial import Delaunay

logger = logging.getLogger(__name__)


def form_network(
    stm: xr.Dataset,
    key_phase: str,
    key_h2ph: str,
    key_Btemp: str,
    key_complex: str = "complex",
    key_xlabel: str = "lon",
    key_ylabel: str = "lat",
    network_method: Literal["redundant", "delaunay"] = "redundant",
    max_length: float = None,
    min_links: int = 12,
    num_partitions: int = 8,
    dphase_method: Literal["conjmult", "subtract"] = "subtract",
) -> xr.Dataset:
    """Generate an STM of arcs from an STM of points.

    Parameters
    ----------
    stm : xr.Dataset
        Space-Time Matrix of scatteres.
    key_phase : str
        Key of the phase values in the STM.
        This phase will be used to compute the differential arc phase.
    key_h2ph : str
        Key of the h2ph values in the STM.
        The arc h2ph will be computed as the average between source and target.
    key_Btemp : str
        Key of the Btemp values in the STM.
    key_complex : str, optional
        Key of the complex values, by default "complex"
    key_xlabel : str, optional
        Key of the x coordinates for calulating arc length, by default "lon"
    key_ylabel : str, optional
        Key of the y coordinates for calulating arc length, by default "lat"
    network_method : Literal["redundant", "delaunay"], optional
        network formation method, by default "redundant"
    max_length : float, optional
        maximum arc length, by default None
    min_links : int, optional
        minimum links per point, by default 12
        only effective when network_method is "redundant"
    num_partitions : int, optional
        number of partitions of searching when forming redundant network, by default 8
        only effective when network_method is "redundant"
    dphase_method : Literal["conjmult", "subtract"], optional
        method of computing phase difference, by default "subtract"
        "subtract" method subtracts the source phase from the target phase;
        "conjmult" method computes the phase difference by conjugate multiplication:
            d_phase = np.angle(complex_target * complex_source.conj())

    Returns
    -------
    xr.Dataset
        Space-Time Matrix of arcs, containing the following variables:
        - d_phase: the arc phase, which is the difference between source and target points
        - Btemp: the temporal baseline, which is the same for all arcs
        - h2ph: the arc h2ph, which is the average between source and target points
    """
    # Generate the network arcs.
    if network_method == "redundant":
        if min_links <= 0:
            logger.error(f"min_links must be strictly positive (currently: {min_links})")
            return
        if num_partitions <= 0:
            logger.error(f"num_partitions must be strictly positive (currently: {num_partitions})")
            return
    elif network_method != "delaunay":
        raise NotImplementedError(f"Unknown network method {network_method}, known are delaunay and redundant")

    # Collect point coordinates.
    indices = [stm[coord] for coord in [key_xlabel, key_ylabel]]
    coordinates = np.column_stack(indices)

    arcs = None

    # Create network arcs as list of tuples of point ids.
    if network_method == "delaunay":
        arcs = _generate_arcs_delaunay(coordinates, max_length)
    elif network_method == "redundant":
        arcs = _generate_arcs_redundant(coordinates, max_length, min_links, num_partitions)

    # Compute the phase difference.
    arcs_unzipped = list(zip(*arcs, strict=False))
    source_idx = list(arcs_unzipped[0])
    target_idx = list(arcs_unzipped[1])
    d_phase = _compute_phase_difference(stm, source_idx, target_idx, key_phase, key_complex, method=dphase_method)

    Btemp = stm[key_Btemp].data

    h2ph = (stm[key_h2ph].isel(space=source_idx).data + stm[key_h2ph].isel(space=target_idx).data) / 2

    stm_arcs = xr.Dataset(
        data_vars={
            "d_phase": (["space", "time"], d_phase),
            "Btemp": (["time"], Btemp),
            "h2ph": (["space", "time"], h2ph),
        },
        coords={"source": (["space"], source_idx), "target": (["space"], target_idx)},
    )

    return stm_arcs


def arc_selection(
    arcs: xr.Dataset,
    threshold: float,
    selection_method: Literal["ens_coh"] = "ens_coh",
    min_n_connection: int = 2,
) -> xr.Dataset:
    """Select aracs based on creteria.

    This function selects arcs in two steps:
    1. It selects arcs based on a threshold value(e.g., ens_coh).
    2. It removes arcs connected to points which have less than a minimum number of connections.

    Parameters
    ----------
    arcs : xr.Dataset
        arcs to select from, in space-time matrix
    threshold : float
        threshold value for selection
    selection_method : Literal["ens_coh"]
        values to use for selection, by default "ens_coh". The available options are:
        - "ens_coh": ensemble coherence, arcs with ens_coh > threshold are selected.
          assumes that arcs have a variable "ens_coh" in the dataset.
    min_n_connection : int, optional
        minimum number of connections, by default 2

    Returns
    -------
    xr.Dataset
        selected arcs in space-time matrix
    """
    # Threshold selection
    match selection_method:
        case "ens_coh":
            mask = np.abs(arcs["ens_coh"]) > threshold  # mask as DataArray
            arcs_selected = arcs.where(mask, drop=True)
        case _:
            raise NotImplementedError

    # Remove arcs which can not be tested
    # These arcs are identified by the points which has <= min_n_connection arcs connected to them
    # All arcs connected to such points are removed
    # An iterative approach is used to remove all arcs connected to such points
    point_ids_all = np.concat(
        [arcs_selected["source"].data, arcs_selected["target"].data]
    )  # all occurrances of point ids
    point_ids_unique, counts = np.unique(point_ids_all, return_counts=True)  # unique point ids and their counts

    while np.any(counts <= min_n_connection):
        # Find points with <=3 arcs connected
        point_ids_to_remove = point_ids_unique[counts <= min_n_connection]
        # Create a mask for arcs to remove
        mask_remove = np.isin(arcs_selected["source"].data, point_ids_to_remove) | np.isin(
            arcs_selected["target"].data, point_ids_to_remove
        )
        idx_select = np.where(~mask_remove)[0]  # indices of arcs to remove
        # Remove these arcs
        arcs_selected = arcs_selected.isel(space=idx_select)

        # Update point ids and counts
        point_ids_all = np.concat([arcs_selected["source"].data, arcs_selected["target"].data])
        point_ids_unique, counts = np.unique(point_ids_all, return_counts=True)

    return arcs_selected


def _get_distance(s, t):
    """Calculate the distance between two points.

    Args:
    ----
        s: the source point.
        t: the target point.

    Returns:
    -------
        The distance between the two points.
    """
    # TODO(tvl) More complex distance functions can be implemented here.
    #
    # For example, network generation gives better results for square Euclidean distances.
    # Non-square coordinate systems (like non-square image coordinates) may distort the circle properties of Delaunay
    # networks into ellipses.
    # Non-Euclidean coordinate systems (like angular lat-lon systems) may distort these same properties depending on
    # the distance to a pole.
    #
    # Coordinate system transformations may be done on the STM before generating the network.
    # However, there may be cases where it is impossible or undesirable to transform the point coordinates in the STM.
    # In such a case, coordinates may be transformed inside this function.

    return math.dist(s, t)


def _generate_arcs_delaunay(coordinates, max_length=None):
    """Create a network using Delaunay triangulation."""
    # Create network and collect neighbors.
    network = Delaunay(coordinates)
    neighbors_ptr, neighbors_idx = network.vertex_neighbor_vertices

    # Convert ptr and idx arrays into list of sorted index pairs.
    arcs = []
    for s in range(len(neighbors_ptr) - 1):
        for t in range(neighbors_ptr[s], neighbors_ptr[s + 1]):
            length = _get_distance(coordinates[int(s)], coordinates[neighbors_idx[t]])
            if max_length is None or length <= max_length:
                arcs.append(tuple(sorted([int(s), int(neighbors_idx[t])])))

    # Remove duplicates and make the list canonical.
    arcs = sorted(list(set(arcs)))

    return arcs


def _generate_arcs_redundant(coordinates, max_length=None, min_links=12, num_partitions=8):
    """Create a network with at least min_links arcs per node.

    Arcs are created ordered by length.
    However, the orientations around the node are split into num_partitions partitions;
    each partition can only get an (x+1)th arc if every other partition either
    already has x arcs connected or already has all allowed arcs connected
    (e.g. there are no more nodes in that partition, or they are all too far away).
    Note that a node may get less than min_links arcs if there are not enough neighbors within max_length.
    """
    arcs = []
    indices = range(len(coordinates))
    for cur_index in indices:
        # Calculate per node, the partition they are in and the distance from the current node.
        partitions = [
            int(math.floor(num_partitions * (0.5 + math.atan2(coordinate[1], coordinate[0]) / math.tau)))
            for coordinate in coordinates - coordinates[cur_index]
        ]
        partitions[cur_index] = num_partitions + 1  # Separate the current node into its own partition.
        distances = [_get_distance(coordinates[cur_index], coordinate) for coordinate in coordinates]

        # Create a list of tuples with the partition, distance, and index, sorted by partition and then distance.
        values = np.array(sorted(list(zip(partitions, distances, indices, strict=False))))

        # Collect the nearest min_links neighbors per partition and discard the partition of the current node.
        partitions_diff = values[1:, 0] - values[:-1, 0]
        separators = np.where(partitions_diff > 0)[0]
        partitions = np.split(values, separators + 1)
        partitions = partitions[: len(partitions) - 1]
        partitions = [partition[:min_links] for partition in partitions]

        # Collect the neighbor 'hierarchies'.
        # Each hierarchy contains the nth nearest neighbor from all partitions.
        neighbor_hierarchies = [[] for _ in range(min_links)]
        count = 0
        for n in range(min_links):
            # Break early if we have gathered enough neighbors.
            if min_links <= count:
                break
            for partition in partitions:
                # Note that we do not break inside this loop,
                # because we want the nth nearest neighbors from all partitions.
                if n < len(partition) and (max_length is None or partition[n][1] <= max_length):
                    neighbor_hierarchies[n].append(partition[n])
                    count = count + 1

        # Sort hierarchies per partition by distance to the current node.
        neighbor_hierarchies = [
            sorted(hierarchy, key=lambda x: x[1]) for hierarchy in neighbor_hierarchies if len(hierarchy) != 0
        ]

        # Add sorted arcs to at least min_links neighbors.
        cur_arcs = [
            tuple(sorted([cur_index, int(neighbor[2])])) for hierarchy in neighbor_hierarchies for neighbor in hierarchy
        ]
        cur_arcs = cur_arcs[:min_links]

        arcs.extend(cur_arcs)

    # Remove duplicates and make the list canonical.
    arcs = sorted(list(set(arcs)))

    return arcs


def _compute_phase_difference(
    stm,
    source_idx,
    target_idx,
    key_phase: str,
    key_complex: str,
    method: Literal["subtract", "conjmult"],
) -> np.ndarray:
    """Calculate the phase difference between two points.

    The method can either be "subtract" or "conjmult".
    """
    if method == "subtract":
        d_phase = stm[key_phase].isel(space=target_idx).data - stm[key_phase].isel(space=source_idx).data
    elif method == "conjmult":
        complex_source = stm[key_complex].isel(space=source_idx).data
        complex_target = stm[key_complex].isel(space=target_idx).data
        d_phase = np.angle(complex_target * complex_source.conj())
    else:
        raise NotImplementedError(f"Unknown difference method {method}, known are subtract and conjmult")
    return d_phase
