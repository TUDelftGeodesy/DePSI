"""Module for creating networks from STM points."""

import logging
import math
from typing import Literal

import networkx as nx
import numpy as np
import sparse
import xarray as xr
from scipy.spatial import Delaunay, KDTree

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

    arcs = xr.Dataset(
        data_vars={
            "d_phase": (["space", "time"], d_phase),
            "Btemp": (["time"], Btemp),
            "h2ph": (["space", "time"], h2ph),
        },
        coords={"source": (["space"], source_idx), "target": (["space"], target_idx)},
    )

    return arcs


def arc_selection(
    arcs: xr.Dataset,
    threshold: float,
    selection_method: Literal["ens_coh"] = "ens_coh",
    min_n_connection: int = 2,
) -> xr.Dataset:
    """Select aracs based on arc quality and connectivity.

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

    # Check if the network has more than one component
    # NetworkX is used. It should have good performance on large datasets.
    G = nx.Graph()
    G.add_edges_from(np.stack((arcs_selected["source"].data, arcs_selected["target"].data)).T)
    if nx.number_connected_components(G) > 1:
        logger.warning(
            "The network has more than one component. Currently, this is not supported by DePSI. "
            "Please adjust the network formation parameters, or decrease the threshold, "
            "to increase the connectivity of the network."
        )

    return arcs_selected


def remove_isolated_stm(stm: xr.Dataset, arcs: xr.Dataset) -> xr.Dataset:
    """Remove isolated points from the STM."""
    # Load source and target indices from arcs
    # these are 1d arrays so should fit in memory
    idx_source = arcs["source"].values
    idx_target = arcs["target"].values

    # Select STM points that are in arcs
    idx_selected = np.sort(np.unique(np.concatenate([idx_source, idx_target])))
    stm_updated = stm.isel(space=idx_selected)

    # The space size of the STM changes, resulting non-contiguous indices in space dimension
    # hence an update in arcs space coordinates is needed
    # Here we use a mapping solution, since the maximum number of network points is usually <100k
    # Map old indices in arcs to new indices
    idx_map = {old_idx: new_idx for new_idx, old_idx in enumerate(idx_selected)}
    # apply the mapping to the source and target indices in arcs
    arcs_updated = arcs.copy()
    arcs_updated["source"] = xr.DataArray(np.vectorize(idx_map.get)(arcs["source"].values), dims="space")
    arcs_updated["target"] = xr.DataArray(np.vectorize(idx_map.get)(arcs["target"].values), dims="space")

    return stm_updated, arcs_updated


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

    # Create a KDTree for fast distance queries.
    # Get all pairs of points within the specified radius.
    tree = KDTree(coordinates)
    if max_length is not None:
        pairs = tree.query_pairs(r=max_length, output_type="ndarray")
    else:
        pairs = tree.query_pairs(r=np.inf, output_type="ndarray")

    # Duplicate pairs with reversed indices to ensure that arcs are undirected.
    # This is necessary because KDTree returns pairs in one direction only.
    pairs = np.concatenate((pairs, np.flip(pairs, axis=1)), axis=0)
    pairs = pairs[np.argsort(pairs[:, 0])]  # Sort pairs by first column (source index).

    # Loop over all indices to collect neighbors.
    # For each index, we will collect the nearest min_links neighbors per partition.
    # The current node is separated into its own partition.
    for cur_index in indices:
        # List of arcs connected to the current node.
        cur_arcs = []

        neighbors = pairs[pairs[:, 0] == cur_index][:, 1].tolist()

        if len(neighbors) == 0:  # skip if there are no neighbors
            continue
        elif len(neighbors) <= min_links:
            # If there are not enough neighbors, we can just connect them all.
            for idx in neighbors:
                arc_to_add = (
                    min(cur_index, idx),
                    max(cur_index, idx),
                )
                cur_arcs.append(arc_to_add)
        else:
            partitions = [
                int(math.floor(num_partitions * (0.5 + math.atan2(coordinate[1], coordinate[0]) / math.tau)))
                for coordinate in coordinates[neighbors] - coordinates[cur_index]
            ]
            distances = [_get_distance(coordinates[cur_index], coordinates[idx]) for idx in neighbors]

            # Create a 3-column array with partition, distance, and index
            # Sort it by partition and then distance
            sorted_arr = np.array(sorted(list(zip(partitions, distances, neighbors, strict=False))))
            list_partitions = sorted_arr[:, 0].astype(int).tolist()  # Convert partitions to int for easier processing
            list_neighbors_idx = sorted_arr[:, 2].astype(int).tolist()  # Convert neighbor indices to int
            list_unique_partitions = np.unique(sorted_arr[:, 0]).astype(int).tolist()  # Unique partitions as integers

            neighbors_candidates = {}
            for (
                partition,
                idx,
            ) in zip(list_partitions, list_neighbors_idx, strict=False):
                if partition not in neighbors_candidates:
                    neighbors_candidates[partition] = []
                neighbors_candidates[partition].append(idx)

            # Loop over the unique partitions and collect the nearest neighbors of that partition.
            # repeat until we have at least min_links neighbors.
            count = 0
            while count < min_links:
                for partition in list_unique_partitions:
                    # make sure the source and target are in ascending order
                    # In each arc, make sure the source is less than the target.
                    # This is done to make the arcs undirected and canonical.

                    if len(neighbors_candidates[partition]) > 0:
                        arc_to_add = (
                            min(cur_index, neighbors_candidates[partition][0]),
                            max(cur_index, neighbors_candidates[partition][0]),
                        )
                        cur_arcs.append(arc_to_add)
                        count += 1

                        if count >= min_links:
                            break

                        # Remove the first element from the partition's neighbors.
                        neighbors_candidates[partition].pop(0)
                    else:
                        continue  # If there are no more neighbors in this partition, skip it

        # Add the current arcs to the list of all arcs.
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


def _network_relation_matirx(idx_source, idx_target, n_points):
    n_arcs = len(idx_source)
    A_sparse_start = sparse.COO(
        (np.arange(n_arcs), idx_source),
        np.full_like(np.arange(n_arcs), -1, dtype=np.int8),
        shape=(n_arcs, n_points),
    )
    A_sparse_end = sparse.COO(
        (np.arange(n_arcs), idx_target),
        np.full_like(np.arange(n_arcs), 1, dtype=np.int8),
        shape=(n_arcs, n_points),
    )
    A_sparse = A_sparse_start + A_sparse_end

    return A_sparse
