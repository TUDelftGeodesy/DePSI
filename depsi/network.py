"""Module for creating networks from STM points."""

import logging
import math
from typing import Literal

import networkx as nx
import numpy as np
import scipy
import sparse
import xarray as xr
from scipy.spatial import Delaunay, KDTree

logger = logging.getLogger(__name__)

# Constants for MHT in network integration
ALPHA0 = 0.1  # Significance level for 1-dimensional test
GAMMA0 = 0.5  # Power of the test


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
    n_links: int = 12,
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
    n_links : int, optional
        target number of links per point, by default 12
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
        if n_links <= 0:
            logger.error(f"n_links must be strictly positive (currently: {n_links})")
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
        arcs = _generate_arcs_redundant(coordinates, max_length, n_links, num_partitions)

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


def _mht_network_adjustment(
    A: scipy.sparse._csr.csr_matrix,
    y: np.ndarray,
    Qyy: np.ndarray,
    k1: float,
    kb_dict: dict,
) -> (int, int):
    """Remove one point/arc from the network to reduce the residual in ambiguity estimation."""
    # Retrive shapes
    N_arcs, N_epochs = y.shape
    N_points = A.shape[1]

    # Inverse of VCM of observations
    if Qyy.ndim == 1:  # Diagonal VCM
        invQy = scipy.sparse.diags(1 / Qyy, 0, shape=(N_arcs, N_arcs))
    elif Qyy.ndim == 2:  # Full VCM
        invQy = np.linalg.inv(Qyy)
    else:
        raise ValueError("Qyy must be either 1D (diagonal VCM) or 2D (full VCM)")

    # Solve ambiguities as float
    _, echeck = _solve_float_ambiguities(A, y, invQy)

    # Post-priori VCM of residuals
    # Qecheck = Qyy - Qycheck = Qyy - A Qxx A'
    Qxx = np.linalg.inv((A.T @ invQy @ A).todense())
    Qecheck = Qyy - (A @ Qxx @ A.T)  # TODO: check how to handle large Qecheck

    # Test statistics TT1 for removing one arc
    Qecheck_diag = np.array(Qecheck.diagonal().flatten()).squeeze()
    w = echeck**2 / np.tile(np.abs(Qecheck_diag), (N_epochs, 1)).T
    TT1 = np.sum(w, axis=1) / k1**2

    # Test statistics for removing one point
    TTq = np.zeros(N_points)
    for pnt_idx in range(N_points):
        arcs_idx = np.where(A[:, pnt_idx].todense() != 0)[0]  # Arcs connected to this point
        arcs_idx = arcs_idx[1:]  # Drop one arc to create basis, see e.g. verhoef97
        echeck_point = echeck[arcs_idx, :]  # Relevant echeck of this point
        Qecheck_point = Qecheck[arcs_idx, :][:, arcs_idx]  # Relevant Qecheck_diag of this point

        # Compute the test statistic for this point
        # TODO: check if this abs is taken correctly
        Tq = np.sum(
            np.abs((echeck_point.T @ np.linalg.inv(Qecheck_point) @ echeck_point).diagonal())
        )  # Before adjust for degree of freedom
        TTq[pnt_idx] = Tq / kb_dict[len(arcs_idx)]

    # Decision one removal strategy
    if max(TT1) > max(TTq):
        idx_removal = np.argmax(TT1)  # index of arc to remove
        flag_removal = 0  # remove arc
    else:
        idx_removal = np.argmax(TTq)  # index of point to remove
        flag_removal = 1  # remove point

    return flag_removal, idx_removal


def _solve_float_ambiguities(A, y, invQy):
    """Solve ambiguities as a float based on Least-Squares."""
    # Solve ambiguities as they are float numbers
    # This solves the equation Ax = y in least square sense
    # With A a sparse matrix
    # And stochastic model Qyy taken into account
    invQyA = invQy @ A  # Avoid repeated computation in vectorized lsmr

    @np.vectorize(signature="(i)->(j)")
    def lsmr(y):
        """Least square iterative solver for sparse data."""
        x, *_ = scipy.sparse.linalg.lsmr(invQyA, y)
        return x

    acheck = lsmr(y.T).T  # float ambiguity estimation
    echeck = y - A @ acheck  # residuals estimation

    return acheck, echeck


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


def remove_isolated_points(stm: xr.Dataset, arcs: xr.Dataset) -> xr.Dataset:
    """Remove isolated points from the network."""
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


def _generate_arcs_delaunay(coordinates, max_length=None):
    """Create a network using Delaunay triangulation."""
    # Create network and collect neighbors.
    network = Delaunay(coordinates)
    neighbors_ptr, neighbors_idx = network.vertex_neighbor_vertices

    # Convert ptr and idx arrays into list of sorted index pairs.
    arcs = []
    for s in range(len(neighbors_ptr) - 1):
        for t in range(neighbors_ptr[s], neighbors_ptr[s + 1]):
            length = math.dist(coordinates[int(s)], coordinates[neighbors_idx[t]])
            if max_length is None or length <= max_length:
                arcs.append(tuple(sorted([int(s), int(neighbors_idx[t])])))

    # Remove duplicates and make the list canonical.
    arcs = sorted(list(set(arcs)))

    return arcs


def _generate_arcs_redundant(coordinates, max_length=None, n_links=12, num_partitions=8):
    """Create a redundant network.

    The redundant network is formed with the following steps:

    1. Create a KDTree and find all pairs of points within the maximum distance.
    2. Loop through each point and find its neighbors within the maximum distance.
    3. Divide neighbors into partitions based on their direction.
    4. Select the nth nearest neighbors from all partitions, starting from n=1.
    5. Sort the selected neighbors by distance, add them to the arcs list. If n_links is not
       exceeded, continue to the n+1th nearest neighbors of all partitions.
    6. Repeat until n_links is reached.
    """
    arcs = []
    indices = range(len(coordinates))

    # Create a KDTree for fast distance queries.
    tree = KDTree(coordinates)
    if max_length is not None:
        pairs = tree.query_pairs(r=max_length, output_type="ndarray")
    else:
        pairs = tree.query_pairs(r=np.inf, output_type="ndarray")

    # Duplicate pairs with reversed indices to ensure that arcs are undirected.
    pairs = np.concatenate((pairs, np.flip(pairs, axis=1)), axis=0)
    pairs = pairs[np.argsort(pairs[:, 0])]  # Sort pairs by first column (source index).

    for cur_index in indices:
        # Get the neighbors of the current node
        neighbors = pairs[pairs[:, 0] == cur_index][:, 1].tolist()

        if len(neighbors) == 0:  # skip if there are no neighbors
            continue
        elif len(neighbors) <= n_links:
            # If there are not enough neighbors, connect them all.
            for idx in neighbors:
                arc_to_add = (min(cur_index, idx), max(cur_index, idx))
                arcs.append(arc_to_add)
        else:
            # Calculate partitions and distances for neighbors
            partitions = [
                int(math.floor(num_partitions * (0.5 + math.atan2(coordinate[1], coordinate[0]) / math.tau)))
                for coordinate in coordinates[neighbors] - coordinates[cur_index]
            ]
            distances = [math.dist(coordinates[cur_index], coordinates[idx]) for idx in neighbors]

            # Create sorted array by partition and then distance
            sorted_arr = np.array(sorted(list(zip(partitions, distances, neighbors, strict=False))))

            # Split into partitions
            partitions_diff = sorted_arr[1:, 0] - sorted_arr[:-1, 0]
            separators = np.where(partitions_diff > 0)[0]
            partitions_split = np.split(sorted_arr, separators + 1)
            partitions_split = [partition[:n_links] for partition in partitions_split]

            # Collect the neighbor 'hierarchies'
            neighbor_hierarchies = [[] for _ in range(n_links)]
            count = 0
            for n in range(n_links):
                # Break early if we have gathered enough neighbors.
                if n_links <= count:
                    break
                for partition in partitions_split:
                    # Note that we do not break inside this loop,
                    # because we want the nth nearest neighbors from all partitions.
                    if n < len(partition) and (max_length is None or partition[n][1] <= max_length):
                        neighbor_hierarchies[n].append(partition[n])
                        count = count + 1

            # Sort hierarchies per partition by distance to the current node
            neighbor_hierarchies = [
                sorted(hierarchy, key=lambda x: x[1]) for hierarchy in neighbor_hierarchies if len(hierarchy) != 0
            ]

            # Add sorted arcs to at least n_links neighbors
            cur_arcs = [
                (min(cur_index, int(neighbor[2])), max(cur_index, int(neighbor[2])))
                for hierarchy in neighbor_hierarchies
                for neighbor in hierarchy
            ]
            cur_arcs = cur_arcs[:n_links]

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


def _network_relation_matrix(idx_source, idx_target, n_points):
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

    return A_sparse.tocsr()  # Convert to csr for efficient arithmetic and matrix vector operations
