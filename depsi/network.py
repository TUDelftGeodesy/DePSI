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

from depsi.mht_utils import pretest

logger = logging.getLogger(__name__)

# Constants for MHT in network integration
ALPHA0 = 0.1  # Significance level for 1-dimensional test
GAMMA0 = 0.5  # Power of the test
OMT_THRES = 1e-7  # Overall Model Test threshold for accepting the network
# In arc/point rejection phase, if OMT < OMT_THRES, stop rejection iteration
# In ambiguity fixing phase, if OMT < OMT_THRES, stop fixing iteration
# In arc/point rejection phase this is hardly triggered
TT1_THRES = 1.0  # Threshold for arc rejection statistics TT1,
# If for all arcs max(TT1) < TT1_THRES, stop rejection iteration
# For most cases this threshold is triggered in rejection phase


def form_network(
    stm: xr.Dataset,
    key_phase: str,
    key_h2ph: str,
    key_Btemp: str,
    key_complex: str = "complex",
    key_xcrds: str = "lon",
    key_ycrds: str = "lat",
    network_method: Literal["redundant", "delaunay"] = "redundant",
    max_length: float = None,
    min_links: int = 16,
    num_partitions: int = 8,
    dphase_method: Literal["conjmult", "subtract"] = "subtract",
) -> xr.Dataset:
    """Generate an STM of arcs from an STM of points.

    Parameters
    ----------
    stm : xr.Dataset
        Space-Time Matrix of scatterers.
    key_phase : str
        Key of the phase values in the STM.
        This phase will be used to compute the differential arc phase.
    key_h2ph : str
        Key of the h2ph values in the STM.
        The arc h2ph will be computed as the average between source and target.
    key_Btemp : str
        Key of the temporal baseline values in the STM.
    key_complex : str, optional
        Key of the complex values, by default "complex"
    key_xcrds  : str, optional
        Key of the x coordinates for calulating arc length, by default "lon"
    key_ycrds  : str, optional
        Key of the y coordinates for calulating arc length, by default "lat"
    network_method : Literal["redundant", "delaunay"], optional
        network formation method, by default "redundant"
    max_length : float, optional
        maximum arc length, by default None
    min_links : int, optional
        minimum links per point, by default 16
        only effective when network_method is "redundant"
    num_partitions : int, optional
        number of partitions of searching when forming redundant network, by default 8
        only effective when network_method is "redundant"
    dphase_method : Literal["conjmult", "subtract"], optional
        method of computing phase difference, by default "subtract"
        "subtract" method subtracts the source phase from the target phase (without re-wrapping);
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
    indices = [stm[coord] for coord in [key_xcrds, key_ycrds]]
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

    # Temporal base line
    Btemp = stm[key_Btemp].data

    # Height to phase factor
    h2ph = (stm[key_h2ph].isel(space=source_idx).data + stm[key_h2ph].isel(space=target_idx).data) / 2

    # Generate a unique identifier of arcs based on source and target for easy indexing
    # This is because when updating network, points can be removed and reindexed
    # Therefore we cannot use 2d index (source, target) as uid
    scale = 10 ** (math.floor(math.log10(stm.sizes["space"])) + 1)  # Scale to ensure no overlap
    uid = scale * (np.array(source_idx) + 1) + (np.array(target_idx) + 1)  # Plus one to avoid zero uid
    uid = uid.astype(np.int64)

    arcs = xr.Dataset(
        data_vars={
            "d_phase": (["space", "time"], d_phase),
            "Btemp": (["time"], Btemp),
            "h2ph": (["space", "time"], h2ph),
        },
        coords={"source": (["space"], source_idx), "target": (["space"], target_idx), "uid": (["space"], uid)},
    )

    return arcs


def _mht_network_adjustment(
    stm_arcs: xr.Dataset,
    stm_pnts: xr.Dataset,
    idx_refpnt: int,
    azimuth_refpnt: int | float,
    range_refpnt: int | float,
    Qyy_diag: np.ndarray,
) -> (xr.Dataset, xr.Dataset):
    # First estimation
    A_sparse = _network_relation_matrix(
        stm_arcs["source"], stm_arcs["target"], stm_pnts.sizes["space"], idx_refpnt
    )  # Network relation matrix A
    y = stm_arcs["ambigs"].data  # Observations y
    invQy = scipy.sparse.diags(
        1 / Qyy_diag, 0, shape=(stm_arcs.sizes["space"], stm_arcs.sizes["space"])
    )  # Stochastic model assuming independent observations
    _, echeck = _solve_float_ambiguities(A_sparse, y, invQy)

    # Setup tests
    kb_dict = {}
    max_con = np.abs(A_sparse).sum(axis=0).max()
    for n_con in range(1, max_con + 1):
        _, k1, kb, _ = pretest(n_con, ALPHA0, GAMMA0)
        kb_dict[n_con] = kb

    # Compute test statistics for Overall Model Test
    OMT = (echeck.T @ invQy @ echeck).diagonal().sum()

    # Initial TT1_max to trigger the while loop
    stm_updated = stm_pnts.copy()
    stm_arcs_updated = stm_arcs.copy()
    TT1max = TT1_THRES + 1.0
    niter = 0
    while (TT1max > TT1_THRES) and (OMT > OMT_THRES) and (niter < stm_arcs.sizes["space"]):
        # In the loop, OMT fail
        # Choose from two Ha: 1) remove an arc; 2) remove a point
        flag_rm, idx_rm, TT1max, TTqmax = _mht_network_adjustment_reject_one(A_sparse, y, Qyy_diag, k1, kb_dict)

        if flag_rm == 0:  # remove arcs
            stm_arcs_updated = stm_arcs_updated.drop_isel(space=idx_rm)  # Remove the arc
        elif flag_rm == 1:  # remove points
            if idx_rm >= idx_refpnt:
                idx_rm += 1  # Adjust index due to removed reference point column in A_sparse

            # Removing points is achieved by removing all arcs connects to the point
            # Later the points will be actually removed when ensuring minimum connections
            # Arc indices connecting to the point to remove
            idx_arcs_selected = np.where(
                ((stm_arcs_updated["source"] != idx_rm) & (stm_arcs_updated["target"] != idx_rm)).data
            )[0]
            # Remove all arcs connects to the point to remove
            stm_arcs_updated = stm_arcs_updated.isel(space=idx_arcs_selected)

        # Ensure all points have at least 3 connections
        previous_size = -1  # Initialize with an impossible value to trigger the while loop
        # Keep iterating until no more points are removed
        while stm_updated.sizes["space"] != previous_size:
            previous_size = stm_updated.sizes["space"]
            # Remove points with <=2 connections
            stm_updated, stm_arcs_updated = remove_network_points_min_connections(
                stm_updated, stm_arcs_updated, min_connections=3
            )

        # Make sure the reference point is still in stm_updated, by checking its azimuth and range
        mask_refpnt = (stm_updated["azimuth"].values == azimuth_refpnt) & (stm_updated["range"].values == range_refpnt)
        if not np.any(mask_refpnt):
            raise ValueError(
                f"Reference point ({azimuth_refpnt}, {range_refpnt}) removed in the MHT process. "
                f"Please choose another reference point."
            )
        idx_refpnt = np.where(mask_refpnt)[0][0]  # Update idx_refpnt

        # Get indices of selected arcs based on uid
        idx_arcs_selected = np.where(stm_arcs_updated["uid"].isin(stm_arcs["uid"]))[0]

        # Update the functional and stochastic model
        Qyy_diag = Qyy_diag[idx_arcs_selected]
        invQy = scipy.sparse.diags(
            1 / Qyy_diag, 0, shape=(stm_arcs_updated.sizes["space"], stm_arcs_updated.sizes["space"])
        )
        A_sparse = _network_relation_matrix(
            stm_arcs_updated["source"], stm_arcs_updated["target"], stm_updated.sizes["space"], idx_refpnt
        )

        y = stm_arcs_updated["ambigs"].data

        # Estimate residual again
        _, echeck = _solve_float_ambiguities(A_sparse, y, invQy)

        OMT = (echeck.T @ invQy @ echeck).diagonal().sum()

        niter += 1

    return stm_arcs_updated, stm_updated


def _mht_network_adjustment_reject_one(
    A: scipy.sparse._csr.csr_matrix,
    y: np.ndarray,
    Qyy_diag: np.ndarray,
    k1: float,
    kb_dict: dict,
) -> (int, int):
    """Remove one point/arc from the network to reduce the residual in ambiguity estimation."""
    # Retrive shapes
    N_arcs, N_epochs = y.shape
    N_points = A.shape[1]

    # Inverse of VCM of observations
    if Qyy_diag.ndim == 1:  # Diagonal VCM
        invQy = scipy.sparse.diags(1 / Qyy_diag, 0, shape=(N_arcs, N_arcs))
        Qyy = scipy.sparse.diags(Qyy_diag, 0, shape=(N_arcs, N_arcs))
    else:
        raise NotImplementedError("Currently only diagonal VCM is supported. Qyy_diag should be an 1d array.")

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
    TT1max = max(TT1)

    # Test statistics for removing one point
    TTq = np.zeros(N_points)
    for pnt_idx in range(N_points):
        arcs_idx = np.where(A[:, pnt_idx].todense() != 0)[0]  # Arcs connected to this point
        arcs_idx = arcs_idx[1:]  # Drop one arc to create basis, see e.g. verhoef97
        echeck_point = echeck[arcs_idx, :]  # Relevant echeck of this point
        Qecheck_point = Qecheck[arcs_idx, :][:, arcs_idx]  # Relevant Qecheck of this point

        # Compute the test statistic for this point
        Tq = np.sum(
            (echeck_point.T @ np.linalg.inv(Qecheck_point) @ echeck_point).diagonal()
        )  # Before adjust for degree of freedom
        TTq[pnt_idx] = Tq / kb_dict[len(arcs_idx)]
    TTqmax = max(TTq)

    # Decision one removal strategy
    if TT1max > TTqmax:
        idx_removal = np.argmax(TT1)  # index of arc to remove
        flag_removal = 0  # remove arc
    else:
        idx_removal = np.argmax(TTq)  # index of point to remove
        flag_removal = 1  # remove point

    return flag_removal, idx_removal, TT1max, TTqmax


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
    min_n_connections: int = 2,
) -> xr.Dataset:
    """Select arcs based on arc quality and connectivity.

    This function selects arcs in two steps:
    1. It selects arcs based on a threshold value (e.g., ens_coh).
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
    min_n_connections : int, optional
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
    # These arcs are identified by the points which have <= min_n_connections arcs connected to them
    # All arcs connected to such points are removed
    # An iterative approach is used to remove all arcs connected to such points
    point_ids_all = np.concat(
        [arcs_selected["source"].data, arcs_selected["target"].data]
    )  # all occurrances of point ids
    point_ids_unique, counts = np.unique(point_ids_all, return_counts=True)  # unique point ids and their counts
    while np.any(counts <= min_n_connections):
        # Find points with <=3 arcs connected
        point_ids_to_remove = point_ids_unique[counts <= min_n_connections]
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


def remove_network_points_min_connections(stm: xr.Dataset, arcs: xr.Dataset, min_connections: int) -> xr.Dataset:
    """Remove points which have less than min_connections arc connections.

    The following steps are performed:

    1. Remove points from stm which have less than min_connections connections in arcs.
    2. Remove arcs which connect to the removed points.
    3. Update the space indices in points/arcs STM accordingly.
       The point indices is always a 0-based continuous array.

    Note that this function does not perform interative removal to assure that all points have
    at least min_connections connections, but only performs one round of removal.
    """
    if min_connections < 1:
        raise ValueError("min_connections must be at least 1")

    # Load source and target indices from arcs
    # these are 1d arrays so should fit in memory
    idx_source = arcs["source"].values
    idx_target = arcs["target"].values

    # Select STM points that are in arcs
    # Only keep points which ids are in arcs, isolated points are removed in idx_selected
    idx_selected, counts = np.unique(np.concatenate([idx_source, idx_target]), return_counts=True)
    idx_selected = idx_selected[counts >= min_connections]  # only keep points with at least min_connections connections

    # If no change, return directly
    if len(idx_selected) == stm.sizes["space"]:
        return stm, arcs

    # Select points
    stm_updated = stm.isel(space=idx_selected)

    # Select arcs that connect selected points (Some arcs may be dropped together with points)
    mask_source = np.isin(idx_source, idx_selected)
    mask_target = np.isin(idx_target, idx_selected)
    mask_arcs = mask_source & mask_target
    arcs = arcs.isel(space=np.where(mask_arcs)[0])

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


def _generate_arcs_redundant(coordinates, max_length, min_links, num_partitions):
    """Create a network with at least min_links arcs per node.

    The redundant network is formed with the following steps:

    1. Create a KDTree and find all pairs of points within the maximum distance.
    2. Loop through each point and find its neighbors within the maximum distance.
    3. Divide neighbors into partitions based on their direction.
    4. Select the nth nearest neighbors from all partitions, starting from n=1.
    5. Sort the selected neighbors by distance, add them to the arcs list. If min_links is not
       exceeded, continue to the n+1th nearest neighbors of all partitions.
    6. Repeat until min_links is reached.
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
        elif len(neighbors) <= min_links:
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
            partitions_split = [partition[:min_links] for partition in partitions_split]

            # Collect the neighbor 'hierarchies'
            neighbor_hierarchies = [[] for _ in range(min_links)]
            count = 0
            for n in range(min_links):
                # Break early if we have gathered enough neighbors.
                if min_links <= count:
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

            # Add sorted arcs to at least min_links neighbors
            cur_arcs = [
                (min(cur_index, int(neighbor[2])), max(cur_index, int(neighbor[2])))
                for hierarchy in neighbor_hierarchies
                for neighbor in hierarchy
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


def _network_relation_matrix(idx_source, idx_target, n_points, idx_refpnt):
    """Create the network relation matrix A as a sparse matrix.

    A network relation matrix has shape (n_arcs, n_points - 1).
    Each row corresponds to an arc, and each column corresponds to a point, excluding the reference point.
    For each arc, the column corresponding to the source point has a value of -1, and the column corresponding
    to the target point has a value of +1.
    All other entries are zero.

    The reference point column removal refers to Eq.4.11 of the following book:
    Kampes, Bert M. Radar interferometry: persistent scatterer technique. Dordrecht: Springer Netherlands, 2006.
    DOI: 10.1007/978-1-4020-4723-7

    Parameters
    ----------
    idx_source : list or np.ndarray
        List of source point indices for each arc.
    idx_target : list or np.ndarray
        List of target point indices for each arc.
    n_points : int
        Total number of points in the network.
    idx_refpnt : int
        Index of the reference point to be excluded from the matrix. This index assumes 0-based indexing of the points.
    """
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

    # Convert to Compressed Sparse Row (CSR) matrix for efficient arithmetic and matrix vector operations
    A_sparse = A_sparse.tocsr()

    # Remove reference point column
    A_sparse = scipy.sparse.hstack([A_sparse[:, :idx_refpnt], A_sparse[:, idx_refpnt + 1 :]])

    return A_sparse
