"""Module for creating networks from STM points."""

import logging
import math

import numpy as np
import xarray as xr
from scipy.spatial import Delaunay

logger = logging.getLogger(__name__)


def generate_arcs(
    stm,
    network_method="redundant",
    x="lon",
    y="lat",
    max_length=None,
    min_links=12,
    num_partitions=8,
    difference="subtract",
) -> xr.Dataset:
    """Get an STM of arcs and phase differences from an STM of points.

    Args:
    ----
        stm: Xarray.Dataset, input Space-Time Matrix.
        network_method: method to form the network; either "delaunay" or "redundant".
        x: str, first coordinate used to describe a point.
        y: str, second coordinate used to describe a point.
        max_length: float, maximum length of any generated arc or None.
        min_links: int, minimum number of arcs per node, limited by max_length. Only used for the redundant method.
        num_partitions: int, number of partitions to split the nodes into based on orientation from the current node.
          Only used for the redundant method.
        difference: str, method for computing the phase difference; either "subtract" or "conjmult".

    Returns:
    -------
        arcs: Xarray.Dataset, STM of arcs, pairs of point indices describing the adjacent nodes and the difference
          between their phases.
          The index pairs are sorted, as is the list of pairs.
          The phase difference depends on the method used:
            either the source phase subtracted from the target phase,
            or the wrapped conjugate multiplication of these phases.

    Raises:
    ------
    NotImplementedError
        Raised when an unknown network or difference method is provided.
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
    indices = [stm[coord] for coord in [x, y]]
    coordinates = np.column_stack(indices)

    arcs = None

    # Create network arcs.
    if network_method == "delaunay":
        arcs = _generate_arcs_delaunay(coordinates, max_length)
    elif network_method == "redundant":
        arcs = _generate_arcs_redundant(coordinates, max_length, min_links, num_partitions)

    # Compute the phase difference.
    arcs_unzipped = list(zip(*arcs, strict=False))
    arcs_unzipped = [list(arcs_unzipped[0]), list(arcs_unzipped[1])]
    d_phase = _compute_phase_difference(stm, arcs_unzipped[0], arcs_unzipped[1], method=difference)

    # Store the phase difference in a DataArray,
    # with source and target coordinates as indices into the points STM.
    d_phase_array = xr.DataArray(
        d_phase,
        name="d_phase",
        dims=("space", "time"),
        coords={
            "source": (["space"], arcs_unzipped[0]),
            "target": (["space"], arcs_unzipped[1]),
            "time": stm.time,
        },
    )

    # Create a dataset to hold the array.
    stm_arcs = xr.Dataset({"d_phase": d_phase_array})

    return stm_arcs


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
    # Create a network with at least min_links arcs per node.
    # Arcs are created ordered by length.
    # However, the orientations around the node are split into num_partitions partitions;
    # each partition can only get an (x+1)th arc if every other partition either
    # already has x arcs connected or already has all allowed arcs connected
    # (e.g. there are no more nodes in that partition, or they are all too far away).
    # Note that a node may get less than min_links arcs if there are not enough neighbors within max_length.

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


def _compute_direct_phase_difference(stm, source_idx, target_idx):
    # Calculate the unwrapped direct phase difference between two points,
    # as the phase of the target minus the phase of the source.
    d_phase = stm.isel(space=target_idx).phase - stm.isel(space=source_idx).phase
    return d_phase


def _compute_wrapped_phase_difference(stm, source_idx, target_idx):
    # Calculate the wrapped phase difference between two points,
    # as the wrapped complex conjugate multiplication of the phases of the target and the source.

    # The original code in `demo_dynamic_estimation.ipynb` used `.sd_complex` (single difference complex),
    # which is computed as the SD (Single (temporal) Difference) phase values between the stm and a mother epoch
    # (`compute_sd` function called from `output_stm.ipynb`; probably imported from `arc_estimation_toolbox`).

    # Extract information of the two points of the arc
    complex_source = stm.isel(space=source_idx).complex
    complex_target = stm.isel(space=target_idx).complex

    # Compute DD phase for the arc
    complex_conj_source = complex_source.conj()
    d_phase = complex_target * complex_conj_source

    # Get the wrapped phase
    d_phase_wrapped = np.angle(d_phase)

    return d_phase_wrapped


def _compute_phase_difference(stm, source_idx, target_idx, method="subtract"):
    # Calculate the phase difference between two points.
    # The method can either be "subtract" or "conjmult".
    if method == "subtract":
        d_phase = _compute_direct_phase_difference(stm, source_idx, target_idx)
    elif method == "conjmult":
        d_phase = _compute_wrapped_phase_difference(stm, source_idx, target_idx)
    else:
        raise NotImplementedError(f"Unknown difference method {method}, known are subtract and conjmult")
    return d_phase
