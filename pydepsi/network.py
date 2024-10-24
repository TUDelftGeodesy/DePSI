"""Module for creating networks from STM points."""

import logging
import math

import numpy as np
from scipy.spatial import Delaunay

logger = logging.getLogger(__name__)


def get_distance(s, t):
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


def generate_arcs(stm_points, method="delaunay", x="lon", y="lat", max_length=None, min_links=12, num_partitions=8):
    """Generate a network from a list of STM points.

    The network is undirected and without self-loops.

    Args:
    ----
        stm_points: Xarray.Dataset, input Space-Time Matrix.
        method: method to form the network; either "delaunay" or "redundant".
        x: str, first coordinate used to describe a point.
        y: str, second coordinate used to describe a point.
        max_length: float, maximum length of any generated arc or None.
        min_links: int, minimum number of arcs per node, limited by max_length. Only used for the redundant method.
        num_partitions: int, number of partitions to split the nodes into based on orientation from the current node.
          Only used for the redundant method.

    Returns:
    -------
        coordinates: list, [x, y] point coordinates extracted from stm_points.
        arcs: list of pairs, point indices describing the adjacent nodes. The pairs are sorted, as is the list.
    """
    if method == "redundant":
        if min_links <= 0:
            logger.warning(f"min_links must be strictly positive (currently: {min_links})")
            return
        if num_partitions <= 0:
            logger.warning(f"num_partitions must be strictly positive (currently: {num_partitions})")
            return

    # Collect point coordinates.
    indexes = [stm_points[coord] for coord in [x, y]]
    coordinates = np.column_stack(indexes)

    arcs = None

    # Create network arcs.
    if method == "delaunay":
        arcs = _generate_arcs_delaunay(coordinates, max_length)
    elif method == "redundant":
        arcs = _generate_arcs_redundant(coordinates, max_length, min_links, num_partitions)

    return coordinates, arcs


def _generate_arcs_delaunay(coordinates, max_length=None):
    # Create network and collect neighbors.
    network = Delaunay(coordinates)
    neighbors_ptr, neighbors_idx = network.vertex_neighbor_vertices

    # Convert ptr and idx arrays into list of sorted index pairs.
    arcs = []
    for s in range(len(neighbors_ptr) - 1):
        for t in range(neighbors_ptr[s], neighbors_ptr[s + 1]):
            length = get_distance(coordinates[int(s)], coordinates[neighbors_idx[t]])
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
        distances = [get_distance(coordinates[cur_index], coordinate) for coordinate in coordinates]

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
