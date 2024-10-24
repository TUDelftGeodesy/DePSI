"""Module for creating networks from STM points."""

import math

import numpy as np
from scipy.spatial import Delaunay


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
    return math.dist(s, t)


def generate_arcs(stm_points, method="delaunay", x="lon", y="lat", max_length=None, max_links=12, num_groups=8):
    """Generate a network from a list of STM points.

    The network is undirected and without self-loops.

    Args:
    ----
        stm_points: Xarray.Dataset, input Space-Time Matrix.
        method: method to form the network; either "delaunay" or "redundant".
        x: str, first coordinate used to describe a point.
        y: str, second coordinate used to describe a point.
        max_length: float, maximum length of any generated arc or None.
        max_links: int, maximum number of arcs per node for the redundant method.
        num_groups: int, number of parts to split the orientations around each node into for the redundant method.

    Returns:
    -------
        arcs: list of pairs, point indices describing the adjacent nodes. The pairs are sorted, as is the list.
    """
    if method == "redundant":
        if max_links <= 0:
            print(f"max_links must be positive (currently: {max_links})")
            return
        if num_groups <= 0:
            print(f"num_groups must be positive (currently: {num_groups})")
            return

    # Collect point coordinates.
    indexes = [stm_points[coord] for coord in [x, y]]
    coordinates = np.column_stack(indexes)

    arcs = None

    # Create network arcs.
    if method == "delaunay":
        arcs = _generate_arcs_delaunay(coordinates, max_length)
    elif method == "redundant":
        arcs = _generate_arcs_redundant(coordinates, max_length, max_links, num_groups)

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


def _generate_arcs_redundant(coordinates, max_length=None, max_links=12, num_groups=8):
    # Create a network with at most max_links arcs per node.
    # Arcs are created ordered by length.
    # However, the orientations around the node are split into num_groups groups;
    # each group can only get an (x+1)th arc if every other group either
    # already has x arcs connected or already has all allowed arcs connected
    # (e.g. there are no more nodes in that group, or they are all too far away).

    arcs = []
    indices = range(len(coordinates))
    for cur_index in indices:
        # Calculate per node, the group they are in and the distance from the current node.
        groups = [
            int(math.floor(num_groups * (0.5 + math.atan2(coordinate[1], coordinate[0]) / math.tau)))
            for coordinate in coordinates - coordinates[cur_index]
        ]
        groups[cur_index] = num_groups + 1  # Separate the current node into its own group.
        distances = [get_distance(coordinates[cur_index], coordinate) for coordinate in coordinates]

        # Create a list of tuples with the group, distance, and index, sorted by group and then distance.
        values = np.array(sorted(list(zip(groups, distances, indices, strict=False))))

        # Collect the nearest max_links neighbors per group and discard the group of the current node.
        groups_diff = values[1:, 0] - values[:-1, 0]
        separators = np.where(groups_diff > 0)[0]
        groups = np.split(values, separators + 1)
        groups = groups[: len(groups) - 1]
        groups = [group[:max_links] for group in groups]

        # Collect the neighbor 'hierarchies'.
        # Each hierarchy contains the nth nearest neighbor from all groups.
        neighbor_hierarchies = [[] for _ in range(max_links)]
        count = 0
        for n in range(max_links):
            # Break early if we have gathered enough neighbors.
            if max_links <= count:
                break
            for group in groups:
                # Note that we do not break inside this loop, because we want the nth nearest neighbors from all groups.
                if n < len(group) and (max_length is None or group[n][1] <= max_length):
                    neighbor_hierarchies[n].append(group[n])
                    count = count + 1
        # Possibly, these loops could be replaced by this.
        # neighbor_hierarchies = [np.stack([el for el in zipped \
        #                         if (el is not None and (max_length is None or el[1] <= max_length))]) \
        #                         for zipped in itertools.zip_longest(*groups)]

        # Sort hierarchies per group by distance to the current node.
        neighbor_hierarchies = [
            sorted(hierarchy, key=lambda x: x[1]) for hierarchy in neighbor_hierarchies if len(hierarchy) != 0
        ]

        # Add sorted arcs to at most max_links neighbors.
        cur_arcs = [
            tuple(sorted([cur_index, int(neighbor[2])])) for hierarchy in neighbor_hierarchies for neighbor in hierarchy
        ]
        cur_arcs = cur_arcs[:max_links]

        arcs.extend(cur_arcs)

    # TODO(tvl): Fix issue with max_links: a later node may add an extra arc to an already 'full' earlier node.
    #           For example, node 1 may be connected to max_links neighbors other than node 2,
    #           before node 2 adds an arc connecting it to node 1 as one of its max_links neighbors.

    # Remove duplicates and make the list canonical.
    arcs = sorted(list(set(arcs)))

    return arcs
