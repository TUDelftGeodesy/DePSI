"""Module for creating networks from STM points."""

import logging
import math
from typing import Literal

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from scipy.spatial import Delaunay, distance_matrix

from depsi.arc_estimation import arc_estimation_control_network
from depsi.utils import get_distance

logger = logging.getLogger(__name__)


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
            logger.error(f"min_links must be strictly positive (currently: {min_links})")
            return
        if num_partitions <= 0:
            logger.error(f"num_partitions must be strictly positive (currently: {num_partitions})")
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
            length = get_distance(coordinates[int(s)], coordinates[neighbors_idx[t]], mode="euclidean")
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
        distances = [get_distance(coordinates[cur_index], coordinate, mode="euclidean") for coordinate in coordinates]

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


def get_ordered_arcs(stm, x_ref_search, y_ref_search, buffer_radius_ref, dist_to_quality, n_max_arcs, nad_max=2):
    """Get a list with ordered arcs based on pnt quality and a search area.

    Args:
    ----
        stm (xr.Dataset): stm containing at least rd_x, rd_y, and slc quality
        x_ref_search (float): x-coordinate of the centre of the search area
        y_ref_search (float): x-coordinate of the centre of the search area
        buffer_radius_ref (float): The buffer (in m) around the centre coordinates where potential arcs are computed
        dist_to_quality (float): parameter that relates arc length to additional sigma
        n_max_arcs (float): The maximum nr of arcs to be outputed
        nad_max (float): The maximum NAD, for the entire time series, for a point to be considered

    Returns:
    -------
        arcs (list): ordered arcs
        arcs_and_quality (list): ordered list with the points and quality
        quality_dict (dictionary): dictionary with the arcs and their quality
    """
    rdx = stm["rd_x"].values
    rdy = stm["rd_y"].values
    slc_quality = stm["slc_quality"].values
    nad_vals = stm["nad_full"].values

    # Find all points within the buffer around the starting location x, y
    idx_pnts_buffer_ref = find_points_within_buffer(rdx, rdy, x_ref_search, y_ref_search, buffer_radius_ref)

    # Get the a-priori quality of all potential arcs that can be made
    arcs_and_quality = _ordered_arcs_all_points(
        rdx, rdy, slc_quality, idx_pnts_buffer_ref, dist_to_quality, nad_vals, nad_max=nad_max
    )

    # We will only work with n_max_arcs, otherwise we need to load an extensive dataset everytime
    arcs_and_quality = arcs_and_quality[0:n_max_arcs]

    # Only store the arcs, and remove the quality values
    arcs = []
    for _, arc in arcs_and_quality:
        arcs.append(arc)

    # Create a dictionary with the quality of the arcs
    quality_dict = {tuple(sorted(arc[1])): arc[0] for arc in arcs_and_quality}

    return arcs, arcs_and_quality, quality_dict


def find_points_within_buffer(x_coords, y_coords, x_pnts, y_pnts, buffer_radius):
    """Find all points located within a specified buffer radius around a given location.

    This function determines which points in a set of coordinates are located within a defined buffer
    distance around one or more specified points.

    The function:
    1. Converts the input x and y coordinates of the search point(s) to arrays if they are not already.
    2. Computes the Euclidean distance between each point in `x_coords` and `y_coords` and the reference points.
    3. Returns the indices of these points for further processing or analysis.

    Args:
    ----
        x_coords (numpy.ndarray): Array of x-coordinates for all points in the dataset.
        y_coords (numpy.ndarray): Array of y-coordinates for all points in the dataset.
        x_pnts (float or array-like): The x-coordinate or array of x-coordinates of the reference point(s).
        y_pnts (float or array-like): The y-coordinate or array of y-coordinates of the reference point(s).
        buffer_radius (float): The radius of the buffer zone around the reference points, specified in meters.

    Returns:
    -------
        indices (numpy.ndarray): Array of indices of the points located within the buffer zone.

    Example:
        indices = _find_points_within_buffer(x_coords, y_coords, x_pnts=10.5, y_pnts=20.3, buffer_radius=5.0)
    """
    # Make sure x_pnts and y_pnts are arrays
    x_pnts = np.atleast_1d(x_pnts)
    y_pnts = np.atleast_1d(y_pnts)

    within_buffer = np.zeros(len(x_coords), dtype=bool)

    # Check for every point, which other points fall within the buffer around it
    for tx, ty in zip(x_pnts, y_pnts, strict=True):
        distances = np.sqrt((x_coords - tx) ** 2 + (y_coords - ty) ** 2)
        within_buffer |= distances <= buffer_radius  # 'OR' operation

    # Find the indices of the points that are within the buffer
    indices = np.where(within_buffer)[0]

    return indices


def _ordered_arcs_all_points(rdx, rdy, slc_quality, idx_pnts_buffer, dist_to_quality, nad_vals, nad_max):
    """Generate a sorted list of unique arcs between all points within buffer, ranked by combined quality and distance.

    This function computes arcs between all points in a given buffer based on their spatial distance and
    quality metrics. The arcs are ranked according to a combination of the Euclidean distance between
    points and the maximum quality value between them.

    The function:
    1. Computes the Euclidean distance matrix for all points within the buffer.
    2. Calculates the arc quality time series for all point pairs based on their SLC quality values.
    3. Extracts the maximum quality for each arc and combines it with the distance between the points.
    4. Uses the lower triangular matrix to avoid duplicate arcs
    5. Sorts the arcs based on the combined distance and quality metric.

    Args:
    ----
        rdx (numpy.ndarray): Array of x-coordinates for all points.
        rdy (numpy.ndarray): Array of y-coordinates for all points.
        slc_quality (numpy.ndarray): Array of SLC quality time series for all points.
        idx_pnts_buffer (numpy.ndarray): Indices of the points within the buffer.
        dist_to_quality (float): Scaling factor for weighting the distance in combination with the quality.
        nad_vals (numpy.ndarray): Array of NAD values of all points.
        nad_max (float): The maximum NAD, for the entire time series, for a point to be considered

    Returns:
    -------
        arcs_and_quality (list): Sorted list of tuples, each tuple contains the combined quality and distance value
                                 and the indices of the two points forming the arc.

    Example:
    -------
        arcs_and_quality = ordered_arcs_all_points(rdx, rdy, slc_quality, idx_pnts_buffer, dist_to_quality)
    """
    # Compute NAD values of the points within the buffer
    nad_buffer_vals = nad_vals[idx_pnts_buffer]

    # Mask points that have an NAD value above the threshold
    mask_nad = nad_buffer_vals < nad_max
    idx_pnts_buffer = idx_pnts_buffer[mask_nad]

    # Compute the distance matrix for all points (below the NAD threshold) within the buffer
    rdx_buffer = rdx[idx_pnts_buffer]
    rdy_buffer = rdy[idx_pnts_buffer]
    coords_points = np.vstack((rdx_buffer, rdy_buffer)).T
    dist_matrix = distance_matrix(coords_points, coords_points)

    # Compute the quality matrix for all potential arcs within the buffer
    slc_quality_buffer = slc_quality[idx_pnts_buffer, :]
    slc_quality_i_buffer = slc_quality_buffer[:, np.newaxis]
    slc_quality_j_buffer = slc_quality_buffer[np.newaxis, :]

    # Compute the arc quality time series
    # quality is based on both sides of the arc
    arc_quality_ts = np.sqrt(slc_quality_i_buffer**2 + slc_quality_j_buffer**2)

    # Get the maximum value in the time dimension per arc
    # This is considered as the 'worst' quality for the entire time period
    arcs_quality_max = np.max(arc_quality_ts, axis=2)

    # Get only the lower triangular matrix
    arcs_quality_max_lower = np.tril(arcs_quality_max)
    arcs_dist = np.tril(dist_matrix)

    # Add additional sigma because of the arc length
    arc_quality_dist_max = arcs_quality_max_lower + arcs_dist * dist_to_quality

    # Order the arcs
    arcs_and_quality = []

    # Loop trough the lower tringular matrix and add values list
    rows, cols = arc_quality_dist_max.shape
    for i in range(rows):
        for j in range(i):  # j < i makes sure that we only get the lower triangular matrix
            arcs_and_quality.append((arc_quality_dist_max[i, j], (idx_pnts_buffer[i], idx_pnts_buffer[j])))

    # Sort the list with all the arcs from best to worse
    arcs_and_quality.sort()

    return arcs_and_quality


def construct_control_network(
    arcs,
    quality_dict_arcs,
    excluded_arcs,
    n_top,
    n_batch,
    deg_threshold,
    min_nodes,
    min_redundancy,
    visualize_network=False,
):
    """Construct a control network of arcs based on quality metrics and structural requirements.

    This function iteratively builds a network from a search area of arcs, filtering out excluded arcs
    (where an excluded arc is an arc we know we don't want to have) and ensuring that the resulting
    network meets specified criteria, such as a minimum number of nodes and redundancy.
    The network is constructed by:
    1. Ranking arcs based on quality metrics provided in `quality_dict_arcs`.
    2. Excluding arcs that are in the `excluded_arcs` list.
    3. Iteratively adding arcs in batches and refining the network to remove nodes with low centrality.
    4. Evaluating the network against requirements such as minimum nodes and average redundancy.

    If the requirements are not met with the given arcs, a warning is issued. Optionally, the network
    construction process can be visualized.

    Args:
    ----
        arcs (list of tuples): List of arcs (pairs of points) within the search area.
        quality_dict_arcs (dict): Dictionary mapping arcs (tuples) to a sigma value.
        excluded_arcs (list of tuples): List of arcs that should be excluded from the network.
        n_top (int): Number of top-ranked arcs to start with.
        n_batch (int): Number of additional arcs to add in each iteration.
        deg_threshold (float): Degree threshold for removing low-centrality nodes from the network.
        min_nodes (int): Minimum number of nodes required in the final network.
        min_redundancy (float): Minimum average redundancy (degree) required in the final network.
        visualize_network (bolean, optional): If set to True, the network construction process is visualized.

    Returns:
    -------
        ref_pnt (int): The reference point used in the final iteration of the network.
        arcs_updated_network_sorted (list of tuples): The sorted list of arcs in the final network
            based on quality metrics.
        ref_pnt_initial (int): The reference point used in the initial network.
        arcs_initial_network (list of tuples): The list of arcs in the network before refinement.

    Example:
    -------
        ref_pnt, final_arcs, initial_ref_pnt, initial_arcs = construct_control_network(
            arcs, quality_dict_arcs, excluded_arcs,
            n_top=50, n_batch=10,
            deg_threshold=1.5, min_nodes=20, min_redundancy=2.0,
            visualize_network=1
        )
    """
    # Sort arcs by quality, removing excluded arcs from the ranked arcs variable
    excluded_arcs_sorted = [tuple(sorted(arc)) for arc in excluded_arcs]
    arcs_without_excluded = [arc for arc in arcs if tuple(sorted(arc)) not in excluded_arcs_sorted]

    # Initialize variables
    current_network = None
    avg_degree = 0
    iteration = 0

    more_arcs_to_test = True
    network_reqs_not_met = True

    # the following while loop will stop when either the network requirements (# nodes and redundancy level) are met,
    # OR when we run out of arcs to text
    while more_arcs_to_test and network_reqs_not_met:
        # Select the arcs for the current iteration
        iteration += 1
        end_idx = n_top + iteration * n_batch
        arcs_to_test = arcs_without_excluded[0:end_idx]

        if not arcs_to_test:  # Break if no more arcs to add
            print("No more arcs to test.")
            more_arcs_to_test = False
            continue

        # Create a new network using the selected arcs
        current_network = nx.Graph()
        current_network, _, ref_pnt = _from_arcs_to_graph(arcs_without_excluded[:end_idx], plot=visualize_network)
        if iteration == 1:
            ref_pnt_initial = ref_pnt
            arcs_initial_network = [tuple(sorted(arc)) for arc in current_network.edges()]

        # Remove nodes with degree equal or lower than the threshold
        current_network, _, ref_pnt = _remove_low_centrality_nodes(
            current_network, deg_threshold=deg_threshold, plot=visualize_network
        )
        arcs_updated_network = [tuple(sorted(arc)) for arc in current_network.edges()]

        # Test the network against requirements
        avg_degree = len(current_network.edges) / len(current_network.nodes) if len(current_network.nodes) > 0 else 0
        print(f"The average degree is {avg_degree:.2f}")

        if len(current_network.nodes) >= min_nodes and avg_degree >= min_redundancy:
            network_reqs_not_met = False
            # Stop if requirements are met

    # Final check if the network meets requirements
    if len(current_network.nodes) < min_nodes or avg_degree < min_redundancy:
        print("Warning: Network could not meet all requirements with the given arcs.")

    # The arcs constructed above are no longer sorted based on their quality
    arcs_updated_network_sorted = sorted(
        arcs_updated_network, key=lambda arc: quality_dict_arcs.get(tuple(sorted(arc)), float("inf"))
    )

    return ref_pnt, arcs_updated_network_sorted, ref_pnt_initial, arcs_initial_network


def _from_arcs_to_graph(arcs, plot=False, save_path="./network.png"):
    """Construct a graph from a set of arcs and identifies key properties of the network.

    Args:
    ----
        arcs : (list of tuple)
            List of arcs (edges) in the graph. Each arc is represented as a tuple
            of two nodes, e.g., [(node1, node2), (node2, node3)].
        plot : (boolean)
            Default is False, If True, the function visualizes the graph structure, including highlighting
            connected components. Default is False (no visualization).
        save_path : (string)
            Default = "./network.png"

    Returns:
    -------
        tuple: A tuple containing:
            - network (networkx.Graph): The constructed graph based on the provided arcs.
            - degree_centrality (dict): A dictionary where keys are nodes and values are their
                                        degree centrality scores.
            - ref_pnt (int): The index of the reference point (point with most connections).
    """
    if not arcs:
        raise ValueError("The input arcs list is empty. Please provide a valid list of arcs.")

    # Compute the network
    network = nx.Graph()
    network.add_edges_from(arcs)

    # Compute the degree centrality
    # Which says something on how connected one particular node is with other nodes
    degree_centrality = nx.degree_centrality(network)

    # Get the reference point (node with the most connections)
    ref_pnt = max(degree_centrality, key=degree_centrality.get)

    if plot:
        # Visualize the network
        pos = nx.spring_layout(network)  # Lay-out for the graph
        nx.draw(
            network,
            pos,
            with_labels=True,
            node_size=500,
            node_color="skyblue",
            font_size=10,
            font_weight="bold",
            edge_color="gray",
        )

        connected_components = list(nx.connected_components(network))
        # Highlight the connected components
        for component in connected_components:
            nx.draw_networkx_nodes(network, pos, nodelist=component, node_color="orange", node_size=700)

        # Highlight the reference point
        nx.draw_networkx_nodes(network, pos, nodelist=[ref_pnt], node_color="red", node_size=800)
        plt.title(f"Network Visualization - {len(connected_components)} Connected Components")
        plt.savefig(save_path, bbox_inches="tight")

    return network, degree_centrality, ref_pnt


def _remove_low_centrality_nodes(network, deg_threshold, plot=False, save_path="./network.png"):
    """Remove nodes with a low degree of centrality.

    Args:
    ----
        network (networkx.Graph): The input graph from which nodes will be removed.
        deg_threshold (int): The degree threshold; nodes with a degree equal to or less than this value will be removed.
        plot (bolean, optional): If True, the function visualizes the updated graph, including its connected components.
                              Default is False (no visualization).
        save_path : (string)
            Default = "./network.png"

    Returns:
    -------
        tuple: A tuple containing:
            - network (networkx.Graph): The updated graph after removing nodes with the specified degree threshold.
            - degree_centrality (dict): A dictionary where keys are nodes and values are their degree centrality scores.
            - ref_pnt (int or None): The node with the highest degree centrality in the updated graph.

    Example:
    -------
        updated_network, degree_centrality, ref_point = remove_low_centrality_nodes(
            network, deg_threshold=1, plot=False
        )
    """
    # Identify nodes to remove
    nodes_to_remove = [node for node, degree in network.degree() if degree <= deg_threshold]
    removed_nodes = set(nodes_to_remove)  # For plot highlighting

    # Remove these nodes from the network
    network.remove_nodes_from(nodes_to_remove)

    # Handle empty graph case
    if not network.nodes:
        return network, {}, None

    # Degree centrality
    degree_centrality = nx.degree_centrality(network)

    # Get the reference point (node with the highest degree centrality)
    ref_pnt = max(degree_centrality, key=degree_centrality.get)

    if plot:
        # Visualize the updated network
        pos = nx.spring_layout(network)
        nx.draw(
            network,
            pos,
            with_labels=True,
            node_size=500,
            node_color="skyblue",
            font_size=10,
            font_weight="bold",
            edge_color="gray",
        )

        connected_components = list(nx.connected_components(network))

        # Highlight connected components
        for component in connected_components:
            nx.draw_networkx_nodes(network, pos, nodelist=component, node_color="orange", node_size=700)

        # Highlight removed nodes in the plot (optional)
        if removed_nodes:
            # Filter valid nodes that exist in the layout
            valid_removed_nodes = [node for node in removed_nodes if node in pos]
            nx.draw_networkx_nodes(network, pos, nodelist=valid_removed_nodes, node_color="red", node_size=800)

        nx.draw_networkx_nodes(network, pos, nodelist=[ref_pnt], node_color="red", node_size=800)

        plt.title(f"Network Visualization after Removing Nodes with Degree ≤ {deg_threshold}")
        plt.savefig(save_path, bbox_inches="tight")

    return network, degree_centrality, ref_pnt


def test_succeeded_arcs_control_network(
    succeeded_arcs, quality_dict_arcs, deg_threshold, min_nodes, min_redundancy, visualize_network=False
):
    """Evaluate a network constructed from a set of arcs to determine if it meets structural requirements.

    This function tests a network built from a set of arcs to ensure it satisfies minimum structural criteria,
    including the number of nodes, redundancy, and centrality. It checks whether any nodes with low degree should
    be removed and then evaluates the network again. The arcs are sorted by quality after refinement.

    Parameters
    ----------
    succeeded_arcs : list of tuples
        A list of arcs (pairs of points) that are used to construct the network.
    quality_dict_arcs : dict
        A dictionary where the keys are the arcs and the values are their respective quality scores.
    deg_threshold : float
        The degree threshold for identifying low-centrality nodes. Nodes with degrees smaller will be removed.
    min_nodes : int
        The minimum number of nodes required for the network to be valid.
    min_redundancy : float
        The minimum average degree (redundancy) required for the network.
    visualize_network : bool, optional
        A flag to control network visualization of the network during evaluation. Default is False (no visualization).

    Returns
    -------
    network_check : int
        A flag indicating whether the network meets the requirements (1 if successful, 0 if not).
    arcs_updated_network_sorted : list of tuples
        The sorted list of arcs based on their quality after refining the network.
    ref_pnt : int
        The reference point used in the network.

    Example
    -------
    network_check, final_arcs, ref_pnt = test_succeeded_arcs_control_network(
        succeeded_arcs, deg_threshold=1.5, min_nodes=20,
        min_redundancy=2.0, visualize_network=1
    )
    """
    arcs_updated_network_sorted = []
    network_check = 0

    # Construct the network based on the input arcs
    current_network = nx.Graph()
    current_network, _, ref_pnt = _from_arcs_to_graph(succeeded_arcs, plot=visualize_network)
    arcs_updated_network = [tuple(sorted(arc)) for arc in current_network.edges()]

    # Check whether there are nodes with degree smaller than one
    nodes_to_remove = [node for node, degree in current_network.degree() if degree <= deg_threshold]

    if not nodes_to_remove:
        # Test the network against requirements
        avg_degree = len(current_network.edges) / len(current_network.nodes) if len(current_network.nodes) > 0 else 0
        print(f"The average degree is {avg_degree:.2f}")

        if len(current_network.nodes) >= min_nodes and avg_degree >= min_redundancy:
            network_check = 1

            # The arcs constructed above are no longer sorted based on their quality
            arcs_updated_network_sorted = sorted(
                arcs_updated_network, key=lambda arc: quality_dict_arcs.get(tuple(sorted(arc)), float("inf"))
            )

    else:
        print("We remove low degree nodes")
        # Remove the nodes with low degree and test the network again
        current_network, _, ref_pnt = _remove_low_centrality_nodes(
            current_network, deg_threshold=deg_threshold, plot=visualize_network
        )
        arcs_updated_network = [tuple(sorted(arc)) for arc in current_network.edges()]

        # Test the network against requirements
        avg_degree = len(current_network.edges) / len(current_network.nodes) if len(current_network.nodes) > 0 else 0
        print(f"The average degree is {avg_degree:.2f}")

        if len(current_network.nodes) >= min_nodes and avg_degree >= min_redundancy:
            network_check = 1

            # The arcs constructed above are no longer sorted based on their quality
            arcs_updated_network_sorted = sorted(
                arcs_updated_network, key=lambda arc: quality_dict_arcs.get(tuple(sorted(arc)), float("inf"))
            )

        else:
            network_check = 0

    return network_check, arcs_updated_network_sorted, ref_pnt


def ordered_arcs_connection_point_and_control_network(
    rdx_connection_point,
    rdy_connection_point,
    slc_quality_connection_point,
    connection_point_idx,
    rdx_control,
    rdy_control,
    slc_quality_control,
    control_idx,
    dist_to_quality,
):
    """Generate a sorted array of unique arcs between a 'connection_point' and the control network.

    This function computes arcs between a connection point and control points, evaluating each arc based on
    a combination of spatial distance and SLC quality. The arcs are sorted from best to worst, where "best"
    is determined by the lowest combined quality and distance value.

    The function:
    1. Computes the Euclidean distance between the connection point and the control points.
    2. Calculates the arc quality time series for each arc.
    3. Combines the distance and maximum quality for each arc to define an overall quality value.
    4. Sorts the arcs by quality values from best (lowest) to worst (highest).
    5. Constructs a sorted array of arcs, ensuring the control point appears first in each arc.

    Args:
    ----
        rdx_connection_point (xarray.DataArray): x-coordinate of the connection point.
        rdy_connection_point (xarray.DataArray): y-coordinate of the connection point.
        slc_quality_connection_point (xarray.DataArray): Array of SLC quality time series for the connection point.
        connection_point_idx (xarray.DataArray): Index of the connection point.
        rdx_control (xarray.DataArray): x-coordinates of the control points.
        rdy_control (xarray.DataArray): y-coordinates of the control points.
        slc_quality_control (xarray.DataArray): Array of SLC quality time series for the control points.
        control_idx (xarray.DataArray): Indices of the control points.
        dist_to_quality (float): Conversion factor to scale distance relative to quality.

    Returns:
    -------
        tuple:
            - arcs (numpy.ndarray): A sorted array of arcs, where each row contains a control point
              followed by the connection point.
            - sorted_quality_values (numpy.ndarray): Sorted quality values corresponding to the arcs.
    """
    # Stack the coordinates of the connection point and the control points
    coords_connection_point = np.vstack((rdx_connection_point, rdy_connection_point)).T
    coords_control = np.vstack((rdx_control, rdy_control)).T

    # Distance matrix between the connection point and the control points
    dist_matrix = distance_matrix(coords_connection_point, coords_control).squeeze()

    # # Compute the quality matrix for connection point (pnt i) and all control points (point j)
    slc_quality_j = slc_quality_control.values  # Quality values of the control points
    slc_quality_i = np.expand_dims(slc_quality_connection_point.values, axis=0)  # Make sure the dimensions match
    slc_quality_i = np.repeat(slc_quality_i, repeats=slc_quality_j.shape[0], axis=0)

    # Compute arc quality time series (which is a function of the slc_quality of point i and point j)
    arc_quality_ts = np.sqrt(slc_quality_i**2 + slc_quality_j**2)
    arc_quality_max = np.max(arc_quality_ts, axis=1)  # Compute the maximum value

    # Combine distance and quality
    arc_quality_dist_max = arc_quality_max + dist_matrix * dist_to_quality

    # Sort the arcs and compute the indices of the control points (since point i, the connection_point, is in all arcs)
    sorted_control = np.argsort(arc_quality_dist_max)
    sorted_control_idx = control_idx.values[sorted_control]  # sort the indices of the control points as well
    sorted_quality_values = np.sort(arc_quality_dist_max)

    # Compute the arcs between the connection_point and the control points.
    # Make sure that the control_points comes first
    arcs = np.zeros((len(sorted_control_idx), 2), dtype=int)
    arcs[:, 0] = sorted_control_idx  # control network points
    arcs[:, 1] = connection_point_idx.values  # connection_point

    return arcs, sorted_quality_values


def construct_control_network_test_arcs(
    stm,
    x_ref_search,
    y_ref_search,
    buffer_radius_ref,
    dist_to_quality,
    N_max_arcs,
    N_top,
    N_batch,
    deg_threshold,
    min_nodes,
    min_redundancy,
    nad_max,
    visualize_network,
    sigma_post_over_sigma_prior,
    nr_max_iter_control,
    bounds,
    m2ph,
    years,
    dates,
    temperature,
    sd_complex,
    slc_quality,
    cr2ph,
    ampl_ts,
    bkps_stm,
    mean_ampl_sd,
    sigma_ampl_sd,
    mad_ampl_sd,
    median_ampl_sd,
    x_coordinates,
    y_coordinates,
    coordinate_type: Literal["euclidean", "geographic"] = "euclidean",
):
    """Test the arcs in the constructed control network.

    Parameters
    ----------
    stm:
        space_time_matrix
    x_ref_search:
        X coordinate of central point to search for reference point
    y_ref_search:
        Y coordinate of central point to search for reference point
    buffer_radius_ref:
        Size of search window for reference point
    dist_to_quality:
        relates arc length to additional sigma
    N_max_arcs:
        Maximum number of arcs
    N_top:
        Number of top-ranked arcs to start with.
    N_batch:
        Number of additional arcs to add in each iteration.
    deg_threshold:
        Degree threshold for removing low-centrality nodes from the network.
    min_nodes:
        Minimum number of nodes required in the final network.
    min_redundancy:
        Minimum average redundancy (degree) required in the final network.
    nad_max:
        Maximum NAD
    visualize_network:
        Boolean whether or not to visualize the network
    sigma_post_over_sigma_prior:
        Upper limit on how much the aposteriori sigma of an arc is allowed to differ from the apriori sigma
    nr_max_iter_control:
        Maximum number of iterations in the network testing
    bounds: list of tuples
        Bounds for parameter estimation in the format (lower_bounds, upper_bounds).
    m2ph: float
        Meters to phase conversion factor
    years: np.ndarray
        Array of decimal years corresponding to the time series epochs.
    dates: np.ndarray
        Array of date indices or timestamps corresponding to the time series.
    temperature: np.ndarray
        Array of temperature values for thermal expansion modeling.
    sd_complex: np.ndarray
        Single-difference complex
    slc_quality: np.ndarray
        Complex-valued standard deviations of the signal for all points.
    cr2ph: np.ndarray
        crossrange-to-phase for each point
    ampl_ts: np.ndarray
        amplitude timeseries
    bkps_stm: np.ndarray
        Breakpoints stm
    mean_ampl_sd: np.ndarray
        Single difference mean amplitude
    sigma_ampl_sd: np.ndarray
        Single difference mean amplitude standard deviation
    mad_ampl_sd: np.ndarray
        Median absolute deviation of the single difference amplitude
    median_ampl_sd: np.ndarray
        Median single difference amplitude
    x_coordinates: np.ndarray
        x coordinates (RD or longitude)
    y_coordinates: np.ndarray
        y coordinates (RD or latitude)
    coordinate_type: Literal["euclidean", "geographic"] = "euclidean"
        Whether the provided coordinates are Euclidean (such as RD) or geographic (such as lon/lat)

    Returns
    -------
    dict
        Results of the accepted network
    list
        Reference point
    list
        Results of the accepted arcs
    """
    # At the start, no arcs are tested yet, so we don't have failed_arcs, and noisy_arcs
    failed_arcs = []
    noisy_arcs = []

    # Get an ordered list of all potential arcs with a buffer area
    arcs_search_area, _, quality_dict_arcs = get_ordered_arcs(
        stm, x_ref_search, y_ref_search, buffer_radius_ref, dist_to_quality, N_max_arcs, nad_max
    )

    network_meets_requirements = False

    while not network_meets_requirements:
        # Create a first network based on the ranked quality and leaving out the failed_arcs and noisy_arcs
        ref_pnt, arcs_updated_network, ref_pnt_initial, arcs_initial_network = construct_control_network(
            arcs_search_area,
            quality_dict_arcs,
            failed_arcs,
            N_top,
            N_batch,
            deg_threshold,
            min_nodes,
            min_redundancy,
            False,
        )

        print("Computing the solutions for the arcs")
        # Test whether solutions for all arcs can be found (sometimes it happens that because of the complex
        # functions no solutions can be found)
        results_initial_control_network_v1 = arc_estimation_control_network(
            arcs_updated_network,
            bounds,
            m2ph,
            nr_max_iter_control,
            years,
            dates,
            temperature,
            sd_complex,
            slc_quality,
            cr2ph,
            ampl_ts,
            bkps_stm,
            mean_ampl_sd,
            sigma_ampl_sd,
            mad_ampl_sd,
            median_ampl_sd,
            x_coordinates,
            y_coordinates,
            coordinate_type,
        )

        # Save the failed arcs to an array (such that they are not taken into account any more)
        succeeded_arcs = [
            arc
            for arc, succeeded in zip(
                arcs_updated_network, results_initial_control_network_v1["succeeded_arcs"], strict=True
            )
            if not np.isnan(succeeded).any()
        ]

        failed_arcs_temp = [
            arc
            for arc, succeeded in zip(
                arcs_updated_network, results_initial_control_network_v1["succeeded_arcs"], strict=True
            )
            if np.isnan(succeeded).any() or not succeeded.any()
        ]
        failed_arcs = list(set(failed_arcs).union(failed_arcs_temp))
        print(f"Failed arcs {failed_arcs}")

        # Remove the failed arcs from the dictionary with all the results (as the estimated parameters etc)
        results_initial_control_network_v2 = {}  # Make a new dictionary where we will not save
        # the results of the failed_arcs

        # Find the indices where the arcs are saved where no solution was found.
        # These arcs have nan values
        valid_indices = ~np.isnan(results_initial_control_network_v1["succeeded_arcs"]).any(axis=1)

        # Filter all variables in the dictionary and only copy the values for the arcs where we found a solution
        for key, value in results_initial_control_network_v1.items():
            if (
                isinstance(value, np.ndarray)
                and value.shape[0] == results_initial_control_network_v1["succeeded_arcs"].shape[0]
            ):
                results_initial_control_network_v2[key] = value[valid_indices]
            else:
                results_initial_control_network_v2[key] = value  # Keep values that are not row-based unchanged

        print("Check if the network with the solved arcs still meets our requirements")

        # Since we have 'failed_arcs', where no solution was found, the new network need to be tested
        # It can happen that we have isolated points,
        # or that the average degree of the network is not high enough anymore
        network_check, arcs_updated_network, ref_pnt = test_succeeded_arcs_control_network(
            succeeded_arcs, quality_dict_arcs, deg_threshold, min_nodes, min_redundancy, visualize_network
        )

        if network_check == 0:
            print("Network fails requirements, starting again")

        if network_check == 1:
            print("Network meets requirements")
            # After the last test, the isoltated arcs are removed (so they are still in
            # 'succeeded_arcs' and in the dictionary)
            # And these arcs need to be removed from the dictionary
            missing_indices = [
                idx for idx, arc in enumerate(succeeded_arcs) if arc not in arcs_updated_network
            ]  # missing indices are the arcs that were isolated and removed. But they are still in the
            # dictionary so there they need to be removed as well

            results_control_network_temp = {}
            for key, value in results_initial_control_network_v2.items():
                if isinstance(value, np.ndarray):
                    if value.ndim in {1, 2}:  # Both 1D and 2D arrays here
                        results_control_network_temp[key] = np.delete(value, missing_indices, axis=0)
                    else:
                        results_control_network_temp[key] = value
                else:
                    results_control_network_temp[key] = value

            # Now we have a network that fullfills requirements but there might be noisy arcs
            # We compute solutions for all arcs and computed RMSE
            print("Calculate whether there are arcs where the solution that we found is noisy")

            est_displ_phase = (
                results_control_network_temp["unwrap_phases_arc"]
                - results_control_network_temp["estimated_cross_range_phase"]
                - results_control_network_temp["estimated_thermal_phase"]
            )
            displ_phase = results_control_network_temp["estimated_displ_phase"]

            residues_per_arc = est_displ_phase - displ_phase
            sigma_post_arc = np.std(residues_per_arc, axis=1)
            mean_sigma_prior_arc = np.mean(results_control_network_temp["sigma_phases_arc"], axis=1)

            idx_bad_arcs = np.where(sigma_post_arc >= sigma_post_over_sigma_prior * mean_sigma_prior_arc)[0]
            idx_good_arcs = np.where(sigma_post_arc < sigma_post_over_sigma_prior * mean_sigma_prior_arc)[0]

            bad_arcs = results_control_network_temp["succeeded_arcs"][idx_bad_arcs]
            good_arcs = results_control_network_temp["succeeded_arcs"][idx_good_arcs].astype(int)
            good_arcs = [tuple(row) for row in good_arcs]  # Change the output to a list

            print(f"The value for sigma_post_over_prior is {sigma_post_over_sigma_prior}")

            print("we removed acs")

            print("The bad arcs are")
            print(idx_bad_arcs)
            print("The good arcs are")
            print(idx_good_arcs)

            # Add the noisy arcs to the list with failed_arcs
            failed_arcs = list(set(failed_arcs).union([tuple(row.astype(int)) for row in bad_arcs]))
            noisy_arcs = list(set(noisy_arcs).union([tuple(row.astype(int)) for row in bad_arcs]))

            # Check whether the network without the noisy arcs still meets requirements
            print("Check whether the network stil meets the requirements, even after the removal of bad arcs")
            network_check_good, arcs_updated_network, ref_pnt = test_succeeded_arcs_control_network(
                good_arcs, quality_dict_arcs, deg_threshold, min_nodes, min_redundancy, visualize_network=False
            )

            if network_check_good == 0:
                print("Network does not meet requirements, start over")

            if network_check_good == 1:
                print("We are happy! The network consisting of the good arcs fulfills the requirements.")

                # If there are noisy arcs, they need to be removed from the final network
                # Convert noisy_arcs to set for quicker lookup
                noisy_arcs_set = {tuple(map(float, arc)) for arc in noisy_arcs}

                # Get the current succeeded_arcs
                succeeded_arcs = results_control_network_temp["succeeded_arcs"]

                # Deterimine which arcs to keep (so NOT in noisy_arcs_set)
                valid_indices = np.array([tuple(row) not in noisy_arcs_set for row in succeeded_arcs])

                # New dict with only valid rows
                results_control_network = {}

                for key, value in results_control_network_temp.items():
                    if isinstance(value, np.ndarray) and value.shape[0] == succeeded_arcs.shape[0]:
                        results_control_network[key] = value[valid_indices]
                    else:
                        results_control_network[key] = value  # Leave unrelated variables unchanged

                network_meets_requirements = True  # Stop the loop

        if not network_meets_requirements:
            # If the network doesn't meet the requirements, start over
            print("Network does not meet the requirements. Recomputing the network with updated failed_arcs...")

    return results_control_network, ref_pnt, arcs_updated_network
