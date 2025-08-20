"""this script is an example of how to connect points to a control network in a distributed manner"""

from dask.delayed import delayed

from depsi.network_adjustment import connect_point_to_control_network


@delayed  # by delaying this function, we can submit a graph to the cluster instead of doing it consecutively
def connect_point_delayed(stm_1_pt, control_network_stm_init, control_network_stm_solved, ref_pnt_idx, dist_to_quality,
                          nr_conn, sigma_post_over_sigma_prior_connection,
                          alpha, bounds,
                          m2ph,
                          max_iter_per_arc,
                          correct_epochs_arc = 1):
    output = connect_point_to_control_network(
        stm_1_pt,
        control_network_stm_init,
        control_network_stm_solved,
        ref_pnt_idx,
        dist_to_quality,
        nr_conn,
        sigma_post_over_sigma_prior_connection,
        alpha, bounds,
        m2ph,
        max_iter_per_arc,
        correct_epochs_arc=correct_epochs_arc
    )
    return output



