"""this script is an example of how to connect points to a control network in a distributed manner"""

import dask
import numpy as np
import xarray as xr
from dask.delayed import delayed

from depsi.network import find_points_within_buffer
from depsi.network_adjustment import connect_point_to_control_network


stm_initial_path = (
    "/Users/sanvandiepen/PycharmProjects/workingEnvironment2/test_zarr/nl_amsterdam_s1_dsc_t037_stm_nosl_view.zarr"
)  # this one is there to speed up the testing
stm_control_path = (
    "/Users/sanvandiepen/PycharmProjects/workingEnvironment2/test_zarr/nl_amsterdam_s1_dsc_t037_stm_CORG.zarr"
)  # This is where the control network will be saved
stm_save_path = (
    "/Users/sanvandiepen/PycharmProjects/workingEnvironment2/test_zarr/nl_amsterdam_s1_dsc_t037_stm_CORG_dens.zarr"
)  # This is where the control network will be saved

# Geolocation
x_crd_label = "rd_x"
y_crd_label = "rd_y"
coordinate_type = "euclidean"
max_point_distance = 50  # meter

# Stochastic model
dist_to_quality = 0.1/1000  # Normal value , rad/m
sigma_post_over_sigma_prior = 4
partition_quality_label = "partition_quality_nmad_2sigma"
alpha = 0.05
min_n_connections = 3

# Input regarding displacement estimates from phases
speed_light = 299792458  # m/s
frequency = 5.404*10**9  # Hz
lam_meters = speed_light * 1 / frequency
ph2m = lam_meters / (4*np.pi)
m2ph = -(4*np.pi)/lam_meters

# Iteration control
nr_max_iter_control = 120
correct_epoch_arcs = 1

# Bounds are for instance defined on the maximum and minimum cross range distance or instantaneous velocity.
A_upper = 10*10**14  # amplitude
A_lower = 10*10**6
a_upper = 5 * -m2ph / 1000  # start point of displacement polynomial  # offset
a_lower = -5 * -m2ph / 1000  # start point of displacement polynomial
b_upper = 5 * -m2ph / 1000  # instanteneous velocity of the displacement polynomial  # snelheid
b_lower = -5 * -m2ph / 1000
c_upper = 1 * -m2ph / 1000  # instanteneous acceleration of the displacement polynomial  # versnelling
c_lower = -1 * -m2ph / 1000

CR_upper = 200 * -m2ph  # Maximum cross range in radians
CR_lower = -200 * -m2ph
exp_upper = 1.5 * -m2ph / 1000
exp_lower = -1.5 * -m2ph / 1000  # Maximum thermal expansion

bounds = (
    A_lower, a_lower, b_lower, c_lower, CR_lower, exp_lower, A_upper, a_upper, b_upper, c_upper, CR_upper, exp_upper
)


@delayed  # by delaying this function, we can submit a graph to the cluster instead of doing it consecutively
def connect_point_delayed(
        stm_1_pt,
        partition_quality_label,
        x_crd_label,
        y_crd_label,
        coordinate_type,
        stm_control_network,
        stm_control_network_solved,
        ref_point_index,
        dist_to_quality,
        min_n_connections,
        sigma_post_over_sigma_prior,
        alpha,
        bounds,
        m2ph,
        nr_max_iter_control,
        correct_epoch_arcs
):
    stm_1_pt = stm_1_pt.compute().squeeze()
    try:
        output = connect_point_to_control_network(
            stm_1_pt,
            partition_quality_label,
            x_crd_label,
            y_crd_label,
            coordinate_type,
            stm_control_network,
            stm_control_network_solved,
            ref_point_index,
            dist_to_quality,
            min_n_connections,
            sigma_post_over_sigma_prior,
            alpha,
            bounds,
            m2ph,
            nr_max_iter_control,
            correct_epochs_arc=correct_epoch_arcs
        )
    except ValueError:
        output = (None, None, None, None, 0)
    return output


# computation
base_stm = xr.open_zarr(stm_initial_path)

control_network_stm = xr.open_zarr(stm_control_path)

control_points = control_network_stm.pnt_idx.values

control_stm = base_stm.sel(space=control_points)

points_within_buffer = find_points_within_buffer(
    x_pnts=np.mean(control_stm[x_crd_label].values),
    y_pnts=np.mean(control_stm[y_crd_label].values),
    x_coords=base_stm[x_crd_label].values,
    y_coords=base_stm[y_crd_label].values,
    coordinate_type=coordinate_type,
    buffer_radius=max_point_distance,
)

# remove the control network from the points within buffer since they already are in the network
points_within_buffer = np.setdiff1d(points_within_buffer, control_points)

print(f"Identified {len(points_within_buffer)} points within {max_point_distance} m of control network")

within_buffer_stm = base_stm.sel(space=points_within_buffer)

stm_control_network = control_stm.compute().squeeze()
stm_control_network_solved = control_network_stm.compute().squeeze()
ref_point_index = stm_control_network_solved.ref_pnt.values[0]

delayed_results = []

for idx in points_within_buffer:
    stm_1_pt = within_buffer_stm.sel(space=idx)  # keep it delayed here, computation will happen in the delayed
    # execution

    delayed_results.append(connect_point_delayed(
        stm_1_pt,
        partition_quality_label,
        x_crd_label,
        y_crd_label,
        coordinate_type,
        stm_control_network,
        stm_control_network_solved,
        ref_point_index,
        dist_to_quality,
        min_n_connections,
        sigma_post_over_sigma_prior,
        alpha,
        bounds,
        m2ph,
        nr_max_iter_control,
        correct_epoch_arcs=correct_epoch_arcs
    ))

# create the full combined STM
combined_stm = stm_control_network_solved.copy()  # first the control network
results = dask.compute(delayed_results, traverse=True)[0]  # then compute all results in a parallellized way
# then merge the results
for res in results:
    if res[-1] == 1:  # estimation successful  (if 0 we ignore it since it failed)
        combined_stm = xr.merge([combined_stm, res[0]])

# finally, save
combined_stm.to_zarr(stm_save_path)
