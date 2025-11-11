"""This script runs through the Confidence-Optimized Robust Geodetic Network (CORG Network) formation"""

import datetime as dt
from pathlib import Path

import numpy as np
import xarray as xr

from depsi.io import read_weather_data
from depsi.network import construct_control_network_test_arcs
from depsi.network_adjustment import adjust_full_corg_control_network
from depsi.utils import npdatetime64_to_datetime
from depsi.viewing_geometry import add_local_viewing_geometry

# the original STM is output from the script script_sidelobe_detection.py
stm_original = (
    "/Users/sanvandiepen/PycharmProjects/workingEnvironment2/test_zarr/nl_amsterdam_s1_dsc_t037_stm_nosl.zarr"
)
stm_viewing_save_path = (
    "/Users/sanvandiepen/PycharmProjects/workingEnvironment2/test_zarr/nl_amsterdam_s1_dsc_t037_stm_nosl_view.zarr"
)  # this one is there to speed up the testing
stm_save_path = (
    "/Users/sanvandiepen/PycharmProjects/workingEnvironment2/test_zarr/nl_amsterdam_s1_dsc_t037_stm_CORG.zarr"
)  # This is where the control network will be saved
knmi_file_path = "/Users/sanvandiepen/PycharmProjects/workingEnvironment2/test_zarr/etmgeg_240.txt"

orbit_cfg_file = Path(
    "/Users/sanvandiepen/PycharmProjects/workingEnvironment2/GitHub_repos/DePSI_group/config/drama/S1_XTI.cfg"
).expanduser()

# User settings
# orbit identifier
orbit_identifier = "s1_dsc_t037"

# Input for the OMT in the network
alpha = 0.05

# Input regarding displacement estimates from phases
speed_light = 299792458  # m/s
frequency = 5.404*10**9  # Hz
lam_meters = speed_light * 1 / frequency
ph2m = lam_meters / (4*np.pi)
m2ph = -(4*np.pi)/lam_meters

# Input variables for the bounds needed in the arc parameter function

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

# Stochastic model
vcm_complex_method = 'mad_median'
dist_to_quality = 0.1/1000  # Normal value , rad/m
sigma_post_over_sigma_prior = 4
partition_quality_label = "partition_quality_nmad_2sigma"

# Geolocation
x_crd_label = "rd_x"
y_crd_label = "rd_y"
coordinate_type = "euclidean"

# Network creation settings
N_max_arcs = 100  # Nr of arcs that we will analyse, otherwise we need to load a very big dataset everytime
N_top = 35  # The Nr of arcs where we start
nad_max = 0.3  # Maximum NAD value for a point to be considered in the control network, otherwise the distance matrices
# to be computed are large
buffer_radius_ref = 700  # Value I work normally with 700

# Initialization
N_batch = 10

# Final limits
min_nodes = 10
min_redundancy = 2.1  # 2.1 is strict value, 1.7 works as well
deg_threshold = 1

# iterative model solution
nr_max_iter_control = 120
visualize_network = False

# orbit settings
orbit_mode = "IWS"
orbit_resolution = 0.01

# Adjustment thresholds
criteria = {
    'cross_range': True,  # Sigma of cross-range height
    'thermal': False,     # Sigma of thermal components
    'displacement': True  # Sigma of displacements
}

criteria_thresholds = {
    'sigma_cross_range': 159, # radians
    'sigma_thermal': 0.005, # m/y/K (ish)
    'sigma_displacement': 0.5  # rad
}

min_points_control_network = 9
min_points_before_estimate = 2
min_degree = 1
max_iter_adjustment = 30
correct_network = True

#
# Calculations
try:
    stm = xr.open_zarr(stm_viewing_save_path)  # if the viewing geometry exists it acts as a checkpoint (for testing)
except IOError:
    # it doesn't exist, so we calculate it
    # open the STM
    stm = xr.open_zarr(stm_original)

    # Add temperature information to it for the thermal component
    timestamps = [npdatetime64_to_datetime(date, tz_aware=False) for date in stm.time.values]
    temperatures = read_weather_data(knmi_file_path, timestamps, requested_data_columns=("TG", ))
    temp = np.array([temperatures[day]["TG"] for day in timestamps])
    stm = stm.assign(temperature=(["time"], temp))

    # add the incidence angle
    print(f"Time: {dt.datetime.now().strftime('%H:%M:%S')}")
    stm = add_local_viewing_geometry(stm, orbit_cfg_file, orbit_resolution, orbit_mode, orbit_identifier)

    print(f"Time: {dt.datetime.now().strftime('%H:%M:%S')}")

    # add the crossrange
    cr2ph = stm["sd_h2ph"] * np.sin(np.radians(stm["local_incidence_angle"].data[:, np.newaxis]))
    stm = stm.assign({"sd_cr2ph": (["space", "time"], cr2ph.data)})

    stm.to_zarr(stm_viewing_save_path, mode='w')

    print(f"Time: {dt.datetime.now().strftime('%H:%M:%S')}")

# Find the center of the AoI
x_min = np.min(stm[x_crd_label].values)
x_max = np.max(stm[x_crd_label].values)
y_min = np.min(stm[y_crd_label].values)
y_max = np.max(stm[y_crd_label].values)
x_center = x_min + (x_max - x_min) * 0.5
y_center = y_min + (y_max - y_min) * 0.5

# construct the CORG control network
results_control_network, ref_pnt, arcs_updated_network = construct_control_network_test_arcs(
    x_center,
    y_center,
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
    stm["years_since_first_img"].values,
    stm["days_since_first_img"].values,
    stm["temperature"].values,
    stm["sd_complex"].values,
    stm[partition_quality_label].values,
    stm["sd_cr2ph"].values,
    stm["amplitude"].values,
    stm["breakpoints"].values,
    stm["partition_sd_amplitude_mean"].values,
    stm["partition_sd_amplitude_sigma"].values,
    stm["partition_sd_mad"].values,
    stm["partition_sd_amplitude_median"].values,
    stm[x_crd_label].values,
    stm[y_crd_label].values,
    stm["full_ts_nad"].values,
    coordinate_type,
)

stm_control_network = adjust_full_corg_control_network(
    stm,
    partition_quality_label,
    results_control_network,
    max_iter_adjustment,
    alpha,
    criteria,
    criteria_thresholds,
    min_points_before_estimate,
    min_points_control_network,
    min_degree,
    ref_pnt,
    m2ph,
    correct_network
)

stm_control_network.to_zarr(stm_save_path)
