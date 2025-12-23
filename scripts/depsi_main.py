"""
This script runs the original PS DePSI workflow as originally implemented in Matlab.
"""
from datetime import datetime

import numpy as np
import sarxarray

from depsi.arc_estimation import periodogram
from depsi.classification import network_stm_selection, ps_selection
from depsi.io import read_slc_stack
from depsi.network import form_network
from depsi.point_quality import detect_side_lobes
from depsi.transformations import radar_to_latlonh
from depsi.utils import add_stm_time_deltas, crop_slc_spacetime, stm_compute_single_time_differences


# ############## INPUT VARIABLES

slc_path = '/Users/sanvandiepen/PycharmProjects/workingEnvironment2/test_zarr/nl_amsterdam_s1_dsc_t037_base.zarr'

# Crop in space
aoi_file = '/Users/sanvandiepen/PycharmProjects/workingEnvironment2/test_zarr/shape/nl_amsterdam_shape.shp'

# metadata
metadata_path = '/Users/sanvandiepen/PycharmProjects/workingEnvironment2/test_zarr/nl_amsterdam_s1_dsc_t037.res'

# STM save path
stm_save_path = '/Users/sanvandiepen/PycharmProjects/workingEnvironment2/test_zarr/nl_amsterdam_s1_dsc_t037_stm.zarr'

# Crop in time
first_date = datetime(2020, 1, 1)
last_date = datetime(2025, 9, 1)

mother_epoch = "auto"

# PS selection method
ps_selection_method = "nmad"
threshold = 0.15
chunks_ps_selection = 1000
start_date_ps_selection = None
end_date_ps_selection = None

# Sidelobes
max_pixel_dist = 2
min_correlation = 0.90

# Network formation
network_crs = "radar"
network_x_crds = "azimuth"
network_y_crds = "range"
min_point_distance = 20  # pixels if radar, otherwise units of the EPSG code in network_crs
max_arc_length = 0.001  # degrees (lat/lon)
network_formation_method="redundant"
min_network_links = 16
network_partition_number = 8
network_dphase_method = "subtract"
min_periodogram_iterations = 10


# 1. Project setup
print("Reading SLC stack...")
slcs = read_slc_stack(slc_path)


# 2. Scatterer selection
print("Start cropping...")
cropped_slcs = crop_slc_spacetime(
    slcs,
    aoi_filename=aoi_file,
    start_date=first_date,
    end_date=last_date
)

# ######## POINT SELECTION WITH THE PARAMETERS ABOVE ############
print("Starting PS selection...")
stm = ps_selection(
    cropped_slcs,
    method=ps_selection_method,
    threshold=threshold,
    ps_selection_start_date=start_date_ps_selection,
    ps_selection_end_date=end_date_ps_selection,
    output_chunks=chunks_ps_selection,
    mem_persist=False
)
print(f"Selected {len(stm.space)} PS.")

print("Computing single differences and temporal baseline...")
stm = stm_compute_single_time_differences(
    stm,
    single_difference_mother=mother_epoch
)
stm = add_stm_time_deltas(stm)
stm = stm.assign({"temporal_baseline": (
    ["time"],
    stm.years_since_first_img.values - stm.years_since_first_img.sel(time=stm.ps_sd_mother).values
)})

# sidelobe removal
print("Removing sidelobes...")
side_lobes_array, _ = detect_side_lobes(stm, max_pixel_dist, min_correlation, "sd_complex", "sd_amplitude_unnormalized")

mask_sidelobes = np.ones(len(stm.space), dtype=bool)
mask_sidelobes[side_lobes_array] = False

stm = stm.isel(space=mask_sidelobes)

print(f"Removed {np.sum(~mask_sidelobes)} sidelobes from the STM, {len(stm.space)} remain.")

# 6a. Geocoding
metadata = sarxarray.read_metadata(metadata_path, driver="doris5")

latlonh = radar_to_latlonh(
    azimuth_coords=stm.azimuth.values,
    range_coords=stm.range.values,
    elevation=stm.h.values,
    metadata=metadata,
)

# 3a. Network construction
mother_epoch_index = np.where(stm.time.values == stm.sel(time=stm.ps_sd_mother).time.values)[0][0]
non_mother = [True] * len(stm.time.values)
non_mother[mother_epoch_index] = False
stm_without_mother_epoch = stm.isel(time=non_mother)

stm_network_pnts = network_stm_selection(
    stm=stm,
    min_dist=min_point_distance,
    crs=network_crs,
    azimuth_spacing=metadata["azimuth_pixel_spacing"],
    range_spacing=metadata["range_pixel_spacing"],
    sortby_var="time_selection_" + ps_selection_method,
    x_var=network_x_crds,
    y_var=network_y_crds,
    include_index=None
)

stm_network_arcs = form_network(
    stm_network_pnts,
    key_phase='sd_phase',
    key_h2ph='sd_h2ph',
    key_Btemp='temporal_baseline',
    max_length=max_arc_length,
    key_xcrds="lon",
    key_ycrds="lat",
    network_method=network_formation_method,
    min_links=min_network_links,
    num_partitions=network_partition_number,
    dphase_method=network_dphase_method,
)

_, ambiguities, _, _, ens_coh = periodogram(
    stm_network_arcs,
    key_dphase='sd_phase',
    key_h2ph='sd_h2ph',
    key_Btemp='temporal_baseline',
    wavelength=metadata["wavelength"],
    min_steps=min_periodogram_iterations,
            )
stm_network_arcs["ambiguities"] = ambiguities
stm_network_arcs["temp_coh"] = ens_coh

stm_network_arcs = stm_network_arcs.compute()


# 4. Atmosphere estimation

# 3b. Network construction

# 5. Densification

# 6b. Geocoding

# 7. STC Calculation

# 8. Scatterer filtering

# 9. Output generation
