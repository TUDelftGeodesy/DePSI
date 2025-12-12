from datetime import datetime

from depsi.io import read_slc_stack
from depsi.utils import crop_slc_spacetime, add_stm_time_deltas, project_stm_coordinates, \
    stm_compute_single_time_differences
from depsi.classification import ps_selection
from depsi.point_quality import detect_outliers_stm, stm_partitioning, stm_add_incremental_recal_nad_nmad

# ############## INPUT VARIABLES

slc_path = '/Users/sanvandiepen/PycharmProjects/workingEnvironment2/test_zarr/nl_amsterdam_s1_dsc_t037_base.zarr'

# Crop in space
aoi_file = '/Users/sanvandiepen/PycharmProjects/workingEnvironment2/test_zarr/shape/nl_amsterdam_shape.shp'

# STM save path
stm_save_path = '/Users/sanvandiepen/PycharmProjects/workingEnvironment2/test_zarr/nl_amsterdam_s1_dsc_t037_stm.zarr'

# Crop in time
first_date = datetime(2014, 1, 1)
last_date = datetime(2025, 9, 1)

# PS Selection based on initialization
start_date_ps_selection = datetime(2017, 1, 1)
initialization_length = 50

# Recalibrated NAD and NMAD settings
increment_mode = "recalibration"
recalibration_jump_size = 10

# PS selection method
ps_selection_method = "nmad"
threshold = 0.25
chunks_ps_selection = 1000

# Input variables for the outlier detection
do_ps_outlier_detection = True
ps_window_size_outliers = 15
ps_outlier_detection_db = True
ps_n_sigma_outliers = 3

# Input variables for the partitioning
do_ps_partitioning = True
ps_partitioning_search_method = 'pelt'
ps_partitioning_cost_function = 'l2'
ps_db_partitioning = False
ps_min_obs_partition = 27

# Compute temporal differences
ps_mother_epoch_sd = '20190806'

# ## FUNCTIONALITY
# ############ LOAD THE SLCS FROM ZARR ###############
print("Reading SLC stack...")
slcs = read_slc_stack(slc_path)

# ############ CROP ##########
print("Start cropping...")
cropped_slcs = crop_slc_spacetime(slcs,
                                  aoi_filename=aoi_file,
                                  start_date=first_date,
                                  end_date=last_date)

# ######## POINT SELECTION WITH THE PARAMETERS ABOVE ############
print("Starting PS selection...")
stm = ps_selection(cropped_slcs,
                   method=ps_selection_method,
                   threshold=threshold,
                   ps_selection_start_date=start_date_ps_selection,
                   ps_selection_end_date=initialization_length,
                   output_chunks=chunks_ps_selection,
                   mem_persist=False)

# Add the incremental or recalibration NAD / NMAD to the STM
print("Add Inc/Recal NAD/NMAD...")
stm = stm_add_incremental_recal_nad_nmad(stm,
                                         mode=increment_mode,
                                         method=ps_selection_method,
                                         recalibration_jump_size=recalibration_jump_size)

# Add RD coordinates to the STM
print("Add RD coordinates...")
stm = project_stm_coordinates(stm, "RD")

# Add time deltas to the STM
print("Add time deltas...")
stm = add_stm_time_deltas(stm)

# Add single differences to the STM
print("Add Single time differences...")
stm = stm_compute_single_time_differences(stm, ps_mother_epoch_sd)

if do_ps_partitioning:
    print("Start normal partitions...")
    stm = stm_partitioning(stm,
                           db_partitioning=ps_db_partitioning,
                           search_method=ps_partitioning_search_method,
                           cost_model=ps_partitioning_cost_function,
                           min_partition_size=ps_min_obs_partition,
                           amplitude_variable_name="amplitude",
                           output_variable_prefix="partition",
                           output_variables=("nmad", "nad", "quality_nmad_2sigma"))

    print("Start Single Difference partitions...")
    stm = stm_partitioning(stm,
                           db_partitioning=ps_db_partitioning,
                           search_method=ps_partitioning_search_method,
                           cost_model=ps_partitioning_cost_function,
                           min_partition_size=ps_min_obs_partition,
                           amplitude_variable_name="sd_amplitude_unnormalized",
                           output_variable_prefix="partition_sd",
                           output_variables=("mad", "amplitude_sigma", "amplitude_mean", "amplitude_median"))

# Rechunk to prevent inconsistent chunks
stm = stm.chunk({"time": -1, "space": chunks_ps_selection})

# Do outlier detection
if do_ps_outlier_detection:
    print("Start outlier detection...")
    stm = detect_outliers_stm(stm, db_outlier_detection=ps_outlier_detection_db, window_size=ps_window_size_outliers,
                              n_sigma=ps_n_sigma_outliers)

# Save
print("Saving...")
stm.to_zarr(stm_save_path, mode='w')
print("Done!")
