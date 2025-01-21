from datetime import datetime

from depsi.io import read_slc_stack
from depsi.utils import crop_slc_spacetime, add_stm_time_deltas, project_stm_coordinates
from depsi.classification import ps_selection
from depsi.point_quality import detect_outliers_stm

# ############## INPUT VARIABLES

slc_path = '/Users/sanvandiepen/PycharmProjects/workingEnvironment2/test_zarr/nl_amsterdam_s1_dsc_t110.zarr'

# Crop in space
aoi_file = '/Users/sanvandiepen/PycharmProjects/workingEnvironment2/test_zarr/shape/nl_amsterdam_shape.shp'

# STM save path
stm_save_path = '/Users/sanvandiepen/PycharmProjects/workingEnvironment2/test_zarr/nl_amsterdam_s1_dsc_t110_stm.zarr'

# Crop in time
last_date = datetime(2025, 9, 1)
first_date = datetime(2014, 1, 1)

# PS Selection based on initialization
start_date_ps_selection = datetime(2017, 1, 1)
initialization_length = 50

# Recalibrated NAD and NMAD settings
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
slcs = read_slc_stack(slc_path)

# ############ CROP ##########
cropped_slcs = crop_slc_spacetime(slcs,
                                  aoi_filename=aoi_file,
                                  start_date=first_date,
                                  end_date=last_date)

# ######## POINT SELECTION WITH THE PARAMETERS ABOVE ############
stm = ps_selection(cropped_slcs,
                   method=ps_selection_method,
                   threshold=threshold,
                   ps_selection_start_date=start_date_ps_selection,
                   ps_selection_end_date=initialization_length,
                   output_chunks=chunks_ps_selection,
                   mem_persist=False,
                   recalibration_jump_size=recalibration_jump_size,
                   do_partitioning=do_ps_partitioning,
                   partitioning_kwargs={"db_partitioning": ps_db_partitioning,
                                        "search_method": ps_partitioning_search_method,
                                        "cost_function": ps_partitioning_cost_function,
                                        "min_obs_partition": ps_min_obs_partition},
                   single_difference_mother=ps_mother_epoch_sd
                   )

# Add RD coordinates to the STM
stm = project_stm_coordinates(stm, "RD")

# Add time deltas to the STM
stm = add_stm_time_deltas(stm)

# Do outlier detection
if do_ps_outlier_detection:
    stm = detect_outliers_stm(stm, db_outlier_detection=ps_outlier_detection_db, window_size=ps_window_size_outliers,
                              n_sigma=ps_n_sigma_outliers)

# SAVE
stm.to_zarr(stm_save_path, mode='w')
