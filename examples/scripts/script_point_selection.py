from datetime import datetime

from depsi.io import read_slc_stack
from depsi.utils import crop_slc_spacetime
from depsi.classification import ps_selection

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
start_date_ps_selection = datetime(2014, 1, 1)
initialization_length = 50

# Recalibrated NAD and NMAD settings
recalibration_jump_size = 10

# PS selection method
nmad_max = 0.25
# nad_max = 0.425 #0.3
ps_selection_method = "nmad"
chunks_ps_selection = 1000

# Input variables for the outlier detection
do_outlier_detection = True
window_size_outliers = 15
outlier_detection_db = True
n_sigma_outliers = 3

# Input variables for the partitioning
do_partitioning = True
search_method = 'pelt'
cost_function = 'l2'
db_partitioning = False
min_obs_partition = 27

# Compute temporal differences
mother_epoch_sd = '20190806'

# ## FUNCTIONALITY
# ############ LOAD THE SLCS FROM ZARR ###############
slcs = read_slc_stack(slc_path)

# ############ CROP ##########
cropped_slcs = crop_slc_spacetime(slcs,
                                  aoi_filename=aoi_file,
                                  start_date=first_date,
                                  end_date=last_date)

# ######## POINT SELECTION WITH THE PARAMETERS ABOVE ############
stm_nmad = ps_selection(cropped_slcs,
                        threshold=nmad_max,
                        method=ps_selection_method,
                        ps_selection_start_date=start_date_ps_selection,
                        ps_selection_end_date=initialization_length,
                        output_chunks=chunks_ps_selection,
                        mem_persist=False,
                        recalibration_jump_size=recalibration_jump_size,
                        do_rd_coordinate_conversion=True,
                        do_partitioning=do_partitioning,
                        partitioning_kwargs={"db_partitioning": db_partitioning,
                                             "search_method": search_method,
                                             "cost_function": cost_function,
                                             "min_obs_partition": min_obs_partition},
                        do_outlier_detection=True,
                        outlier_detection_kwargs={"db_outlier_detection": outlier_detection_db,
                                                  "window_size": window_size_outliers,
                                                  "n_sigma": n_sigma_outliers},
                        single_difference_mother=mother_epoch_sd
                        )

# SAVE
stm_nmad.to_zarr(stm_save_path, mode='w')
