from depsi.config import ConfigGenerateSTM
import socket
from dask.distributed import Client
from dask_jobqueue import SLURMCluster

from depsi.io import read_slc_stack
from depsi.utils import add_stm_time_deltas, project_stm_coordinates, \
    stm_compute_single_time_differences
from depsi.classification import ps_selection
from depsi.point_quality import detect_outliers_stm, stm_partitioning, stm_add_incremental_recal_nad_nmad

# Get config path either from the environment variable or default config
# diretory or from current directory
config_path = "./generate_stm_extended_s1_asc_t088.yml"

# Read the config file
cfg = ConfigGenerateSTM.from_yaml(config_path)
print(f"Using config file {config_path}.")
print(cfg)
stop_here

## FUNCTIONALITY
# Start cluster

def get_free_port():
    """Get a non-occupied port number."""
    sock = socket.socket()
    sock.bind(("", 0))  # Bind a port, it will be busy now
    freesock = sock.getsockname()[1]  # get the port number
    sock.close()  # Free the port, so it can be used later
    return freesock

N_WORKERS = 10 # Manual input: number of workers to spin-up
FREE_SOCKET = get_free_port() # Get a free port
cluster = SLURMCluster(
    name="dask-worker",  # Name of the Slurm job
    queue="normal", # Name of the node partition on your SLURM system
    cores=4, # Number of cores per worker
    memory="30 GB",  # Total amount of memory per worker
    processes=1,  # Number of Python processes per worker
    walltime="4-00:00:00",  # Reserve each worker for X hour
    scheduler_options={"dashboard_address": f":{FREE_SOCKET}"},  # Host Dashboard in a free socket
)

cluster.scale(jobs=N_WORKERS)
client = Client(cluster)


# ############ LOAD THE SLCS FROM ZARR ###############
slcs = read_slc_stack(cfg.paths['slc_file'])
print("Finished reading.")

# ######## POINT SELECTION WITH THE PARAMETERS ABOVE ############
stm = ps_selection(
    slcs,
    method=cfg.ps_selection.method,
    threshold=cfg.ps_selection.threshold,
    output_chunks=cfg.ps_selection.output_chunks,
)
print(f"Selected {stm.sizes['space']} PS points")

# Add the incremental or recalibration NAD / NMAD to the STM
stm = stm_add_incremental_recal_nad_nmad(
    stm,
    mode=cfg.stm_add_incremental.mode,
    method=cfg.stm_add_incremental.method,
    recalibration_jump_size=cfg.stm_add_incremental.recalibration_jump_size,
)
print("Added incremental/recal NAD/NMAD")

# Add RD coordinates to the STM
stm = project_stm_coordinates(stm, cfg.project_stm_coordinates_projection)
print("Projected coordinates")

# Add time deltas to the STM
stm = add_stm_time_deltas(stm)
print("Calculated time deltas")

# Add single differences to the STM
stm = stm_compute_single_time_differences(
    stm,
    cfg.stm_single_difference_mother,
)
print("Calculated single differences")

if cfg.do_stm_partitioning:
    partitioning_cfg = cfg.stm_partitioning
    partitioning_normal = cfg.stm_partitioning_normal
    stm = stm_partitioning(
        stm,
        db_partitioning=partitioning_cfg.db_partitioning,
        search_method=partitioning_cfg.search_method,
        cost_model=partitioning_cfg.cost_model,
        min_partition_size=partitioning_cfg.min_partition_size,
        amplitude_variable_name=partitioning_normal.amplitude_variable_name,
        output_variable_prefix=partitioning_normal.output_variable_prefix,
        output_variables=partitioning_normal.output_variables,
    )
    print("Finished normal partitions")
    partitioning_single_difference = cfg.stm_partitioning_single_difference
    stm = stm_partitioning(
        stm,
        db_partitioning=partitioning_cfg.db_partitioning,
        search_method=partitioning_cfg.search_method,
        cost_model=partitioning_cfg.cost_model,
        min_partition_size=partitioning_cfg.min_partition_size,
        amplitude_variable_name=partitioning_single_difference.amplitude_variable_name,
        output_variable_prefix=partitioning_single_difference.output_variable_prefix,
        output_variables=partitioning_single_difference.output_variables,
    )
    print("Finished single difference partitions")

# Rechunk to prevent inconsistent chunks
stm = stm.chunk({"time": -1, "space": cfg.ps_selection.output_chunks})

# Do outlier detection
if cfg.do_detect_outliers_stm:
    stm = detect_outliers_stm(
        stm,
        db_outlier_detection=cfg.detect_outliers_stm.db_outlier_detection,
        window_size=cfg.detect_outliers_stm.window_size,
        n_sigma=cfg.detect_outliers_stm.n_sigma,
    )
    print("Finished outlier detection")

# Save
stm.to_zarr(cfg.paths['stm_file'], mode='w')
print("Done!")

client.close()
