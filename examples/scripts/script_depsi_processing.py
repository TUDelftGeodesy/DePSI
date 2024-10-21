"""Example script for selecting PS from SLCs.

This .py script is designed to be executed with a Dask SLURMCluster on a SLURM managed HPC system.
It should be executed through a SLURM script by `sbatch` command.
Please do not run this script by "python xxx.py" on a login node.
"""

import logging
import os
import socket
import xarray as xr
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
from dask.distributed import Client
from dask_jobqueue import SLURMCluster
import sarxarray
import stmtools

from pydepsi.io import read_metadata
from pydepsi.classification import ps_selection, network_stm_seletcion

# Make a logger to log the stages of processing
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()  # create console handler
ch.setLevel(logging.INFO)
logger.addHandler(ch)


def get_free_port():
    """Get a non-occupied port number."""
    sock = socket.socket()
    sock.bind(("", 0))  # Bind a port, it will be busy now
    freesock = sock.getsockname()[1]  # get the port number
    sock.close()  # Free the port, so it can be used later
    return freesock


# ---- Config 1: Human Input ----
# Input data paths
path_slc_zarr = Path("/project/caroline/slc_file.zarr")  # Zarr file of all SLCs
path_metadata = Path("/project/caroline/metadata.res")  # Metadata file

# Parameters PS selection
ps_selection_method = 'nmad'  # Method for PS selection
ps_selection_threshold = 0.45  # Threshold for PS selection

# Parameters network selection
network_stm_quality_metric = 'nmad' # Quality metric for network selection
network_stm_quality_threshold = 0.45 # Quality threshold for network selection
dist_thres = 200  # Distance threshold for network selection, in meters
include_index = [101] # Force including the 101th point, use None if no point need to be included


# Output config
overwrite_zarr = False  # Flag for zarr overwrite
chunk_space = 10000  # Output chunk size in space dimension
path_figure = Path("./figure")  # Output path for figure
path_ps_zarr = Path("./ps.zarr") # output file for selected PS

path_figure.mkdir(exist_ok=True)    # Make figure directory if not exists

# ---- Config 2: Dask configuration ----

# Option 1: Initiate a new SLURMCluster
# Uncomment the following part to setup a new Dask SLURMCluster
# N_WORKERS = 4 # Manual input: number of workers to spin-up
# FREE_SOCKET = get_free_port() # Get a free port
# cluster = SLURMCluster(
#     name="dask-worker",  # Name of the Slurm job
#     queue="normal", # Name of the node partition on your SLURM system
#     cores=4, # Number of cores per worker
#     memory="32 GB",  # Total amount of memory per worker
#     processes=1,  # Number of Python processes per worker
#     walltime="3:00:00",  # Reserve each worker for X hour
#     scheduler_options={"dashboard_address": f":{FREE_SOCKET}"},  # Host Dashboard in a free socket
# )
# logger.info(f"Dask dashboard hosted at port: {FREE_SOCKET}.")
# logger.info(
#     f"If you are forwarding Jupyter Server to a local port 8889, \
#     you can access it at: localhost:8889/proxy/{FREE_SOCKET}/status"
# )

# Option 2: Use an existing SLURMCluster by giving the schedular address 
# Uncomment the following part to use an existing Dask SLURMCluster
ADDRESS = "tcp://XX.X.X.XX:12345" # Manual input: Dask schedular address
SOCKET = 12345 # Manual input: port number. It should be the number after ":" of ADDRESS
cluster = None  # Keep this None, needed for an if statement
logger.info(f"Dask dashboard hosted at port: {SOCKET}.")
logger.info(
    f"If you are forwarding Jupyter Server to a local port 8889, \
    you can access it at: localhost:8889/proxy/{SOCKET}/status"
)

if __name__ == "__main__":
    # ---- Processing Stage 0: Initialization ----
    logger.info("Initializing ...")

    # Initiate a Dask client
    if cluster is None:
        # Use existing cluster
        client = Client(ADDRESS)
    else:
        # Scale a certain number workers
        # each worker will appear as a Slurm job
        cluster.scale(jobs=N_WORKERS)
        client = Client(cluster)

    # Load metadata
    metadata = read_metadata(path_metadata)

    # ---- Processing Stage 1: Pixel Classification ----
    # Load the SLC data
    logger.info("Processing Stage 1: Pixel Classification")
    logger.info("Loading SLC data ...")
    ds = xr.open_zarr(path_slc_zarr) # Load the zarr file as a xr.Dataset
    # Construct SLCs from xr.Dataset
    # this construct three datavariables: complex, amplitude, and phase 
    slcs = sarxarray.from_dataset(slcs)

    # A rechunk might be needed to make a optimal usage of the resources
    # Uncomment the following line to apply a rechunk after loading the data
    # slcs = slcs.chunk({"azimuth":1000, "range":1000, "time":-1})

    # Select PS
    logger.info("PS Selection ...")
    stm_ps = ps_selection(method, threshold, method=ps_selection_method, output_chunks=chunk_space)

    # Re-order the PS to make the spatially adjacent PS in the same chunk
    logger.info("Reorder selected scatterers ...")
    stm_ps_reordered = stm_ps.stm.reorder(xlabel='lon', ylabel='lat')

    # Save the PS to zarr
    logger.info("Writting seleced pixels to Zarr ...")
    if overwrite_zarr:
        stm_ps_reordered.to_zarr(path_ps_zarr, mode="w")
    else:
        stm_ps_reordered.to_zarr(path_ps_zarr)

    # ---- Processing Stage 2: Network Processing ----
    # Uncomment the following line to load the PS data from zarr
    # stm_ps_reordered = xr.open_zarr(path_ps_zarr)

    # Select network points
    logger.info("Select network scatterers ...")
    # Apply a pre-filter
    stm_network_candidates = xr.where(stm_ps_reordered[network_stm_quality_metric]<network_stm_quality_threshold)
    # Select based on sparsity and quality
    stm_network = network_stm_seletcion(stm_network_candidates, 
                                        dist_thres,
                                        include_index=include_index,
                                        sortby_var=network_stm_quality_metric,
                                        azimuth_spacing=metadata['azimuth_spacing'],
                                        range_spacing=metadata['range_spacing'])

    # Close the client when finishing
    client.close()
