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

from pydepsi.classification import ps_selection

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


# Parameters
method = 'nmad'  # Method for selection
threshold = 0.45  # Threshold for selection

# Input data paths
path_slc_zarr = Path("/project/caroline/slc_file.zarr")  # Zarr file of all SLCs

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
    logger.info("Initializing ...")

    if cluster is None:
        # Use existing cluster
        client = Client(ADDRESS)
    else:
        # Scale a certain number workers
        # each worker will appear as a Slurm job
        cluster.scale(jobs=N_WORKERS)
        client = Client(cluster)

    # Load the SLC data
    logger.info("Loading data ...")
    ds = xr.open_zarr(path_slc_zarr) # Load the zarr file as a xr.Dataset
    # Construct SLCs from xr.Dataset
    # this construct three datavariables: complex, amplitude, and phase 
    slcs = sarxarray.from_ds(slcs)

    # A rechunk might be needed to make a optimal usage of the resources
    # Uncomment the following line to apply a rechunk after loading the data
    # slcs = slcs.chunk({"azimuth":1000, "range":1000, "time":-1})

    # Select PS
    stm_ps = ps_selection(method, threshold, method='nmad', output_chunks=chunk_space)

    # Re-order the PS to make the spatially adjacent PS in the same chunk
    stm_ps_reordered = stm_ps.stm.reorder(xlabel='lon', ylabel='lat')

    # Save the PS to zarr
    if overwrite_zarr:
        stm_ps_reordered.to_zarr(path_ps_zarr, mode="w")
    else:
        stm_ps_reordered.to_zarr(path_ps_zarr)

    # Close the client when finishing
    client.close()
