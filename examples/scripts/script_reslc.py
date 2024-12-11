"""Example script for calculating SLCs from interferograms.

This .py script is designed to be executed with a Dask SLURMCluster on a SLURM managed HPC system.
It should be executed through a SLURM script by `sbatch` command.
Please do not run this script by "python xxx.py" on a login node.
"""

import logging
import os
import socket
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
from dask.distributed import Client
from dask_jobqueue import SLURMCluster
import sarxarray

from depsi.slc import ifg_to_slc
from depsi.io import read_metadata

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


# ---- Config 1: Path Configuration ----
# Paths and files
stack_dir = Path("/project/caroline/Share/stacks/nl_veenweiden_s1_dsc_t037/stack")  # ifg stack dir
mother_dir = stack_dir / "20200328"  # Mother image dir
reading_chunks = (2000, 2000)  # Reading chunks (azimuth, range) from binary

# Output config
overwrite_zarr = False  # Flag for zarr overwrite
writing_chunks = (2000, 2000)  # Writing chunks to zarr, (azimuth, range)
path_figure = Path("./figure")  # Output path for figure
path_figure.mkdir(exist_ok=True)    # Make figure directory if not exists


# ---- Config 2: Dask configuration ----

# Option 1: Intiate a new SLURMCluster
# Uncomment the following part to setup a new Dask SLURMCluster
# N_WORKERS = 4 # Manual input: number of workers to spin-up
# FREE_SOCKET = get_free_port() # Get a free port
# cluster = SLURMCluster(
#     name="dask-worker",  # Name of the Slurm job
#     queue="normal", # Name of the node partition on your SLURM system
#     cores=4, # Number of cores per worker
#     memory="30 GB",  # Total amount of memory per worker
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

    logger.info("Loading data ...")
    # Metadata
    f_mother_res = mother_dir / "slave.res"
    metadata = read_metadata(f_mother_res)

    # Coordinates
    f_lam = mother_dir / "lam.raw"
    f_phi = mother_dir / "phi.raw"

    # Mother SLC
    f_mother_slc = mother_dir / "slave_rsmp_reramped.raw"

    # List of SLC
    f_ifgs = list(sorted(stack_dir.rglob("2*/cint_srd.raw")))
    f_ifgs = f_ifgs[:3]

    shape = (metadata["n_lines"], metadata["n_pixels"])
    dtype_slc_ifg = np.dtype([("re", np.float32), ("im", np.float32)])
    dtype_lam_phi = np.float32

    # Lazy loading mother SLC and ifg stack
    mother = sarxarray.from_binary([f_mother_slc], shape, dtype=dtype_slc_ifg, chunks=reading_chunks)
    ifgs = sarxarray.from_binary(f_ifgs, shape, dtype=dtype_slc_ifg, chunks=reading_chunks)

    # Generate reconstructed SLCs
    slc_recon = ifg_to_slc(mother, ifgs)

    # Extract real and image part. remove other fields. convert to float16
    slc_recon_output = slc_recon.copy()
    slc_recon_output = slc_recon_output.assign(
        {
            "real": slc_recon_output["complex"].real.astype(np.float16),
            "imag": slc_recon_output["complex"].imag.astype(np.float16),
        }
    )
    slc_recon_output = slc_recon_output.drop_vars(["complex", "amplitude", "phase"])

    # Rechunk and write as zarr
    slc_recon_output = slc_recon_output.chunk({"azimuth": writing_chunks[0], "range": writing_chunks[1]})
    if overwrite_zarr:
        slc_recon.to_zarr("nl_veenweiden_s1_dsc_t037.zarr", mode="w")
    else:
        slc_recon.to_zarr("nl_veenweiden_s1_dsc_t037.zarr")

    # Close the client when finishing
    client.close()
