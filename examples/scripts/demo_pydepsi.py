import logging
import os
import socket
from pathlib import Path

import geopandas as gpd
import numpy as np
import sarxarray
import xarray as xr
from dask.distributed import Client
from dask_jobqueue import SLURMCluster
from matplotlib import pyplot as plt


def get_free_port():
    """Get a non occupied port number."""
    sock = socket.socket()
    sock.bind(("", 0))  # Bind a port, it will be busy now
    freesock = sock.getsockname()[1]  # get the port number
    sock.close()  # Free the port, so it can be used later
    return freesock


## Setup processing
# Make a logger to log the stages of processing
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()  # create console handler
ch.setLevel(logging.INFO)
logger.addHandler(ch)

# Paths and files
cwd = Path(os.getcwd())
path_slc = Path(cwd / "data/nl_amsterdam_s1_asc_t088")  # SLC stack processed by Doris V5
f_slc = "cint_srd.raw"  # Data file in each date folder under path_slc
f_lat = [path_slc / "lat.raw"]  # Geo referenced coordinates, lat
f_lon = [path_slc / "lon.raw"]  # Geo referenced coordinates, lon
overwrite_zarr = True  # Flag for zarr overwrite
path_stm = Path("./stm.zarr")  # Zarr output storage for STM
path_figure = Path("./figure")  # Output path for figure
path_polygon = Path(cwd / "data/bag_light_AMS_WGS84.gpkg")  # Path to the BRP polygon of NL.

# Metadata of the SLC stack
shape = (2000, 4000)  # Shape per SLC image
dtype = np.dtype(np.float64)  # Data type per image
reading_chunks = (500, 500)  # Reading chunk size

# Size of subset slice. this demoe will only process the slice
azimuth_subset = range(0, 2000)  # Subset slice, azimuth direction
range_subset = range(0, 4000)  # Subset slice, range direction

# Dask setup
n_workers = 16  # number of workers
# Config SlurmCluster
freesock = get_free_port()
cluster = SLURMCluster(
    name="dask-worker",  # Name of the Slurm job
    queue="normal",
    cores=4,
    memory="30 GB",  # Total amount of memory per job
    processes=1,  # Number of Python processes per job
    walltime="1:00:00",  # reserve each worker for 1 hour
    scheduler_options={"dashboard_address": f":{freesock}"},  # Host Dashboard in a free socket
)
logger.info(f"Dask dashboard hosted at port: {freesock}.")
logger.info(
    f"If you are forwarding Jupyter Server to a local port 8889, \
    you can access it at: localhost:8889/proxy/{freesock}/status"
)

if __name__ == "__main__":
    ## Step0: Setup environment
    logger.info("Initializing ...")
    # Scale a certain number workers, each worker will appear as a Slurm job
    cluster.scale(jobs=n_workers)
    client = Client(cluster)

    # Make figure directory if not exists
    path_figure.mkdir(exist_ok=True)

    ## Step1: Data loading
    logger.info("Loading data ...")
    # Build slcs lists
    list_slcs = [p for p in path_slc.rglob("*_cint_srd.raw")]
    list_slcs.sort()

    # Load complex data
    stack = sarxarray.from_binary(list_slcs, shape, dtype=np.complex64, chunks=reading_chunks)

    # Load coordinates
    lat = sarxarray.from_binary(f_lat, shape, vlabel="lat", dtype=np.float32, chunks=reading_chunks)
    lon = sarxarray.from_binary(f_lon, shape, vlabel="lon", dtype=np.float32, chunks=reading_chunks)
    stack = stack.assign_coords(
        lat=(("azimuth", "range"), lat.squeeze().lat.data), lon=(("azimuth", "range"), lon.squeeze().lon.data)
    )

    ## Step2: Make a spatial subset
    logger.info("Slicing SLC stack ...")
    stack_subset = stack.sel(azimuth=azimuth_subset, range=range_subset)

    ## Step3: Make mean reflection map (MRM)
    logger.info("Computing MRM ...")
    mrm = stack_subset.slcstack.mrm()
    mrm = mrm.compute()

    fig, ax = plt.subplots()
    ax.imshow(mrm)
    ax.set_aspect(2)
    im = mrm.plot(ax=ax, cmap="gray")
    im.set_clim([0, 40000])
    fig.savefig(path_figure / "mrm.png")

    ## Step4: Point selection
    logger.info("Point selection ...")
    stmat = stack_subset.slcstack.point_selection(threshold=4, method="amplitude_dispersion", chunks=5000)

    fig, ax = plt.subplots()
    plt.scatter(stmat.lon.data, stmat.lat.data, s=0.005)
    fig.savefig(path_figure / "selected_points.png")

    # Export point selection to Zarr
    if overwrite_zarr:
        stmat.to_zarr(path_stm, mode="w")
    else:
        if not path_stm.exists():
            stmat.to_zarr(path_stm)

    ## Step5: STM enrichment from Polygon file
    logger.info("STM enrichment ...")
    # Load SpaceTime Matrix from Zarr
    stm_demo = xr.open_zarr(path_stm)

    # Compute the bounding box
    xmin, ymin, xmax, ymax = [
        stm_demo["lon"].data.min().compute(),
        stm_demo["lat"].data.min().compute(),
        stm_demo["lon"].data.max().compute(),
        stm_demo["lat"].data.max().compute(),
    ]
    polygons = gpd.read_file(path_polygon, bbox=(xmin, ymin, xmax, ymax))
    polygons.plot()

    # Data enrichment
    fields_to_query = ["bouwjaar"]
    stm_demo = stm_demo.stm.enrich_from_polygon(polygons, fields_to_query)

    # Subset by Polygons
    stm_demo_subset = stm_demo.stm.subset(method="polygon", polygon=path_polygon)
    bouwjaar = stm_demo_subset["bouwjaar"].compute()

    # Visualize the classes
    import matplotlib.cm as cm

    colormap = cm.jet
    fig, ax = plt.subplots()
    plt.title("Construction year, PS")
    plt.scatter(stm_demo_subset.lon.data, stm_demo_subset.lat.data, c=bouwjaar, s=0.002, cmap=colormap)
    plt.clim([1900, 2023])
    plt.colorbar()
    fig.savefig(path_figure / "construction_year.png")

    ## Close the client when finishing
    client.close()
