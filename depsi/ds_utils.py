import json
import os
import re
import xml.etree.ElementTree as ET  # noqa: N817
from datetime import datetime
from pathlib import Path

import dask.array as da
import geopandas as gpd
import numpy as np
import sarxarray
import xarray as xr
from slc import *  # noqa: F403


def extract_mother_date(xml_file):
    """Function to extract mother date as an integer from doris_input.xml.

    Parameters
    ----------
    xml_file : file
        An xml containing the setup for doris input, inc. the mother date.

    Returns
    -------
    int
        Mother date as an integer with yyyyMMdd format.

    Raises
    ------
    ValueError
        Raise when mother date element not found in the XML file.
    """  # noqa: D401
    # Parse the XML file
    tree = ET.parse(xml_file)
    root = tree.getroot()

    # Find the mother element
    mother_date_element = root.find(".//master_date")

    if mother_date_element is not None:
        # Extract the date string
        mother_date_str = mother_date_element.text

        # Convert the date string to a datetime object
        mother_date_obj = datetime.strptime(mother_date_str, "%Y-%m-%d")

        # Format the date as an integer in yyyyMMdd format
        mother_date_int = int(mother_date_obj.strftime("%Y%m%d"))

        return mother_date_int
    else:
        raise ValueError("master_date element not found in the XML file.")


def identify_stacks(settings):
    """Identify stack data information and store it to a metadata file.
    
    This function works for coregistered SAR data using DORIS.

    Parameters
    ----------
    settings : dict
        Variable containing initial parameters for analysis.

    Returns
    -------
    dict
        Two dictionaries: settings and stack_meta.
        Stack meta contains information of one track
    """
    ## Detect number of mother images / tracks
    dirlist = os.listdir(settings["stack_root_dir"])
    pattern = re.compile(rf'^{settings["stack_prefix"]}_?s1_[ad]sc_t\d{{3}}$')
    stack_dirs = sorted([os.path.join(settings["stack_root_dir"], i) for i in dirlist if pattern.match(i)])
    num_tracks_init = len(stack_dirs)

    print(f"Found {num_tracks_init} DORIS v5 stacks in the specified location")

    assert num_tracks_init >= 1, "Could not find SAR data in the specified directory"

    ## Check if stack intersects the region of interest
    num_tracks = 0
    for i in range(num_tracks_init):
        gdf_stack = gpd.read_file(os.path.join(stack_dirs[i], "stackburst_coverage.shp"))
        gdf_aoi = gpd.read_file(settings["aoi_shapefile"])

        ## Ensure aoi shapefiles are in the WGS84 coordinate system
        if gdf_aoi.crs != "EPSG:4326":
            gdf_aoi = gdf_aoi.to_crs(gdf_stack.crs)

        ## Dissolve all features in each shapefile into single geometries
        boundary_stack = gdf_stack.unary_union
        boundary_aoi = gdf_aoi.unary_union

        ## Check for intersection
        if boundary_stack.intersects(boundary_aoi):
            num_tracks += 1
        else:
            print(f"Discarding stack {stack_id[i]}")  # noqa: F405
            stack_dirs.remove(stack_dirs[i])

    assert num_tracks >= 1, "Could not find SAR data in the specified directory"

    ## Add parameter to settings and save it
    settings["num_tracks"] = num_tracks
    with open(os.path.join(settings["proj_dir"], settings["settings_filename"]), "w") as file:
        json.dump(settings, file, indent=2)

    ## Extract the metadata of each stack
    ## Read variables
    start_date = settings["start_date"]
    end_date = settings["end_date"]

    stack_meta_list = []
    for i in range(num_tracks):
        mother_date = extract_mother_date(os.path.join(stack_dirs[i], "doris_input.xml"))

        assert (
            mother_date << start_date or mother_date >> end_date
        ), "Mother image is outside specified date range, add it later"

        if settings["processor"] == "flinsar":
            mother_dir = os.path.join(stack_dirs[i], str(mother_date))
            full_dates = np.loadtxt(os.path.join(stack_dirs[i], "dates.txt"), dtype=int)
        if settings["processor"] == "caroline":
            mother_dir = os.path.join(stack_dirs[i], "stack", str(mother_date))
            full_dates = np.loadtxt(os.path.join(stack_dirs[i], "stack", "dir.txt"), dtype=int)
        sid = np.argwhere(full_dates >= start_date)[0][0]
        eid = np.argwhere(full_dates <= end_date)[-1][0]
        slc_dates = full_dates[sid : eid + 1].tolist()
        if mother_date < start_date:
            slc_dates.insert(0, mother_date)
        if mother_date > end_date:
            slc_dates.append(mother_date)
        nslcs = len(slc_dates)

        ## Extract the name of the stack folder
        stack_name = stack_dirs[i].split("/")[-1]
        stack_id = "_".join(stack_name.split("_")[-3:])

        ## Create list of dates without the mother
        ifg_dates = slc_dates.copy()
        ifg_dates.remove(mother_date)
        nifgs = len(ifg_dates)

        print(f"Track: {stack_id}  Mother: {mother_date}  nifgs: {nifgs}")

        ## Create a list of stack paths
        stack_dir = os.path.dirname(mother_dir)
        mother_idx = slc_dates.index(mother_date)

        if settings["reslc"] == "yes":
            stack_type = "cint"
            filelist = sorted(
                [
                    fp
                    for fp in Path(stack_dir).glob("**/cint_srd.raw")
                    if "swath" not in str(fp) or "burst" not in str(fp)
                ]
            )
            filelist = [str(path) for path in filelist]
            slc_paths = [path for path in filelist if extract_date(path) in ifg_dates]
            if settings["processor"] == "caroline":
                slc_paths.insert(mother_idx, os.path.join(mother_dir, "slave_rsmp_reramped.raw"))
            if settings["processor"] == "flinsar":
                slc_paths.insert(mother_idx, os.path.join(mother_dir, "slc_srd.raw"))

        if settings["reslc"] == "no":
            stack_type = "slc"
            if settings["processor"] == "caroline":
                filelist = sorted(Path(stack_dir).glob("**/slave_rsmp_reramped.raw"))
            if settings["processor"] == "flinsar":
                filelist = sorted(Path(stack_dir).glob("**/slc_srd.raw"))
            slc_paths = [path for path in filelist if extract_date(path) in slc_dates]

        ## Read data from the master res file
        if settings["processor"] == "caroline":
            with open(os.path.join(mother_dir, "master.res")) as f:
                lines = f.readlines()
                swath = lines[37].strip().split()[-1]
                mode = lines[38].strip().split()[-1]
                r_px_spacing = float(lines[46].strip().split()[-1])
                az_px_spacing = float(lines[47].strip().split()[-1])
                npixels_res = int(lines[98].strip().split()[-1])
                nlines_res = int(lines[99].strip().split()[-1])
        elif settings["processor"] == "flinsar":
            with open(os.path.join(stack_dirs[i], "nlines_crp.txt")) as f:
                lines = f.readlines()
                nlines_res = int(lines[0])
            with open(os.path.join(stack_dirs[i], "npixels_crp.txt")) as f:
                lines = f.readlines()
                npixels_res = int(lines[0])
            r_px_spacing = "n/a"
            az_px_spacing = "n/a"

        ## Store variables to the metadata file
        stack_meta = {
            "processor": settings["processor"],
            "aoi_path": settings["aoi_shapefile"],
            "stack_prefix": settings["stack_prefix"],
            "meta_dir": os.path.join(settings["meta_dir"], stack_id),
            "do_reslc": settings["reslc"],
            "mother_date": mother_date,
            "stack_dir": stack_dir,
            "mother_dir": mother_dir,
            "slc_paths": slc_paths,
            "stack_type": stack_type,
            "stack_id": stack_id,
            "nslcs": nslcs,
            "nifgs": nifgs,
            "slc_dates": slc_dates,
            "ifg_dates": ifg_dates,
            "r_px_spacing": r_px_spacing,
            "az_px_spacing": az_px_spacing,
            "npixels_res": npixels_res,
            "nlines_res": nlines_res,
        }
        if not os.path.exists(stack_meta["meta_dir"]):
            os.makedirs(stack_meta["meta_dir"])
        with open(os.path.join(stack_meta["meta_dir"], "stack_meta_" + stack_id + ".json"), "w") as file:
            json.dump(settings, file, indent=2)

        stack_meta_list.append(stack_meta)

    return settings, stack_meta_list


def create_processing_folders(settings):
    """Create directories if not already there.

    Parameters
    ----------
    settings : dict
        Variable containing parameters for analysis.

    Returns
    -------
    dict
        An updated variable settings.
    """
    ## Do not modify
    settings["run_dir"] = os.path.join(settings["proj_dir"], settings["run_name"])
    settings["meta_dir"] = os.path.join(settings["run_dir"], "metadata/")
    settings["phase_est_dir"] = os.path.join(settings["run_dir"], "phase_estimation/")
    settings["stm_dir"] = os.path.join(settings["run_dir"], "stm/")

    if not os.path.exists(settings["run_dir"]):
        os.makedirs(settings["run_dir"])

    if not os.path.exists(settings["meta_dir"]):
        os.makedirs(settings["meta_dir"])

    if not os.path.exists(settings["phase_est_dir"]):
        os.makedirs(settings["phase_est_dir"])

    if not os.path.exists(settings["stm_dir"]):
        os.makedirs(settings["stm_dir"])

    return settings


def extract_date(path):  # noqa: D103
    match = re.search(r"/(\d{8})/", path)
    return int(match.group(1)) if match else None


def load_slc_stack(stack_meta, chunks=(500, 500)):
    """Load stack (SLCs or interferograms), assign coordinates, subset to the AoI and do re-slc if desired.

    Parameters
    ----------
    stack_meta : dict
        Contains information of a stack
    chunks : tuple, optional
        By default (500,500)

    Returns
    -------
    xr.Dataset
        SLC stack with three variables: (complex, amplitude, phase)
        and two coordinates: space (lat, lon, azimuth, range) and time.
    dict
        An updated stack metadata.
    """
    ## Read variables
    mother_dir = stack_meta["mother_dir"]
    mother_date = stack_meta["mother_date"]
    filelist = stack_meta["slc_paths"]
    npixels = stack_meta["npixels_res"]
    nlines = stack_meta["nlines_res"]
    slc_dates = stack_meta["slc_dates"]

    ## Load coordinates
    lat = sarxarray.from_binary(
        [os.path.join(mother_dir, "phi.raw")], shape=(nlines, npixels), vlabel="lat", dtype=np.float32, chunks=chunks
    )
    lon = sarxarray.from_binary(
        [os.path.join(mother_dir, "lam.raw")], shape=(nlines, npixels), vlabel="lon", dtype=np.float32, chunks=chunks
    )

    ## Load SLCs and crop to aoi
    if stack_meta["do_reslc"] == "no":
        slc_stack = sarxarray.from_binary(filelist, shape=(nlines, npixels), dtype=np.complex64, chunks=chunks)

        ## Drop amplitude and phase, keep complex only
        slc_stack = slc_stack.drop_vars(["amplitude", "phase"])

        ## Assign coordinates to the stack
        slc_stack = slc_stack.assign_coords(
            lat=(("azimuth", "range"), lat.squeeze().lat.data), lon=(("azimuth", "range"), lon.squeeze().lon.data)
        )
        # slc_stack = slc_stack.assign({'lon':lon, 'lat':lat})

        ## Assign datetime as time coordinates
        slc_stack["time"] = [datetime.strptime(str(date_int), "%Y%m%d") for date_int in slc_dates]

        ## Extract aoi indices for cropping stack to the area of interest
        l0, lN, p0, pN = extract_stack_aoi_indices(slc_stack, stack_meta) 

        ## Crop stack
        slc_stack_subset = slc_stack.sel(azimuth=slice(l0, lN), range=slice(p0, pN))

        ## Add amplitude and phase as attributes
        slc_stack_subset = get_amplitude(slc_stack_subset)
        slc_stack_subset = get_phase(slc_stack_subset)

    ## Recompute SLCs from IFGs
    if stack_meta["do_reslc"] == "yes":
        ## Load IFGs
        ifg_stack = sarxarray.from_binary(filelist, shape=(nlines, npixels), dtype=np.complex64, chunks=chunks)

        ## Drop amplitude and phase, keep complex only
        ifg_stack = ifg_stack.drop_vars(["amplitude", "phase"])

        ## Assign coordinates to the stack
        ifg_stack = ifg_stack.assign_coords(
            lat=(("azimuth", "range"), lat.squeeze().lat.data), lon=(("azimuth", "range"), lon.squeeze().lon.data)
        )

        ## Assign datetime as time coordinates
        ifg_stack["time"] = [datetime.strptime(str(date_int), "%Y%m%d") for date_int in slc_dates]

        ## Load mother SLC
        mother_idx = slc_dates.index(mother_date)
        slc_mother = ifg_stack.isel(time=slice(mother_idx, mother_idx + 1))

        ## Extract aoi indices for cropping stack to the area of interest
        l0, lN, p0, pN = extract_stack_aoi_indices(ifg_stack, stack_meta)

        ## Crop stack
        ifg_stack_subset = ifg_stack.sel(azimuth=slice(l0, lN), range=slice(p0, pN))
        slc_mother_subset = slc_mother.sel(azimuth=slice(l0, lN), range=slice(p0, pN))

        ## Re-SLC
        slc_stack_out = ifg_to_slc(slc_mother_subset, ifg_stack_subset)

        ## Insert mother slc to the re-slc stack
        ## TODO: check the mother slc
        slc_stack_subset = (
            xr.concat([slc_stack_out, slc_mother_subset], dim="time")
            .drop_duplicates(dim="time", keep="last")
            .sortby("time")
        )

        ## Add amplitude and phase as attributes
        slc_stack_subset = get_amplitude(slc_stack_subset)
        slc_stack_subset = get_phase(slc_stack_subset)

    ## Add the cropped shape to stack_meta
    stack_meta["nlines"] = int(lN - l0 + 1)
    stack_meta["npixels"] = int(pN - p0 + 1)

    ## Save the updated stack metadata
    with open(os.path.join(stack_meta["meta_dir"], "stack_meta_" + stack_meta["stack_id"] + ".json"), "w") as file:
        json.dump(stack_meta, file, indent=2)

    return slc_stack_subset, stack_meta


def extract_stack_aoi_indices(stack, stack_meta):
    """Find the start indices (l0, p0) and end indices (lN, pN) of the aoi on the stack.

    Parameters
    ----------
    stack : xr.Dataset
        SLC stack with three variables: (complex, amplitude, phase)
        and two coordinates: space (lat, lon, azimuth, range) and time.
    stack_meta : dict
        Contains information of the corresponding stack

    Returns
    -------
    int (four variables)
        Start and end indices of azimuth and range.
    """
    ## Load area of interest as geodataframe
    gdf_aoi = gpd.read_file(stack_meta["aoi_path"])
    assert gdf_aoi.crs == "EPSG:4326", "Area of interest is not using WGS84 coordinate reference system"
    lon_min, lat_min, lon_max, lat_max = gdf_aoi.total_bounds

    ## Create a mask for cropping stack to the area of interest
    mask = (stack["lat"] > lat_min) & (stack["lat"] < lat_max) & (stack["lon"] > lon_min) & (stack["lon"] < lon_max)
    mask_idx = np.argwhere(mask.values)
    # mask = mask.compute()

    ## Extract the start and the end of lines/azimuth and pixel/range
    l0, lN = min(mask_idx[:, 0]), max(mask_idx[:, 0])
    p0, pN = min(mask_idx[:, 1]), max(mask_idx[:, 1])

    return l0, lN, p0, pN


################################################################################################################
#    Source:                                                                                                   #
#    https://github.com/TUDelftGeodesy/sarxarray/blob/main/sarxarray/stack.py                                  #
#                                                                                                              #
################################################################################################################
def get_amplitude(slc):  # noqa: D103
    slc_out = slc.copy()
    meta_arr = np.array((), dtype=np.float32)
    amplitude = da.apply_gufunc(_compute_amp, "()->()", slc["complex"], meta=meta_arr)
    slc_out = slc_out.assign({"amplitude": (("azimuth", "range", "time"), amplitude)})
    return slc_out


def get_phase(slc):  # noqa: D103
    slc_out = slc.copy()
    meta_arr = np.array((), dtype=np.float32)
    phase = da.apply_gufunc(_compute_phase, "()->()", slc["complex"], meta=meta_arr)
    slc_out = slc_out.assign({"phase": (("azimuth", "range", "time"), phase)})
    return slc_out


def _compute_amp(complex):
    return np.abs(complex)


def _compute_phase(complex):
    return np.angle(complex)


################################################################################################################
