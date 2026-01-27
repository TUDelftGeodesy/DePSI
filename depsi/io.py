"""io methods."""

import json
import os
import re
import warnings
from datetime import datetime
from glob import glob
from io import StringIO
from typing import Literal

import geopandas as gpd
import numpy as np
import pandas as pd
import sarxarray
import scipy.spatial as scs
import xarray as xr
from shapely.geometry import Point, Polygon

from depsi.constants import SPEED_OF_LIGHT
from depsi.utils import _orbit_fit, npdatetime64_to_datetime

# Define constants
SC_N_PATTERN = r"\s+([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)"
ALLOWED_KNMI_DATA_COLUMNS = ["TG", "TN", "TX", "RH", "RXH", "EV24"]
SHAPEFILE_PROJECTIONS = {
    "RD": {
        "EPSG_code": "EPSG:28992",
        "x_crd_layer": "x_euclidean_proj_epsg28992",
        "y_crd_layer": "y_euclidean_proj_epsg28992",
    },
    "WGS84": {
        "EPSG_code": "EPSG:4326",
        "x_crd_layer": "lon",
        "y_crd_layer": "lat",
    },
}
SHAPEFILE_FIELD_NAMES = {
    "geometry": "Point",
    "properties": [
        "ID",
        "X (RD) [m]",
        "Y (RD) [m]",
        "H [m-NAP]",
        "Lat (WGS84) [deg]",
        "Lon (WGS84) [deg]",
        "h (WGS84) [m]",
        "Azimuth",
        "Range",
        "FUNC_INSERTS_MODEL_PARAMS_HERE",
        "Std linear [mm/y]",
        "STC [mm]",
        "Coherence [0-1]",
        "Std [mm]",
    ],
}
CSV_FIELD_NAMES = [
    "ID",
    "X (RD) [m]",
    "Y (RD) [m]",
    "H [m-NAP]",
    "Lat (WGS84) [deg]",
    "Lon (WGS84) [deg]",
    "h (WGS84) [m]",
    "Azimuth",
    "Range",
    "FUNC_INSERTS_MODEL_PARAMS_HERE",
    "Std linear [mm/y]",
    "STC [mm]",
    "Coherence [0-1]",
    "Std [mm]",
    "FUNC_INSERTS_TIMESERIES_HERE",
    "FUNC_INSERTS_AMP_HERE",
]
PORTAL_CSV_FIELD_NAMES = [
    "pnt_id",
    "pnt_lat",
    "pnt_lon",
    "pnt_demheight",
    "pnt_height",
    "pnt_azimuth",
    "pnt_range",
    "pnt_quality",
    "pnt_linear",
    "FUNC_INSERTS_MODEL_PARAMS_HERE",
    "FUNC_INSERTS_TIMESERIES_HERE",
    "FUNC_INSERTS_AMP_HERE",
]


def read_metadata(resfile, mode="raw", **kwargs):
    """Read metadata from a DORIS v5 resfile.

    Modified from the original functions in:
    https://github.com/Pbaz98/Caroline-Radar-Coding-Toolbox/blob/main/gecoris/dorisUtils.py.
    """
    # Deprecation warning for this function
    warning_msg = (
        "The depsi.io.read_metadata function is deprecated and will be removed in a "
        "future release. For reading metadata of coregistered SLC stacks, "
        "please use the sarxarray.read_metadata function instead."
    )
    warnings.warn(warning_msg, DeprecationWarning, stacklevel=2)

    # check crop_flag
    if mode == "coreg" and "crop" in kwargs:
        crop_flag = kwargs["crop"]
    else:
        crop_flag = 0

    # Open the file
    with open(resfile) as file:
        content = file.read()

    # +++   - Satellite ID
    pattern = r"Product type specifier:\s+(.*?(?=\n))"
    match = re.search(pattern, content)
    sat_id = match.group(1).upper()

    # ++++1 - Geometry [DESCENDING or ASCENDING]
    pattern = r"PASS:\s+(.*?(?=\n))"
    match = re.search(pattern, content)
    geometry = match.group(1).upper()

    # ++++ 2 - Acquisition Date [dictionary with __datetime__ and str, in the format 'yyyy-mm-dd hh:mm:ss']
    pattern = r"First_pixel_azimuth_time \(UTC\):\s+(\d+-\w+-\d+\s+)(\d+:\d+:\d+)"
    match = re.search(pattern, content)

    # --- extract the datetime string
    datetime_toconvert = match.group(1) + match.group(2)
    # Parse the original datetime string
    acq_date = datetime.strptime(datetime_toconvert, "%Y-%b-%d %H:%M:%S")
    acq_date.strftime("%Y%m%d")

    # ++++ 3 - azimuth0time

    # convert from time format to seconds of the day
    pattern = r"(\d+):(\d+):(\d+.\d+)"
    match = re.search(pattern, content)
    azimuth0time = int(match.group(1)) * 3600 + int(match.group(2)) * 60 + float(match.group(3))

    # ++++ 4 - range0time
    pattern = r"Range_time_to_first_pixel \(2way\) \(ms\):" + SC_N_PATTERN
    match = re.search(pattern, content)
    # devide by 2 to balance the two way travel
    range0time = float(match.group(1)) * 1e-3 / 2

    # ++++ 5 - prf
    pattern = r"Pulse_Repetition_Frequency \(computed, Hz\):" + SC_N_PATTERN
    match = re.search(pattern, content)
    prf = float(match.group(1))

    # ++++ 6 - rsr
    pattern = r"Range_sampling_rate \(computed, MHz\):" + SC_N_PATTERN
    match = re.search(pattern, content)
    rsr = float(match.group(1)) * 1e6 * 2

    # ++++ 7 - wavelength
    pattern = r"Radar_wavelength \(m\):" + SC_N_PATTERN
    match = re.search(pattern, content)
    wavelength = float(match.group(1))

    # ++++ 8 - orbit_fit

    # Define the regular expression pattern to match the table rows
    pattern = r"(\d+)\s+([-+]?\d+\.\d+(?:\.\d+)?)\s+([-+]?\d+\.\d+(?:\.\d+)?)\s+([-+]?\d+\.\d+(?:\.\d+)?)"

    # extract the table rows
    table_rows = re.findall(pattern, content)

    orbit = np.ones((len(table_rows), 4))

    for i in range(len(table_rows)):
        for j in range(4):
            orbit[i][j] = float(table_rows[i][j])

    # Generate the orbfit dictionary
    orbfit = _orbit_fit(orbit, verbose=0)

    # ++++ 9 - range_spacing
    pattern = r"rangePixelSpacing:" + SC_N_PATTERN
    match = re.search(pattern, content)
    range_spacing = float(match.group(1))

    # ++++ 10 - azimuth_spacing
    pattern = r"azimuthPixelSpacing:" + SC_N_PATTERN
    match = re.search(pattern, content)
    azimuth_spacing = float(match.group(1))

    # ++++ 11 - center_lon
    pattern = r"Scene_centre_longitude:" + SC_N_PATTERN
    match = re.search(pattern, content)
    center_lon = float(match.group(1))

    # ++++ 12 - center_lat
    pattern = r"Scene_centre_latitude:" + SC_N_PATTERN
    match = re.search(pattern, content)
    center_lat = float(match.group(1))

    # ++++ 13 - center_h
    pattern = r"Scene_center_heading:" + SC_N_PATTERN
    match = re.search(pattern, content)
    center_h = float(match.group(1))

    # ++++ 14 - n_azimuth
    pattern = r"Number_of_lines_original:" + SC_N_PATTERN
    match = re.search(pattern, content)
    n_azimuth = int(match.group(1))

    # ++++ 15 - n_range
    pattern = r"Number_of_pixels_original:" + SC_N_PATTERN
    match = re.search(pattern, content)
    n_range = int(match.group(1))

    # ++++ 16 - swath
    pattern = r"SWATH:\s+IW(\d+)"
    match = re.search(pattern, content)
    swath = int(match.group(1))

    # ++++ 17 - center_azimuth
    center_azimuth = np.round(n_azimuth / 2)

    # ++++ 18 - beta0, rank, chirprate
    beta0 = 237
    if swath == 1:
        rank = 9
        chirp_rate = 1078230321255.894
    elif swath == 2:
        rank = 8
        chirp_rate = 779281727512.0481
    elif swath == 3:
        rank = 10
        chirp_rate = 801450949070.5804

    # resolutions [from s1 annual performance reports]
    az_resolutions = np.array([21.76, 21.89, 21.71])
    sr_resolutions = np.array([2.63, 3.09, 3.51])  # slant range resolution
    azimuth_resolution = az_resolutions[swath - 1]

    # ++++ 20 - range_resolution
    pattern = r"Total_range_band_width \(MHz\):" + SC_N_PATTERN
    match = re.search(pattern, content)
    range_resolution = SPEED_OF_LIGHT / (2 * float(match.group(1)) * 1e6)

    # ++++ 21 - nBursts
    burst_n = None

    # ++++ 23 - steering_rate
    pattern = r"Azimuth_steering_rate \(deg/s\):" + SC_N_PATTERN
    match = re.search(pattern, content)
    steering_rate = float(match.group(1)) * np.pi / 180

    # ++++ 24 and 25 - azFmRateArray and dcPolyArray
    # Are skipped because the io.datetimeToMJD function is missing

    # ++++ 26 - pri
    pattern = r"Pulse_Repetition_Frequency_raw_data\(TOPSAR\):" + SC_N_PATTERN
    match = re.search(pattern, content)
    pri = 1 / float(match.group(1))

    # ++++ 27 - rank
    # See Beta0 section

    # ++++ 28 - chirp_rate

    # ++++ 29 - n_azimuth
    if crop_flag:
        crop_file = "/".join(str(resfile).split("/")[0:-2]) + "/nlines_crp.txt"
        with open(crop_file) as file:
            content = file.readlines()
            n_lines, first_line, last_line = (
                int(content[0].strip()),
                int(content[1].strip()),
                int(content[2].strip()),
            )

    else:
        # Extract first
        pattern = r"First_line \(w.r.t. original_image\):" + SC_N_PATTERN
        match = re.search(pattern, content)
        first_line = int(match.group(1))
        # Extract last
        pattern = r"Last_line \(w.r.t. original_image\):" + SC_N_PATTERN
        match = re.search(pattern, content)
        last_line = int(match.group(1))
        # difference
        n_lines = last_line - first_line + 1

    # ++++ 30 - n_range
    if crop_flag:
        crop_file = "/".join(str(resfile).split("/")[0:-2]) + "/npixels_crp.txt"
        with open(crop_file) as file:
            content = file.readlines()
            n_pixels, first_pixel, last_pixel = (
                int(content[0].strip()),
                int(content[1].strip()),
                int(content[2].strip()),
            )
    else:
        # Extract first
        pattern = r"First_pixel \(w.r.t. original_image\):" + SC_N_PATTERN
        match = re.search(pattern, content)
        first_pixel = int(match.group(1))
        # Extract last
        pattern = r"Last_pixel \(w.r.t. original_image\):" + SC_N_PATTERN
        match = re.search(pattern, content)
        last_pixel = int(match.group(1))
        # difference
        n_pixels = last_pixel - first_pixel + 1

    # ----------------------------------------

    # Fill the dictionary
    datewise_metadata = {
        "sat_id": sat_id,
        "orbit": geometry,
        "acq_date": acq_date,
        "azimuth0time": azimuth0time,
        "range0time": range0time,
        "prf": prf,
        "rsr": rsr,
        "wavelength": wavelength,
        "orbit_fit": orbfit,
        "range_spacing": range_spacing,
        "azimuth_spacing": azimuth_spacing,
        "center_lon": center_lon,
        "center_lat": center_lat,
        "center_h": center_h,
        "n_azimuth": n_azimuth,
        "n_range": n_range,
        "1stAzimuth": first_line,
        "1stRange": first_pixel,
        "swath": swath,
        "center_azimuth": center_azimuth,
        "beta0": beta0,
        "azimuth_resolution": azimuth_resolution,
        "range_resolution": range_resolution,
        "slant_range_resolution": sr_resolutions,
        "nBursts": 1,
        "burstInfo": burst_n,
        "steering_rate": steering_rate,
        "pri": pri,
        "rank": rank,
        "chirp_rate": chirp_rate,
        "n_lines": n_lines,
        "n_pixels": n_pixels,
        # -------------------------------------------------------------------
    }

    return datewise_metadata


def read_weather_data(filename: str, dates: list, requested_data_columns: tuple = ("TG", "RH")) -> dict:
    """Read columns of a KNMI weather data file at specific dates into a dictionary.

    The weather file is downloadable from https://www.knmi.nl/nederland-nu/klimatologie/daggegevens . Values that are
    not available in the file are replaced by `np.nan`

    Parameters
    ----------
    filename : str
        absolute filepath to the KNMI weather data file
    dates : list
        list of datetime.datetime objects of days at which the data is to be returned
    requested_data_columns : tuple
        tuple of strings of the column names. Currently implemented:
        - "TG": average daily temperature [deg C]
        - "TN": minimum temperature [deg C]
        - "TX": maximum temperature [deg C]
        - "RH": total daily precipitation [mm]
        - "RXH": maximum hourly precipitation [mm]
        - "EV24": reference evapotranspiration following Makkink [mm]

    Returns
    -------
    dict
        Dictionary with as keys the requested dates, as argument a dictionary with as keys requested columns, as
    argument the value
    """
    # check if the input is valid
    assert os.path.exists(filename), f"The requested file {filename} does not exist!"
    assert np.all([isinstance(date, datetime) for date in dates]), "Not all dates are of type datetime.datetime!"

    assert np.all(
        [requested_data_column in ALLOWED_KNMI_DATA_COLUMNS for requested_data_column in requested_data_columns]
    ), (
        f"Invalid requested data column detected in"
        f"{requested_data_columns}, allowed are "
        f"{ALLOWED_KNMI_DATA_COLUMNS}. See documentation"
        f" for explanation on abbreviations."
    )

    # read the file, and add a datetime column to the Pandas Dataframe
    weather_data = pd.read_csv(filename, sep=",", skiprows=51, skipinitialspace=True)
    weather_data["DateTime"] = weather_data.iloc[:, 1].apply(lambda x: pd.to_datetime(str(x), format="%Y%m%d"))

    # generate the output by looping over the requested dates and columns
    datewise_data = {}
    for date in dates:
        datewise_data[date] = {}
        # get the data at that date
        for data_column in requested_data_columns:
            value = weather_data[weather_data["DateTime"] == date][data_column].values[0]
            # convert the value to the correct units
            if np.isnan(value):
                datewise_data[date][data_column] = np.nan
            else:
                match data_column:
                    case "TG" | "TN" | "TX" | "EV24":  # provided in 0.1 deg C or 0.1 mm
                        datewise_data[date][data_column] = int(value) / 10
                    case "RH" | "RXH":  # provided in 0.1 mm, where -1 indicates < 0.05 mm
                        match value:
                            case -1:  # turn values smaller than 0.05 mm to 0
                                datewise_data[date][data_column] = 0.0
                            case _:  # otherwise, divide by 10
                                datewise_data[date][data_column] = int(value) / 10
                    case _:
                        raise NotImplementedError(f"Unknown requested column {data_column}!")

    return datewise_data


def read_slc_stack(
    filename: str, engine: str = "zarr", nlines_file: str = None, npixels_file: str = None, chunks=(500, 500)
) -> xr.Dataset:
    """Read a stack of SLCs into an xarray.Dataset.

    Supports different engines for reading:
    - zarr: reads a zarr archive (default)
    - doris: reads using the doris engine

    Parameters
    ----------
    filename : str
        Absolute filepath to the data archive (zarr folder or doris stack folder).
    engine : str, optional
        Engine to use for reading the data. Defaults to 'zarr'.
    nlines_file : str, optional
        Required for doris engine. Path to the file containing number of lines in the stack.
    npixels_file : str, optional
        Required for doris engine. Path to the file containing number of pixels in the stack.
    chunks : tuple, optional
        Tuple specifying the chunk size for loading doris stacks (default is (500, 500)).

    Returns
    -------
    xarray.Dataset
        Lazily loaded dataset with:
        - coordinates azimuth, range, lat, lon, time
        - variables h2ph, complex, amplitude, phase
    """
    assert os.path.exists(filename), f"The requested file/folder {filename} does not exist!"

    if engine.lower() == "zarr":
        # Load the zarr file as a xr.Dataset
        dataset = xr.open_zarr(filename)
        # Add complex, amplitude, and phase to the dataset
        slcs = sarxarray.from_dataset(dataset)
        return slcs

    elif engine.lower() == "doris":
        if nlines_file is None or npixels_file is None:
            raise ValueError(
                "For doris engine, 'nlines_file' and 'npixels_file' must be provided. Recommended 500x500."
            )

        # Collect file paths of the SLC stack
        stack_list = glob(os.path.join(filename, "*", "slc_srd.raw"))
        if not stack_list:
            raise FileNotFoundError(f"No SLC files found in {filename} matching pattern */slc_srd.raw")

        # Read the number of lines and pixels from the configuration files
        try:
            with open(nlines_file) as f:
                nlines = int(f.readline().strip())
            with open(npixels_file) as f:
                npixels = int(f.readline().strip())
        except Exception as e:
            raise RuntimeError("Failed to read number of lines or pixels. ") from e

        # Load the SLC stack using sarxarray
        try:
            slc_stack = sarxarray.from_binary(stack_list, shape=(nlines, npixels), dtype=np.complex64, chunks=chunks)
            return slc_stack
        except Exception as e:
            raise RuntimeError("Failed to load the SLC stack. ") from e

    else:
        raise ValueError(f"Unsupported engine '{engine}'. Use 'zarr' or 'doris'.")


def read_rcs_csv(file_path):
    r"""Load an STM-like CSV resulting from the RadarCoding Toolbox.

    Load an STM-like CSV resulting from the RadarCoding Toolbox, filter out metadata or header marked by
    '*****' markers, and extract dates from the header.

    Parameters
    ----------
    file_path : str
        Path to the CSV file to load.
    skip_rows : int, optional
        Number of rows to skip before reading the actual data. Default is 5.
    delimiter : str, optional
        The delimiter used in the CSV file. Default is tab (`\t`).

    Returns
    -------
    pd.DataFrame
        The cleaned DataFrame containing the STM with 0/1 flags indicating the existence of the targets' data.
    list of str
        A list of date strings extracted from the header.
    """
    try:
        # Read the file to find metadata and filter lines
        with open(file_path) as file:
            lines = file.readlines()

        # Find the indices of lines containing '*****'
        start_idx = None
        end_idx = None
        for i, line in enumerate(lines):
            if "*****" in line:
                if start_idx is None:
                    start_idx = i  # First occurrence of '*****'
                else:
                    end_idx = i  # Second occurrence of '*****'
                    break

        # Filter out the content between the '*****' markers
        filtered_lines = lines[:start_idx] + lines[end_idx + 1 :]

        # Join the filtered lines into a single string for pandas to read
        filtered_data = "".join(filtered_lines)

        # Load the cleaned data into a DataFrame
        df = pd.read_csv(StringIO(filtered_data))

        # Extract the header row (first row of the DataFrame)
        header_row = df.iloc[0]

        # Get the columns starting from the 8th index onward (i.e., index 7 corresponds to 20150430)
        dates = header_row.index[7:].tolist()

        # Print success message
        print(f"Radar Coding (RC) Toolbox output file '{file_path}' successfully loaded.")

        return df, dates

    except Exception as e:
        # If an error occurs, print the error message
        raise RuntimeError(f"Error loading the Radar Coding (RC) Toolbox output file '{file_path}") from e


def get_targets_from_slc(slc_stack, targets):
    """Extract target-matched data from a SLC stack.

    This function matches given azimuth/range coordinates to the closest points
    in the SLC stack, extracts all data variables, and aligns detection_flag
    time series from the targets dataset with the SLC time dimension. Then return
    a Space-Time Matrix (STM) containing the aligned data.

    Parameters
    ----------
    slc_stack : xarray.Dataset
        Dataset containing SLC data with dimensions (azimuth, range, time),
        as well as variables like lat, lon, azimuth, range, etc.
    targets : xarray.Dataset
        A Space-Time Matrix containing target information with coordinates (space, time)
        and variables such as azimuth, range, target, and detection_flag.

    Returns
    -------
    matching_scatterers : xarray.Dataset
        Dataset containing the extracted SLC data for matching targets.

    Notes
    -----
        - detection_flag values are aligned to the SLC timestamps.
    """
    # Get bounds for SLC stack
    # Select targets within the bound
    targets_in_bound = targets.where(
        (targets["azimuth"] >= slc_stack["azimuth"].min())
        & (targets["azimuth"] <= slc_stack["azimuth"].max())
        & (targets["range"] >= slc_stack["range"].min())
        & (targets["range"] <= slc_stack["range"].max()),
        drop=True,
    )

    # Using nearest neighbor to select slc pixels matching targets
    matching_scatterers = slc_stack.sel(
        azimuth=targets_in_bound["azimuth"], range=targets_in_bound["range"], method="nearest"
    )

    # Compute detection flag masks
    # First linear interpolate in time dimension
    # This only take into account the two neighbors in time
    detection_flag = targets_in_bound["detection_flag"].interp(time=matching_scatterers["time"], method="linear")
    # Epochs outside the original 1 periods will have values <1, Set them to 0
    detection_flag = detection_flag.where(detection_flag >= 1.0, 0)

    # Insert detection_flag to matching_scatterers
    matching_scatterers["detection_flag"] = xr.DataArray(detection_flag.data, dims=("space", "time"))

    return matching_scatterers


def export_to_csv(
    stm: xr.Dataset,
    save_path: str,
    model_parameter_layer_names: tuple | list,
    ts_proj: Literal["los", "vertical"],
    point_annotation_label: str,
) -> None:
    """Export an STM to CSV-format.

    Parameters
    ----------
    stm: xr.Dataset
        The STM to export
    save_path: str
        Full path to where to save the CSV
    model_parameter_layer_names: tuple | list
        Tuple or list with the layer names of the model parameters in the order that they will be stored in the csv
    ts_proj: Literal["los", "vertical"]
        Whether the saved time series is Line-of-Sight or projected onto the vertical
    point_annotation_label: str
        An extra annotation given to the point IDs in the CSV (`point_annotation_label`_az########r########)

    Raises
    ------
    AssertionError
        - if the provided path does not end in .csv
    ValueError
        - If unknown or unmodeled parameters have been added to `CSV_FIELD_NAMES`
    """
    assert save_path.split(".")[-1] == "csv", f"Provided path {save_path} is not a csv!"

    # Store everything into memory
    for layer in ["azimuth", "range", "lat", "lon", "h", "amplitude", "stc"]:
        stm[layer].data = stm[layer].values
    for layer in ["unwrapped_phase"]:  # these ones are projectable onto the vertical
        if ts_proj == "los":
            stm[layer].data = stm[layer].values
        else:
            stm[f"{layer}_pov"].data = stm[f"{layer}_pov"].values
    for layer in model_parameter_layer_names:
        if ts_proj == "los":
            stm[layer].data = stm[layer].values
        else:
            stm[f"{layer}_pov"].data = stm[f"{layer}_pov"].values
    for layer in ["x_euclidean_proj_epsg28992", "y_euclidean_proj_epsg28992"]:  # these ones might not exist
        if layer in stm.variables.keys():
            stm[layer].data = stm[layer].values

    fmt_dates = [npdatetime64_to_datetime(date).strftime("%Y%m%d") for date in stm["time"].values]
    headers = []
    for name in CSV_FIELD_NAMES:
        if "FUNC_INSERTS" not in name:
            headers.append(name)
        else:
            match name:
                case "FUNC_INSERTS_MODEL_PARAMS_HERE":
                    for param in model_parameter_layer_names:
                        headers.append(param)
                case "FUNC_INSERTS_TIMESERIES_HERE":
                    for fmt_date in fmt_dates:
                        headers.append(fmt_date)
                case "FUNC_INSERTS_AMP_HERE":
                    for fmt_date in fmt_dates:
                        headers.append(f"a_{fmt_date}")
                case _:
                    raise ValueError(f"Function insert {name} requested but not defined!")

    f = open(save_path, "w")
    f.write(",".join(headers) + "\n")
    for point in stm["space"].values:
        point_values = []
        for value in headers:
            match value:
                case "ID":
                    az = int(stm.azimuth.sel(space=point).values)
                    r = int(stm.range.sel(space=point).values)
                    point_values.append(f"{point_annotation_label}_az{az:0>8d}r{r:0>8d}")
                case "X (RD) [m]":
                    if "x_euclidean_proj_epsg28992" in stm.variables.keys():
                        point_values.append(round(float(stm.x_euclidean_proj_epsg28992.sel(space=point).values), 2))
                    else:
                        point_values.append("NULL")
                case "Y (RD) [m]":
                    if "y_euclidean_proj_epsg28992" in stm.variables.keys():
                        point_values.append(round(float(stm.y_euclidean_proj_epsg28992.sel(space=point).values), 2))
                    else:
                        point_values.append("NULL")
                case "H [m-NAP]":
                    if "rd_h" in stm.variables.keys():
                        # point_values.append(round(float(stm.rd_h.sel(space=point).values), 4))
                        point_values.append("NULL")
                    else:
                        point_values.append("NULL")
                case "Lat (WGS84) [deg]":
                    point_values.append(round(float(stm.lat.sel(space=point).values), 8))
                case "Lon (WGS84) [deg]":
                    point_values.append(round(float(stm.lon.sel(space=point).values), 8))
                case "h (WGS84) [m]":
                    point_values.append(round(float(stm.h.sel(space=point).values), 3))
                case "Azimuth":
                    point_values.append(int(stm.azimuth.sel(space=point).values))
                case "Range":
                    point_values.append(int(stm.range.sel(space=point).values))
                case "Std linear [mm/y]":
                    # point_values.append(round(float(stm.linear_std.sel(space=point).values), 3))
                    point_values.append("NULL")
                case "STC [mm]":
                    point_values.append(round(float(stm.stc.sel(space=point).values), 3))
                case "Coherence [0-1]":
                    # point_values.append(round(float(stm.coherence.sel(space=point).values), 4))
                    point_values.append("NULL")
                case "Std [mm]":
                    # point_values.append(round(float(stm.ts_std.sel(space=point).values), 3))
                    point_values.append("NULL")
                case _:
                    if value in model_parameter_layer_names:
                        if ts_proj == "vertical":
                            point_values.append(round(float(stm[f"{value}_pov"].sel(space=point).values), 5))
                        elif ts_proj == "los":
                            point_values.append(round(float(stm[value].sel(space=point).values), 5))
                    elif value in fmt_dates:
                        if ts_proj == "vertical":
                            point_values.append(
                                round(
                                    float(stm.unwraped_phase_pov.sel(space=point).isel(time=fmt_dates.index(value))), 5
                                )
                            )
                        elif ts_proj == "los":
                            point_values.append(
                                round(float(stm.unwrapped_phase.sel(space=point).isel(time=fmt_dates.index(value))), 5)
                            )
                    elif value[:2] == "a_" and value[2:] in fmt_dates:
                        point_values.append(
                            round(float(stm.amplitude.sel(space=point).isel(time=fmt_dates.index(value[2:]))), 3)
                        )
                    else:
                        raise ValueError(f"Requested header {value} but this is undefined!")

        f.write(",".join([str(p) for p in point_values]) + "\n")
    f.close()


def export_to_skygeo_portal(
    stm: xr.Dataset,
    save_path: str,
    ts_proj: Literal["los", "vertical"],
    point_annotation_label: str,
    satellite: str,
    asc_dsc: Literal["asc", "dsc"],
    azimuth_spacing: float,
    range_spacing: float,
    model_names: list[str],
    model_parameter_layer_names: list | tuple,
) -> None:
    """Export an STM to CSV-format.

    Parameters
    ----------
    stm: xr.Dataset
        The STM to export
    save_path: str
        Full path to where to save the CSV
    ts_proj: Literal["los", "vertical"]
        Whether the saved time series is Line-of-Sight or projected onto the vertical
    point_annotation_label: str
        An extra annotation given to the point IDs in the CSV (`point_annotation_label`_az########r########)
    satellite: str
        Name of the satellite that acquired the imagery
    asc_dsc: Literal["asc", "dsc"]
        Whether the viewing geometry is ascending or descending
    azimuth_spacing: float
        The pixel spacing in azimuth direction
    range_spacing: float
        The pixel spacing in range direction
    model_names: list[str]
        List containing the names of all models applied in the parameter estimation
    model_parameter_layer_names: tuple | list
        Tuple or list with the layer names of the model parameters in the order that they will be stored in the csv

    Raises
    ------
    AssertionError
        - if the provided path does not end in .csv
    ValueError
        - If unknown or unmodeled parameters have been added to `CSV_FIELD_NAMES`
    """
    assert save_path.split(".")[-1] == "csv", f"Provided path {save_path} is not a csv!"

    # Store everything into memory
    for layer in ["azimuth", "range", "lat", "lon", "h", "amplitude", "local_incidence_angle"]:
        stm[layer].data = stm[layer].values
    for layer in ["pnt_velocity", "unwrapped_phase"]:  # these ones are projectable onto the vertical
        if ts_proj == "los":
            stm[layer].data = stm[layer].values
        else:
            stm[f"{layer}_pov"].data = stm[f"{layer}_pov"].values
    for layer in model_parameter_layer_names:
        if ts_proj == "los":
            stm[layer].data = stm[layer].values
        else:
            stm[f"{layer}_pov"].data = stm[f"{layer}_pov"].values

    # first generate the CSV file

    if save_path.split(".")[-2].split("_")[-1] != "portal":
        save_path = ".".join(save_path.split(".")[:-1]) + "_portal.csv"

    fmt_dates = [npdatetime64_to_datetime(date).strftime("%Y%m%d") for date in stm["time"].values]
    headers = []
    for name in PORTAL_CSV_FIELD_NAMES:
        if "FUNC_INSERTS" not in name:
            headers.append(name)
        else:
            match name:
                case "FUNC_INSERTS_TIMESERIES_HERE":
                    for fmt_date in fmt_dates:
                        headers.append(f"d_{fmt_date}")
                case "FUNC_INSERTS_AMP_HERE":
                    for fmt_date in fmt_dates:
                        headers.append(f"a_{fmt_date}")
                case "FUNC_INSERTS_MODEL_PARAMS_HERE":
                    for param in model_parameter_layer_names:
                        if param != "pnt_height":  # this one already gets added in a different place
                            headers.append(param)
                case _:
                    raise ValueError(f"Function insert {name} requested but not defined!")

    f = open(save_path, "w")
    f.write(",".join(headers) + "\n")
    for point in stm["space"].values:
        point_values = []
        for value in headers:
            match value:
                case "pnt_id":
                    az = int(stm.azimuth.sel(space=point).values)
                    r = int(stm.range.sel(space=point).values)
                    point_values.append(f"{point_annotation_label}_az{az:0>8d}r{r:0>8d}")
                case "pnt_lat":
                    point_values.append(round(float(stm.lat.sel(space=point).values), 8))
                case "pnt_lon":
                    point_values.append(round(float(stm.lon.sel(space=point).values), 8))
                case "pnt_demheight":
                    point_values.append(round(float(stm.h.sel(space=point).values), 3))
                case "pnt_azimuth":
                    point_values.append(int(stm.azimuth.sel(space=point).values))
                case "pnt_range":
                    point_values.append(int(stm.range.sel(space=point).values))
                case "pnt_quality":
                    # point_values.append(round(float(stm.ens_coh_local.sel(space=point).values), 3))
                    point_values.append(1)
                case "pnt_linear":
                    if ts_proj == "vertical":
                        point_values.append(round(float(stm.pnt_velocity_pov.sel(space=point).values), 5))
                    elif ts_proj == "los":
                        point_values.append(round(float(stm.pnt_velocity.sel(space=point).values), 5))
                case _:
                    if value[:2] == "d_" and value[2:] in fmt_dates:
                        if ts_proj == "vertical":
                            point_values.append(
                                round(
                                    float(
                                        stm.unwrapped_phase_pov.sel(space=point).isel(time=fmt_dates.index(value[2:]))
                                    ),
                                    5,
                                )
                            )
                        elif ts_proj == "los":
                            point_values.append(
                                round(
                                    float(stm.unwrapped_phase.sel(space=point).isel(time=fmt_dates.index(value[2:]))), 5
                                )
                            )
                    elif value[:2] == "a_" and value[2:] in fmt_dates:
                        point_values.append(
                            round(float(stm.amplitude.sel(space=point).isel(time=fmt_dates.index(value[2:]))), 3)
                        )
                    elif value in model_parameter_layer_names:
                        if ts_proj == "vertical":
                            point_values.append(round(float(stm[f"{value}_pov"].sel(space=point).values), 5))
                        elif ts_proj == "los":
                            point_values.append(round(float(stm[value].sel(space=point).values), 5))
                    else:
                        raise ValueError(f"Requested header {value} but this is undefined!")

        f.write(",".join([str(p) for p in point_values]) + "\n")
    f.close()

    ref_point_idx = [stm.idx_refpnt]
    ref_pt_coords = [
        f"{round(float(stm.lat.sel(space=point).values), 8)}, {round(float(stm.lon.sel(space=point).values), 8)}"
        for point in ref_point_idx
    ]
    ref_pts = " ; ".join(ref_pt_coords)

    # then generate the JSON file
    json_dict = {
        "acquisition_period": f"{min(fmt_dates)} - {max(fmt_dates)}",
        "number_of_observations_in_time": len(fmt_dates),
        "number_of_measurements_in_AoI": len(list(stm["space"].values)),
        "resolution": f"{round(azimuth_spacing, 1)} x {round(range_spacing, 1)} m",
        "deformation_direction": "Line of Sight" if ts_proj == "los" else "Projected onto Vertical",
        "DEM": "SRTM",
        "reference_point_location": ref_pts,
        "satellite_name": satellite,
        "satellite_incidence_angle": round(np.mean(stm.local_incidence_angle.values)[0], 1),
        "satellite_pass_direction": asc_dsc,
        "processing_id": point_annotation_label,
        "DePSI_version": "DePSI_group",
        "description": "",
        "estimated_models_for_time_series": ", ".join(model_names),
    }
    json_filename = ".".join(save_path.split(".")[:-1]) + ".json"
    with open(json_filename, "w") as f:
        json.dump(json_dict, f, ensure_ascii=False, indent=4)


def export_to_shapefile(
    stm: xr.Dataset,
    save_path: str,
    projection: Literal["RD", "WGS84"],
    model_parameter_layer_names: tuple | list,
    point_annotation_label: str,
) -> None:
    """Export an STM to a shapefile.

    Parameters
    ----------
    stm: xr.Dataset
        The STM to export
    save_path: str
        Full path to where to save the shapefile
    projection: Literal["RD", "WGS84"]
        Whether to output the shapefile in RD or in WGS84 (properties will contain both if available regardless, this
        only affects the coordinate system of the shapefile itself)
    model_parameter_layer_names: tuple | list
        Tuple or list with the layer names of the model parameters in the order that they will be stored in the
        shapefile
    point_annotation_label: str
        An extra annotation given to the point IDs in the CSV (`point_annotation_label`_az########r########)

    Raises
    ------
    AssertionError
        - when an unknown projection is passed
        - when the save path does not end in .shp
    ValueError
        - When an undefined property is requested from SHAPEFILE_SCHEMA
    """
    assert projection in SHAPEFILE_PROJECTIONS.keys(), f"Unknown requested projection {projection}!"
    assert save_path.split(".")[-1] == "shp", f"Provided path {save_path} is not a shapefile!"

    # Store coordinates into memory
    for layer in [SHAPEFILE_PROJECTIONS[projection]["x_crd_layer"], SHAPEFILE_PROJECTIONS[projection]["y_crd_layer"]]:
        stm[layer].data = stm[layer].values

    schema = []
    for name in SHAPEFILE_FIELD_NAMES["properties"]:
        if "FUNC_INSERTS" not in name:
            schema.append(name)
        else:
            match name:
                case "FUNC_INSERTS_MODEL_PARAMS_HERE":
                    for param in model_parameter_layer_names:
                        schema.append(param)
                case _:
                    raise ValueError(f"Function insert {name} requested but not defined!")

    properties = {}
    geometry = []
    for point in stm["space"].values:
        geometry.append(
            Point(
                float(stm[SHAPEFILE_PROJECTIONS[projection]["x_crd_layer"]].sel(space=point).values),
                float(stm[SHAPEFILE_PROJECTIONS[projection]["y_crd_layer"]].sel(space=point).values),
            )
        )
    for value in schema:
        match value:
            case "ID":
                properties[value] = []
                for point in stm["space"].values:
                    az = int(stm.azimuth.sel(space=point).values)
                    r = int(stm.range.sel(space=point).values)
                    properties[value].append(f"{point_annotation_label}_az{az:0>8d}r{r:0>8d}")
            case "X (RD) [m]":
                if "x_euclidean_proj_epsg28992" in stm.variables.keys():
                    properties[value] = [round(float(val), 2) for val in stm.x_euclidean_proj_epsg28992.values]
                else:
                    properties[value] = [np.nan for _ in stm["space"].values]
            case "Y (RD) [m]":
                if "y_euclidean_proj_epsg28992" in stm.variables.keys():
                    properties[value] = [round(float(val), 2) for val in stm.y_euclidean_proj_epsg28992.values]
                else:
                    properties[value] = [np.nan for _ in stm["space"].values]
            case "H [m-NAP]":
                if "rd_h" in stm.variables.keys():
                    # properties[value] = [round(float(val), 3) for val in stm.rd_h.values]
                    properties[value] = [np.nan for _ in stm["space"].values]
                else:
                    properties[value] = [np.nan for _ in stm["space"].values]
            case "Lat (WGS84) [deg]":
                properties[value] = [round(float(val), 8) for val in stm.lat.values]
            case "Lon (WGS84) [deg]":
                properties[value] = [round(float(val), 8) for val in stm.lon.values]
            case "h (WGS84) [m]":
                properties[value] = [round(float(val), 3) for val in stm.h.values]
            case "Azimuth":
                properties[value] = [int(val) for val in stm.azimuth.values]
            case "Range":
                properties[value] = [int(val) for val in stm.range.values]
            case "Std linear [mm/y]":
                # properties[value] = [round(float(val), 3) for val in stm.linear_std.values]
                properties[value] = [np.nan for _ in stm["space"].values]
            case "STC [mm]":
                properties[value] = [round(float(val), 3) for val in stm.stc.values]
            case "Coherence [0-1]":
                # properties[value] = [round(float(val), 4) for val in stm.coherence.values]
                properties[value] = [np.nan for _ in stm["space"].values]
            case "Std [mm]":
                # properties[value] = [round(float(val), 3) for val in stm.ts_std.values]
                properties[value] = [np.nan for _ in stm["space"].values]
            case _:
                if value in model_parameter_layer_names:
                    properties[value] = [round(float(val), 5) for val in stm[value].values]
                else:
                    raise ValueError(f"Requested header {value} but this is undefined!")

    properties["geometry"] = geometry
    dataframe = gpd.GeoDataFrame(properties, crs=SHAPEFILE_PROJECTIONS[projection]["EPSG_code"])
    dataframe.to_file(save_path)


def export_convex_hull_to_shapefile(stm: xr.Dataset, save_path: str, projection: Literal["RD", "WGS84"]) -> None:
    """Export the convex hull of an STM to a shapefile.

    Parameters
    ----------
    stm: xr.Dataset
        The STM to export
    save_path: str
        Full path to where to save the shapefile, ending in .shp
    projection: Literal["RD", "WGS84"]
        Whether to output the convex hull in RD or in WGS84

    Raises
    ------
    AssertionError
        When an unknown projection is provided
        When the provided save_path does not end in .shp
    """
    assert projection in SHAPEFILE_PROJECTIONS.keys(), f"Unknown requested projection {projection}!"
    assert save_path.split(".")[-1] == "shp", f"Provided path {save_path} is not a shapefile!"

    point_coords_x = stm[SHAPEFILE_PROJECTIONS[projection]["x_crd_layer"]].values.flatten()
    point_coords_y = stm[SHAPEFILE_PROJECTIONS[projection]["y_crd_layer"]].values.flatten()
    point_coords = np.vstack([point_coords_x, point_coords_y]).T

    hull = scs.ConvexHull(point_coords)

    hull_vertices = hull.points[hull.vertices]
    listified_hull = [list(vertex) for vertex in hull_vertices]
    listified_hull.append(listified_hull[0])  # to make it a closed hull

    hull_gdf = gpd.GeoDataFrame(
        {"geometry": [Polygon(listified_hull)]}, crs=SHAPEFILE_PROJECTIONS[projection]["EPSG_code"]
    )

    hull_gdf.to_file(save_path)
