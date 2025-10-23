"""io methods."""

import os
import re
from datetime import datetime
from glob import glob
from io import StringIO

import numpy as np
import pandas as pd
import sarxarray
import xarray as xr

from depsi.utils import _orbit_fit

# Define constants
SC_N_PATTERN = r"\s+([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)"
SPEED_OF_LIGHT = 299792458.0  # m/s
ALLOWED_KNMI_DATA_COLUMNS = ["TG", "TN", "TX", "RH", "RXH", "EV24"]


def read_metadata(resfile, mode="raw", **kwargs):
    """Read metadata from a DORIS v5 resfile.

    Modified from the original functions in:
    https://github.com/Pbaz98/Caroline-Radar-Coding-Toolbox/blob/main/gecoris/dorisUtils.py.
    """
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


def export_to_csv(stm: xr.Dataset, save_path: str) -> None:
    """Export an STM to CSV-format.

    Parameters
    ----------
    stm: xr.Dataset
        The STM to export
    save_path: str
        Full path to where to save the CSV
    """
    pass


def export_to_skygeo_portal(stm: xr.Dataset, save_path: str) -> None:
    """Export an STM to the files necessary for uploading to the SkyGeo portal.

    This function produces both a CSV and a JSON, which together can be uploaded to the SkyGeo portal.

    Parameters
    ----------
    stm: xr.Dataset
        The STM to export
    save_path: str
        Full path to where to save the CSV file. The JSON file will be saved in the same directory.
    """
    pass


def export_to_shapefile(stm: xr.Dataset, save_path: str) -> None:
    """Export an STM to a shapefile.

    Parameters
    ----------
    stm: xr.Dataset
        The STM to export
    save_path: str
        Full path to where to save the shapefile
    """
    pass


def export_convex_hull_to_shapefile(stm: xr.Dataset, save_path: str) -> None:
    """Export the convex hull of an STM to a shapefile.

    Parameters
    ----------
    stm: xr.Dataset
        The STM to export
    save_path: str
        Full path to where to save the shapefile
    """
    pass
