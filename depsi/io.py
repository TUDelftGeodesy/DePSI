"""io methods."""

import os
import re
import subprocess
from datetime import datetime
from glob import glob
from io import StringIO

import numpy as np
import pandas as pd
import sarxarray
import xarray as xr
import yaml

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
    range0time = float(match.group(1)) * 1e-3 / 2  # devide by 2 to balance the two way travel

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


def read_slc_stack(filename: str) -> xr.Dataset:
    """Read a zarr stack of SLCs into a xarray.Dataset.

    Reads a zarr archive, and converts it to an xarray dataset compatible with the point selection functions.

    Parameters
    ----------
    filename : str
        absolute filepath to the zarr archive. The zarr archive should contain:
        - coordinates azimuth, range, lat, lon, time
        - variables h2ph, imag, real

    Returns
    -------
    xarray.Dataset
        Lazily loaded dataset with:
        - coordinates azimuth, range, lat, lon, time
        - variables h2ph, complex, amplitude, phase
    """
    assert os.path.exists(filename), f"The requested file {filename} does not exist!"

    # Load the zarr file as a xr.Dataset
    dataset = xr.open_zarr(filename)
    # Add complex, amplitude, and phase to the dataset
    slcs = sarxarray.from_dataset(dataset)

    return slcs


def collect_srdraw_file_paths(folder_path):
    """Collect all paths of `slc_srd.raw` files in subfolders of the doris stack specified by the folder_path.

    Parameters
    ----------
    folder_path : str
        The path to the main folder containing subfolders.

    Returns
    -------
    list of str
        A list of full paths to the `slc_srd.raw` files.
    """
    # Use glob to search for slc_srd.raw files in all subfolders
    file_paths = glob(os.path.join(folder_path, "*", "slc_srd.raw"))
    return file_paths


def read_first_line(file_path):
    """Read the first line from the specified text file.

    Reads the first line from the specified text file.
    Meant to read number of pixels and lines from stack stiching crp.txt files

    Parameters
    ----------
    file_path : str
        Path to the text file.

    Returns
    -------
    str
        The first line of the file, stripped of leading and trailing whitespace.
    """
    with open(file_path) as file:
        first_line = file.readline().strip()
    return first_line


def load_stm_rcscsv(file_path):
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
        print(f"Error loading the Radar Coding (RC) Toolbox output file '{file_path}': {e}")
        return None, None


def run_script_subprocess(script_path, args):
    """Run a Python script with arguments using `subprocess`.

    Meant to run the RadarCoding Toolbox within Python scripts.

    Parameters
    ----------
    script_path : str
        Path to the script to execute.
    args : list of str
        Arguments to pass to the script.

    Raises
    ------
    subprocess.CalledProcessError
        If the script execution fails.
    """
    try:
        subprocess.run(["python", script_path] + args, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error while running the RCS Toolbox: {e}")


def doris_sar_stack_to_xarray(stack_folder, nlines_file, npixels_file, chunks=(500, 500)):
    """Load a stack of SLCs in XArray format from a given folder and configuration files.

    Parameters
    ----------
    stack_folder : str
        Path to the folder containing the stack of SLCs.
    nlines_file : str
        Path to the file containing the number of lines in the stack (nlines_crp.txt).
    npixels_file : str
        Path to the file containing the number of pixels in the stack (npixels_crp.txt).
    chunks : tuple, optional
        Tuple specifying the chunk size to use when loading the stack (default is (500, 500)).

    Returns
    -------
    xarray.DataArray or None
        Returns the loaded SLC stack as an xarray DataArray if successful, otherwise None.
    """
    # Collect file paths of the SLC stack
    stack_list = collect_srdraw_file_paths(stack_folder)

    # Read the number of lines and pixels from the configuration files
    try:
        nlines = int(read_first_line(nlines_file))
        npixels = int(read_first_line(npixels_file))
    except Exception as e:
        print(f"Warning: Failed to read the number of lines or pixels. Error: {e}")
        return None

    # Attempt to load the SLC stack using sarxarray
    try:
        slc_stack = sarxarray.from_binary(stack_list, shape=(nlines, npixels), dtype=np.complex64, chunks=chunks)
        print("Successfully loaded the SLC stack.")
        return slc_stack
    except Exception as e:
        print(f"Warning: Failed to load the SLC stack. Error: {e}")
        return None


def ensure_rc_csv_exists(rcsoutput_folder, rcsanalysis, arguments):
    """Ensure that an output of the RadarCoding Toolbox (*_RC.csv file) exists in the specified folder.

    Ensure that an output of the RadarCoding Toolbox (*_RC.csv file) exists
    in the specified folder, if not, run the RadarCoding tool.

    Parameters
    ----------
    rcsoutput_folder : str
        Path to the folder where *_RC.csv file is expected to be found.
    rcsanalysis : str
        Path to the RadarCoding tool script to execute.
    arguments : list of str
        Arguments to pass to the RCS tool script.

    Returns
    -------
    str or None
        The path of the *_RC.csv file if found or generated, otherwise None.

    Notes
    -----
    This function checks if a file matching the pattern `*_RC.csv` exists in the
    specified output folder. If no such file is found, it runs the RCS tool
    (`detectDesignatedTargets.py`) with the provided arguments and checks if
    the output file is generated. It will print appropriate messages based
    on the file’s existence or failure to generate.
    """
    # Check for any file that matches the pattern *_RC.csv
    output_file = None
    for file_name in os.listdir(rcsoutput_folder):
        if file_name.endswith("_RC.csv"):
            output_file = os.path.join(rcsoutput_folder, file_name)
            break

    # Check if the file exists
    if not output_file:
        print("No *_RC.csv file found. Running the RCS tool...")
        run_script_subprocess(rcsanalysis, arguments)

        # Check again if any *_RC.csv file was generated
        output_file = None
        for file_name in os.listdir(rcsoutput_folder):
            if file_name.endswith("_RC.csv"):
                output_file = os.path.join(rcsoutput_folder, file_name)
                break

        if output_file:
            print(f"Found the output file: {output_file}")
        else:
            print("Failed to generate *_RC.csv file. Please check the RCS tool output and parameters.")
    else:
        print(f"Found the output file: {output_file}. Proceed with the identifying the targets in the stack.")

    return output_file


def extract_dttarget_data_from_slc(slc_stack, matching_indices, targets):
    """Extract matching data from the slc_stack of a designated target.

    Description
    -----------
    Extract matching data from the slc_stack for the given matching azimuth and range indices
    and target information of a designated target.

    Parameters
    ----------
    slc_stack : xarray.Dataset
        The dataset containing the slc data with azimuth, range, and other variables.

    matching_indices : list of tuples
        A list of tuples, where each tuple contains the indices (azimuth_idx, range_idx) of the matching targets
        in the `slc_stack`.

    targets : xarray.Dataset
        A dataset containing target information, including azimuth and range coordinates.

    Returns
    -------
    dict
        A dictionary containing the extracted data for the matching targets, with variable names as keys
        and the corresponding values as lists of data arrays.

    list
        A list of the target names corresponding to the matching indices.

    list
        A list of latitudes for the matching targets.

    list
        A list of longitudes for the matching targets.
    """
    # Initialize dictionary to store the extracted data
    matched_slc_targets_dict = {}

    # Lists to store latitudes, longitudes, and target names
    lat_vals = []
    lon_vals = []
    target_names = []

    # Iterate over each matching coordinate
    for az_idx, rg_idx in matching_indices:
        # Slice slc_stack for each matching (azimuth, range) pair
        matching_point_slc_data = slc_stack.isel(azimuth=az_idx, range=rg_idx)

        # Add the sliced data to the dictionary with 'target' as the dimension
        for var in matching_point_slc_data.data_vars:
            if var not in matched_slc_targets_dict:
                matched_slc_targets_dict[var] = []
            matched_slc_targets_dict[var].append(matching_point_slc_data[var].values)

        # Append the lat and lon values for each matching target
        lat_vals.append(matching_point_slc_data["lat"].values)  # lat is scalar for each target
        lon_vals.append(matching_point_slc_data["lon"].values)  # lon is scalar for each target

        # Append the corresponding target name (from the `targets` dataset)
        target_idx = np.where(targets["azimuth"].data == slc_stack["azimuth"].data[az_idx])[0][0]
        target_name = targets["target"].data[target_idx]
        target_names.append(target_name)

    return matched_slc_targets_dict, lat_vals, lon_vals, target_names


def load_ymlparams(yml_file):
    """Load input parameters from a YAML file.

    Parameters
    ----------
    yml_file : str
        Path to the YAML file to load.

    Returns
    -------
    dict
        Dictionary containing the loaded parameters.
    """
    try:
        with open(yml_file) as file:
            params = yaml.safe_load(file)
        print(f"Parameters successfully loaded from '{yml_file}'.")
        return params
    except Exception as e:
        print(f"Error loading YAML file '{yml_file}': {e}")
        return None
