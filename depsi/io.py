"""io methods."""

import re
from datetime import datetime

import numpy as np

from depsi.utils import _orbit_fit

# Define constants
SC_N_PATTERN = r"\s+([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)"
SPEED_OF_LIGHT = 299792458.0  # m/s


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
