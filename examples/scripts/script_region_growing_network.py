"""This script runs through the region growing first order network formation"""
"""Confidence-Optimized Robust Geodetic Network (CORG Network) is de nieuwe naam"""

import xarray as xr

from depsi.io import read_weather_data
from depsi.utils import npdatetime64_to_datetime

# the original STM is output from the script script_sidelobe_detection.py
stm_original = (
    "/Users/sanvandiepen/PycharmProjects/workingEnvironment2/test_zarr/nl_amsterdam_s1_dsc_t110_stm_nosl.zarr"
)
stm_save_path = (
    "/Users/sanvandiepen/PycharmProjects/workingEnvironment2/test_zarr/nl_amsterdam_s1_dsc_t110_stm_nosl.zarr"
)
knmi_file_path = "/Users/sanvandiepen/PycharmProjects/workingEnvironment2/test_zarr/etmgeg_240.txt"

stm = xr.open_zarr(stm_original)
timestamps = [npdatetime64_to_datetime(date, tz_aware=False) for date in stm.time.values]
read_weather_data(knmi_file_path, timestamps, requested_data_columns=("TG", ))
