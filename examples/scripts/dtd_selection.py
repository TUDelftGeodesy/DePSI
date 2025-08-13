"""
Executable that identifies designated targets within a corregistered SLC stack
by loading the stack, loading the RadarCoding Toolbox output, querrying for targets in the stack

An example of the input parameter file (dtd_params.yml), 
an SLC stack in zarr frmat (nl_groningen_s1_dsc_t037_haren), 
the output of the RadarCoding Toolbox (s1_dsc037_RC.csv) can be found at:  
https://figshare.com/account/items/28218506/edit

Please place these files in ./examples/scripts/data

Created on Fri Jan 10 11:08:57 2025

@author: Alex Lapadat
"""

# Load Modules
import os
import sys
sys.path.append('/home/parallels/Sprint_Mobyle/DePSI_group/depsi/')

import sarxarray
import subprocess
import numpy as np
import xarray as xr
import pandas as pd
from depsi import io


## STEP 0 - LOAD INPUT PARAMETERS
# Specify folders and files to be used:
print("STEP 0 - LOAD INPUT PARAMS")
# Load the YAML parameters containing folders and files to be used:
params = io.load_ymlparams("/home/parallels/Sprint_Mobyle/DePSI_group/examples/scripts/data/dtd_params.yml")

if params is None:
    raise ValueError("Failed to load parameters from ''.")

# Access paths from the YAML file
doris_stack_folder = params['paths']['doris_stack_folder']
nlines_file        = params['paths']['nlines_file']
npixels_file       = params['paths']['npixels_file']

rcsoutput_folder   = params['paths']['rcs_output_folder']
rcsAnalizer        = params['paths']['rcsAnalizer']
rcs_AnalizerArgs   = [params['paths']['rcs_AnalizerArgs']]

rcsoutput_folder   = "/home/parallels/CarolineRadarCodingToolbox1/Caroline-Radar-Coding-Toolbox/example/nl_groningen_s1_dsc_t037_test/RadarCoordinates/"
rcsAnalizer        = "/home/parallels/CarolineRadarCodingToolbox1/Caroline-Radar-Coding-Toolbox/detectDesignatedTargets.py"
rcs_AnalizerArgs   = ["/home/parallels/CarolineRadarCodingToolbox1/Caroline-Radar-Coding-Toolbox/example/aoi_groningen/nl_groningen_s1_dsct037.parms"]

print("")

## STEP 1 - LOAD COREISTERED STACK OF SLC IN XARRAY FORMAT
print("STEP 1 - LOAD CORREGISTRED STACK")
# Check if the folder is a .zarr directory
if doris_stack_folder.endswith(".zarr"):
    # Load the stack using io.read_slc_stack
    slc_stack   = io.read_slc_stack(doris_stack_folder)
    print(f"Loaded stack from Zarr format: {doris_stack_folder}")
else:
    # Ensure nlines_file and npixels_file paths are provided
    if not (nlines_file and npixels_file):
        raise ValueError(
            "nlines_file and npixels_file paths must be specified for non-Zarr stack folders."
        )
    # Load the stack using io.doris_sar_stack_to_xarray
    slc_stack   = io.doris_sar_stack_to_xarray(doris_stack_folder, nlines_file, npixels_file)
    print(f"Loaded stack from folder with nlines and npixels files: {doris_stack_folder}")

print("")


## STEP 2 - Run RADAR CODING TOOLBOX, get targets radard coordinates. 
# In case results are already produced by RC Toolbox read them in only.
# Note: Please configure your .parms and stacksRC.json beforehand 
print("STEP 2 - RUN RADAR CODING TOOLBOX, LOAD RESULTS")
# Run the script with rcs_AnalizerArgs if output does not exist already
rcsoutput_file   = io.ensure_rc_csv_exists(rcsoutput_folder, rcsAnalizer, rcs_AnalizerArgs)
print("")

## STEP 3 - LOAD THE RC TOOLBOX OUTPUT CSV FILE IN XARRAY FORMAT
print("STEP 3 - LOAD RADAR CODING TOOLBOX RESULTS")
dt_df, dt_dates  = io.load_stm_rcscsv(rcsoutput_file) # Extract information from RC Toolbox output csv

# Convert to xarray.Dataset
targets          = xr.Dataset(
                        {
                            "existing_flag": (["target", "time"], dt_df.iloc[:, 7:].values),
                        },
                        coords={
                            "target": dt_df["ID"].values,
                            "range": ("range", dt_df["Range"].values),
                            "azimuth": ("azimuth", dt_df["Azimuth"].values),
                            "lat": ("latitude", dt_df["Lat"].values),
                            "lon": ("longitude", dt_df["Lon"].values),
                            "height": ("height", dt_df["Height"].values),
                            "time": pd.to_datetime(dt_dates, format="%Y%m%d"),
                        },
                    )
print("")

## STEP 4 - EXTRACT THE TARGETS FROM THE STACK 

print("STEP 4 - EXTRACT DESIGNATED TARGETS FROM STACK")
# Create target coordinates set
target_coords    = set(zip(
                    targets['azimuth'].data.flatten(),
                    targets['range'].data.flatten()
                    ))

# Create meshgrid for azimuth and range in slc_stack
slc_azimuth, slc_range = np.meshgrid(
                            slc_stack['azimuth'].data,
                            slc_stack['range'].data,
                            indexing='ij'
                            )

# Find matching coordinates in slc_stack
matching_indices = [
                    (az_idx, rg_idx) 
                    for az_idx, az in enumerate(slc_stack['azimuth'].data) 
                    for rg_idx, rg in enumerate(slc_stack['range'].data)
                    if (az, rg) in target_coords
                    ]

if not matching_indices:
    print("No matches found.")
else:
    print(f"Found {len(matching_indices)} matching targets.")

## In case one is interested in querring for 1 target only, use this commented out part.

# # Step 1: Identify the index of the matching azimuth and range in the slc_stack
# azimuth_idx = matching_indices[0][0]
# range_idx = matching_indices[0][1]

# # Step 2: Slice the slc_stack to extract all time-related data for the matching point
# matching_point_slc_data = slc_stack.isel(azimuth=azimuth_idx, range=range_idx)

## For querring more targets
# Call the function to extract the data for the matching targets
matched_slc_targets_dict, lat_vals, lon_vals, target_names = io.extract_dttarget_data_from_slc(slc_stack, matching_indices, targets)

# Convert lists to numpy arrays
lat_vals        = np.array(lat_vals)
lon_vals        = np.array(lon_vals)
target_names    = np.array(target_names)

# Add lat/lon and target_names to the dataset dictionary, without the time dimension
matched_slc_targets_dict['lat']      = lat_vals
matched_slc_targets_dict['lon']      = lon_vals
matched_slc_targets_dict['target_name'] = target_names  # New coordinate for target names

# Convert the dictionary into a new xarray.Dataset
matching_scatterer_slc_data_f        = {var: (('target', 'time'), np.array(matched_slc_targets_dict[var])) 
                                    for var in matched_slc_targets_dict if var not in ['lat', 'lon', 'target_name']}

# Add lat and lon as 1D coordinates for the target dimension
matching_scatterer_slc_data_f['lat'] = (['target'], lat_vals)
matching_scatterer_slc_data_f['lon'] = (['target'], lon_vals)

# Create the target_name dictionary for the attribute
target_name_dict = {
    name: (lon, lat) for name, lon, lat in zip(target_names, lon_vals, lat_vals)
}

# Create the final xarray.Dataset with the target dimension
xar_matching_scatterer_slc_data_f = xr.Dataset(matching_scatterer_slc_data_f)

# Add the target_name dictionary as an attribute
xar_matching_scatterer_slc_data_f.attrs['target_name'] = target_name_dict

# Rename the 'target' dimension to 'space' in the xarray Dataset
xar_matching_scatterer_slc_data_f = xar_matching_scatterer_slc_data_f.rename({'target': 'space'})

print("")
print(xar_matching_scatterer_slc_data_f)