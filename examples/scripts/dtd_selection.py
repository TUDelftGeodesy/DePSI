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
import yaml
import sarxarray
import subprocess
import numpy as np
import xarray as xr
import pandas as pd
from depsi import io
sys.path.append('/home/parallels/Sprint_Mobyle/DePSI_group/depsi/')


## FUNCTIONS
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

if __name__ == "__main__":
    ## STEP 0 - LOAD INPUT PARAMETERS
    # Specify folders and files to be used:
    print("STEP 0 - LOAD INPUT PARAMS")
    # Load the YAML parameters containing folders and files to be used:
    params = load_ymlparams("/home/parallels/Sprint_Mobyle/DePSI_group/examples/scripts/data/dtd_params.yml")
    
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
        slc_stack = io.read_slc_stack(filename=doris_stack_folder, engine="doris", nlines_file=nlines_file, npixels_file=npixels_file, chunks=(500, 500))
        print(f"Loaded stack from folder with nlines and npixels files: {doris_stack_folder}")
    
    print("")
    
    
    ## STEP 2 - Run RADAR CODING TOOLBOX, get targets radard coordinates. 
    # In case results are already produced by RC Toolbox read them in only.
    # Note: Please configure your .parms and stacksRC.json beforehand 
    print("STEP 2 - RUN RADAR CODING TOOLBOX, LOAD RESULTS")
    # Run the script with rcs_AnalizerArgs if output does not exist already
    rcsoutput_file   = ensure_rc_csv_exists(rcsoutput_folder, rcsAnalizer, rcs_AnalizerArgs)
    print("")
    
    ## STEP 3 - LOAD THE RC TOOLBOX OUTPUT CSV FILE IN XARRAY FORMAT
    print("STEP 3 - LOAD RADAR CODING TOOLBOX RESULTS")
    dt_df, dt_dates  = io.read_rcs_csv(rcsoutput_file) # Extract information from RC Toolbox output csv
    
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
                                "detection_flag": (["target", "time"], dt_df.iloc[:, 7:].values),
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
    
    # # Find matching coordinates in slc_stack
    # matching_indices = [
    #                     (az_idx, rg_idx) 
    #                     for az_idx, az in enumerate(slc_stack['azimuth'].data) 
    #                     for rg_idx, rg in enumerate(slc_stack['range'].data)
    #                     if (az, rg) in target_coords
    #                     ]
    
    # if not matching_indices:
    #     print("No matches found.")
    # else:
    #     print(f"Found {len(matching_indices)} matching targets.")
    
    # Find matching coordinates in slc_stack (value-based instead of index-based)
    matching_coords = [
        (az, rg)
        for az in slc_stack['azimuth'].data
        for rg in slc_stack['range'].data
        if (az, rg) in target_coords
    ]
    
    if not matching_coords:
        print("No matches found.")
    else:
        print(f"Found {len(matching_coords)} matching targets.")
    
    
    ## In case one is interested in querring for 1 target only, use this commented out part.
    
    # # Step 1: Identify the index of the matching azimuth and range in the slc_stack
    # azimuth_idx = matching_indices[0][0]
    # range_idx = matching_indices[0][1]
    
    # # Step 2: Slice the slc_stack to extract all time-related data for the matching point
    # matching_point_slc_data = slc_stack.isel(azimuth=azimuth_idx, range=range_idx)
    
    ## For querring more targets
    # Call the function to extract the data for the matching targets
    matched_slc_targets_dict, lat_vals, lon_vals, target_names, detection_flag_space_time, target_space_indices = io.extract_dttarget_data_from_slc(
        slc_stack, matching_coords, targets, verbose=True
    )
    
    # Convert lists to numpy arrays
    lat_vals     = np.array(lat_vals)
    lon_vals     = np.array(lon_vals)
    target_names = np.array(target_names)
    
    # Add 1D coordinate variables (target coordinates) to the dictionary
    matched_slc_targets_dict['lat']     = lat_vals
    matched_slc_targets_dict['lon']     = lon_vals
    matched_slc_targets_dict['azimuth'] = np.array(matched_slc_targets_dict['azimuth'])
    matched_slc_targets_dict['range']   = np.array(matched_slc_targets_dict['range'])
    
    # Build dataset from 2D variables (time series)
    matching_scatterer_slc_data_f = {
        var: (('target', 'time'), np.array(matched_slc_targets_dict[var]))
        for var in matched_slc_targets_dict
        if np.array(matched_slc_targets_dict[var]).size > 0 
           and np.array(matched_slc_targets_dict[var]).ndim == 2
    }
    
    # Add 1D variables (target coordinates)
    for var in ['lat', 'lon', 'azimuth', 'range']:
        matching_scatterer_slc_data_f[var] = (['target'], np.array(matched_slc_targets_dict[var]))
    
    # Create target_name dictionary for attributes
    target_name_dict = {
        name: (lon, lat) for name, lon, lat in zip(target_names, lon_vals, lat_vals)
    }
    
    # Create Dataset
    xar_matching_scatterer_slc_data_f = xr.Dataset(matching_scatterer_slc_data_f)
    
    # Add target_name dictionary as an attribute (not a variable)
    xar_matching_scatterer_slc_data_f.attrs['target_name'] = target_name_dict
    
    # Rename 'target' dimension to 'space'
    xar_matching_scatterer_slc_data_f = xar_matching_scatterer_slc_data_f.rename({'target': 'space'})
    
    # Add the 'time' coordinate
    xar_matching_scatterer_slc_data_f = xar_matching_scatterer_slc_data_f.assign_coords(
        time=slc_stack['time'].data
    )

    print(xar_matching_scatterer_slc_data_f)


## Note: Key variables
# xar_matching_scatterer_slc_data_f  -> xarray Dataset containing SLC values for the matching targets
#                                       along with 0/1 detection flags per target.
# detection_flag_space_time          -> xarray DataArray of shape (azimuth, range, time) containing
#                                       0/1 detection flags for all pixels in the original SLC stack,
#                                       indicating whether a target was detected at each time step.
