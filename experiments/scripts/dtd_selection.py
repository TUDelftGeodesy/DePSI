"""
Executable that identifies designated targets within a corregistered SLC stack
by loading the stack, loading the RadarCoding Toolbox output, querrying for targets in the stack

An example of the input parameter file (dtd_params.yml),
an SLC stack in zarr frmat (nl_groningen_s1_dsc_t037_haren),
the output of the RadarCoding Toolbox (s1_dsc037_RC.csv) can be found at:
https://figshare.com/ndownloader/files/51712790

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

sys.path.append("/home/parallels/Sprint_Mobyle/DePSI_group/depsi/")


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
    params = load_ymlparams("/home/parallels/Sprint_Mobyle/DePSI_group/experiments/scripts/data/dtd_params.yml")

    if params is None:
        raise ValueError("Failed to load parameters from ''.")

    # Access paths from the YAML file
    doris_stack_folder = params["paths"]["doris_stack_folder"]
    nlines_file = params["paths"]["nlines_file"]
    npixels_file = params["paths"]["npixels_file"]

    rcsoutput_folder = params["paths"]["rcs_output_folder"]
    rcsAnalizer = params["paths"]["rcsAnalizer"]
    rcs_AnalizerArgs = [params["paths"]["rcs_AnalizerArgs"]]

    print("")

    ## STEP 1 - LOAD COREISTERED STACK OF SLC IN XARRAY FORMAT
    print("STEP 1 - LOAD CORREGISTRED STACK")
    # Check if the folder is a .zarr directory
    if doris_stack_folder.endswith(".zarr"):
        # Load the stack using io.read_slc_stack
        slc_stack = io.read_slc_stack(doris_stack_folder)
        print(f"Loaded stack from Zarr format: {doris_stack_folder}")
    else:
        # Ensure nlines_file and npixels_file paths are provided
        if not (nlines_file and npixels_file):
            raise ValueError("nlines_file and npixels_file paths must be specified for non-Zarr stack folders.")
        # Load the stack using io.doris_sar_stack_to_xarray
        slc_stack = io.read_slc_stack(
            filename=doris_stack_folder,
            engine="doris",
            nlines_file=nlines_file,
            npixels_file=npixels_file,
            chunks=(500, 500),
        )
        print(f"Loaded stack from folder with nlines and npixels files: {doris_stack_folder}")

    print("")

    ## STEP 2 - Run RADAR CODING TOOLBOX, get targets radard coordinates.
    # In case results are already produced by RC Toolbox read them in only.
    # Note: Please configure your .parms and stacksRC.json beforehand
    print("STEP 2 - RUN RADAR CODING TOOLBOX, LOAD RESULTS")
    # Run the script with rcs_AnalizerArgs if output does not exist already
    rcsoutput_file = ensure_rc_csv_exists(rcsoutput_folder, rcsAnalizer, rcs_AnalizerArgs)
    print("")

    ## STEP 3 - LOAD THE RC TOOLBOX OUTPUT CSV FILE IN XARRAY FORMAT
    print("STEP 3 - LOAD RADAR CODING TOOLBOX RESULTS")
    dt_df, dt_dates = io.read_rcs_csv(rcsoutput_file)  # Extract information from RC Toolbox output csv

    # Convert to xarray.Dataset
    targets = xr.Dataset(
        {
            "existing_flag": (["space", "time"], dt_df.iloc[:, 7:].values),
        },
        coords={
            "target": ("space", dt_df["ID"].to_numpy(dtype=object)),
            "range": ("space", dt_df["Range"].values),
            "azimuth": ("space", dt_df["Azimuth"].values),
            "lat": ("space", dt_df["Lat"].values),
            "lon": ("space", dt_df["Lon"].values),
            "height": ("space", dt_df["Height"].values),
            "time": pd.to_datetime(dt_dates, format="%Y%m%d"),
            "detection_flag": (["space", "time"], dt_df.iloc[:, 7:].values),
        },
    )
    ## STEP 4 - EXTRACT THE TARGETS FROM THE STACK
    matching_scatterers = io.get_targets_from_slc(slc_stack, targets)

    # Re-organize the target names to the attributes, since they are currently string coords
    target_names = dict()
    for name in matching_scatterers.coords["target"].values:
        tgt = matching_scatterers.where(matching_scatterers["target"] == name, drop=True)
        target_names[name] = {
            "lon": tgt["lon"].values.item(),
            "lat": tgt["lat"].values.item(),
            "azimuth": tgt["azimuth"].values.item(),
            "range": tgt["range"].values.item(),
        }

    matching_scatterers.attrs = target_names

    # Print the number of identified targets and their names
    print(f"Identified {len(matching_scatterers.attrs)} targets in the SLC stack:")
    for name, info in matching_scatterers.attrs.items():
        print(f" - {name}: (lon={info['lon']}, lat={info['lat']}, az={info['azimuth']}, rg={info['range']})")

    matching_scatterers = matching_scatterers.drop(["target"])


## Note: Key variables
# matching_scatterers  -> xarray Dataset containing SLC values for the matching targets
#                                       along with 0/1 detection flags per target.
