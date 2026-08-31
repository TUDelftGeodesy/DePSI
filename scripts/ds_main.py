import os
import sys
from datetime import datetime
from pathlib import Path

import geopandas as gpd
import h5py
import numpy as np
import pandas as pd
import sarxarray
import xarray as xr

from depsi.arc_estimation import periodogram
from depsi.atmosphere_estimation import estimate_atmosphere_phase
from depsi.classification import network_stm_selection, ps_selection
from depsi.densification import densification
from depsi.ds import assign_parcel_id, ds_phase_estimation, _open_datatree_compat, select_common_fop_ref, ps_ds_arc
from depsi.io import (
    export_to_csv,
    export_to_skygeo_portal,
    export_convex_hull_to_shapefile,
    export_to_shapefile,
    read_knmi_txt,
    read_slc_stack
)
from depsi.model_estimation import MODEL_NAMES_PARAMS, estimate_model_params
from depsi.network import form_network, spatial_integration
from depsi.point_quality import compute_spatiotemporal_consistency, detect_side_lobes
from depsi.postprocessing import stm_point_filter
from depsi.transformations import radar_to_latlonh
from depsi.utils import (
    add_stm_time_deltas,
    convert_geographic_coords_to_euclidean,
    crop_slc_spacetime,
    stm_compute_single_time_differences,
)
from depsi.viewing_geometry import add_local_viewing_geometry


#########################################
#                                       #
#           INPUT VARIABLES             #
#                                       #
#########################################
ROOT_DIR = ""
proj_dir = ROOT_DIR + "/projects/nieuwolda/insar/"
run_name = "runp_test"
log_filename = None #"runp_test_1sel.txt"  # Set to None to disable logging

# - Contextual data path
aoi_shp_filepath = ROOT_DIR + "/projects/nieuwolda/contextual_data/aoi/nieuwolda_aoi.shp"
parcel_shp_filepath = ROOT_DIR + "/projects/nieuwolda/contextual_data/parcels/nieuwolda_attributes_for_depsi_test.shp"
meteo_dir = ROOT_DIR + "/projects/nieuwolda/contextual_data/knmi/"

# - Stack
stack_root_dir = ROOT_DIR + "/projects/nieuwolda/stacks/"
stack_prefix = "nl_nieuwolda"
mission = "s1"
wavelength = 0.055465763 # m
stack_ids = [
    "s1_asc_t088",
    "s1_dsc_t037",
]
metadata_paths = [os.path.join(stack_root_dir, f"{stack_prefix}_{stack_id}", "master.res") for stack_id in stack_ids]
mother_epochs = [
    datetime(2020, 3, 25),
    datetime(2020, 3, 28),
]
start_date = datetime(2020, 1, 1)
end_date = datetime(2020, 12, 31)

# - PS selection
ps_selection_method = "nad"
threshold = 0.35
chunks_ps_selection = 5000
start_date_ps_selection = None
end_date_ps_selection = None

# - Sidelobes
do_sidelobe_detection = False
max_pixel_dist = 2
min_correlation = 0.90

# - Network
network_crs = "radar"
network_x_crds = "azimuth"
network_y_crds = "range"
min_point_distance = 500  # meters
max_arc_length = 0.045  # degrees (lat/lon)
network_formation_method = "redundant"
min_network_links = 16
network_partition_number = 8
network_dphase_method = "subtract"
min_periodogram_iterations = 10
arc_quality_threshold = 0.5  # ensemble coherence
model_types = ["offset", "velocity", "height"]
reference_point_index = None

# - APS variables
euclidean_epsg_code_number = 28992
atmo_unmodeled_displacement_filter_length = 0.5
atmo_unmodeled_displacement_sampling_rate = 1000
atmo_unmodeled_displacement_filter_type = "gaussian"
atmo_kriging_n_nearest_neighbours = None
atmo_kriging_backend = "vectorized"
atmo_empirical_variogram_method = "standard"
atmo_empirical_variogram_nlags = 50
atmo_empirical_variogram_cutoff_distance = 10000.0
atmo_variogram_model = "gaussian"
atmo_variogram_drift_terms = "regional_linear"

# - Densification
n_densification_connections = 1

# - Spatio-temporal consistency
stc_min_dist = 30  # meters
stc_max_dist = 100  # meters

# - Viewing geometry
orbit_mode = "IWS"
orbit_resolution = 0.01
orbit_file = "../config/drama/S1_XTI.cfg"

# - Output settings
filter_dict = {"h": [-2000, 2000]}
output_types = ["csv_web_portal", "csv", "shapefile", "convex_hull"]

# - DS analysis
ds_id_list = None # Options: None, list of ds_ids. Specify only to reestimate esm phase for a subset of parcels
ds_multilooking_window = "polygon"
ds_min_cells = 40
ds_shp_test = None  # Options: None, "ks_test"
ds_min_group_size = 4
ds_coh_threshold = 0.2
ds_min_seg_len = 10
igrs_codes = None
igrs_locs = None

# - zarr export & checkpoints
# - Example: ps_stm_s1_asc_t088_4atmo.zarr
# -          ds_stm_s1_asc_t088.zarr
ps_stm_save_name = "ps_stm_{}_{}.zarr"
ds_dtree_save_name = "ds_dtree_{}.zarr"
arc_dtree_save_name = "arc_dtree_{}.zarr"
########################################


# - Create and add directory lists
run_dir = os.path.join(proj_dir, run_name)
phase_est_dir = os.path.join(run_dir, "phase_estimation/")
stm_dir = os.path.join(run_dir, "stm/")

if not os.path.exists(run_dir):
    os.makedirs(run_dir)
if not os.path.exists(phase_est_dir):
    os.makedirs(phase_est_dir)
if not os.path.exists(stm_dir):
    os.makedirs(stm_dir)

# - Log file
if log_filename is not None:
    logdir = os.path.join(run_dir, "logs/")
    if not os.path.exists(logdir):
        os.makedirs(logdir)
    sys.stdout = open(os.path.join(logdir, log_filename), "w")

# - Load parcel shapefile
gdf_parcel = gpd.read_file(parcel_shp_filepath)
if "id"in gdf_parcel.columns:
    ds_ids = gdf_parcel["id"].unique()
elif "int_id" in gdf_parcel.columns:
    ds_ids = gdf_parcel["int_id"].unique()
else:
    raise ValueError("Parcel shapefile must contain 'id' or 'int_id' attribute.")

# - Load meteorological data
filelist = sorted(Path(meteo_dir).glob("**/etmgeg_*.txt"))
f1 = [None] * len(filelist)
for i, file in enumerate(filelist):
    f1[i] = read_knmi_txt(file)
df_meteo = pd.concat(f1)


# ========================================= #
#           1. PS DS Selection              #
# ========================================= #
print("========== 1. PS and DS selection ...")
for i, stack_id in enumerate(stack_ids):
    print(f"Processing stack {stack_id} ...")
    metadata = sarxarray.read_metadata(metadata_paths[i], driver="doris5")

    pe_dir = os.path.join(phase_est_dir, stack_id)
    if not os.path.exists(pe_dir):
        os.makedirs(pe_dir)

    print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Load SLC stack ...")
    slc_filename = f"{stack_prefix}_{stack_id}"
    slc_path = os.path.join(stack_root_dir, slc_filename, slc_filename + ".zarr")
    slc_stack = read_slc_stack(slc_path)
    print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} SLC stack are loaded.")

    # - PS selection
    print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} PS selection ...")
    ps_filepath_1sel = os.path.join(stm_dir, ps_stm_save_name.format(stack_id, "1sel"))
    if os.path.isdir(ps_filepath_1sel):
        print("PS stm file already exist. Skipping ...")
    else:
        ps_stm = ps_selection(
                slcs=slc_stack,
                method=ps_selection_method,
                threshold=threshold,
                output_chunks=chunks_ps_selection,
            )

        print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Computing single differences and temporal baseline...")
        ps_stm = stm_compute_single_time_differences(
            ps_stm,
            single_difference_mother=mother_epochs[i]
        )
        ps_stm = add_stm_time_deltas(ps_stm)
        ps_stm = ps_stm.assign({"temporal_baseline": (
            ["time"],
            ps_stm.years_since_first_img.values - ps_stm.years_since_first_img.sel(time=ps_stm.ps_sd_mother).values
        )})
        
        # - Sidelobe removal
        print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Removing sidelobes...")
        if do_sidelobe_detection:
            side_lobes_array, _ = detect_side_lobes(ps_stm, max_pixel_dist, min_correlation, "sd_complex", "sd_amplitude_unnormalized")
            mask_sidelobes = np.ones(len(ps_stm.space), dtype=bool)
            mask_sidelobes[side_lobes_array] = False
        
        else:
            print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Skipping sidelobe detection...")
            mask_sidelobes = np.ones(len(ps_stm.space), dtype=bool)
        
        ps_stm = ps_stm.isel(space=mask_sidelobes)
        
        print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Removed {np.sum(~mask_sidelobes)} sidelobes from the STM, "
              f"{len(ps_stm.space)} remain.")
        
        # - Geocoding
        print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Starting geocoding...")
        ps_stm = ps_stm.assign_attrs({"wavelength": metadata["wavelength"]})
        
        latlonh = radar_to_latlonh(
            azimuth_coords=ps_stm.azimuth.values,
            range_coords=ps_stm.range.values,
            elevation=ps_stm.h.values,
            metadata=metadata,
        )
        
        ps_stm["lat"].data = latlonh[0].flatten()
        ps_stm["lon"].data = latlonh[1].flatten()
        ps_stm["h"].data = latlonh[2].flatten()
        
        # - Add Euclidean coords in preparation for atmosphere estimation (here so it is also present in the first order network)
        print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Adding Euclidean coordinates...")
        ps_crd_x, ps_crd_y = convert_geographic_coords_to_euclidean(
            ps_stm["lon"].values,
            ps_stm["lat"].values,
            target_crs=f"EPSG:{euclidean_epsg_code_number}"
        )
        ps_stm = ps_stm.assign_coords({
            f"x_euclidean_proj_epsg{euclidean_epsg_code_number}": (["space"], ps_crd_x),
            f"y_euclidean_proj_epsg{euclidean_epsg_code_number}": (["space"], ps_crd_y),
        })

        ps_stm = ps_stm.chunk({"time": -1, "space": "auto"})
        ps_stm.to_zarr(ps_filepath_1sel, mode="w", zarr_format=2)
    print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} PS selection done.")

    # - DS selection
    print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} DS selection ...")
    ds_filepath = os.path.join(stm_dir, ds_dtree_save_name.format(stack_id))
    ds_mask_filepath = os.path.join(pe_dir, "id_pixel_" + stack_id + ".h5")

    if os.path.isfile(ds_mask_filepath):
        print("DS mask h5 file already exist. Loading ...")
        with h5py.File(ds_mask_filepath, "r") as f:
            ds_mask = f["ds_mask"][()]
            ds_ids = f["ds_id"][()]
            centroid = f["centroid"][()]
            az_centroid = f["az_centroid"][()]
            rg_centroid = f["rg_centroid"][()]
    else:
        print("Assigning radar pixel_id to the corresponding parcel ...")
        ds_mask, ds_ids, centroid, az_centroid, rg_centroid = assign_parcel_id(
            slc_stack,
            parcel_shp_filepath,
            ds_min_cells,
            ps_stm,
        )

        print("Saving pixel_id into an HDF file ...")
        dataset_name=["ds_mask", "ds_ids", "centroid", "az_centroid", "rg_centroid"]
        dataset=[ds_mask, ds_ids, centroid, az_centroid, rg_centroid]
        data_dict = {}
        for j, dset_name in enumerate(dataset_name):
            data_dict.update({dset_name: dataset[j]})
        with h5py.File(ds_mask_filepath, "w") as f:
            for dset_name in data_dict:
                f.create_dataset(dset_name, data=data_dict[dset_name])

    # - DS phase estimation
    if os.path.isdir(ds_filepath) and ds_id_list is None:
        print("DS datatree file already exist. Skipping ...")
    else:
        print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} DS phase estimation ...")
        ds_dtree = ds_phase_estimation(
            slc_stack,
            metadata,
            stack_id,
            wavelength,
            centroid,
            az_centroid,
            rg_centroid,
            gdf_parcel,
            ds_mask,
            ds_ids,
            ds_id_list=ds_id_list,
            shp_test=ds_shp_test,
            min_seg_len=ds_min_seg_len,
            coh_threshold=ds_coh_threshold,
            ds_filepath=ds_filepath,
        )
        
        # - Add local incident angle
        print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Add local incidence angle...")
        ds_stm = ds_dtree["ds_stm"].to_dataset()
        if mission == "s1" and ds_stm["local_incident_angle"].isnull().all() and orbit_file is not None:
            ds_stm = add_local_viewing_geometry(
                stm=ds_stm,
                orbit_config_file=orbit_file,
                orbit_res=orbit_resolution,
                orbit_mode=orbit_mode,
                orbit=stack_id,
            )
            z2ph = (-4 * np.pi) / (wavelength) * np.cos(np.radians(ds_stm["local_incident_angle"].values))
            ds_stm = ds_stm.assign({"z2ph": (["space"], z2ph)})
        
        # - Add projected coordinates for atmosphere estimation
        print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Adding Euclidean coordinates...")
        ds_crd_x, ds_crd_y = convert_geographic_coords_to_euclidean(
            ds_stm["lon"].values,
            ds_stm["lat"].values,
            target_crs=f"EPSG:{euclidean_epsg_code_number}"
        )
        ds_stm = ds_stm.assign_coords({
            f"x_euclidean_proj_epsg{euclidean_epsg_code_number}": (["space"], ds_crd_x),
            f"y_euclidean_proj_epsg{euclidean_epsg_code_number}": (["space"], ds_crd_y),
        })
        ds_dtree["ds_stm"] = ds_stm
        ds_dtree.to_zarr(ds_filepath, mode="w")

    print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} DS selection done.")
    print(f"PS and DS selection for stack {stack_id} done.")
print("========== 1. PS and DS selection done.")


# ========================================= #
#        2. Atmospheric Phase Screen        #
# ========================================= #
print("========== 2. Atmospheric phase screen ...")
for i, stack_id in enumerate(stack_ids):
    print(f"Processing stack {stack_id} ...")
    metadata = sarxarray.read_metadata(metadata_paths[i], driver="doris5")

    ps_filepath_2atmo = os.path.join(stm_dir, ps_stm_save_name.format(stack_id, "2atmo"))
    if os.path.isdir(ps_filepath_2atmo):
        print("PS stm atmospheric phase screen file already exist. Skipping ...")
    else:
        ps_filepath_1sel = os.path.join(stm_dir, ps_stm_save_name.format(stack_id, "1sel"))
        if os.path.isdir(ps_filepath_1sel):
            ps_stm = xr.open_zarr(ps_filepath_1sel, consolidated=True)
        else:
            raise FileNotFoundError(f"PS stm file {ps_filepath_1sel} not found. Please run PS selection first.")
        
        # - Network construction
        print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Removing mother from network STM...")
        mother_epoch_index = np.where(ps_stm.time.values == ps_stm.sel(time=ps_stm.ps_sd_mother).time.values)[0][0]
        non_mother = [True] * len(ps_stm.time.values)
        non_mother[mother_epoch_index] = False
        ps_stm_without_mother_epoch = ps_stm.isel(time=non_mother)
        
        print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Selecting first-order network points...")
        stm_network_pnts = network_stm_selection(
            stm=ps_stm_without_mother_epoch,
            min_dist=min_point_distance,
            crs=network_crs,
            azimuth_spacing=metadata["azimuth_pixel_spacing"],
            range_spacing=metadata["range_pixel_spacing"],
            sortby_var="time_selection_" + ps_selection_method,
            x_var=network_x_crds,
            y_var=network_y_crds,
            include_index=None
        )
        
        print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Forming first-order network...")
        stm_network_arcs = form_network(
            stm_network_pnts,
            key_phase='sd_phase',
            key_h2ph='h2ph',
            key_Btemporal='temporal_baseline',
            max_length=max_arc_length,
            key_xcrds="lon",
            key_ycrds="lat",
            network_method=network_formation_method,
            min_links=min_network_links,
            num_partitions=network_partition_number,
            dphase_method=network_dphase_method,
        )
        print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Formed first-order network with {len(stm_network_pnts.space)} "
              f"points and {len(stm_network_arcs.space)} arcs.")
        
        print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Estimating ambiguities through periodogram...")
        _, ambiguities, _, _, ens_coh = periodogram(
            stm_network_arcs,
            key_dphase='d_phase',
            key_h2ph='h2ph',
            key_Btemporal='Btemp',
            min_steps=min_periodogram_iterations,
                    )
        stm_network_arcs["ambiguities"] = ambiguities
        stm_network_arcs["temp_coh"] = ens_coh

        print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Computing periodogram output...")
        stm_network_pnts = stm_network_pnts.compute()
        stm_network_arcs = stm_network_arcs.compute()
                
        print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Integrating the spatial network...")
        _, stm_pnts_output = spatial_integration(
            stm_network_pnts,
            stm_network_arcs,
            key_arc_quality="temp_coh",
            threshold_arc_quality=arc_quality_threshold,
            idx_refpnt=reference_point_index,
        )

        # - Atmosphere estimation
        print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Estimating the mother atmosphere...")
        stm_pnts_output, _ = estimate_model_params(
            stm=stm_pnts_output,
            models=model_types,
            key_observations="unwrapped_phase",
            key_h2ph="sd_h2ph",
            key_time="temporal_baseline"
        )

        stm_pnts_output = stm_pnts_output.rename_vars({
            "pnt_offset": "mother_atmosphere",
        })
        # - Rename the space dimension since it will otherwise clash with the space dimension of the predicted coordinates
        stm_pnts_output = stm_pnts_output.rename_dims({"space": "space_fon"})

        print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Creating the atmosphere prediction coordinate dataset...")
        # - PS coordinates for the prediction of the atmospheric phase screens
        print("PS coordinates...")
        ps_atmosphere_prediction_coords = xr.Dataset(
            coords={
                f"x_euclidean_proj_epsg{euclidean_epsg_code_number}": (
                    "space", ps_stm[f"x_euclidean_proj_epsg{euclidean_epsg_code_number}"].values
                ),
                f"y_euclidean_proj_epsg{euclidean_epsg_code_number}": (
                    "space", ps_stm[f"y_euclidean_proj_epsg{euclidean_epsg_code_number}"].values
                ),
            }
        )
        # - DS coordinates here so that it will be included in the evaluation of the atmospheric phase screens
        ds_filepath = os.path.join(stm_dir, ds_dtree_save_name.format(stack_id))
        if os.path.isdir(ds_filepath):
            ds_dtree = _open_datatree_compat(ds_filepath)
            ds_stm = ds_dtree["ds_stm"].to_dataset()
        else:
            raise FileNotFoundError(f"DS datatree file {ds_filepath} not found. Please run DS selection first.")
        print("DS coordinates...")
        ds_atmosphere_prediction_coords = xr.Dataset(
            coords={
                f"x_euclidean_proj_epsg{euclidean_epsg_code_number}": (
                    "space", ds_stm[f"x_euclidean_proj_epsg{euclidean_epsg_code_number}"].values
                ),
                f"y_euclidean_proj_epsg{euclidean_epsg_code_number}": (
                    "space", ds_stm[f"y_euclidean_proj_epsg{euclidean_epsg_code_number}"].values
                ),
            }
        )

        # - Merge them so that we have a single prediction coordinates including PS and DS
        # - Extract the number of PS and DS points for later splitting of the atmospheric phase screens result
        print("Merging the atmosphere prediction coordinate dataset...")
        nps, nds = len(ps_stm.space), len(ds_stm.space)
        atmosphere_prediction_coords = xr.concat([ps_atmosphere_prediction_coords, ds_atmosphere_prediction_coords], dim="space")

        print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Estimating the atmospheric phase screens...")
        atmospheric_phase_screens = estimate_atmosphere_phase(
            stm=stm_pnts_output,
            prediction_coords=atmosphere_prediction_coords,
            psc_phase_residuals="phase_residuals",
            atmosphere_mother="mother_atmosphere",
            key_Btemporal="temporal_baseline",
            unmodeled_displacement_args={
                "filter_length": atmo_unmodeled_displacement_filter_length,
                "sampling_rate": atmo_unmodeled_displacement_sampling_rate,
                "filter_type": atmo_unmodeled_displacement_filter_type,
            },
            kriging_args={
                "n_nearest_neighbors": atmo_kriging_n_nearest_neighbours,
                "kriging_backend": atmo_kriging_backend,
                "empirical_variogram_args": {
                    "method": atmo_empirical_variogram_method,
                    "nlags": atmo_empirical_variogram_nlags,
                    "cutoff": atmo_empirical_variogram_cutoff_distance,
                },
                "variogram_args": {
                    "variogram_model": atmo_variogram_model,
                    "variogram_parameters": None,  # to be estimated from the empirical variogram
                    "drift_terms": atmo_variogram_drift_terms,
                },
            },
            stm_coords_metadata={
                "mode": "euclidean",
                "x_label": f"x_euclidean_proj_epsg{euclidean_epsg_code_number}",
                "y_label": f"y_euclidean_proj_epsg{euclidean_epsg_code_number}",
                "space_dim_name": "space_fon",
            },
            prediction_coords_metadata={
                "mode": "euclidean",
                "x_label": f"x_euclidean_proj_epsg{euclidean_epsg_code_number}",
                "y_label": f"y_euclidean_proj_epsg{euclidean_epsg_code_number}",
                "space_dim_name": "space"
            },
        )

        # - Split the atmosphere for PS and DS points
        print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Splitting the atmospheric phase screens for PS and DS points...")
        ps_aps = atmospheric_phase_screens.isel(space=slice(0, nps))
        ds_aps = atmospheric_phase_screens.isel(space=slice(nps, nps + nds))

        # - Create the zero layer for the mother atmosphere to be appended to the interferometric atmospheric phase screens
        ps_mother_atmo_stm = xr.Dataset(
            data_vars={
                "atmosphere_predicted": ("space", np.zeros((len(ps_stm.space), ))),
                "atmosphere_sigmasq": ("space", np.zeros((len(ps_stm.space), ))), },
            coords=ps_stm["phase"].sel(time=ps_stm.ps_sd_mother).coords
        )

        ds_mother_atmo_stm = xr.Dataset(
            data_vars={
                "atmosphere_predicted": ("space", np.zeros((len(ds_stm.space), ))),
                "atmosphere_sigmasq": ("space", np.zeros((len(ds_stm.space), ))), },
            coords=ds_stm["ds_phi_esm_full"].sel(time=mother_epochs[i]).coords
        )

        ps_atmosphere = xr.concat([ps_aps, ps_mother_atmo_stm], dim="time").sortby("time")
        ds_atmosphere = xr.concat([ds_aps, ds_mother_atmo_stm], dim="time").sortby("time")

        # - Add the atmosphere into the original STM
        print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Adding the atmospheric phase screens into the STM...")
        ps_stm["atmosphere_predicted"] = ps_atmosphere["atmosphere_predicted"].transpose("space", "time")
        ps_stm["atmosphere_sigmasq"] = ps_atmosphere["atmosphere_sigmasq"].transpose("space", "time")
        ds_stm["atmosphere_predicted"] = ds_atmosphere["atmosphere_predicted"].transpose("space", "time")
        ds_stm["atmosphere_sigmasq"] = ds_atmosphere["atmosphere_sigmasq"].transpose("space", "time")

        # - Remove the atmosphere from the data
        print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Removing the atmospheric phase screens...")
        ps_stm["sd_phase_minus_atmo"] = (ps_stm["sd_phase"] - ps_stm["atmosphere_predicted"] + np.pi) % (2 * np.pi) - np.pi
        ds_stm["ds_phi_aps_full"] = (ds_stm["ds_phi_esm_full"] - ds_stm["atmosphere_predicted"] + np.pi) % (2 * np.pi) - np.pi

        print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Atmospheric phase screen done. Saving to zarr...")
        ps_stm = ps_stm.chunk({"time": -1, "space": "auto"})
        ps_stm.to_zarr(ps_filepath_2atmo, mode="w", zarr_format=2)
        ds_dtree["ds_stm"] = ds_stm
        ds_dtree.to_zarr(ds_filepath, mode="w")
        print("Saved!")
        print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} APS for stack {stack_id} done.")
print("========== 2. Atmospheric phase screen done.")


# ========================================= #
#          3a. PS First Order Points        #
# ========================================= #
print("========== 3a. PS first order points selection ...")
stm_network_pnts_list = []
for i, stack_id in enumerate(stack_ids):
    print(f"Processing stack {stack_id} ...")
    metadata = sarxarray.read_metadata(metadata_paths[i], driver="doris5")

    ps_filepath_3fop = os.path.join(stm_dir, ps_stm_save_name.format(stack_id, "3fop"))
    if os.path.isdir(ps_filepath_3fop):
        print("PS stm first order points files already exist. Loading ...")
        stm_network_pnts = xr.open_zarr(ps_filepath_3fop)
        stm_network_pnts_list.append(stm_network_pnts)
    else:
        ps_filepath_2atmo = os.path.join(stm_dir, ps_stm_save_name.format(stack_id, "2atmo"))
        if os.path.isdir(ps_filepath_2atmo):
            print("PS stm atmospheric phase screen file already exist. Loading ...")
            ps_stm = xr.open_zarr(ps_filepath_2atmo)
        else:
            raise FileNotFoundError(f"PS stm atmospheric phase screen file {ps_filepath_2atmo} not found. Please run APS estimation first.")

        # - First order points selection
        print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Removing mother from network STM...")
        mother_epoch_index = np.where(ps_stm.time.values == ps_stm.sel(time=ps_stm.ps_sd_mother).time.values)[0][0]
        non_mother = [True] * len(ps_stm.time.values)
        non_mother[mother_epoch_index] = False
        stm_atmo_corr_without_mother_epoch = ps_stm.isel(time=non_mother)
        stm_atmo_corr_without_mother_epoch = stm_atmo_corr_without_mother_epoch.chunk({"time": -1})
                
        print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Selecting first-order network points...")
        stm_network_pnts = network_stm_selection(
            stm=stm_atmo_corr_without_mother_epoch,
            min_dist=min_point_distance,
            crs=network_crs,
            azimuth_spacing=metadata["azimuth_pixel_spacing"],
            range_spacing=metadata["range_pixel_spacing"],
            sortby_var="time_selection_" + ps_selection_method,
            x_var=network_x_crds,
            y_var=network_y_crds,
            include_index=None
        )

        print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Append first-order network points to the list...")
        stm_network_pnts = stm_network_pnts.chunk({"time": -1, "space": "auto"})
        stm_network_pnts.to_zarr(ps_filepath_3fop, mode="w", zarr_format=2)
        stm_network_pnts_list.append(stm_network_pnts)

        print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Selecting first order PS points for stack {stack_id} done.")
print("========== 3a. PS first order points selection done.")


# ========================================= #
#        3b. Reference Point Selection      #
# ========================================= #
print("========== 3b. Common first order points and reference point selection ...")
print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Selecting common points and a reference point...")
ps_ref_idxs, common_fop_idxs = select_common_fop_ref(
    ps_stm_list=stm_network_pnts_list,
    proj_crs=euclidean_epsg_code_number,
    dist_ub=5.0,
    quality_var="full_ts_" + ps_selection_method,
)

print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Assigning the selected reference point index to stm...")
for i, stack_id in enumerate(stack_ids):
    print(f"Checking if the selected reference point index is present in the first order network points for stack {stack_id}...")
    ps_filepath_3fopn = os.path.join(stm_dir, ps_stm_save_name.format(stack_id, "3fopn"))
    if os.path.isdir(ps_filepath_3fopn):
        print(f"Reference point index has been assigned for stack {stack_id}. Skipping assignment.")
    else:
        stm_network_pnts = stm_network_pnts_list[i]
        stm_network_pnts_space = stm_network_pnts["space"].values
        ref_pnt_space = stm_network_pnts_space[ps_ref_idxs[i]]
        stm_network_pnts = stm_network_pnts.isel(space=common_fop_idxs[i])
        ps_ref_idx = np.where(stm_network_pnts["space"].values == ref_pnt_space)[0][0]
        stm_network_pnts = stm_network_pnts.assign_attrs({"ps_ref_idx": ps_ref_idx})
        stm_network_pnts = stm_network_pnts.chunk({"time": -1, "space": "auto"})
        stm_network_pnts.to_zarr(ps_filepath_3fopn, mode="w", zarr_format=2)

print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Reference point selection done.")
print("========== 3b. Common first order points and reference point selection done.")


# ========================================= #
#         3c. PS First Order Network        #
# ========================================= #
print("========== 3c. PS first order network spatial integration ...")
for i, stack_id in enumerate(stack_ids):
    print(f"Processing stack {stack_id} ...")
    metadata = sarxarray.read_metadata(metadata_paths[i], driver="doris5")

    ps_filepath_3fon = os.path.join(stm_dir, ps_stm_save_name.format(stack_id, "3fon"))
    if os.path.isdir(ps_filepath_3fon):
        print("PS stm first order network file already exist. Skipping...")
    else:
        ps_filepath_3fopn = os.path.join(stm_dir, ps_stm_save_name.format(stack_id, "3fopn"))
        if os.path.isdir(ps_filepath_3fopn):
            stm_network_pnts = xr.open_zarr(ps_filepath_3fopn)
        else:
            raise FileNotFoundError(f"PS stm first order points file {ps_filepath_3fopn} not found. Please run PS first order points selection first.")
        
        # - Network construction
        print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Forming first-order network...")
        stm_network_arcs = form_network(
            stm_network_pnts,
            key_phase='sd_phase_minus_atmo',
            key_h2ph='h2ph',
            key_Btemporal='temporal_baseline',
            max_length=max_arc_length,
            key_xcrds="lon",
            key_ycrds="lat",
            network_method=network_formation_method,
            min_links=min_network_links,
            num_partitions=network_partition_number,
            dphase_method=network_dphase_method,
        )
        print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Formed first-order network with {len(stm_network_pnts.space)} "
              f"points and {len(stm_network_arcs.space)} arcs.")
                        
        print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Estimating ambiguities through periodogram...")
        _, ambiguities, _, _, ens_coh = periodogram(
            stm_network_arcs,
            key_dphase='d_phase',
            key_h2ph='h2ph',
            key_Btemporal='Btemp',
            min_steps=min_periodogram_iterations,
        )
        stm_network_arcs["ambiguities"] = ambiguities
        stm_network_arcs["temp_coh"] = ens_coh
        
        print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Computing periodogram output...")
        stm_network_pnts = stm_network_pnts.compute()
        stm_network_arcs = stm_network_arcs.compute()

        print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Integrating the spatial network...")
        _, stm_pnts_output = spatial_integration(
            stm_network_pnts,
            stm_network_arcs,
            key_sdphase="sd_phase_minus_atmo",
            key_arc_quality="temp_coh",
            threshold_arc_quality=arc_quality_threshold,
            idx_refpnt=stm_network_pnts.attrs["ps_ref_idx"],
        )

        print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Spatial integration done. Saving to zarr...")
        stm_pnts_output = stm_pnts_output.chunk({"time": -1, "space": "auto"})
        stm_pnts_output.to_zarr(ps_filepath_3fon, mode="w", zarr_format=2)
        print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} First order network spatial integration for stack {stack_id} done.")
print("========== 3c. PS first order network spatial integration done.")


# ========================================= #
#            4. PS Densification            #
# ========================================= #
print("========== 4. PS densification ...")
for i, stack_id in enumerate(stack_ids):
    print(f"Processing stack {stack_id} ...")
    metadata = sarxarray.read_metadata(metadata_paths[i], driver="doris5")

    ps_filepath_4dens = os.path.join(stm_dir, ps_stm_save_name.format(stack_id, "4dens"))
    if os.path.isdir(ps_filepath_4dens):
        print("PS stm densification file already exist. Skipping...")
    else:
        ps_filepath_3fon = os.path.join(stm_dir, ps_stm_save_name.format(stack_id, "3fon"))
        if os.path.isdir(ps_filepath_3fon):
            stm_pnts_output = xr.open_zarr(ps_filepath_3fon)
        else:
            raise FileNotFoundError(f"PS stm first order network file {ps_filepath_3fon} not found. Please run PS first order network formation first.")

        ps_filepath_2atmo = os.path.join(stm_dir, ps_stm_save_name.format(stack_id, "2atmo"))
        if os.path.isdir(ps_filepath_2atmo):
            ps_stm = xr.open_zarr(ps_filepath_2atmo)
        else:
            raise FileNotFoundError(f"PS stm atmospheric phase screen file {ps_filepath_2atmo} not found. Please run APS estimation first.")

        # - Prepare input for densification
        mother_epoch_index = np.where(ps_stm.time.values == ps_stm.sel(time=ps_stm.ps_sd_mother).time.values)[0][0]
        non_mother = [True] * len(ps_stm.time.values)
        non_mother[mother_epoch_index] = False
        stm_atmo_corr_without_mother_epoch = ps_stm.isel(time=non_mother)
        stm_atmo_corr_without_mother_epoch = stm_atmo_corr_without_mother_epoch.chunk({"time": -1})

        print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Starting densification...")
        stm_densified = densification(
            stm_atmo_corr_without_mother_epoch,
            stm_pnts_output,
            idx_refpnt=stm_network_pnts.attrs["ps_ref_idx"],
            n_connections=n_densification_connections,
            key_sdphase="sd_phase_minus_atmo",
            key_h2ph="sd_h2ph",
            key_Btemporal="temporal_baseline"
        )

        print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Estimating the model parameters...")
        stm_densified, model_parameter_layer_names = estimate_model_params(
            stm=stm_densified,
            models=model_types,
            key_observations="unwrapped_phase",
            key_h2ph="sd_h2ph",
            key_time="temporal_baseline"
        )

        # - Geocoding
        print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Starting geocoding...")
        latlonh = radar_to_latlonh(
            azimuth_coords=stm_densified.azimuth.values,
            range_coords=stm_densified.range.values,
            elevation=stm_densified.pnt_height.values,
            metadata=metadata,
        )

        stm_densified["lat"].data = latlonh[0].flatten()
        stm_densified["lon"].data = latlonh[1].flatten()
        stm_densified["h"].data = latlonh[2].flatten()

        print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Calculating the unwrapped phase timeseries...")
        stm_densified["unwrapped_phase"].data = stm_densified["unwrapped_phase"].values

        print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Calculating spatiotemporal consistency...")
        stm_densified = compute_spatiotemporal_consistency(
            stm_densified,
            min_dist=stc_min_dist,
            max_dist=stc_max_dist,
            x_crd_layer_name="lon",
            y_crd_layer_name="lat",
            coordinate_type="geographic"
        )

        print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Filtering scatterers...")
        for layer in filter_dict.keys():
            stm_densified = stm_point_filter(
                stm=stm_densified,
                layer_to_filter=layer,
                vmin=filter_dict[layer][0],
                vmax=filter_dict[layer][1],
                return_removed=False
            )

        print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Adding viewing geometry...")
        stm_densified = add_local_viewing_geometry(
            stm=stm_densified,
            orbit_config_file=orbit_file,
            orbit_res=orbit_resolution,
            orbit_mode=orbit_mode,
            orbit=stack_id
        )

        print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Projecting phases and models onto the vertical...")
        stm_densified["unwrapped_phase_pov"] = stm_densified["unwrapped_phase"] / np.cos(
            np.radians(stm_densified["local_incidence_angle"])
        )
        for param in model_parameter_layer_names:
            stm_densified[param].data = stm_densified[param].values
            stm_densified[f"{param}_pov"] = stm_densified[param] / np.cos(np.radians(stm_densified["local_incidence_angle"]))

        print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Saving to zarr...")
        stm_densified = stm_densified.chunk({"time": 100, "space": "auto"})
        stm_densified.to_zarr(ps_filepath_4dens, mode="w", zarr_format=2)

    print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} PS densification for stack {stack_id} done.")
print("========== 4. PS densification done.")


# ========================================= #
#          5. PS-DS Arc Formation           #
# ========================================= #
print("========== 5. PS-DS arc formation ...")
for i, stack_id in enumerate(stack_ids):
    print(f"Processing stack {stack_id} ...")
    metadata = sarxarray.read_metadata(metadata_paths[i], driver="doris5")

    arc_filepath = os.path.join(stm_dir, arc_dtree_save_name.format(stack_id))
    if os.path.isdir(arc_filepath):
        print("PS-DS arc datatree file already exist. Skipping ...")
    else:
        print("Loading PS first order network file and DS datatree...")
        ps_filepath_3fon = os.path.join(stm_dir, ps_stm_save_name.format(stack_id, "3fon"))
        if os.path.isdir(ps_filepath_3fon):
            stm_pnts_output = xr.open_zarr(ps_filepath_3fon)
        else:
            raise FileNotFoundError(f"PS stm first order network file {ps_filepath_3fon} not found. Please run PS first order network formation first.")

        ds_filepath = os.path.join(stm_dir, ds_dtree_save_name.format(stack_id))
        if os.path.isdir(ds_filepath):
            ds_dtree = _open_datatree_compat(ds_filepath)
        else:
            raise FileNotFoundError(f"DS datatree file {ds_filepath} not found. Please run PS DS selection first.")
    
        print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} PS-DS arc formation...")
        arc_dtree = ps_ds_arc(
            stm_pnts_output, 
            ds_dtree, 
            n_connections=1, 
            key_xcoord=f"x_euclidean_proj_epsg{euclidean_epsg_code_number}", 
            key_ycoord=f"y_euclidean_proj_epsg{euclidean_epsg_code_number}",
            key_h2ph="h2ph",
            key_sd_phase_ps="sd_phase_minus_atmo",
            key_sd_phase_ds="ds_phi_aps_full",
        )

        arc_dtree.to_zarr(arc_filepath, mode="w")
        print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} PS-DS arc formation for {stack_id} done.")
print("========== 5. PS-DS arc formation done.")


# ========================================= #
#     6. Displacement Model Parameters      #
# ========================================= #
print("========== 6. Displacement Model Parameters ...")
# - TODO: Estimate displacement model parameters


if log_filename is not None:
    sys.stdout.flush()
