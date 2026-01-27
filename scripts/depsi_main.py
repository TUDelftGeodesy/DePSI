"""
This script runs the original PS DePSI workflow as originally implemented in Matlab.
"""
from datetime import datetime

import numpy as np
import sarxarray
import xarray

from depsi.arc_estimation import periodogram
from depsi.atmosphere_estimation import estimate_atmosphere_phase
from depsi.classification import network_stm_selection, ps_selection
from depsi.densification import densification
from depsi.io import (
    export_to_csv,
    export_to_skygeo_portal,
    export_convex_hull_to_shapefile,
    export_to_shapefile,
    read_slc_stack
)
from depsi.model_estimation import estimate_model_params
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


# ############## INPUT VARIABLES

slc_path = '/Users/sanvandiepen/PycharmProjects/workingEnvironment2/test_zarr/nl_amsterdam_s1_dsc_t037_base.zarr'

# Crop in space
aoi_file = '/Users/sanvandiepen/PycharmProjects/workingEnvironment2/test_zarr/shape/nl_amsterdam_shape.shp'

# metadata
metadata_path = '/Users/sanvandiepen/PycharmProjects/workingEnvironment2/test_zarr/nl_amsterdam_s1_dsc_t037.res'
satellite = "s1"
track = 37
direction = "dsc"

# Crop in time
first_date = datetime(2020, 1, 1)
last_date = datetime(2025, 9, 1)

mother_epoch = "auto"
reference_point_index = None

# PS selection method
ps_selection_method = "nmad"
threshold = 0.15
chunks_ps_selection = 1000
start_date_ps_selection = None
end_date_ps_selection = None

# Sidelobes
max_pixel_dist = 2
min_correlation = 0.90

# Network formation
network_crs = "radar"
network_x_crds = "azimuth"
network_y_crds = "range"
min_point_distance = 20  # pixels if radar, otherwise units of the EPSG code in network_crs
max_arc_length = 0.001  # degrees (lat/lon)
network_formation_method="redundant"
min_network_links = 16
network_partition_number = 8
network_dphase_method = "subtract"
min_periodogram_iterations = 10
arc_quality_threshold = 0.5  # ensemble coherence
model_types = ["offset", "velocity", "height"]


# atmosphere
euclidean_epsg_code_number = 28992  # 28992 is Dutch RD, for other AoIs visit https://epsg.io/ , units must be meters
atmo_unmodeled_displacement_filter_length = 1
atmo_unmodeled_displacement_sampling_rate = 1000
atmo_unmodeled_displacement_filter_type = "gaussian"

atmo_kriging_n_nearest_neighbours = None
atmo_kriging_backend = "vectorized"

atmo_empirical_variogram_method = "standard"
atmo_empirical_variogram_nlags = 50
atmo_empirical_variogram_cutoff_distance = 10000.0

atmo_variogram_model = "gaussian"
atmo_variogram_drift_terms = "regional_linear"

# Densification
n_densification_connections = 1

# Spatio-temporal consistency
stc_min_dist = 50  # meters
stc_max_dist = 200  # meters

# Viewing geometry
orbit_mode = "IWS"
orbit_resolution = 0.01
orbit_file = "/Users/sanvandiepen/PycharmProjects/workingEnvironment2/GitHub_repos/DePSI_group/config/drama/S1_XTI.cfg"

# output settings
filter_dict = {"h": [-2000, 2000]}
output_types = ["csv_web_portal", "csv", "shapefile", "convex_hull", "zarr"]
# csv, csv_web_portal, shapefile, convex_hull, zarr

# csv export
csv_save_path = "/Users/sanvandiepen/PycharmProjects/workingEnvironment2/test_zarr/nl_amsterdam_s1_dsc_t037.csv"
csv_ts_proj = "los"
csv_point_annotation_label = f"nl_amsterdam_{satellite}_{direction}_t{track:0>3d}"

# csv web portal export
csv_web_save_path = (
    "/Users/sanvandiepen/PycharmProjects/workingEnvironment2/test_zarr/nl_amsterdam_s1_dsc_t037_portal.csv"
)
csv_web_ts_proj = "los"
csv_web_point_annotation_label = f"nl_amsterdam_{satellite}_{direction}_t{track:0>3d}"

# shapefile export
shape_save_path = "/Users/sanvandiepen/PycharmProjects/workingEnvironment2/test_zarr/nl_amsterdam_s1_dsc_t037.shp"
shape_projection = "RD"  # RD or WGS84
shape_point_annotation_label = f"nl_amsterdam_{satellite}_{direction}_t{track:0>3d}"

# convex hull export
chull_save_path = "/Users/sanvandiepen/PycharmProjects/workingEnvironment2/test_zarr/nl_amsterdam_s1_dsc_t037_chull.shp"
chull_projection = "RD"  # RD or WGS84

# zarr export
stm_save_path = '/Users/sanvandiepen/PycharmProjects/workingEnvironment2/test_zarr/nl_amsterdam_s1_dsc_t037_result.zarr'

# 1. Project setup
print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Reading SLC stack...")
slcs = read_slc_stack(slc_path)


# 2. Scatterer selection
print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Start cropping...")
cropped_slcs = crop_slc_spacetime(
    slcs,
    aoi_filename=aoi_file,
    start_date=first_date,
    end_date=last_date
)

# ######## POINT SELECTION WITH THE PARAMETERS ABOVE ############
print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Starting PS selection...")
stm = ps_selection(
    cropped_slcs,
    method=ps_selection_method,
    threshold=threshold,
    ps_selection_start_date=start_date_ps_selection,
    ps_selection_end_date=end_date_ps_selection,
    output_chunks=chunks_ps_selection,
    mem_persist=False
)
print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Selected {len(stm.space)} PS.")

print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Computing single differences and temporal baseline...")
stm = stm_compute_single_time_differences(
    stm,
    single_difference_mother=mother_epoch
)
stm = add_stm_time_deltas(stm)
stm = stm.assign({"temporal_baseline": (
    ["time"],
    stm.years_since_first_img.values - stm.years_since_first_img.sel(time=stm.ps_sd_mother).values
)})

# sidelobe removal
print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Removing sidelobes...")
side_lobes_array, _ = detect_side_lobes(stm, max_pixel_dist, min_correlation, "sd_complex", "sd_amplitude_unnormalized")

mask_sidelobes = np.ones(len(stm.space), dtype=bool)
mask_sidelobes[side_lobes_array] = False

stm = stm.isel(space=mask_sidelobes)

print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Removed {np.sum(~mask_sidelobes)} sidelobes from the STM, "
      f"{len(stm.space)} remain.")

# 6a. Geocoding
print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Starting geocoding...")
metadata = sarxarray.read_metadata(metadata_path, driver="doris5")
stm = stm.assign_attrs({"wavelength": metadata["wavelength"]})

latlonh = radar_to_latlonh(
    azimuth_coords=stm.azimuth.values,
    range_coords=stm.range.values,
    elevation=stm.h.values,
    metadata=metadata,
)

stm["lat"].data = latlonh[0].flatten()
stm["lon"].data = latlonh[1].flatten()
stm["h"].data = latlonh[2].flatten()

# Add Euclidean coords in preparation for atmosphere estimation (here so it is also present in the first order network)
print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Adding Euclidean coordinates...")
crd_x, crd_y = convert_geographic_coords_to_euclidean(
    stm["lon"].values,
    stm["lat"].values,
    target_crs=f"EPSG:{euclidean_epsg_code_number}"
)
stm = stm.assign_coords({
    f"x_euclidean_proj_epsg{euclidean_epsg_code_number}": (["space"], crd_x),
    f"y_euclidean_proj_epsg{euclidean_epsg_code_number}": (["space"], crd_y),
})

# 3a. Network construction
print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Removing mother from network STM...")
mother_epoch_index = np.where(stm.time.values == stm.sel(time=stm.ps_sd_mother).time.values)[0][0]
non_mother = [True] * len(stm.time.values)
non_mother[mother_epoch_index] = False
stm_without_mother_epoch = stm.isel(time=non_mother)

print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Selecting first-order network points...")
stm_network_pnts = network_stm_selection(
    stm=stm_without_mother_epoch,
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

stm_network_arcs = stm_network_arcs.compute()

print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Integrating the spatial network...")
_, stm_pnts_output = spatial_integration(
    stm_network_pnts,
    stm_network_arcs,
    key_arc_quality="temp_coh",
    threshold_arc_quality=arc_quality_threshold,
    idx_refpnt=reference_point_index,
)

# 4. Atmosphere estimation
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
# Rename the space dimension since it will otherwise clash with the space dimension of the predicted coordinates
stm_pnts_output = stm_pnts_output.rename_dims({"space": "space_fon"})

print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Creating the atmosphere prediction coordinate dataset...")
atmosphere_prediction_coords = xarray.Dataset(
    coords={
        f"x_euclidean_proj_epsg{euclidean_epsg_code_number}": (
            "space", stm[f"x_euclidean_proj_epsg{euclidean_epsg_code_number}"].values
        ),
        f"y_euclidean_proj_epsg{euclidean_epsg_code_number}": (
            "space", stm[f"y_euclidean_proj_epsg{euclidean_epsg_code_number}"].values
        ),
    }
)

print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Estimating the atmospheric phase screens...")
# import pdb; pdb.set_trace()
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

# create the zero layer for the mother atmosphere to be appended to the interferometric atmospheric phase screens
mother_atmo_stm = xarray.Dataset(
    data_vars={
        "atmosphere_predicted": ("space", np.zeros((len(stm.space), ))),
        "atmosphere_sigmasq": ("space", np.zeros((len(stm.space), ))), },
    coords=stm["phase"].sel(time=stm.ps_sd_mother).coords
)

atmosphere = xarray.concat([atmospheric_phase_screens, mother_atmo_stm], dim="time").sortby("time")

# Add the atmosphere into the original STM
stm["atmosphere_predicted"] = atmosphere["atmosphere_predicted"].transpose("space", "time")
stm["atmosphere_sigmasq"] = atmosphere["atmosphere_sigmasq"].transpose("space", "time")

print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Removing the atmospheric phase screens...")
stm["phase_minus_atmo"] = (stm["phase"] - stm["atmosphere_predicted"] + np.pi) % (2 * np.pi) - np.pi
stm["sd_phase_minus_atmo"] = \
    (stm["phase_minus_atmo"] - stm["phase_minus_atmo"].sel(time=stm.ps_sd_mother) + np.pi) % (2 * np.pi) - np.pi

# 3b. Network construction
stm_atmo_corr_without_mother_epoch = stm.isel(time=non_mother)
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

stm_network_arcs = stm_network_arcs.compute()

print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Integrating the spatial network...")
_, stm_firstordernetwork = spatial_integration(
    stm_network_pnts,
    stm_network_arcs,
    key_arc_quality="temp_coh",
    threshold_arc_quality=arc_quality_threshold,
    idx_refpnt=reference_point_index,
)

# 5. Densification
print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Starting densification...")
stm_densified = densification(
    stm_atmo_corr_without_mother_epoch,
    stm_firstordernetwork,
    idx_refpnt=reference_point_index,
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

# 6b. Geocoding
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

# 7. STC Calculation
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

# 8. Scatterer filtering
print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Filtering scatterers...")
for layer in filter_dict.keys():
    stm_densified = stm_point_filter(
        stm=stm_densified,
        layer_to_filter=layer,
        vmin=filter_dict[layer][0],
        vmax=filter_dict[layer][1],
        return_removed=False
    )

# 9. Output generation
print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Adding viewing geometry...")
stm_densified = add_local_viewing_geometry(
    stm=stm_densified,
    orbit_config_file=orbit_file,
    orbit_res=orbit_resolution,
    orbit_mode=orbit_mode,
    orbit=f"{satellite}_{direction}_t{track:0>3d}"
)

print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Projecting phases and models onto the vertical...")
stm_densified["unwrapped_phase_pov"] = stm_densified["unwrapped_phase"] / np.cos(
    np.radians(stm_densified["local_incidence_angle"])
)
for param in model_parameter_layer_names:
    stm_densified[param].data = stm_densified[param].values
    stm_densified[f"{param}_pov"] = stm_densified[param] / np.cos(np.radians(stm_densified["local_incidence_angle"]))

print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Starting export...")
if "csv" in output_types:
    print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Exporting to CSV...")
    export_to_csv(
        stm=stm_densified,
        save_path=csv_save_path,
        model_parameter_layer_names=model_parameter_layer_names,
        ts_proj=csv_ts_proj,
        point_annotation_label=csv_point_annotation_label,
    )

if "csv_web_portal" in output_types:
    print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Exporting to Portal CSV...")
    export_to_skygeo_portal(
        stm=stm_densified,
        save_path=csv_web_save_path,
        ts_proj=csv_web_ts_proj,
        point_annotation_label=csv_web_point_annotation_label,
        satellite=satellite,
        asc_dsc=direction,
        azimuth_spacing=metadata["azimuth_pixel_spacing"],
        range_spacing=metadata["range_pixel_spacing"],
        model_names=model_types,
        model_parameter_layer_names=model_parameter_layer_names,
    )

if "shapefile" in output_types:
    print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Exporting to shapefile...")
    export_to_shapefile(
        stm=stm_densified,
        save_path=shape_save_path,
        projection=shape_projection,
        model_parameter_layer_names=model_parameter_layer_names,
        point_annotation_label=shape_point_annotation_label,
    )

if "convex_hull" in output_types:
    print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Exporting to convex hull...")
    export_convex_hull_to_shapefile(
        stm=stm_densified,
        save_path=chull_save_path,
        projection=chull_projection,
    )

if "zarr" in output_types:
    print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Exporting to zarr...")
    stm_densified.to_zarr(stm_save_path, mode="w")

print(f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} Finished!")
