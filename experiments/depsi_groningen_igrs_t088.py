"""
Run a PS + optional regional-PS + designated-target DePSI workflow.

The PS workflow is the main workflow. DT processing is optional:
- if a DT Radar Coding CSV is available, DT points are selected from the SLC
  stack and used as additional atmosphere-prediction locations;
- if no DT CSV is available, the script runs as a normal PS-only workflow.

When DT is available, PS and DT are forced onto the same time axis before any
selection is done. This prevents mismatches when the SLC stack has been updated
with newer acquisitions but the DT Radar Coding CSV has not.

Regional PS processing is also optional. It applies a separate, potentially
looser PS threshold inside either a point buffer or a supplied polygon.
The selected focus PS is reserved for regional densification. Other regional
points may enter the global first-order network when they also satisfy the
global PS threshold; those overlaps are removed from the regional product
after the final network has been formed. Global PS, regional PS, and DT results
remain separate through filtering and export.
"""

# AML changes applied
import math
# AML changes applied end
import csv
import inspect
import json
import os
import warnings
from datetime import datetime
from pathlib import Path

import geopandas as gpd
# AML changes applied
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
# AML changes applied end
# AML changes applied
import networkx as nx
# AML changes applied end
# <>
import drama.mission.timeline.no_drama as nd
from drama.io import cfg
from drama.performance.sar import SARModeFromCfg
# <|>
import numpy as np
import drama.geo.geometry as drama_geometry
import pandas as pd
import sarxarray
import xarray as xr
# <>
from shapely.geometry import Point, Polygon, box
from shapely.ops import unary_union
# <|>

from depsi.arc_estimation import periodogram
from depsi.atmosphere_estimation import estimate_atmosphere_phase
from depsi.classification import (
    designated_target_selection,
    network_stm_selection,
    ps_selection,
)
# <>
from depsi.io import (
    export_convex_hull_to_shapefile,
    export_to_csv as _export_to_csv,
    export_to_shapefile,
    export_to_skygeo_portal as _export_to_skygeo_portal,
    read_rcs_csv,
    read_slc_stack,
)
# <|>
from depsi.model_estimation import estimate_model_params
import depsi.network as depsi_network
from depsi.network import form_network, spatial_integration
from depsi.point_quality import detect_side_lobes
from depsi.transformations import radar_to_latlonh
from depsi.densification import _determine_dens_connections
from depsi.model_definition import MODEL_NAMES_PARAMS, construct_design_matrix

from depsi.utils import (
    add_stm_time_deltas,
    compute_phase_difference,
    convert_geographic_coords_to_euclidean,
    crop_slc_spacetime,
    generate_pnt_uids,
    stm_compute_single_time_differences,
)

from depsi.point_quality import compute_spatiotemporal_consistency
from depsi.postprocessing import stm_point_filter
from depsi.viewing_geometry import add_local_viewing_geometry
# <>
from depsi.utils import identify_s1_orbits_in_aoi
# <|>

# ############## INPUT VARIABLES

# AML changes applied
# Start a fresh run because the spatial crop and Step 4 densification scope
# differ from the existing small-area checkpoints.
start_at_checkpoint = 0
# 0 = run full workflow
# 1 = start from PS/DT selection checkpoints
# 2 = start from atmosphere-corrected checkpoints
# 3 = start from integrated first-order PS network
# 4 = start from densified/model-estimated PS/DT checkpoints
# 5 = start from geocoded/viewing-geometry result checkpoints
# 6 = export only

proj_dir = "/project/caroline/Share/users/caroline-alapadat/software/py_depsi"
run_name = "test_dtd_groningen_60x60km_asct088_v2"

stack_id = "s1_asc_t088"
satellite = "s1"
track = 88
direction = "asc"
wavelength = 0.055465763

# DT diagnostic quantity: "phase" for radians or "displacement" for LOS mm.
dt_plot_quantity = "displacement"

# - Original full SLC zarr
#   This is the full stack before the notebook/script crop is applied.
slc_path = Path(
    "/project/caroline/Share/stacks_zarr/nl_groningen_s1_asc_t088/nl_groningen_s1_asc_t088.zarr"
)

# - Metadata
#   Used for PS geocoding and radar-coordinate network spacing.
metadata_path = (
    "/project/caroline/Share/stacks/nl_groningen_s1_asc_t088/stack/20200325/master.res"
)

# - Crop in time
#   If last_date is None, the crop runs until the last available SLC epoch.
first_date = datetime(2018, 1, 1)
last_date = None

# - Crop in space
#   Bounding box is given in lon/lat coordinates and covers the approximately
#   60 x 60 km Groningen processing region including all IGRS.
lon_min = 6.287104
lon_max = 7.252834
lat_min = 52.982511
lat_max = 53.545538
# AML changes applied end

# - Optional DT Radar Coding CSV
#   Coordinates must be in the same radar-coordinate convention as the SLC stack.
#   Set to None, or to a non-existing path, to run PS-only.
dt_rcs_csv_path = (
    "/project/caroline/Share/users/caroline-alapadat/software/radarcoding_toolbox/"
    "DelftRadarCodingToolbox/results/nl_groningen_s1_asc_t088_workflow_ovs6/"
    "RadarCoordinates/s1_asc088_RC.csv"
)

# - Reference epoch
#   Derive the mother epoch from the dated directory containing master.res so
#   the processing date cannot drift from the selected metadata file.
mother_epoch = datetime.strptime(Path(metadata_path).parent.name, "%Y%m%d")
reference_point_index = None

# - PS selection
ps_selection_method = "nad"
threshold = 0.5
chunks_ps_selection = 5000
start_date_ps_selection = None
end_date_ps_selection = None

# - Sidelobe removal
do_sidelobe_detection = False
max_pixel_dist = 2
min_correlation = 0.90

# - Network formation
#   The first-order PS network is used to estimate the atmospheric phase screen.
network_crs = "radar"
network_x_crds = "azimuth"
network_y_crds = "range"
# AML changes applied
# Use a 1 km spacing to keep the first-order network computationally manageable
# over the approximately 60 x 60 km AOI. Arc length remains in degrees because
# form_network currently uses longitude and latitude coordinates.
min_point_distance = 1000  # metres
max_arc_length = 0.06  # degrees (longitude/latitude)
network_formation_method = "redundant"
min_network_links = 16
network_partition_number = 8
network_dphase_method = "subtract"
min_periodogram_iterations = 10
# <>
arc_quality_threshold = 0.6
# <|>
min_arc_connections = 3
# AML changes accepted
# Keep the largest connected component automatically when it contains at least
# <>
# 50% of the points remaining after quality and minimum-connection filtering.
largest_component_ratio = 0.5
skip_MHT_spatial_integration = True
# <|>
# AML changes accepted end
# AML changes applied end
model_types = ["offset", "velocity", "height"]

# - Atmosphere estimation
#   Euclidean coordinates are used by kriging, so lon/lat is projected first.
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

# Densification
n_densification_connections = 1

# AML changes applied
# - Optional regional PS selection and densification
#   None        = do not select or densify additional regional PS points
#   "buffer"    = create a circular region around the configured focus point
#   "shapefile" = use the polygon supplied by regional_ps_region_file
#
# The focus point is configured in longitude/latitude for convenience and is
# transformed to EPSG:28992 before metric distance calculations.
regional_ps_mode = "buffer"
# regional_ps_mode = "buffer"
# regional_ps_mode = "shapefile"

regional_focus_lon = 6.933417019109445 #None
regional_focus_lat = 53.326531768919104 #None
regional_focus_point_crs = "EPSG:4326"
regional_focus_search_radius_m = 50.0

# Used only in buffer mode.
regional_buffer_radius_m = 500.0

# Used only in shapefile mode. EPSG:28992 is recommended, although any file
# with a valid CRS is accepted and reprojected internally.
regional_ps_region_file = None
# regional_ps_region_file = "/path/to/tide_gauge_region.shp"

# Regional PS selection may be looser than the full-area PS selection.
# For NAD, a larger threshold is looser because selection uses NAD < threshold.
regional_ps_selection_method = "nad"
regional_ps_threshold = 0.7
# AML changes applied end

# AML changes accepted
# Spatio-temporal consistency and filtering can be configured independently
# for the global PS and regional PS products. DT points are not STC-filtered.
stc_min_dist = 30  # meters
stc_max_dist = 100  # meters
ps_compute_stc = True
ps_filter_dict = {"h": [-2000, 2000]}
regional_ps_compute_stc = True
regional_ps_filter_dict = {}
# AML changes accepted end

# Viewing geometry
orbit_mode = "IWS"
orbit_resolution = 0.01
orbit_file = "/project/caroline/Share/users/caroline-alapadat/software/py_depsi/DePSI/config/drama/S1_XTI.cfg"

# - Output/checkpoint paths
#   Checkpoints allow restarting from selected points or atmosphere-corrected STMs.
run_dir = os.path.join(proj_dir, run_name)
crop_dir = os.path.join(run_dir, "input_crops")
stm_dir = os.path.join(run_dir, "stm")

aoi_file = os.path.join(crop_dir, "groningen_60x60km_bbox.shp")
# cropped_slc_path = os.path.join(
#     crop_dir,
#     "nl_groningen_s1_asc_t088_60x60km_20180101.zarr",
# )
cropped_slc_path = os.path.join(
    proj_dir,
    "test_dtd_groningen_60x60km_asct088",  # original run folder
    "input_crops",
    "nl_groningen_s1_asc_t088_60x60km_20180101.zarr",
)

# AML changes applied
# Normalized regional geometries written by the script. Metric operations use
# RD New, while crop_slc_spacetime receives the WGS84 copy.
regional_ps_region_rd_path = os.path.join(
    crop_dir,
    "regional_ps_region_epsg28992.shp",
)
regional_ps_region_wgs84_path = os.path.join(
    crop_dir,
    "regional_ps_region_epsg4326.shp",
)
# AML changes applied end

# AML changes applied
ps_stm_1sel_path = os.path.join(stm_dir, f"ps_stm_{stack_id}_1sel.zarr")
dt_stm_1sel_path = os.path.join(stm_dir, f"dt_stm_{stack_id}_1sel.zarr")
regional_ps_stm_1sel_path = os.path.join(
    stm_dir,
    f"regional_ps_stm_{stack_id}_1sel.zarr",
)
ps_stm_4atmo_path = os.path.join(stm_dir, f"ps_stm_{stack_id}_4atmo.zarr")
dt_stm_4atmo_path = os.path.join(stm_dir, f"dt_stm_{stack_id}_4atmo.zarr")
regional_ps_stm_4atmo_path = os.path.join(
    stm_dir,
    f"regional_ps_stm_{stack_id}_4atmo.zarr",
)
# AML changes applied end

# <>
# APS periodogram checkpoints store the selected network points and the
# computed arc ambiguities/coherences. Keeping both allows spatial integration
# to be retried without repeating network formation or the periodogram.
aps_periodogram_points_path = os.path.join(
    stm_dir,
    f"ps_stm_{stack_id}_2aps_periodogram_points.zarr",
)
aps_periodogram_arcs_path = os.path.join(
    stm_dir,
    f"ps_stm_{stack_id}_2aps_periodogram_arcs.zarr",
)
# <|>

ps_stm_5fnet_path = os.path.join(stm_dir, f"ps_stm_{stack_id}_5fnet.zarr")
ps_stm_5dens_path = os.path.join(stm_dir, f"ps_stm_{stack_id}_5dens.zarr")
# AML changes accepted
# Regional PS remain separate after densification so their quality filtering
# and export can be controlled independently from the global PS product.
regional_ps_stm_5dens_path = os.path.join(
    stm_dir,
    f"regional_ps_stm_{stack_id}_5dens.zarr",
)
# AML changes accepted end
dt_stm_5dens_path = os.path.join(stm_dir, f"dt_stm_{stack_id}_5dens.zarr")

ps_stm_result_path = os.path.join(stm_dir, f"ps_stm_{stack_id}_result.zarr")
# AML changes accepted
regional_ps_stm_result_path = os.path.join(
    stm_dir,
    f"regional_ps_stm_{stack_id}_result.zarr",
)
# AML changes accepted end
dt_stm_result_path = os.path.join(stm_dir, f"dt_stm_{stack_id}_result.zarr")

os.makedirs(crop_dir, exist_ok=True)
os.makedirs(stm_dir, exist_ok=True)

# - Output settings
output_types = ["csv_web_portal", "csv", "shapefile", "convex_hull"]
# csv, csv_web_portal, shapefile, convex_hull, zarr

output_dir = os.path.join(run_dir, "outputs")
os.makedirs(output_dir, exist_ok=True)

output_label = f"groningen_{satellite}_{direction}_t{track:0>3d}"

# - PS CSV export
ps_csv_save_path = os.path.join(output_dir,f"ps_{stack_id}.csv")
csv_ts_proj = "los"
csv_point_annotation_label = output_label

# - PS CSV web portal export
ps_csv_web_save_path = os.path.join(output_dir,f"ps_{stack_id}_portal.csv")
csv_web_ts_proj = "los"
csv_web_point_annotation_label = output_label

# - PS shapefile export
ps_shape_save_path = os.path.join(output_dir,f"ps_{stack_id}.shp")
shape_projection = "RD"  # RD or WGS84
shape_point_annotation_label = output_label

# - PS convex hull export
ps_chull_save_path = os.path.join(output_dir,f"ps_{stack_id}_chull.shp")
chull_projection = "RD"  # RD or WGS84

# AML changes accepted
# - Regional PS exports use their own files and annotation label.
regional_ps_output_label = (
    f"groningen_regional_{satellite}_{direction}_t{track:0>3d}"
)
regional_ps_csv_save_path = os.path.join(
    output_dir,
    f"regional_ps_{stack_id}.csv",
)
regional_ps_csv_web_save_path = os.path.join(
    output_dir,
    f"regional_ps_{stack_id}_portal.csv",
)
regional_ps_shape_save_path = os.path.join(
    output_dir,
    f"regional_ps_{stack_id}.shp",
)
regional_ps_chull_save_path = os.path.join(
    output_dir,
    f"regional_ps_{stack_id}_chull.shp",
)
# AML changes accepted end

# - DT CSV export
dt_csv_save_path = os.path.join(output_dir,f"dt_{stack_id}.csv")

# - DT CSV web portal export
dt_csv_web_save_path = os.path.join(output_dir,f"dt_{stack_id}_portal.csv")

# - DT shapefile export
dt_shape_save_path = os.path.join(output_dir,f"dt_{stack_id}.shp")

# - DT convex hull export
dt_chull_save_path = os.path.join(output_dir,f"dt_{stack_id}_chull.shp")

# AML changes applied
# - DT diagnostic time-series figure
dt_timeseries_plot_path = os.path.join(
    output_dir,
    f"dt_{stack_id}_{dt_plot_quantity}_amplitude.png",
)
# AML changes applied end

# =============================================================================
# Local workflow functions
# =============================================================================


# -----------------------------------------------------------------------------
# General workflow utilities
# -----------------------------------------------------------------------------


def timestamp():
    """Return a timestamp for progress messages.

    Returns
    -------
    str
        Local date and time formatted as ``YYYY-MM-DDTHH:MM:SS``.

    Notes
    -----
    This is a script utility and does not replace a DePSI function.
    """
    return datetime.now().strftime("%Y-%m-%dT%H:%M:%S")


def dt_csv_exists():
    """Check whether optional DT processing has a readable Radar Coding CSV.

    Returns
    -------
    bool
        ``True`` when a DT CSV path is configured and exists.

    Notes
    -----
    This is a script utility and does not replace a DePSI function.
    """
    return dt_rcs_csv_path is not None and os.path.isfile(dt_rcs_csv_path)


def format_epochs(times):
    """Format acquisition epochs for readable mismatch diagnostics.

    Parameters
    ----------
    times : array-like
        Datetime-like acquisition epochs.

    Returns
    -------
    list[str]
        Epochs formatted as ``YYYY-MM-DD``.

    Notes
    -----
    This is a script utility and does not replace a DePSI function.
    """
    return [pd.Timestamp(t).strftime("%Y-%m-%d") for t in times]


def validate_mother_epoch(times, label):
    """Verify that the configured mother epoch exists on a time axis.

    Parameters
    ----------
    times : array-like
        Processing epochs to validate.
    label : str
        Human-readable dataset name used in the error message.

    Raises
    ------
    ValueError
        If the configured mother epoch is absent.

    Notes
    -----
    This precondition protects the DePSI single-difference operations after
    SLC and DT time-axis alignment; it does not replace an upstream function.
    """
    epochs = pd.DatetimeIndex(times)
    mother = pd.Timestamp(mother_epoch)
    if mother not in epochs:
        raise ValueError(
            f"Mother epoch {mother:%Y-%m-%d} is absent from the {label} "
            "time axis. Check metadata_path and the PS/DT common epochs."
        )



# -----------------------------------------------------------------------------
# Storage and serialization
# -----------------------------------------------------------------------------


def save_zarr_v2(dataset, path, spatial_chunk=4000):
    """Write a Zarr v2 checkpoint using compatible, aligned chunks.

    Parameters
    ----------
    dataset : xarray.Dataset
        Dataset to save.
    path : str or path-like
        Destination Zarr store.
    spatial_chunk : int, default 4000
        Maximum chunk length for azimuth and range dimensions.

    Returns
    -------
    None

    Notes
    -----
    This helper replaces direct ``Dataset.to_zarr`` calls in this script. It
    resolves Zarr v2/v3 codec errors and Dask/Zarr chunk-alignment failures by
    rechunking rasters, selecting the installed xarray version's format
    keyword, and enabling chunk alignment when supported.
    """
    if "azimuth" in dataset.dims and "range" in dataset.dims:
        chunks = {
            "azimuth": min(spatial_chunk, dataset.sizes["azimuth"]),
            "range": min(spatial_chunk, dataset.sizes["range"]),
        }
        if "time" in dataset.dims:
            chunks["time"] = 1
        dataset = dataset.chunk(chunks)

    to_zarr_parameters = inspect.signature(dataset.to_zarr).parameters

    write_options = {
        "mode": "w",
    }

    if "zarr_format" in to_zarr_parameters:
        write_options["zarr_format"] = 2
    else:
        write_options["zarr_version"] = 2

    if "align_chunks" in to_zarr_parameters:
        write_options["align_chunks"] = True

    dataset.to_zarr(path, **write_options)


def make_dt_stm_zarr_safe(dt_stm):
    """Convert DT string metadata to a Zarr-safe fixed-width Unicode dtype.

    Parameters
    ----------
    dt_stm : xarray.Dataset or None
        DT space-time matrix to prepare.

    Returns
    -------
    xarray.Dataset or None
        Dataset with a fixed-width Unicode ``target`` field.

    Notes
    -----
    This prepares input for ``Dataset.to_zarr`` rather than replacing it. It
    resolves the ``expected unicode string, found 1`` VLenUTF8 error caused by
    object-dtype target identifiers.
    """
    if dt_stm is None:
        return None

    if "target" in dt_stm:
        target_values = dt_stm["target"].values.astype(str)

        if target_values.size == 0:
            max_target_len = 1
        else:
            max_target_len = max(1, max(len(value) for value in target_values))

        dt_stm["target"] = xr.DataArray(
            target_values.astype(f"<U{max_target_len}"),
            dims=dt_stm["target"].dims,
            coords=dt_stm["target"].coords,
            attrs=dt_stm["target"].attrs,
        )

    return dt_stm



# -----------------------------------------------------------------------------
# STM preparation and overlap handling
# -----------------------------------------------------------------------------


def add_single_differences_and_time(stm, mother_epoch):
    """Add single differences, elapsed time, and temporal baselines to an STM.

    Parameters
    ----------
    stm : xarray.Dataset
        Input space-time matrix.
    mother_epoch : datetime-like
        Epoch used as the single-difference reference.

    Returns
    -------
    xarray.Dataset
        STM containing single-difference and temporal-baseline variables.

    Notes
    -----
    This convenience wrapper combines
    ``stm_compute_single_time_differences`` and ``add_stm_time_deltas``; it
    does not change their numerical behavior.
    """
    stm = stm_compute_single_time_differences(
        stm,
        single_difference_mother=mother_epoch,
    )
    stm = add_stm_time_deltas(stm)
    stm = stm.assign({
        "temporal_baseline": (
            ["time"],
            stm.years_since_first_img.values
            - stm.years_since_first_img.sel(time=stm.ps_sd_mother).values,
        )
    })
    return stm


def add_euclidean_coordinates(stm):
    """Add projected coordinates used by atmospheric kriging.

    Parameters
    ----------
    stm : xarray.Dataset
        STM containing longitude and latitude.

    Returns
    -------
    xarray.Dataset
        STM with projected x/y coordinates.

    Notes
    -----
    This wraps ``convert_geographic_coords_to_euclidean`` using the configured
    EPSG code so PS and DT atmosphere predictions use the same metric system.
    """
    x, y = convert_geographic_coords_to_euclidean(
        stm["lon"].values,
        stm["lat"].values,
        target_crs=f"EPSG:{euclidean_epsg_code_number}",
    )

    return stm.assign_coords({
        f"x_euclidean_proj_epsg{euclidean_epsg_code_number}": (["space"], x),
        f"y_euclidean_proj_epsg{euclidean_epsg_code_number}": (["space"], y),
    })


def remove_radar_coordinate_overlaps(stm, exclusion_stm):
    """Remove points that share an integer radar sample with exclusion points.

    Parameters
    ----------
    stm : xarray.Dataset
        PS candidates from which overlaps are removed.
    exclusion_stm : xarray.Dataset or None
        DT or other STM whose azimuth/range pairs take precedence.

    Returns
    -------
    tuple[xarray.Dataset, int]
        Filtered STM and number of removed points.

    Notes
    -----
    This is new PS+DT workflow logic and has no ``depsi_main.py`` equivalent.
    It prevents one SLC sample from appearing both as a PS and as a DT, which
    would duplicate its measurement and potentially give it double influence.
    """
    if exclusion_stm is None or exclusion_stm.sizes.get("space", 0) == 0:
        return stm, 0

    stm_coords = pd.MultiIndex.from_arrays([
        stm["azimuth"].values,
        stm["range"].values,
    ])
    exclusion_coords = pd.MultiIndex.from_arrays([
        exclusion_stm["azimuth"].values,
        exclusion_stm["range"].values,
    ])
    overlap = stm_coords.isin(exclusion_coords)
    return stm.isel(space=~overlap), int(overlap.sum())



# -----------------------------------------------------------------------------
# Optional regional PS selection
# -----------------------------------------------------------------------------


def regional_ps_enabled():
    """Check whether optional regional PS processing is configured.

    Returns
    -------
    bool
        ``True`` for buffer or shapefile regional modes.

    Notes
    -----
    This is a custom workflow helper and does not replace a DePSI function.
    """
    return regional_ps_mode in {"buffer", "shapefile"}


def build_regional_ps_geometry():
    """Build the regional selection geometry and projected focus point (PS of interest).

    Returns
    -------
    tuple[geopandas.GeoDataFrame or None, shapely.geometry.Point or None]
        Regional geometry and focus point in the configured metric CRS.

    Raises
    ------
    ValueError
        If required coordinates, mode, CRS, or shapefile input are invalid.

    Notes
    -----
    This implements optional buffer/shapefile selection absent from the
    original ``depsi_main.py`` workflow.
    """
    if not regional_ps_enabled():
        return None, None

    if regional_focus_lon is None or regional_focus_lat is None:
        raise ValueError(
            "regional_focus_lon and regional_focus_lat must be configured "
            "when regional PS processing is enabled."
        )

    metric_crs = f"EPSG:{euclidean_epsg_code_number}"
    focus = gpd.GeoDataFrame(
        {"id": [1]},
        geometry=[Point(regional_focus_lon, regional_focus_lat)],
        crs=regional_focus_point_crs,
    ).to_crs(metric_crs)

    if regional_ps_mode == "buffer":
        region_rd = gpd.GeoDataFrame(
            {"id": [1]},
            geometry=[focus.geometry.iloc[0].buffer(regional_buffer_radius_m)],
            crs=metric_crs,
        )

    else:
        if regional_ps_region_file is None:
            raise ValueError(
                "regional_ps_region_file must be configured in shapefile mode."
            )

        region_path = Path(regional_ps_region_file)
        if not region_path.is_file():
            raise FileNotFoundError(
                f"Regional PS polygon does not exist: {region_path}"
            )

        region_input = gpd.read_file(region_path)
        if region_input.empty:
            raise ValueError(
                f"Regional PS polygon contains no geometries: {region_path}"
            )
        if region_input.crs is None:
            raise ValueError(
                f"Regional PS polygon has no CRS: {region_path}"
            )

        region_rd = region_input.to_crs(metric_crs)

    region_geometry = region_rd.geometry.unary_union
    focus_geometry = focus.geometry.iloc[0]

    if not region_geometry.covers(focus_geometry):
        raise ValueError(
            "The configured regional focus point is outside the regional PS polygon."
        )

    # Save one normalized metric copy for inspection and a WGS84 copy for the
    # SLC crop helper used below.
    region_rd_normalized = gpd.GeoDataFrame(
        {"id": [1]},
        geometry=[region_geometry],
        crs=metric_crs,
    )
    region_rd_normalized.to_file(regional_ps_region_rd_path)
    region_rd_normalized.to_crs("EPSG:4326").to_file(
        regional_ps_region_wgs84_path
    )

    return region_rd_normalized, focus


def select_ps_inside_region(stm, region_rd):
    """Select STM points located inside the regional PS geometry.

    Parameters
    ----------
    stm : xarray.Dataset
        Candidate PS STM with projected coordinates.
    region_rd : geopandas.GeoDataFrame
        Regional selection geometry in the metric workflow CRS.

    Returns
    -------
    xarray.Dataset
        Spatially selected PS candidates.

    Notes
    -----
    This is custom regional-PS logic and does not replace a DePSI function.
    """
    x_label = f"x_euclidean_proj_epsg{euclidean_epsg_code_number}"
    y_label = f"y_euclidean_proj_epsg{euclidean_epsg_code_number}"
    ps_points = gpd.GeoSeries(
        gpd.points_from_xy(stm[x_label].values, stm[y_label].values),
        crs=region_rd.crs,
    )
    region_geometry = region_rd.geometry.unary_union
    selected = ps_points.within(region_geometry) | ps_points.touches(
        region_geometry
    )

    # <>
    # Materialize the one-dimensional polygon mask before applying isel().
    # This prevents a lazy boolean selection from propagating an unknown
    # Dask space length (represented as NaN) into regional densification,
    # where xarray would otherwise raise:
    #   ValueError: conflicting sizes for dimension 'space': length nan
    selected_values = np.asarray(
        selected.to_numpy(),
        dtype=bool,
    )
    selected_indices = np.flatnonzero(selected_values)
    selected_stm = stm.isel(space=selected_indices)
    # <|>

    return selected_stm


# <>
def validate_known_space_chunks(stm, label):
    """Validate concrete Dask dimensions before densification.

    This catches the unknown-size Dask state that previously caused xarray
    ``isel`` to fail with ``conflicting sizes for dimension 'space': length
    nan`` during regional PS densification.
    """
    if "space" not in stm.dims:
        raise ValueError(
            f"{label} does not contain a space dimension."
        )

    for name, variable in stm.variables.items():
        chunks = getattr(variable.data, "chunks", None)

        if chunks is None:
            continue

        for axis_chunks in chunks:
            if any(
                chunk_size is None
                or not np.isfinite(chunk_size)
                for chunk_size in axis_chunks
            ):
                raise ValueError(
                    f"{label} variable '{name}' has unknown Dask "
                    "chunk sizes. Save and reopen the STM before "
                    "densification."
                )

    print(
        f"{timestamp()} {label}: confirmed known space dimension "
        f"with {stm.sizes['space']} point(s)."
    )


def open_densification_checkpoint(path, label):
    """Open a Zarr checkpoint with concrete dimensions for Step 4.

    Reopening replaces inherited lazy selection graphs with dimensions read
    from Zarr metadata. This bypasses the xarray ``length nan`` failure seen
    when regional PS products continued directly from an in-memory Dask graph.
    """
    if not os.path.isdir(path):
        raise FileNotFoundError(
            f"{label} checkpoint does not exist: {path}"
        )

    stm = xr.open_zarr(path)
    validate_known_space_chunks(stm, label)

    return stm
# <|>


def mark_regional_focus_ps(regional_stm, focus_rd):
    """Mark the regional PS nearest to the configured focus location.

    Parameters
    ----------
    regional_stm : xarray.Dataset
        Regional PS candidates.
    focus_rd : shapely.geometry.Point
        Focus location in the metric workflow CRS.

    Returns
    -------
    xarray.Dataset
        STM containing the boolean ``is_focus_ps`` marker.

    Raises
    ------
    ValueError
        If no candidate lies within the configured search radius.

    Notes
    -----
    This custom helper reserves the focus PS for regional densification; the
    original ``depsi_main.py`` workflow has no equivalent concept.
    """
    x_label = f"x_euclidean_proj_epsg{euclidean_epsg_code_number}"
    y_label = f"y_euclidean_proj_epsg{euclidean_epsg_code_number}"
    focus_x = focus_rd.geometry.iloc[0].x
    focus_y = focus_rd.geometry.iloc[0].y

    distances = np.hypot(
        regional_stm[x_label].values - focus_x,
        regional_stm[y_label].values - focus_y,
    )
    nearby = np.flatnonzero(distances <= regional_focus_search_radius_m)
    is_focus_ps = np.zeros(regional_stm.sizes["space"], dtype=bool)

    if nearby.size == 0:
        message = (
            "No non-DT regional PS satisfying "
            f"{regional_ps_selection_method.upper()} < {regional_ps_threshold} "
            f"was found within {regional_focus_search_radius_m:.1f} m of the "
            "configured focus point."
        )
        warnings.warn(message, stacklevel=2)
        print(f"{timestamp()} WARNING: {message}")
    else:
        quality_key = f"time_selection_{regional_ps_selection_method}"
        nearby_quality = regional_stm[quality_key].isel(space=nearby).values
        focus_index = nearby[np.nanargmin(nearby_quality)]
        is_focus_ps[focus_index] = True

        focus_quality = float(
            np.asarray(
                regional_stm[quality_key].isel(space=focus_index).values
            ).item()
        )
        focus_azimuth = regional_stm["azimuth"].isel(
            space=focus_index
        ).values.item()
        focus_range = regional_stm["range"].isel(
            space=focus_index
        ).values.item()

        print(
            f"{timestamp()} Focus PS selected at "
            f"{distances[focus_index]:.2f} m from the configured point; "
            f"{regional_ps_selection_method.upper()}="
            f"{focus_quality:.4f}, radar sample=({focus_azimuth}, "
            f"{focus_range})."
        )

    regional_stm["is_focus_ps"] = (("space",), is_focus_ps)
    return regional_stm



# -----------------------------------------------------------------------------
# Network diagnostics and compatibility patch
# -----------------------------------------------------------------------------


def report_network_connectivity(
    stm_pnts,
    stm_arcs,
    label,
    quality_key=None,
    quality_threshold=None,
):
    """Report connected components before spatial integration.

    Parameters
    ----------
    stm_pnts : xarray.Dataset
        Network points.
    stm_arcs : xarray.Dataset
        Network arcs containing source and target indices.
    label : str
        Description printed with the report.
    quality_key : str, optional
        Arc-quality variable used for optional filtering.
    quality_threshold : float, optional
        Minimum absolute arc quality retained in the report.

    Returns
    -------
    list[int]
        Connected-component sizes in descending order.

    Notes
    -----
    This diagnostic supplements ``spatial_integration`` rather than replacing
    it. It exposes disconnected subnetworks before DePSI raises a less
    informative largest-component error.
    """
    source = np.asarray(stm_arcs["source"].values)
    target = np.asarray(stm_arcs["target"].values)

    if quality_key is not None:
        quality = np.asarray(stm_arcs[quality_key].values)
        keep = np.abs(quality) > quality_threshold
        source = source[keep]
        target = target[keep]

    graph = nx.Graph()
    graph.add_nodes_from(range(stm_pnts.sizes["space"]))
    graph.add_edges_from(zip(source, target, strict=True))

    component_sizes = sorted(
        (len(component) for component in nx.connected_components(graph)),
        reverse=True,
    )
    n_points = stm_pnts.sizes["space"]
    largest_fraction = component_sizes[0] / n_points if n_points else 0.0
    component_preview = component_sizes[:10]
    preview_suffix = " ..." if len(component_sizes) > 10 else ""

    print(f"{timestamp()} Connectivity report: {label}")
    print(f"  points: {n_points}")
    print(f"  retained arcs: {len(source)}")
    print(f"  connected components: {len(component_sizes)}")
    print(f"  largest-component fraction: {largest_fraction:.3f}")
    print(f"  largest component sizes: {component_preview}{preview_suffix}")

    return component_sizes


def validate_arc_quality(stm_arcs, label):
    """Validate that periodogram coherence filtering retains usable arcs.

    Parameters
    ----------
    stm_arcs : xarray.Dataset
        Arcs containing ``temp_coh``.
    label : str
        Network description used in diagnostics.

    Raises
    ------
    RuntimeError
        If no finite arc exceeds the configured quality threshold.

    Notes
    -----
    This guard runs before ``spatial_integration`` and replaces an indirect
    empty-array failure with a clear periodogram/network diagnostic.
    """
    quality = np.asarray(stm_arcs["temp_coh"].values)
    finite = np.isfinite(quality)
    retained = finite & (np.abs(quality) > arc_quality_threshold)

    print(
        f"{timestamp()} {label} arc-quality summary: "
        f"{finite.sum()}/{quality.size} finite; "
        f"{retained.sum()} above {arc_quality_threshold}."
    )

    if not np.any(retained):
        raise RuntimeError(
            f"No {label} arcs survived coherence filtering. Check the mother "
            "epoch, h2ph values, temporal baselines, and periodogram output."
        )


def _ambiguity_adaptation_with_iteration_limit(
    stm_arcs: xr.Dataset,
    stm_pnts: xr.Dataset,
    idx_refpnt: int,
    sparse_mode: bool,
    arc_estimation_method: str,
) -> tuple[xr.Dataset, xr.Dataset, int]:
    """Adapt network ambiguities with a bounded number of attempts per epoch.

    Parameters
    ----------
    stm_arcs : xarray.Dataset
        Integrated-network arc observations.
    stm_pnts : xarray.Dataset
        First-order network points.
    idx_refpnt : int
        Positional reference-point index.
    sparse_mode : bool
        Whether to construct sparse network matrices.
    arc_estimation_method : str
        Arc-estimation method; currently ``periodogram``.

    Returns
    -------
    tuple[xarray.Dataset, xarray.Dataset, int]
        Updated arcs, updated points, and reference index.

    Notes
    -----
    This runtime patch replaces ``depsi.network._ambiguity_adaptation`` with
    the safeguard from DePSI PR #143. The original loop could continue when
    the overall model test never passed; attempts are now capped by arc count.
    """
    A = depsi_network._network_relation_matrix(
        stm_arcs["source"],
        stm_arcs["target"],
        stm_pnts.sizes["space"],
        idx_refpnt,
        sparse_mode,
    )

    if arc_estimation_method == "periodogram":
        Qyy_diag = 1 - stm_arcs["temp_coh"].values
    else:
        raise NotImplementedError(
            f"arc_estimation_method '{arc_estimation_method}' "
            "is not supported in form_network."
        )

    invQy = np.diag(1 / Qyy_diag)

    # The reference point is excluded from the estimated point ambiguities.
    acheck = np.zeros(
        (
            stm_pnts.sizes["space"] - 1,
            stm_pnts.sizes["time"],
        )
    )

    stm_arcs_updated = stm_arcs.copy()
    stm_pnts_updated = stm_pnts.copy()

    for epoch in range(stm_pnts.sizes["time"]):
        depsi_network.logger.debug(
            f"Adjusting ambiguities for epoch {epoch}"
        )

        y = stm_arcs["ambiguities"].isel(time=epoch).data

        acheck_ifg, echeck_ifg, _ = (
            depsi_network._solve_float_ambiguities(
                A,
                y,
                invQy,
            )
        )

        OMT = echeck_ifg.T @ invQy @ echeck_ifg
        idx_previous_arc_fix = -1

        # PR #143 change 1/3: reset the counter for every epoch.
        counter = 0

        # PR #143 change 2/3: cap attempts at the network arc count.
        while (
            OMT >= depsi_network.OMT_THRES
            and counter < stm_arcs_updated.sizes["space"]
        ):
            idx_sort = np.argsort(np.abs(echeck_ifg))[::-1]
            idx_max_echeck = idx_sort[0]

            if idx_max_echeck == idx_previous_arc_fix:
                idx_max_echeck = idx_sort[1]

            if np.round(abs(echeck_ifg[idx_max_echeck])) >= 1:
                y[idx_max_echeck] -= np.round(
                    echeck_ifg[idx_max_echeck]
                )
            elif echeck_ifg[idx_max_echeck] > 0:
                y[idx_max_echeck] -= 1.0
            elif echeck_ifg[idx_max_echeck] < 0:
                y[idx_max_echeck] += 1.0

            idx_previous_arc_fix = idx_max_echeck

            acheck_ifg, echeck_ifg, _ = (
                depsi_network._solve_float_ambiguities(
                    A,
                    y,
                    invQy,
                )
            )

            OMT = echeck_ifg.T @ invQy @ echeck_ifg

            depsi_network.logger.debug(
                f"Fixing arc index {idx_max_echeck}, "
                f"new OMT={OMT:.2e}"
            )

            # PR #143 change 3/3: count this adaptation attempt.
            counter += 1

        stm_arcs_updated["ambiguities"][:, epoch] = y
        acheck[:, epoch] = acheck_ifg

    acheck = np.round(acheck).astype(np.int16)

    # Assign zero ambiguities to the reference point.
    acheck_full = np.zeros(
        (
            stm_pnts.sizes["space"],
            stm_pnts.sizes["time"],
        ),
        dtype=np.int16,
    )
    acheck_full[
        np.arange(stm_pnts.sizes["space"]) != idx_refpnt,
        :,
    ] = acheck

    stm_pnts_updated["ambiguities"] = (
        ["space", "time"],
        acheck_full,
    )

    return (
        stm_arcs_updated,
        stm_pnts_updated,
        idx_refpnt,
    )

# spatial_integration() resolves this module attribute at runtime.
depsi_network._ambiguity_adaptation = (
    _ambiguity_adaptation_with_iteration_limit
)
print(
    "Applied local PR #143 ambiguity-adaptation "
    "iteration-limit patch."
)



# -----------------------------------------------------------------------------
# Densification and masked DT model estimation
# -----------------------------------------------------------------------------


def densification_with_optional_mask(
    stm_dens_pnts: xr.Dataset,
    stm_network_pnts: xr.Dataset,
    idx_refpnt: int | None = None,
    n_connections: int = 1,
    key_xcoord: str = "azimuth",
    key_ycoord: str = "range",
    key_h2ph: str = "h2ph_values",
    key_sdphase: str = "sd_phase",
    key_Btemporal: str = "years",
    **kwargs_arc_estimation,
) -> xr.Dataset:
    # AML changes accepted
    """Densify DT or regional PS points against an integrated PS network.

    Parameters
    ----------
    stm_dens_pnts : xarray.Dataset
        Points requiring densification.
    stm_network_pnts : xarray.Dataset
        Integrated first-order PS network.
    idx_refpnt : int, optional
        Positional network reference index.
    n_connections : int, default 1
        Number of network connections per densified point.
    key_xcoord, key_ycoord : str
        Spatial coordinate variable names.
    key_h2ph, key_sdphase, key_Btemporal : str
        Observation-variable names used by periodogram estimation.
    **kwargs_arc_estimation
        Additional arguments forwarded to ``periodogram``.

    Returns
    -------
    xarray.Dataset
        Newly densified points without the attached PS network.

    Notes
    -----
    This adapts ``depsi.densification.densification``. It forwards an optional
    per-DT coherence mask, materializes UID overlap indices before ``isel`` to
    avoid unknown-size Dask indexing, and keeps PS, regional PS, and DT output
    schemas separate instead of concatenating the first-order network.
    """
    # AML changes accepted end
    if "pnt_uid" not in stm_dens_pnts:
        stm_dens_pnts = generate_pnt_uids(stm_dens_pnts)
    if "pnt_uid" not in stm_network_pnts:
        stm_network_pnts = generate_pnt_uids(stm_network_pnts)

    if "wavelength" not in stm_network_pnts.attrs:
        raise ValueError(
            "Wavelength is not provided and not found in the network STM attributes."
        )

    wavelength = stm_network_pnts.attrs["wavelength"]

    if idx_refpnt is None:
        if "idx_refpnt" in stm_network_pnts.attrs:
            idx_refpnt = stm_network_pnts.attrs["idx_refpnt"]
        else:
            raise ValueError(
                "The attribute 'idx_refpnt' was not found in stm_network_pnts. "
                "Please provide idx_refpnt explicitly."
            )

    if n_connections % 2 == 0:
        raise ValueError("n_connections should be an odd number.")
    if n_connections < 1:
        raise ValueError("n_connections should be at least 1.")
    if n_connections > 1:
        raise NotImplementedError("Currently only n_connections=1 is supported.")

    if not np.array_equal(stm_dens_pnts["time"].data, stm_network_pnts["time"].data):
        raise ValueError("stm_dens_pnts and stm_network_pnts must have the same time coordinates.")

    # - Remove points that are already in the first-order network.
    densification_uids = np.asarray(
        stm_dens_pnts["pnt_uid"].values
    )
    network_uids = np.asarray(
        stm_network_pnts["pnt_uid"].values
    )

    mask_network_overlap = np.isin(
        densification_uids,
        network_uids,
    )

    stm_dens_pnts = stm_dens_pnts.isel(
        space=np.flatnonzero(~mask_network_overlap)
    )

    # - Make a copy of the points that still need densification.
    stm_dens_pnts_output = stm_dens_pnts.copy()

    # - Connect each densification point to its nearest first-order network point.
    idx_dens_pnts, idx_network_pnts = _determine_dens_connections(
        stm_dens_pnts,
        stm_network_pnts,
        n_connections,
        key_xcoord,
        key_ycoord,
    )

    # - Form densification arcs.
    h2ph = (
        stm_dens_pnts[key_h2ph].isel(space=idx_dens_pnts).data
        + stm_network_pnts.isel(space=idx_network_pnts)[key_h2ph].data
    ) / 2

    dd_phase = compute_phase_difference(
        stm_network_pnts.isel(space=idx_network_pnts)[key_sdphase].data,
        stm_dens_pnts[key_sdphase].isel(space=idx_dens_pnts).data,
        "subtract",
    )

    Btemp = stm_dens_pnts[key_Btemporal].data

    stm_densification_arcs = xr.Dataset(
        coords={
            "idx_dens": (("space",), idx_dens_pnts),
            "idx_network": (("space",), idx_network_pnts),
            "Btemp": (("time",), Btemp),
        },
        data_vars={
            "h2ph": (("space", "time"), h2ph),
            "dd_phase": (("space", "time"), dd_phase),
        },
        attrs={"wavelength": wavelength},
    )

    # - Attach optional periodogram coherence mask.
    #   For DT points this is created from existing_flag before calling this
    #   function. The mask is copied onto the arc STM so periodogram can use it
    #   while estimating DT-network arc ambiguities.
    key_coh_mask = kwargs_arc_estimation.get("key_coh_mask")

    if key_coh_mask is not None:
        if key_coh_mask not in stm_dens_pnts:
            raise ValueError(
                f"Coherence mask variable '{key_coh_mask}' is not found "
                "in the densification points."
            )

        stm_densification_arcs[key_coh_mask] = (
            ("space", "time"),
            stm_dens_pnts[key_coh_mask].isel(space=idx_dens_pnts).data,
        )

    # - Estimate densification arc ambiguities and temporal coherence.
    _, ambiguities, _, _, temporal_coh_arc = periodogram(
        stm_densification_arcs,
        key_dphase="dd_phase",
        key_h2ph="h2ph",
        key_Btemporal="Btemp",
        **kwargs_arc_estimation,
    )

    stm_densification_arcs["ambiguities"] = ambiguities

    if n_connections == 1:
        stm_densification_arcs = stm_densification_arcs.drop_vars("idx_dens")

    # - Transfer network ambiguities to densification points.
    estimated_ambiguities = (
        stm_network_pnts.isel(space=stm_densification_arcs["idx_network"])["ambiguities"].data
        + stm_densification_arcs["ambiguities"].data
    )
    stm_dens_pnts_output["ambiguities"] = (("space", "time"), estimated_ambiguities)

    # - Compute densified unwrapped phase relative to the same reference point.
    ref_phase = stm_network_pnts[key_sdphase].isel(space=idx_refpnt).data
    unwrapped_phase_densification = (
        stm_dens_pnts_output[key_sdphase].data
        + estimated_ambiguities * 2 * np.pi
        - np.tile(ref_phase, (stm_dens_pnts_output.sizes["space"], 1))
    )
    stm_dens_pnts_output["unwrapped_phase"] = (
        ("space", "time"),
        unwrapped_phase_densification,
    )

    stm_dens_pnts_output["local_temp_coh"] = (
        ("space",),
        temporal_coh_arc.data,
    )

    # AML changes accepted
    # - Return only the newly densified points. The caller keeps the integrated
    #   first-order network as a separate global PS product.
    return stm_dens_pnts_output


def estimate_model_params_masked(
    stm: xr.Dataset,
    models: list[str] = None,
    key_observations: str = "unwrapped_phase",
    key_h2ph: str = "sd_h2ph",
    key_time: str = "temporal_baseline",
    key_mask: str = "dt_coh_mask",
    key_cross_range: str = "cr2ph",
    key_temperature: str = "temperature",
    wavelength: float | None = None,
) -> tuple[xr.Dataset, list]:
    """Estimate model parameters using each point's valid epochs only.

    Parameters
    ----------
    stm : xarray.Dataset
        Densified DT STM.
    models : list[str], optional
        Requested deformation-model components.
    key_observations, key_h2ph, key_time, key_mask : str
        Observation, height sensitivity, time, and validity-mask variables.
    key_cross_range, key_temperature : str
        Optional model-input variable names.
    wavelength : float, optional
        Radar wavelength in metres.

    Returns
    -------
    tuple[xarray.Dataset, list[str]]
        STM with fitted parameters/predictions/residuals and parameter names.

    Notes
    -----
    This adapts ``depsi.model_estimation.estimate_model_params``. The original
    all-epoch inversion propagates DT ``NaN`` values; this version builds a
    point-specific valid-epoch design matrix and restores full-length outputs
    with ``NaN`` at epochs where the DT does not exist.
    """
    if models is None:
        models = ["offset", "velocity", "height"]
    else:
        if len(models) == 0:
            raise ValueError("At least one model must be specified.")

        for model in models:
            if model not in MODEL_NAMES_PARAMS:
                raise NotImplementedError(
                    f"Model '{model}' is not supported. "
                    f"Available models are: {list(MODEL_NAMES_PARAMS.keys())}"
                )

        ordered_models = []
        for model_name in MODEL_NAMES_PARAMS:
            if model_name in models:
                ordered_models.append(model_name)
        models = ordered_models

    if wavelength is None:
        if "wavelength" in stm.attrs:
            wavelength = stm.attrs["wavelength"]
        else:
            raise ValueError("Wavelength must be provided either as an argument or in the dataset attributes.")

    m2ph = -4 * np.pi / wavelength

    param_names = []
    for model in models:
        param_names.extend(MODEL_NAMES_PARAMS[model])

    observations = stm[key_observations]

    #*--#
    # - Load the per-point/per-epoch validity mask.
    #   For DT points this is dt_coh_mask, which is built from existing_flag.
    #   Only epochs where this mask is True are allowed into the inversion.
    mask = stm[key_mask].astype(bool)
    #*--*#

    st_args = []
    st_args_keys = []

    if "height" in models:
        st_args.append(stm[key_h2ph])
        st_args_keys.append("h2ph")

    if "cross_range" in models:
        st_args.append(stm[key_cross_range])
        st_args_keys.append("cross_range")

    kwargs = {
        "models": models,
        "m2ph": m2ph,
        "time": stm[key_time].values,
        "st_args_keys": st_args_keys,
        "param_count": len(param_names),
        "n_time": stm.sizes["time"],
    }

    if "temperature" in models:
        kwargs["temperature"] = stm[key_temperature].values

    #*--#
    # - Apply the one-point model estimation with both observations and mask.
    #   The returned y_hat and residual arrays keep the full original time
    #   length; invalid epochs are left as NaN.
    params, yhat, ehat = xr.apply_ufunc(
        _estimate_model_params_one_point_masked,
        observations,
        mask,
        *st_args,
        input_core_dims=[["time"], ["time"]] + [["time"]] * len(st_args),
        output_core_dims=[["params"], ["time"], ["time"]],
        vectorize=True,
        dask="parallelized",
        output_dtypes=[float, float, float],
        dask_gufunc_kwargs={
            "output_sizes": {
                "params": len(param_names),
                "time": stm.sizes["time"],
            },
        },
        kwargs=kwargs,
    )
    #*--*#

    stm_out = stm.copy()

    for i, param_name in enumerate(param_names):
        stm_out = stm_out.assign({
            param_name: (["space"], params.isel(params=i).data)
        })

    stm_out = stm_out.assign({
        "expected_phases_yhat": (["space", "time"], yhat.data),
        "phase_residuals": (["space", "time"], ehat.data),
    })

    return stm_out, param_names


def _estimate_model_params_one_point_masked(
    obs: np.ndarray,
    mask: np.ndarray,
    *st_args,
    models: list[str],
    m2ph: float,
    time: np.ndarray,
    st_args_keys: list[str],
    param_count: int,
    n_time: int,
    temperature: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Fit model parameters for one point using finite, mask-valid epochs.

    Parameters
    ----------
    observations : numpy.ndarray
        One point's observation time series.
    mask : numpy.ndarray
        Boolean validity mask for the same epochs.
    models : list[str]
        Requested model components.
    m2ph : float
        Metre-to-phase conversion factor.
    time : numpy.ndarray
        Temporal baseline.
    st_args : tuple
        Additional spatial/temporal model inputs.
    st_args_keys : list[str]
        Names corresponding to ``st_args``.
    param_count : int
        Expected number of model parameters.

    Returns
    -------
    tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]
        Parameters, full-length predictions, and full-length residuals.

    Notes
    -----
    This is the point-level masked replacement used by
    ``estimate_model_params_masked``; invalid epochs remain ``NaN``.
    """
    obs = np.asarray(obs, dtype=float)

    #*--#
    # - Prepare full-length output arrays.
    #   y_hat and e_hat have the same length as the original time series.
    #   Epochs where the DT does not exist remain NaN.
    mask = np.asarray(mask, dtype=bool)

    params = np.full(param_count, np.nan, dtype=float)
    y_hat = np.full(n_time, np.nan, dtype=float)
    e_hat = np.full(n_time, np.nan, dtype=float)

    valid = mask & np.isfinite(obs)
    #*--*#

    model_inputs_full = dict(zip(st_args_keys, st_args, strict=True))

    #*--#
    # - Keep only epochs where all model inputs are finite.
    for key, value in model_inputs_full.items():
        value = np.asarray(value, dtype=float)
        valid &= np.isfinite(value)

    if temperature is not None:
        temperature = np.asarray(temperature, dtype=float)
        valid &= np.isfinite(temperature)

    n_valid = int(np.count_nonzero(valid))

    if n_valid < param_count:
        return params, y_hat, e_hat
    #*--*#

    #*--#
    # - Build the design matrix only for valid epochs.
    model_inputs_valid = {
        "n_epochs": n_valid,
        "time": time[valid],
    }

    if "h2ph" in model_inputs_full:
        model_inputs_valid["h2ph"] = np.asarray(model_inputs_full["h2ph"])[valid]

    if "cross_range" in model_inputs_full:
        model_inputs_valid["cross_range"] = np.asarray(model_inputs_full["cross_range"])[valid]

    if temperature is not None:
        model_inputs_valid["temperature"] = temperature[valid]

    A = construct_design_matrix(models, m2ph, **model_inputs_valid)
    obs_valid = obs[valid]
    #*--*#

    #*--#
    # - Estimate model parameters using only valid epochs.
    #   Write predictions and residuals back only at valid epochs.
    try:
        params_valid, _, _, _ = np.linalg.lstsq(A, obs_valid, rcond=None)
    except np.linalg.LinAlgError:
        return params, y_hat, e_hat

    y_hat_valid = A @ params_valid
    e_hat_valid = obs_valid - y_hat_valid.reshape(obs_valid.shape)

    params[:] = params_valid
    y_hat[valid] = y_hat_valid
    e_hat[valid] = e_hat_valid
    #*--*#

    return params, y_hat, e_hat



# -----------------------------------------------------------------------------
# Viewing geometry compatibility
# -----------------------------------------------------------------------------


def r_eci_to_ecef(epoch):
    # SciPy 1.17 requires vectorized single-axis Euler angles to have shape (N, 1).
    """Build ECI-to-ECEF rotations compatible with vectorized SciPy Euler input.

    Parameters
    ----------
    epoch : array-like
        Epochs used to calculate Greenwich apparent sidereal time.

    Returns
    -------
    scipy.spatial.transform.Rotation
        Rotation objects for the requested epochs.

    Notes
    -----
    This runtime patch replaces ``drama.geo.geometry.r_eci_to_ecef``. Newer
    SciPy requires vectorized single-axis Euler angles to have shape ``(N, 1)``;
    the original DRaMa input shape ``(N,)`` raised a ``ValueError``.
    """
    gm_st = (
        drama_geometry.Time(epoch)
        .sidereal_time("apparent", "greenwich")
        .radian
    )

    angles = np.expand_dims(
        -np.asarray(gm_st),
        axis=-1,
    )

    return drama_geometry.R.from_euler(
        "z",
        angles,
        degrees=False,
    )

drama_geometry.r_eci_to_ecef = r_eci_to_ecef


def add_local_viewing_geometry(
    stm: xr.Dataset,
    orbit_config_file: str,
    orbit_res: float,
    orbit_mode: str,
    orbit: str,
) -> xr.Dataset:
    """Add local incidence and azimuth angles for the requested Sentinel-1 track.

    Parameters
    ----------
    stm : xarray.Dataset
        Geocoded STM points.
    orbit_config_file : str
        DRaMa orbit configuration path.
    orbit_res : float
        Orbit sampling resolution.
    orbit_mode : str
        SAR acquisition mode; currently ``IWS``.
    orbit : str
        Requested labelled orbit, for example ``s1_dsc_t037``.

    Returns
    -------
    xarray.Dataset
        STM with local incidence and azimuth angles.

    Notes
    -----
    This local definition replaces
    ``depsi.viewing_geometry.add_local_viewing_geometry``. For each point, ASF
    provides the labelled Sentinel-1 tracks whose footprints cover that point,
    whereas DRaMa returns possible incidence/azimuth solutions without track
    labels.

    When the number of finite DRaMa solutions equals the number of covering
    same-direction tracks, both sets are ordered consistently: ASF tracks are
    ordered west-to-east and DRaMa solutions by incidence angle. The solution
    having the same rank as the requested track is then selected.

    At footprint-overlap boundaries, ASF may report multiple covering tracks
    while DRaMa returns only one finite geometry. In that case, the sole finite
    DRaMa solution is assigned to the requested track because no alternative
    finite solution exists at that point; the other nominal footprint matches
    do not produce a physically valid DRaMa geometry there.

    If the counts differ and DRaMa returns more than one finite solution, the
    function raises an error instead of guessing the track-to-solution mapping.
    """
    if orbit_mode != "IWS":
        raise ValueError(f"Orbit mode {orbit_mode} is not implemented.")

    np.NaN = np.nan

    longitudes = np.asarray(stm["lon"].values)
    latitudes = np.asarray(stm["lat"].values)
    orbits, footprints = identify_s1_orbits_in_aoi(longitudes, latitudes)

    if orbit not in orbits:
        raise ValueError(
            f"Requested orbit {orbit} was not found. Available orbits are {orbits}."
        )

    mode = SARModeFromCfg(cfg.ConfigFile(orbit_config_file), orbit_mode)
    _, _, asc_inc, desc_inc, asc_alpha, desc_alpha = nd.viewing_geometry(
        orbit_config_file,
        mode,
        orbit_res,
        latitudes,
        longitudes,
    )

    orbit_direction = orbit.split("_")[1]
    same_direction_orbits = [
        candidate
        for candidate in orbits
        if candidate.split("_")[1] == orbit_direction
    ]
    track_geometries = {
        candidate: unary_union(
            [Polygon(footprint) for footprint in footprints[candidate]]
        )
        for candidate in same_direction_orbits
    }

    if orbit_direction == "dsc":
        incidence_candidates = desc_inc
        azimuth_candidates = desc_alpha
    elif orbit_direction == "asc":
        incidence_candidates = asc_inc
        azimuth_candidates = asc_alpha
    else:
        raise ValueError(f"Cannot determine orbit direction from {orbit}.")

    incidence_angles = np.full(longitudes.shape, np.nan, dtype=float)
    azimuth_angles = np.full(longitudes.shape, np.nan, dtype=float)

    for point_index, (longitude, latitude) in enumerate(
        zip(longitudes, latitudes)
    ):
        point = Point(float(longitude), float(latitude))
        covering_orbits = [
            candidate
            for candidate in same_direction_orbits
            if track_geometries[candidate].covers(point)
        ]
        covering_orbits.sort(
            key=lambda candidate: track_geometries[candidate].bounds[0]
        )

        if orbit not in covering_orbits:
            raise ValueError(
                f"Point {point_index} at ({longitude:.6f}, {latitude:.6f}) "
                f"is not covered by requested orbit {orbit}. Covering "
                f"same-direction orbits are {covering_orbits}."
            )

        incidence_row = incidence_candidates[point_index]
        azimuth_row = azimuth_candidates[point_index]
        finite_columns = np.flatnonzero(
            np.isfinite(incidence_row) & np.isfinite(azimuth_row)
        )

        if len(finite_columns) == 1:
            selected_column = finite_columns[0]
        else:
            if len(finite_columns) != len(covering_orbits):
                raise ValueError(
                    f"Cannot match viewing geometries at point {point_index}: "
                    f"DRaMa returned {len(finite_columns)} finite solution(s), but "
                    f"ASF identified {len(covering_orbits)} covering "
                    f"{orbit_direction} track(s): {covering_orbits}."
                )

            ordered_columns = finite_columns[
                np.argsort(incidence_row[finite_columns])
            ]
            if orbit_direction == "asc":
                ordered_columns = ordered_columns[::-1]

            target_rank = covering_orbits.index(orbit)
            selected_column = ordered_columns[target_rank]
        incidence_angles[point_index] = incidence_row[selected_column]
        azimuth_angles[point_index] = azimuth_row[selected_column]

    return stm.assign(
        local_incidence_angle=(["space"], np.degrees(incidence_angles)),
        local_azimuth_angle=(["space"], np.degrees(azimuth_angles)),
    )



# -----------------------------------------------------------------------------
# Reference handling and export compatibility
# -----------------------------------------------------------------------------


def filter_stm_preserve_reference(
    stm,
    layer_to_filter,
    vmin=None,
    vmax=None,
):
    """Filter global PS points while preserving and re-indexing the reference PS.

    Parameters
    ----------
    stm : xarray.Dataset
        Global PS STM.
    layer_to_filter : str
        One-dimensional point variable used for filtering.
    vmin, vmax : float, optional
        Inclusive lower and upper limits.

    Returns
    -------
    xarray.Dataset
        Filtered STM with an updated positional ``idx_refpnt``.

    Notes
    -----
    This replaces ``stm_point_filter`` only for global PS filtering. The
    original function can remove the reference or leave its positional index
    stale after preceding rows are removed; this version retains and reindexes
    the phase datum explicitly.
    """
    if "idx_refpnt" not in stm.attrs:
        raise ValueError("The PS STM has no positional 'idx_refpnt' attribute.")

    reference_index = int(stm.attrs["idx_refpnt"])
    if not 0 <= reference_index < stm.sizes["space"]:
        raise IndexError(
            f"Reference index {reference_index} is outside the PS STM with "
            f"{stm.sizes['space']} points."
        )

    values = np.asarray(stm[layer_to_filter].values)
    if values.ndim != 1 or values.size != stm.sizes["space"]:
        raise ValueError(
            f"Filter layer '{layer_to_filter}' must have dimensions ('space',)."
        )

    mask = np.ones(stm.sizes["space"], dtype=bool)
    if vmin is not None:
        mask &= values >= vmin
    if vmax is not None:
        mask &= values <= vmax

    if not mask[reference_index]:
        print(
            f"{timestamp()} Reference PS falls outside the '{layer_to_filter}' "
            "filter limits and will be retained."
        )
    mask[reference_index] = True

    retained_positions = np.flatnonzero(mask)
    filtered = stm.isel(space=retained_positions)
    filtered.attrs = dict(stm.attrs)
    filtered.attrs["idx_refpnt"] = int(
        np.flatnonzero(retained_positions == reference_index)[0]
    )
    return filtered


def get_reference_coordinates(stm):
    """Return verified coordinates for the positional global PS reference.

    Parameters
    ----------
    stm : xarray.Dataset
        Global PS result containing latitude, longitude, and ``idx_refpnt``.

    Returns
    -------
    tuple[float, float]
        Reference latitude and longitude.

    Notes
    -----
    This replaces the lookup embedded in ``export_to_skygeo_portal``, which
    incorrectly used ``sel(space=idx_refpnt)`` although ``idx_refpnt`` is a
    position. Using ``isel`` resolves failures such as ``KeyError: 145`` and
    avoids silently selecting an unrelated coordinate label.
    """
    if "reference_lat" in stm.attrs and "reference_lon" in stm.attrs:
        return float(stm.attrs["reference_lat"]), float(stm.attrs["reference_lon"])

    reference_index = int(stm.attrs["idx_refpnt"])
    if not 0 <= reference_index < stm.sizes["space"]:
        raise IndexError(
            f"Reference index {reference_index} is outside the PS STM with "
            f"{stm.sizes['space']} points."
        )

    reference_lat = float(stm["lat"].isel(space=reference_index).values)
    reference_lon = float(stm["lon"].isel(space=reference_index).values)
    return reference_lat, reference_lon


# <>
def add_displacement_quality_metrics(stm, wavelength):
    """Add displacement-domain time series and export quality metrics.

    Parameters
    ----------
    stm : xarray.Dataset
        Model-estimated STM containing ``unwrapped_phase`` in radians and,
        when available, ``phase_residuals`` and ``local_temp_coh``.
    wavelength : float
        Radar wavelength in metres.

    Returns
    -------
    xarray.Dataset
        A shallow copy with LOS displacement in metres, LOS residuals and
        their standard deviation in millimetres, and temporal coherence.

    Notes
    -----
    This workflow helper complements ``depsi.model_estimation``. The DePSI
    model-estimation functions correctly retain ``unwrapped_phase`` and
    ``phase_residuals`` in radians, but STC expects a displacement time series
    in metres and the standard CSV expects displacement plus quality metrics.
    Invalid DT epochs remain NaN, so ``residual_std_mm`` uses only epochs for
    which the DT was processed. The population standard deviation (ddof=0)
    is used, matching NumPy/xarray's default definition.
    """
    if not np.isfinite(wavelength) or wavelength <= 0:
        raise ValueError("Wavelength must be a positive finite value.")
    if "unwrapped_phase" not in stm:
        raise KeyError("STM has no 'unwrapped_phase' variable.")

    output = stm.copy(deep=False)
    phase_to_m = -float(wavelength) / (4.0 * np.pi)

    output["displacement_los_m"] = phase_to_m * output["unwrapped_phase"]
    output["displacement_los_m"].attrs = {
        "long_name": "line-of-sight displacement",
        "units": "m",
    }

    if "phase_residuals" in output:
        output["residual_los_mm"] = (
            1000.0 * phase_to_m * output["phase_residuals"]
        )
        output["residual_los_mm"].attrs = {
            "long_name": "line-of-sight model residual",
            "units": "mm",
        }
        output["residual_std_mm"] = output["residual_los_mm"].std(
            dim="time",
            skipna=True,
            ddof=0,
        )
        output["residual_std_mm"].attrs = {
            "long_name": "standard deviation of LOS model residuals",
            "units": "mm",
        }

    if "local_temp_coh" in output:
        output["coherence"] = output["local_temp_coh"].astype(float)
        output["coherence"].attrs = {
            "long_name": "local temporal coherence",
            "units": "1",
        }

    return output


def prepare_stm_for_standard_exports(stm, product_name):
    """Prepare existing result variables for standard DePSI exporters.

    Parameters
    ----------
    stm : xarray.Dataset
        Product to export. Scientific displacement and residual quantities
        must already have been calculated by
        ``add_displacement_quality_metrics``.
    product_name : str
        Name used in informational messages.

    Returns
    -------
    xarray.Dataset
        Export-only shallow copy containing all scalar quality fields expected
        by the local CSV exporter.

    Notes
    -----
    This retains the previous workaround for ``depsi.io.export_to_csv`` and
    ``depsi.io.export_to_shapefile``, which always access ``stc`` and caused
    ``KeyError: 'stc'`` for DT products. It now also supplies NaN fallbacks
    when coherence or residual standard deviation is genuinely unavailable.
    No scientific unit conversion is performed here.
    """
    export_stm = stm.copy(deep=False)
    quality_variables = {
        "stc": "STC estimate",
        "coherence": "temporal coherence estimate",
        "residual_std_mm": "LOS residual standard deviation",
    }

    for variable, description in quality_variables.items():
        if variable in export_stm:
            continue

        print(
            f"{timestamp()} {product_name} has no {description}; "
            f"exporting {variable} as NaN."
        )
        export_stm[variable] = xr.DataArray(
            np.full(stm.sizes["space"], np.nan, dtype=float),
            dims=("space",),
            coords={"space": stm["space"]},
        )

    return export_stm


def _prepare_displacement_export_stm(stm, ts_proj):
    """Map explicit displacement variables to the legacy exporter fields."""
    if "displacement_los_m" not in stm:
        raise KeyError(
            "STM has no 'displacement_los_m'. Run "
            "add_displacement_quality_metrics() before export."
        )

    export_stm = stm.copy(deep=False)

    if ts_proj == "los":
        export_stm["unwrapped_phase"] = stm["displacement_los_m"]
    elif ts_proj == "vertical":
        if "local_incidence_angle" not in stm:
            raise KeyError(
                "Vertical displacement export requires "
                "'local_incidence_angle'."
            )
        export_stm["unwrapped_phase_pov"] = (
            stm["displacement_los_m"]
            / np.cos(np.radians(stm["local_incidence_angle"]))
        )
    else:
        raise ValueError(
            f"Unknown time-series projection '{ts_proj}'. "
            "Expected 'los' or 'vertical'."
        )

    return export_stm


def export_to_csv(
    stm,
    save_path,
    model_parameter_layer_names,
    ts_proj,
    point_annotation_label,
):
    """Export displacement and quality metrics to the standard DePSI CSV.

    Parameters
    ----------
    stm : xarray.Dataset
        Prepared result STM.
    save_path : str
        Destination CSV path.
    model_parameter_layer_names : tuple or list
        Model parameter variables written by the DePSI exporter.
    ts_proj : {"los", "vertical"}
        Requested displacement projection.
    point_annotation_label : str
        Prefix used for exported point identifiers.

    Returns
    -------
    None

    Notes
    -----
    This local function takes over ``depsi.io.export_to_csv`` while preserving
    its signature and output layout. The original writes radian-valued
    ``unwrapped_phase`` into the date columns and hardcodes ``Coherence [0-1]``
    and ``Std [mm]`` as NULL. This wrapper passes displacement in metres to the
    legacy writer, then fills those two columns from ``coherence`` and
    ``residual_std_mm``. It avoids changing the checked-out DePSI package.
    """
    export_stm = _prepare_displacement_export_stm(stm, ts_proj)
    _export_to_csv(
        stm=export_stm,
        save_path=save_path,
        model_parameter_layer_names=model_parameter_layer_names,
        ts_proj=ts_proj,
        point_annotation_label=point_annotation_label,
    )

    with open(save_path, newline="", encoding="utf-8") as file:
        rows = list(csv.reader(file))

    if not rows:
        raise ValueError(f"CSV exporter produced an empty file: {save_path}")

    headers = rows[0]
    coherence_column = headers.index("Coherence [0-1]")
    std_column = headers.index("Std [mm]")
    coherence_values = np.asarray(stm["coherence"].values).reshape(-1)
    std_values = np.asarray(stm["residual_std_mm"].values).reshape(-1)

    if len(rows) - 1 != stm.sizes["space"]:
        raise ValueError(
            "CSV row count does not match the STM space dimension: "
            f"{len(rows) - 1} versus {stm.sizes['space']}."
        )

    for row, coherence, residual_std in zip(
        rows[1:],
        coherence_values,
        std_values,
    ):
        row[coherence_column] = (
            f"{float(coherence):.4f}" if np.isfinite(coherence) else "NULL"
        )
        row[std_column] = (
            f"{float(residual_std):.3f}"
            if np.isfinite(residual_std)
            else "NULL"
        )

    with open(save_path, "w", newline="", encoding="utf-8") as file:
        csv.writer(file).writerows(rows)
# <|>


def export_to_skygeo_portal_with_reference(
    *,
    stm,
    reference_lat,
    reference_lon,
    save_path,
    **kwargs,
):
    """Export portal files with the verified global PS reference coordinates.

    Parameters
    ----------
    stm : xarray.Dataset
        PS, regional PS, or DT product to export.
    reference_lat, reference_lon : float
        Global PS reference coordinates.
    save_path : str or path-like
        Requested portal CSV path.
    **kwargs
        Remaining arguments for ``depsi.io.export_to_skygeo_portal``.

    Returns
    -------
    None

    Notes
    -----
    This wraps ``depsi.io.export_to_skygeo_portal``. The original treats a
    positional reference index as a ``space`` label and assumes every product
    contains the reference. The wrapper prevents that lookup error and writes
    the actual global PS reference into each portal JSON file. It also maps
    ``displacement_los_m`` to the legacy phase field so portal date columns
    contain displacement in metres instead of unwrapped phase in radians.
    """
    if stm.sizes.get("space", 0) == 0:
        raise ValueError("Cannot export an empty STM to the SkyGeo portal format.")

    export_stm = stm.copy(deep=False)
    export_stm.attrs = dict(stm.attrs)

    # <>
    # The original portal exporter also reads phase variables directly.
    # Supply displacement in metres while retaining its reference workaround.
    export_stm = _prepare_displacement_export_stm(
        export_stm,
        kwargs["ts_proj"],
    )
    # <|>
    first_space_label = export_stm["space"].values[0]
    if isinstance(first_space_label, np.generic):
        first_space_label = first_space_label.item()
    export_stm.attrs["idx_refpnt"] = first_space_label

    _export_to_skygeo_portal(
        stm=export_stm,
        save_path=save_path,
        **kwargs,
    )

    portal_csv_path = Path(save_path)
    if portal_csv_path.stem.split("_")[-1] != "portal":
        portal_csv_path = portal_csv_path.with_name(
            f"{portal_csv_path.stem}_portal.csv"
        )
    portal_json_path = portal_csv_path.with_suffix(".json")

    with portal_json_path.open("r", encoding="utf-8") as file:
        portal_metadata = json.load(file)
    portal_metadata["reference_point_location"] = (
        f"{reference_lat:.8f}, {reference_lon:.8f}"
    )
    with portal_json_path.open("w", encoding="utf-8") as file:
        json.dump(portal_metadata, file, ensure_ascii=False, indent=4)

print(f"{timestamp()} Reading metadata...")
metadata = sarxarray.read_metadata(metadata_path, driver="doris5")


# ========================================= #
#             0. Crop SLC Stack             #
# ========================================= #

if start_at_checkpoint < 1:
    print("========== Step 0: crop SLC stack ...")

    # - Load existing cropped SLC stack
    #   This avoids recropping the full original stack on repeated runs.
    if os.path.isdir(cropped_slc_path):
        print("Cropped SLC zarr already exists. Loading ...")
        cropped_slcs = xr.open_zarr(cropped_slc_path)

    else:
        # - Create AOI shapefile
        #   crop_slc_spacetime expects the spatial crop as a shapefile.
        print(f"{timestamp()} Creating AOI shapefile...")
        aoi_gdf = gpd.GeoDataFrame(
            {"id": [1]},
            geometry=[box(lon_min, lat_min, lon_max, lat_max)],
            crs="EPSG:4326",
        )
        aoi_gdf.to_file(aoi_file)

        # - Read original full SLC stack
        print(f"{timestamp()} Reading original SLC stack...")
        slcs = read_slc_stack(slc_path)

        # - Determine temporal crop end
        if last_date is None:
            last_date_crop = pd.Timestamp(slcs["time"].max().item()).to_pydatetime()
        else:
            last_date_crop = last_date

        # - Crop original SLC stack
        #   Original azimuth/range coordinate labels are preserved.
        print(f"{timestamp()} Cropping SLC stack...")
        cropped_slcs = crop_slc_spacetime(
            slcs,
            aoi_filename=aoi_file,
            start_date=first_date,
            end_date=last_date_crop,
        )

        # - Save cropped SLC stack
        #   The saved crop already contains complex/amplitude/phase.
        print(f"{timestamp()} Saving cropped SLC stack...")
        save_zarr_v2(cropped_slcs, cropped_slc_path)
        cropped_slcs = xr.open_zarr(cropped_slc_path)

    print(cropped_slcs)
    print(f"{timestamp()} Step 0 done.")

if start_at_checkpoint < 1:
    print("========== Step 1: PS and optional DT selection ...")

    # - Default to PS-only processing
    #   If DT input is found below, this is replaced by the common SLC/DT time axis.
    dt_stm = None
    slcs_for_processing = cropped_slcs
    validate_mother_epoch(
        slcs_for_processing["time"].values,
        "cropped SLC",
    )

    # - Optional DT selection
    #   DT is only created when the Radar Coding CSV exists.
    if dt_csv_exists():
        print(f"{timestamp()} Reading DT Radar Coding CSV...")
        dt = read_rcs_csv(dt_rcs_csv_path)

        # - Compare SLC and DT acquisition dates
        #   The cropped SLC stack may include newer epochs than the DT CSV.
        slc_times = pd.DatetimeIndex(cropped_slcs["time"].values)
        dt_times = pd.DatetimeIndex(dt["time"].values)

        epochs_only_in_slc = slc_times.difference(dt_times)
        epochs_only_in_dt = dt_times.difference(slc_times)

        if len(epochs_only_in_slc) > 0:
            print(
                f"{timestamp()} Epochs present in cropped SLC but missing from DT CSV "
                f"and left out of PS/DT processing: {format_epochs(epochs_only_in_slc)}"
            )

        if len(epochs_only_in_dt) > 0:
            print(
                f"{timestamp()} Epochs present in DT CSV but missing from cropped SLC "
                f"and left out of PS/DT processing: {format_epochs(epochs_only_in_dt)}"
            )

        # - Define common processing time axis
        #   Both PS and DT will use this time axis when DT is available.
        common_times = slc_times.intersection(dt_times)
        common_times = slc_times[slc_times.isin(common_times)]

        if len(common_times) == 0:
            print(
                f"{timestamp()} No common acquisition dates found between cropped SLC "
                "stack and DT CSV. Continuing with PS-only workflow."
            )
        else:
            print(
                f"{timestamp()} Using {len(common_times)} common epoch(s) for both "
                "PS and DT processing."
            )

            slcs_for_processing = cropped_slcs.sel(time=common_times)
            dt = dt.sel(time=common_times)
            validate_mother_epoch(
                slcs_for_processing["time"].values,
                "common PS/DT",
            )

            # - Select DT points from the cropped SLC stack
            #   designated_target_selection samples the SLC stack at the nearest
            #   radar pixel to each DT subpixel coordinate.
            print(f"{timestamp()} Selecting DT points...")
            dt_stm = designated_target_selection(slcs_for_processing, dt)

            # - Add DT single differences and time variables
            dt_stm = add_single_differences_and_time(dt_stm, mother_epoch)

            # - Add DT projected coordinates for atmosphere prediction
            dt_stm = add_euclidean_coordinates(dt_stm)

            # AML changes applied
            # - Add stable radar-pixel identities used to exclude DT overlaps
            #   from both global and regional PS selections.
            dt_stm = generate_pnt_uids(dt_stm)
            # AML changes applied end

            # - Save DT selection checkpoint
            dt_stm = make_dt_stm_zarr_safe(dt_stm)
            dt_stm = dt_stm.chunk({
                "time": -1,
                "space": dt_stm.sizes["space"],
            })
            save_zarr_v2(dt_stm, dt_stm_1sel_path)

            print(f"{timestamp()} Selected and saved {dt_stm.sizes['space']} DT point(s).")
    else:
        print(
            f"{timestamp()} No DT Radar Coding CSV found. "
            "Continuing with PS-only workflow."
        )

    # - Select PS points
    #   If DT exists, PS uses the same common time axis as DT.
    #   If DT does not exist, PS uses the full cropped SLC time axis.
    print(f"{timestamp()} Selecting PS points...")
    ps_stm = ps_selection(
        slcs_for_processing,
        method=ps_selection_method,
        threshold=threshold,
        ps_selection_start_date=start_date_ps_selection,
        ps_selection_end_date=end_date_ps_selection,
        output_chunks=chunks_ps_selection,
        mem_persist=False,
    )
    print(f"{timestamp()} Selected {len(ps_stm.space)} PS.")

    # - Add PS single differences and time variables
    ps_stm = add_single_differences_and_time(ps_stm, mother_epoch)

    # - Optional sidelobe removal
    print(f"{timestamp()} Removing sidelobes...")
    if do_sidelobe_detection:
        side_lobes_array, _ = detect_side_lobes(
            ps_stm,
            max_pixel_dist,
            min_correlation,
            "sd_complex",
            "sd_amplitude_unnormalized",
        )
        mask_sidelobes = np.ones(len(ps_stm.space), dtype=bool)
        mask_sidelobes[side_lobes_array] = False
    else:
        print(f"{timestamp()} Skipping sidelobe detection...")
        mask_sidelobes = np.ones(len(ps_stm.space), dtype=bool)

    ps_stm = ps_stm.isel(space=mask_sidelobes)

    # - Optional PS-DT overlap removal
    #   A DT should not also be treated as an ordinary PS point.
    if dt_stm is not None:
        ps_stm, n_global_dt_overlap = remove_radar_coordinate_overlaps(
            ps_stm,
            dt_stm,
        )

        print(
            f"{timestamp()} Removed {n_global_dt_overlap} global PS points "
            "overlapping DT radar pixels."
        )
    else:
        print(f"{timestamp()} No DT points available; skipping PS-DT overlap removal.")

    # - Geocode PS points
    #   PS coordinates are selected in radar geometry and geocoded here.
    print(f"{timestamp()} Geocoding PS points...")
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

    # - Add PS projected coordinates for atmosphere estimation
    ps_stm = add_euclidean_coordinates(ps_stm)

    # AML changes applied
    # - Add stable radar-pixel identities. Regional/global overlap remains
    #   temporarily allowed so these high-quality global points can still be
    #   selected into the post-APS first-order network.
    ps_stm = generate_pnt_uids(ps_stm)

    # - Optionally select additional PS candidates inside a metric regional
    #   buffer or user-supplied polygon using a separate, looser threshold.
    regional_ps_stm = None

    if regional_ps_enabled():
        print(f"{timestamp()} Building optional regional PS geometry...")
        regional_region_rd, regional_focus_rd = build_regional_ps_geometry()

        print(f"{timestamp()} Cropping SLC stack to the regional PS area...")
        regional_slcs = crop_slc_spacetime(
            slcs_for_processing,
            aoi_filename=regional_ps_region_wgs84_path,
            start_date=pd.Timestamp(
                slcs_for_processing["time"].min().item()
            ).to_pydatetime(),
            end_date=pd.Timestamp(
                slcs_for_processing["time"].max().item()
            ).to_pydatetime(),
        )

        print(
            f"{timestamp()} Selecting regional PS with "
            f"{regional_ps_selection_method.upper()} < "
            f"{regional_ps_threshold}..."
        )
        regional_ps_stm = ps_selection(
            regional_slcs,
            method=regional_ps_selection_method,
            threshold=regional_ps_threshold,
            ps_selection_start_date=start_date_ps_selection,
            ps_selection_end_date=end_date_ps_selection,
            output_chunks=chunks_ps_selection,
            mem_persist=False,
        )
        regional_ps_stm = add_single_differences_and_time(
            regional_ps_stm,
            mother_epoch,
        )

        # - Apply the same optional sidelobe filtering used for global PS.
        if do_sidelobe_detection:
            regional_sidelobes, _ = detect_side_lobes(
                regional_ps_stm,
                max_pixel_dist,
                min_correlation,
                "sd_complex",
                "sd_amplitude_unnormalized",
            )
            regional_keep = np.ones(
                regional_ps_stm.sizes["space"],
                dtype=bool,
            )
            regional_keep[regional_sidelobes] = False
            regional_ps_stm = regional_ps_stm.isel(space=regional_keep)

        # - Geocode regional candidates, add metric coordinates, and enforce
        #   the exact polygon rather than only its rectangular SLC crop.
        regional_ps_stm = regional_ps_stm.assign_attrs({
            "wavelength": metadata["wavelength"]
        })
        regional_latlonh = radar_to_latlonh(
            azimuth_coords=regional_ps_stm.azimuth.values,
            range_coords=regional_ps_stm.range.values,
            elevation=regional_ps_stm.h.values,
            metadata=metadata,
        )
        regional_ps_stm["lat"].data = regional_latlonh[0].flatten()
        regional_ps_stm["lon"].data = regional_latlonh[1].flatten()
        regional_ps_stm["h"].data = regional_latlonh[2].flatten()
        regional_ps_stm = add_euclidean_coordinates(regional_ps_stm)
        regional_ps_stm = select_ps_inside_region(
            regional_ps_stm,
            regional_region_rd,
        )

        # - DT has the highest classification priority. A radar pixel occupied
        #   by a DT must not also be processed as a regional PS.
        regional_ps_stm, n_regional_dt_overlap = (
            remove_radar_coordinate_overlaps(regional_ps_stm, dt_stm)
        )
        print(
            f"{timestamp()} Removed {n_regional_dt_overlap} regional PS "
            "points overlapping DT radar pixels."
        )

        regional_ps_stm = generate_pnt_uids(regional_ps_stm)
        regional_ps_stm = mark_regional_focus_ps(
            regional_ps_stm,
            regional_focus_rd,
        )

        print(
            f"{timestamp()} Selected {regional_ps_stm.sizes['space']} "
            "regional PS candidate(s) after polygon and DT filtering."
        )

        if regional_ps_stm.sizes["space"] == 0:
            regional_ps_stm = None
        else:
            regional_ps_stm = regional_ps_stm.chunk({
                "time": -1,
                "space": "auto",
            })
            save_zarr_v2(regional_ps_stm, regional_ps_stm_1sel_path)
    else:
        print(f"{timestamp()} Regional PS processing is disabled.")
    # AML changes applied end

    # - Save PS selection checkpoint
    ps_stm = ps_stm.chunk({"time": -1, "space": "auto"})
    save_zarr_v2(ps_stm, ps_stm_1sel_path)

    print(f"{timestamp()} Step 1 done.")

else:
    # - Load Step 1 checkpoints
    #   This restarts from the selected PS/DT points without recropping or reselection.
    print(f"{timestamp()} Loading Step 1 checkpoints...")
    ps_stm = xr.open_zarr(ps_stm_1sel_path)

    if dt_csv_exists() and os.path.isdir(dt_stm_1sel_path):
        dt_stm = xr.open_zarr(dt_stm_1sel_path)
        dt_stm = make_dt_stm_zarr_safe(dt_stm)
    else:
        dt_stm = None
        print(f"{timestamp()} No DT Step 1 checkpoint loaded; continuing PS-only.")

    # AML changes applied
    if regional_ps_enabled() and os.path.isdir(regional_ps_stm_1sel_path):
        regional_ps_stm = xr.open_zarr(regional_ps_stm_1sel_path)
    else:
        regional_ps_stm = None
        if regional_ps_enabled():
            print(
                f"{timestamp()} No regional PS Step 1 checkpoint found; "
                "regional PS processing is unavailable for this restart."
            )
    # AML changes applied end


# AML changes applied
# - Report the point populations available after Step 1.
#   This block is outside the fresh-selection and checkpoint-loading branches,
#   so the same summary is printed for both execution paths.
n_primary_ps = ps_stm.sizes["space"]
n_dt = dt_stm.sizes["space"] if dt_stm is not None else 0
n_regional_ps = (
    regional_ps_stm.sizes["space"]
    if regional_ps_stm is not None
    else 0
)

print("========== Step 1 selection summary ==========")
print(f"Primary/global PS selected: {n_primary_ps}")
print(f"DT selected:                {n_dt}")
print(f"Regional PS selected:       {n_regional_ps}")
print(
    "Note: global and regional PS may overlap until the post-APS "
    "first-order network is removed from the regional densification input."
)
# AML changes applied end


# ========================================= #
#        2. Atmospheric Phase Screen        #
# ========================================= #

if start_at_checkpoint < 2:
    print("========== Step 2: atmospheric phase screen estimation and removal ...")

    # - Verify PS/DT time-axis consistency
    #   Step 1 should already enforce this when DT is available.
    if dt_stm is not None and not np.array_equal(ps_stm["time"].values, dt_stm["time"].values):
        ps_times = pd.DatetimeIndex(ps_stm["time"].values)
        dt_times = pd.DatetimeIndex(dt_stm["time"].values)

        print(
            f"{timestamp()} Epochs present in PS STM but missing from DT STM: "
            f"{format_epochs(ps_times.difference(dt_times))}"
        )
        print(
            f"{timestamp()} Epochs present in DT STM but missing from PS STM: "
            f"{format_epochs(dt_times.difference(ps_times))}"
        )

        raise ValueError(
            "PS STM and DT STM do not have the same time axis. "
            "Recreate Step 1 checkpoints."
        )

    # AML changes applied
    # Regional PS are selected from the same processing SLC stack and must use
    # the same epochs before APS prediction and removal.
    if regional_ps_stm is not None and not np.array_equal(
        ps_stm["time"].values,
        regional_ps_stm["time"].values,
    ):
        ps_times = pd.DatetimeIndex(ps_stm["time"].values)
        regional_times = pd.DatetimeIndex(regional_ps_stm["time"].values)
        print(
            f"{timestamp()} Epochs present in global PS but missing from "
            f"regional PS: {format_epochs(ps_times.difference(regional_times))}"
        )
        print(
            f"{timestamp()} Epochs present in regional PS but missing from "
            f"global PS: {format_epochs(regional_times.difference(ps_times))}"
        )
        raise ValueError(
            "Global and regional PS STMs do not have the same time axis. "
            "Recreate Step 1 checkpoints."
        )
    # AML changes applied end

    # <>
    # - Reuse the expensive APS periodogram when both checkpoint stores exist.
    #   On the first retry these stores do not exist, so the network and
    #   periodogram are computed once and saved before spatial integration.
    aps_periodogram_checkpoint_exists = (
        os.path.isdir(aps_periodogram_points_path)
        and os.path.isdir(aps_periodogram_arcs_path)
    )

    if aps_periodogram_checkpoint_exists:
        print(f"{timestamp()} Loading APS periodogram checkpoints...")
        stm_network_pnts = xr.open_zarr(aps_periodogram_points_path).compute()
        stm_network_arcs = xr.open_zarr(aps_periodogram_arcs_path).compute()
        print(f"{timestamp()} APS periodogram checkpoints loaded.")
    else:
        # - Remove mother epoch from PS network STM
        #   The first-order network is formed on single-difference epochs only.
        print(f"{timestamp()} Removing mother epoch from PS network STM...")
        mother_epoch_index = np.where(
            ps_stm.time.values == ps_stm.sel(time=ps_stm.ps_sd_mother).time.values
        )[0][0]

        non_mother = [True] * len(ps_stm.time.values)
        non_mother[mother_epoch_index] = False
        ps_stm_without_mother_epoch = ps_stm.isel(time=non_mother)

        # - Select first-order PS network points
        #   These PS points are used to form the network for APS estimation.
        print(f"{timestamp()} Selecting first-order PS network points...")
        stm_network_pnts = network_stm_selection(
            stm=ps_stm_without_mother_epoch,
            min_dist=min_point_distance,
            crs=network_crs,
            azimuth_spacing=metadata["azimuth_pixel_spacing"],
            range_spacing=metadata["range_pixel_spacing"],
            sortby_var="time_selection_" + ps_selection_method,
            x_var=network_x_crds,
            y_var=network_y_crds,
            include_index=None,
        )

        # - Form first-order PS network arcs
        print(f"{timestamp()} Forming first-order PS network...")
        stm_network_arcs = form_network(
            stm_network_pnts,
            key_phase="sd_phase",
            key_h2ph="h2ph",
            key_Btemporal="temporal_baseline",
            max_length=max_arc_length,
            key_xcrds="lon",
            key_ycrds="lat",
            network_method=network_formation_method,
            min_links=min_network_links,
            num_partitions=network_partition_number,
            dphase_method=network_dphase_method,
        )

        # - Report connectivity before periodogram quality filtering.
        report_network_connectivity(
            stm_network_pnts,
            stm_network_arcs,
            label="APS first-order network before coherence filtering",
        )

        # - Estimate arc ambiguities and temporal coherence
        print(f"{timestamp()} Estimating ambiguities through periodogram...")
        _, ambiguities, _, _, ens_coh = periodogram(
            stm_network_arcs,
            key_dphase="d_phase",
            key_h2ph="h2ph",
            key_Btemporal="Btemp",
            min_steps=min_periodogram_iterations,
        )

        stm_network_arcs["ambiguities"] = ambiguities
        stm_network_arcs["temp_coh"] = ens_coh

        # - Compute both lazy periodogram outputs in one Dask computation.
        print(f"{timestamp()} Computing APS periodogram outputs...")
        stm_network_arcs = stm_network_arcs.compute()
        stm_network_pnts = stm_network_pnts.compute()
        print(f"{timestamp()} APS periodogram outputs computed.")

        print(f"{timestamp()} Saving APS periodogram checkpoints...")
        save_zarr_v2(stm_network_pnts, aps_periodogram_points_path)
        save_zarr_v2(stm_network_arcs, aps_periodogram_arcs_path)
        print(f"{timestamp()} APS periodogram checkpoints saved.")
    # <|>

    # AML changes applied
    # - Report the components that would survive the coherence threshold.
    report_network_connectivity(
        stm_network_pnts,
        stm_network_arcs,
        label="APS first-order network after coherence filtering",
        quality_key="temp_coh",
        quality_threshold=arc_quality_threshold,
    )
    validate_arc_quality(stm_network_arcs, "APS first-order network")
    # AML changes applied end

    # - Integrate PS network
    #   This produces unwrapped PS network point phases used to estimate APS.
    print(f"{timestamp()} Integrating the spatial network...")
    _, stm_pnts_output = spatial_integration(
        stm_network_pnts,
        stm_network_arcs,
        key_arc_quality="temp_coh",
        threshold_arc_quality=arc_quality_threshold,
        min_arc_connections=min_arc_connections,
        # AML changes applied
        largest_component_ratio=largest_component_ratio,
        skip_network_adaptation = skip_MHT_spatial_integration,
        # AML changes applied end
        idx_refpnt=reference_point_index,
    )

    # - Estimate atmosphere model from PS network residuals
    print(f"{timestamp()} Estimating mother atmosphere...")
    stm_pnts_output, _ = estimate_model_params(
        stm=stm_pnts_output,
        models=model_types,
        key_observations="unwrapped_phase",
        key_h2ph="sd_h2ph",
        key_time="temporal_baseline",
    )

    stm_pnts_output = stm_pnts_output.rename_vars({"pnt_offset": "mother_atmosphere"})

    # - Rename network dimension
    #   This avoids a dimension-name clash with prediction coordinates.
    stm_pnts_output = stm_pnts_output.rename_dims({"space": "space_fon"})

    # - Build atmosphere prediction coordinates
    #   APS is estimated from the global first-order PS network, but queried at
    #   all global PS, optional regional PS, and optional DT points. All global
    #   PS remain required because Step 3 performs a fresh post-APS selection.
    print(f"{timestamp()} Creating atmosphere prediction coordinates...")
    x_label = f"x_euclidean_proj_epsg{euclidean_epsg_code_number}"
    y_label = f"y_euclidean_proj_epsg{euclidean_epsg_code_number}"

    ps_prediction_coords = xr.Dataset(
        coords={
            x_label: ("space", ps_stm[x_label].values),
            y_label: ("space", ps_stm[y_label].values),
        }
    )

    # AML changes applied
    nps = ps_stm.sizes["space"]
    nregional = (
        regional_ps_stm.sizes["space"]
        if regional_ps_stm is not None
        else 0
    )
    ndt = dt_stm.sizes["space"] if dt_stm is not None else 0

    prediction_coord_parts = [ps_prediction_coords]

    if regional_ps_stm is not None:
        regional_prediction_coords = xr.Dataset(
            coords={
                x_label: ("space", regional_ps_stm[x_label].values),
                y_label: ("space", regional_ps_stm[y_label].values),
            }
        )
        prediction_coord_parts.append(regional_prediction_coords)

    if dt_stm is not None:
        dt_prediction_coords = xr.Dataset(
            coords={
                x_label: ("space", dt_stm[x_label].values),
                y_label: ("space", dt_stm[y_label].values),
            }
        )
        prediction_coord_parts.append(dt_prediction_coords)

    atmosphere_prediction_coords = xr.concat(
        prediction_coord_parts,
        dim="space",
    )

    print(
        f"{timestamp()} Atmosphere will be predicted at {nps} global PS, "
        f"{nregional} regional PS, and {ndt} DT point(s)."
    )
    # AML changes applied end

    # - Estimate atmospheric phase screens
    print(f"{timestamp()} Estimating atmospheric phase screens...")
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
                "variogram_parameters": None,
                "drift_terms": atmo_variogram_drift_terms,
            },
        },
        stm_coords_metadata={
            "mode": "euclidean",
            "x_label": x_label,
            "y_label": y_label,
            "space_dim_name": "space_fon",
        },
        prediction_coords_metadata={
            "mode": "euclidean",
            "x_label": x_label,
            "y_label": y_label,
            "space_dim_name": "space",
        },
    )

    # - Add zero atmosphere for the mother epoch
    #   In single-difference convention, the mother epoch atmosphere is zero.
    n_prediction_points = atmosphere_prediction_coords.sizes["space"]
    mother_atmo_stm = xr.Dataset(
        data_vars={
            "atmosphere_predicted": ("space", np.zeros(n_prediction_points)),
            "atmosphere_sigmasq": ("space", np.zeros(n_prediction_points)),
        },
        coords={"time": ps_stm.sel(time=ps_stm.ps_sd_mother).time.values},
    )

    atmosphere = xr.concat(
        [atmospheric_phase_screens, mother_atmo_stm],
        dim="time",
    ).sortby("time")

    # - Split atmosphere prediction back into global PS, optional regional PS,
    #   and optional DT products.
    ps_atmosphere = atmosphere.isel(space=slice(0, nps)).assign_coords(space=ps_stm.space)

    ps_stm["atmosphere_predicted"] = ps_atmosphere[
        "atmosphere_predicted"
    ].transpose("space", "time")
    ps_stm["atmosphere_sigmasq"] = ps_atmosphere[
        "atmosphere_sigmasq"
    ].transpose("space", "time")

    # AML changes applied
    regional_start = nps
    regional_stop = regional_start + nregional
    dt_start = regional_stop
    dt_stop = dt_start + ndt

    if regional_ps_stm is not None:
        regional_atmosphere = atmosphere.isel(
            space=slice(regional_start, regional_stop)
        ).assign_coords(space=regional_ps_stm.space)

        regional_ps_stm["atmosphere_predicted"] = regional_atmosphere[
            "atmosphere_predicted"
        ].transpose("space", "time")
        regional_ps_stm["atmosphere_sigmasq"] = regional_atmosphere[
            "atmosphere_sigmasq"
        ].transpose("space", "time")

    if dt_stm is not None:
        dt_atmosphere = atmosphere.isel(
            space=slice(dt_start, dt_stop)
        ).assign_coords(space=dt_stm.space)

        dt_stm["atmosphere_predicted"] = dt_atmosphere[
            "atmosphere_predicted"
        ].transpose("space", "time")
        dt_stm["atmosphere_sigmasq"] = dt_atmosphere[
            "atmosphere_sigmasq"
        ].transpose("space", "time")
    # AML changes applied end

    # - Remove atmospheric phase
    #   PS is always corrected. DT is corrected only when DT points are available.
    print(f"{timestamp()} Removing atmospheric phase screens...")
    ps_stm["sd_phase_minus_atmo"] = (
        ps_stm["sd_phase"] - ps_stm["atmosphere_predicted"] + np.pi
    ) % (2 * np.pi) - np.pi

    # AML changes applied
    if regional_ps_stm is not None:
        regional_ps_stm["sd_phase_minus_atmo"] = (
            regional_ps_stm["sd_phase"]
            - regional_ps_stm["atmosphere_predicted"]
            + np.pi
        ) % (2 * np.pi) - np.pi
    # AML changes applied end

    if dt_stm is not None:
        dt_stm["sd_phase_minus_atmo"] = (
            dt_stm["sd_phase"] - dt_stm["atmosphere_predicted"] + np.pi
        ) % (2 * np.pi) - np.pi

    # - Save atmosphere-corrected PS STM
    print(f"{timestamp()} Saving atmosphere-corrected PS STM...")
    ps_stm = ps_stm.chunk({"time": -1, "space": "auto"})
    save_zarr_v2(ps_stm, ps_stm_4atmo_path)

    # - Save atmosphere-corrected DT STM
    #   Only written when DT processing actually happened.
    if dt_stm is not None:
        print(f"{timestamp()} Saving atmosphere-corrected DT STM...")
        dt_stm = make_dt_stm_zarr_safe(dt_stm)
        dt_stm = dt_stm.chunk({
            "time": -1,
            "space": dt_stm.sizes["space"],
        })
        save_zarr_v2(dt_stm, dt_stm_4atmo_path)

    # AML changes applied
    # - Save atmosphere-corrected regional PS STM when enabled.
    if regional_ps_stm is not None:
        print(f"{timestamp()} Saving atmosphere-corrected regional PS STM...")
        regional_ps_stm = regional_ps_stm.chunk({
            "time": -1,
            "space": "auto",
        })
        save_zarr_v2(regional_ps_stm, regional_ps_stm_4atmo_path)
    # AML changes applied end

    print(
        f"{timestamp()} Finished global PS + optional regional PS/DT "
        "atmosphere correction."
    )

else:
    # - Load Step 2 checkpoints
    #   This restarts from the atmosphere-corrected products.
    print(f"{timestamp()} Loading atmosphere-corrected checkpoints...")
    ps_stm = xr.open_zarr(ps_stm_4atmo_path)

    if dt_csv_exists() and os.path.isdir(dt_stm_4atmo_path):
        dt_stm = xr.open_zarr(dt_stm_4atmo_path)
        dt_stm = make_dt_stm_zarr_safe(dt_stm)
    else:
        dt_stm = None
        print(f"{timestamp()} No DT atmosphere-corrected checkpoint loaded; continuing PS-only.")

    # AML changes applied
    if regional_ps_enabled() and os.path.isdir(regional_ps_stm_4atmo_path):
        regional_ps_stm = xr.open_zarr(regional_ps_stm_4atmo_path)
    else:
        regional_ps_stm = None
        if regional_ps_enabled():
            print(
                f"{timestamp()} No atmosphere-corrected regional PS checkpoint "
                "found; regional PS densification is unavailable for this restart."
            )
    # AML changes applied end
        
# ========================================= #
#  3. First-order PS network after APS      #
# ========================================= #

if start_at_checkpoint < 3:
    print("========== Step 3: first-order PS network after APS removal ...")

    # - Use the atmosphere-corrected PS STM.
    #   Step 2 leaves ps_stm in memory when it is run, and the checkpoint branch
    #   also loads ps_stm from ps_stm_4atmo_path. DT points are intentionally
    #   not used here because they should not participate in the first-order
    #   reference network.
    if "sd_phase_minus_atmo" not in ps_stm:
        raise ValueError(
            "Variable 'sd_phase_minus_atmo' is missing from the PS STM. "
            "Run Step 2 before Step 3."
        )

    # - Remove the mother epoch.
    #   Single-difference phase is zero at the mother epoch, so it is excluded
    #   from arc ambiguity estimation and spatial integration.
    print(f"{timestamp()} Removing mother epoch from PS STM...")

    mother_epoch_index = np.where(
        ps_stm.time.values == ps_stm.sel(time=ps_stm.ps_sd_mother).time.values
    )[0][0]

    non_mother = [True] * len(ps_stm.time.values)
    non_mother[mother_epoch_index] = False

    ps_stm_without_mother_epoch = ps_stm.isel(time=non_mother)
    ps_stm_without_mother_epoch = ps_stm_without_mother_epoch.chunk({"time": -1})

    # AML changes accepted
    # - Reserve the selected focus PS for regional densification.
    #   It remains in regional_ps_stm, but is removed from the global network
    #   candidate pool so it cannot disappear into the first-order PS product.
    network_candidate_stm = ps_stm_without_mother_epoch

    if regional_ps_stm is not None:
        if "is_focus_ps" not in regional_ps_stm:
            raise ValueError(
                "Regional PS processing is enabled, but 'is_focus_ps' is "
                "missing. Re-run Step 1 before forming the post-APS network."
            )

        focus_mask = regional_ps_stm["is_focus_ps"].values.astype(bool)
        if not focus_mask.any():
            raise ValueError(
                "No focus PS was selected within the configured search radius. "
                "Adjust the focus location, search radius, or regional threshold."
            )

        if "pnt_uid" not in network_candidate_stm:
            network_candidate_stm = generate_pnt_uids(network_candidate_stm)
        if "pnt_uid" not in regional_ps_stm:
            regional_ps_stm = generate_pnt_uids(regional_ps_stm)

        focus_uids = regional_ps_stm["pnt_uid"].values[focus_mask]
        keep_network_candidate = ~np.isin(
            network_candidate_stm["pnt_uid"].values,
            focus_uids,
        )
        network_candidate_stm = network_candidate_stm.isel(
            space=np.flatnonzero(keep_network_candidate)
        )

        print(
            f"{timestamp()} Reserved {len(focus_uids)} focus PS point(s) "
            "for regional densification."
        )
    # AML changes accepted end

    # - Select first-order PS network points.
    #   These are a spatially thinned subset of PS points used to form a stable
    #   reference network. The selection is PS-only.
    print(f"{timestamp()} Selecting first-order PS network points...")

    stm_network_pnts = network_stm_selection(
        stm=network_candidate_stm,
        min_dist=min_point_distance,
        crs=network_crs,
        azimuth_spacing=metadata["azimuth_pixel_spacing"],
        range_spacing=metadata["range_pixel_spacing"],
        sortby_var="time_selection_" + ps_selection_method,
        x_var=network_x_crds,
        y_var=network_y_crds,
        include_index=None,
    )

    # - Form first-order PS-PS network arcs.
    #   The arc phase is now formed from sd_phase_minus_atmo, not sd_phase,
    #   because atmospheric phase screens have already been removed.
    print(f"{timestamp()} Forming first-order PS network...")

    stm_network_arcs = form_network(
        stm_network_pnts,
        key_phase="sd_phase_minus_atmo",
        key_h2ph="h2ph",
        key_Btemporal="temporal_baseline",
        max_length=max_arc_length,
        key_xcrds="lon",
        key_ycrds="lat",
        network_method=network_formation_method,
        min_links=min_network_links,
        num_partitions=network_partition_number,
        dphase_method=network_dphase_method,
    )

    # AML changes applied
    # - Report connectivity before periodogram quality filtering.
    report_network_connectivity(
        stm_network_pnts,
        stm_network_arcs,
        label="post-APS first-order network before coherence filtering",
    )
    # AML changes applied end

    print(
        f"{timestamp()} Formed first-order PS network with "
        f"{len(stm_network_pnts.space)} points and "
        f"{len(stm_network_arcs.space)} arcs."
    )

    # - Estimate PS-PS arc ambiguities and temporal coherence.
    #   No DT coherence mask is used here. This first-order network contains
    #   only PS points, and all selected PS epochs are treated as valid.
    print(f"{timestamp()} Estimating ambiguities through periodogram...")

    _, ambiguities, _, _, ens_coh = periodogram(
        stm_network_arcs,
        key_dphase="d_phase",
        key_h2ph="h2ph",
        key_Btemporal="Btemp",
        min_steps=min_periodogram_iterations,
    )

    stm_network_arcs["ambiguities"] = ambiguities
    stm_network_arcs["temp_coh"] = ens_coh

    # AML changes accepted
    # - Materialize both post-APS periodogram outputs in one Dask computation.
    #   The following connectivity report then reads in-memory arrays and does
    #   not start a second partial periodogram calculation.
    print(f"{timestamp()} Computing post-APS periodogram outputs...")
    stm_network_arcs = stm_network_arcs.compute()
    stm_network_pnts = stm_network_pnts.compute()
    print(f"{timestamp()} Post-APS periodogram outputs computed.")
    # AML changes accepted end

    # AML changes applied
    # - Report the components that would survive the coherence threshold.
    report_network_connectivity(
        stm_network_pnts,
        stm_network_arcs,
        label="post-APS first-order network after coherence filtering",
        quality_key="temp_coh",
        quality_threshold=arc_quality_threshold,
    )
    validate_arc_quality(stm_network_arcs, "post-APS first-order network")
    # AML changes applied end

    # - Spatially integrate the first-order PS network.
    #   Periodogram provides unwrapped PS-PS arc phase differences.
    #   Spatial integration converts those arc differences into unwrapped
    #   single-difference phase time series for the selected PS network points.
    print(f"{timestamp()} Integrating the spatial network...")

    _, stm_firstordernetwork = spatial_integration(
        stm_network_pnts,
        stm_network_arcs,
        key_arc_quality="temp_coh",
        threshold_arc_quality=arc_quality_threshold,
        min_arc_connections=min_arc_connections,
        # AML changes applied
        largest_component_ratio=largest_component_ratio,
        # AML changes applied end
        idx_refpnt=reference_point_index,
    )

    # AML changes applied
    # - Report the surviving post-APS first-order network size separately from
    #   the full primary/global PS candidate population reported after Step 1.
    print(
        f"{timestamp()} Post-APS first-order network contains "
        f"{stm_firstordernetwork.sizes['space']} PS point(s)."
    )
    # AML changes applied end

    # - Save the first-order PS network checkpoint.
    #   This checkpoint is the input network for the next densification step.
    print(f"{timestamp()} Saving first-order PS network checkpoint...")

    stm_firstordernetwork = stm_firstordernetwork.chunk({
        "time": -1,
        "space": "auto",
    })
    save_zarr_v2(stm_firstordernetwork, ps_stm_5fnet_path)

    print(f"{timestamp()} Step 3 done.")

else:
    # - Load Step 3 checkpoint.
    #   This allows restarting directly from the integrated first-order network
    #   before continuing to densification.
    print(f"{timestamp()} Loading first-order PS network checkpoint...")

    stm_firstordernetwork = xr.open_zarr(ps_stm_5fnet_path)

    print(f"{timestamp()} Step 3 checkpoint loaded.")
    # AML changes applied
    print(
        f"{timestamp()} Post-APS first-order network contains "
        f"{stm_firstordernetwork.sizes['space']} PS point(s)."
    )
    # AML changes applied end
    
# ========================================= #
#  4. Densification and model estimation    #
# ========================================= #

if start_at_checkpoint < 4:
    print("========== Step 4: PS/regional-PS/DT densification and model estimation ...")

    # <>
    # Reopen the Step 2 and Step 3 products before densification. Zarr stores
    # concrete dimension sizes, so this removes unknown-size Dask metadata
    # inherited from lazy regional PS selection. It prevents the previous:
    #   ValueError: conflicting sizes for dimension 'space': length nan
    print(
        f"{timestamp()} Reopening Step 4 densification inputs "
        "from saved checkpoints..."
    )

    stm_firstordernetwork = open_densification_checkpoint(
        path=ps_stm_5fnet_path,
        label="Post-APS first-order PS network",
    )

    if regional_ps_enabled() and os.path.isdir(
        regional_ps_stm_4atmo_path
    ):
        regional_ps_stm = open_densification_checkpoint(
            path=regional_ps_stm_4atmo_path,
            label="Atmosphere-corrected regional PS STM",
        )
    else:
        regional_ps_stm = None

        if regional_ps_enabled():
            print(
                f"{timestamp()} No regional PS atmosphere-corrected "
                "checkpoint found; regional densification is skipped."
            )

    if dt_csv_exists() and os.path.isdir(dt_stm_4atmo_path):
        dt_stm = open_densification_checkpoint(
            path=dt_stm_4atmo_path,
            label="Atmosphere-corrected DT STM",
        )
        dt_stm = make_dt_stm_zarr_safe(dt_stm)
    else:
        dt_stm = None
        print(
            f"{timestamp()} No DT atmosphere-corrected checkpoint "
            "found; DT densification is skipped."
        )

    print(
        f"{timestamp()} Step 4 densification inputs reopened."
    )
    # <|>

    # - Build the common non-mother time mask for regional PS and DT points.
    #   The first-order network checkpoint already has the mother epoch removed.
    #   Densification and model estimation must therefore use the same non-mother
    #   time axis.
    print(f"{timestamp()} Preparing densification time axis...")

    mother_epoch_index = np.where(
        ps_stm.time.values == ps_stm.sel(time=ps_stm.ps_sd_mother).time.values
    )[0][0]

    non_mother = [True] * len(ps_stm.time.values)
    non_mother[mother_epoch_index] = False

    # AML changes accepted
    # - Keep the integrated post-APS network as the global PS product.
    #   Regional PS and DT points are densified and saved independently below,
    #   so their later STC/filter choices cannot alter the global PS output.
    ps_stm_densified = stm_firstordernetwork.copy()

    print(f"{timestamp()} Estimating global PS model parameters...")

    ps_stm_densified, model_parameter_layer_names = estimate_model_params(
        stm=ps_stm_densified,
        models=model_types,
        key_observations="unwrapped_phase",
        key_h2ph="sd_h2ph",
        key_time="temporal_baseline",
        wavelength=wavelength,
    )

    # <>
    # Model estimation returns phase-domain observations and residuals. Add
    # explicit displacement-domain variables before saving the checkpoint.
    ps_stm_densified = add_displacement_quality_metrics(
        ps_stm_densified,
        wavelength,
    )
    # <|>

    print(f"{timestamp()} Saving global PS checkpoint...")

    ps_stm_densified = ps_stm_densified.chunk({
        "time": -1,
        "space": "auto",
    })
    save_zarr_v2(ps_stm_densified, ps_stm_5dens_path)

    # - Densify regional PS points without attaching the first-order network.
    #   Network points are removed after the final network is known. The focus
    #   PS was excluded from that network in Step 3 and must remain here.
    regional_ps_stm_densified = None

    if regional_ps_stm is not None:
        print(f"{timestamp()} Preparing regional PS densification points...")

        regional_ps_without_mother_epoch = regional_ps_stm.isel(time=non_mother)
        regional_ps_without_mother_epoch = (
            regional_ps_without_mother_epoch.chunk({"time": -1})
        )

        if "pnt_uid" not in stm_firstordernetwork:
            stm_firstordernetwork = generate_pnt_uids(stm_firstordernetwork)
        if "pnt_uid" not in regional_ps_without_mother_epoch:
            regional_ps_without_mother_epoch = generate_pnt_uids(
                regional_ps_without_mother_epoch
            )

        if "is_focus_ps" not in regional_ps_without_mother_epoch:
            raise ValueError(
                "Regional PS STM has no 'is_focus_ps' marker. Re-run Step 1."
            )

        focus_mask = regional_ps_without_mother_epoch[
            "is_focus_ps"
        ].values.astype(bool)
        if not focus_mask.any():
            raise ValueError(
                "No focus PS was selected. Adjust the focus configuration and "
                "re-run from Step 1."
            )

        focus_uids = regional_ps_without_mother_epoch["pnt_uid"].values[
            focus_mask
        ]
        network_overlap = np.isin(
            regional_ps_without_mother_epoch["pnt_uid"].values,
            stm_firstordernetwork["pnt_uid"].values,
        )
        regional_ps_to_densify = regional_ps_without_mother_epoch.isel(
            space=np.flatnonzero(~network_overlap)
        )

        print(
            f"{timestamp()} Removed {int(network_overlap.sum())} regional PS "
            "point(s) already present in the post-APS first-order network."
        )

        focus_is_available = np.isin(
            focus_uids,
            regional_ps_to_densify["pnt_uid"].values,
        ).all()
        if not focus_is_available:
            raise ValueError(
                "The selected focus PS entered the post-APS network and was "
                "removed from regional densification. Re-run Step 3 with the "
                "focus-reservation change enabled."
            )

        print(
            f"{timestamp()} Focus PS retained in the regional densification "
            "product."
        )
        print(f"{timestamp()} Starting regional PS densification...")

        regional_ps_stm_densified = densification_with_optional_mask(
            regional_ps_to_densify,
            stm_firstordernetwork,
            idx_refpnt=reference_point_index,
            n_connections=n_densification_connections,
            key_sdphase="sd_phase_minus_atmo",
            key_h2ph="sd_h2ph",
            key_Btemporal="temporal_baseline",
            min_steps=min_periodogram_iterations,
        )

        print(f"{timestamp()} Estimating regional PS model parameters...")
        regional_ps_stm_densified, _ = estimate_model_params(
            stm=regional_ps_stm_densified,
            models=model_types,
            key_observations="unwrapped_phase",
            key_h2ph="sd_h2ph",
            key_time="temporal_baseline",
            wavelength=wavelength,
        )

        # <>
        regional_ps_stm_densified = add_displacement_quality_metrics(
            regional_ps_stm_densified,
            wavelength,
        )
        # <|>

        print(f"{timestamp()} Saving regional PS checkpoint...")
        regional_ps_stm_densified = regional_ps_stm_densified.chunk({
            "time": -1,
            "space": "auto",
        })
        save_zarr_v2(
            regional_ps_stm_densified,
            regional_ps_stm_5dens_path,
        )
    else:
        print(f"{timestamp()} No regional PS STM available; skipping it.")
    # AML changes accepted end

    # - Densify and model DT points only when DT processing is available.
    #   The helper returns only DT points and therefore keeps this product
    #   separate from the global and regional PS products.
    if dt_stm is not None:
        print(f"{timestamp()} Preparing DT densification points...")

        dt_stm_without_mother_epoch = dt_stm.isel(time=non_mother)
        dt_stm_without_mother_epoch = dt_stm_without_mother_epoch.chunk({"time": -1})

        # - Build the DT coherence mask.
        #   True means this DT exists at that epoch and can be used by
        #   periodogram. False epochs are ignored and remain NaN in the
        #   ambiguity and unwrapped-phase output.
        dt_stm_without_mother_epoch["dt_coh_mask"] = (
            dt_stm_without_mother_epoch["existing_flag"].astype(bool)
        )

        print(f"{timestamp()} Starting DT densification...")

        dt_stm_densified = densification_with_optional_mask(
            dt_stm_without_mother_epoch,
            stm_firstordernetwork,
            idx_refpnt=reference_point_index,
            n_connections=n_densification_connections,
            key_sdphase="sd_phase_minus_atmo",
            key_h2ph="sd_h2ph",
            key_Btemporal="temporal_baseline",
            key_coh_mask="dt_coh_mask",
            min_steps=min_periodogram_iterations,
        )

        # - Estimate model parameters for DT points.
        #   DT unwrapped phase can contain NaNs where existing_flag == 0, so
        #   model estimation uses only epochs where dt_coh_mask is True and
        #   the required values are finite.
        print(f"{timestamp()} Estimating DT model parameters...")

        dt_stm_densified, dt_model_parameter_layer_names = estimate_model_params_masked(
            stm=dt_stm_densified,
            models=model_types,
            key_observations="unwrapped_phase",
            key_h2ph="sd_h2ph",
            key_time="temporal_baseline",
            key_mask="dt_coh_mask",
            wavelength=wavelength,
        )

        # <>
        # DT residual statistics ignore non-existing epochs because the masked
        # estimator leaves those phase residuals as NaN.
        dt_stm_densified = add_displacement_quality_metrics(
            dt_stm_densified,
            wavelength,
        )
        # <|>

        # - Save DT densification/model checkpoint.
        print(f"{timestamp()} Saving DT densification checkpoint...")

        dt_stm_densified = make_dt_stm_zarr_safe(dt_stm_densified)
        dt_stm_densified = dt_stm_densified.chunk({
            "time": -1,
            "space": "auto",
        })
        save_zarr_v2(dt_stm_densified, dt_stm_5dens_path)
    else:
        dt_stm_densified = None
        dt_model_parameter_layer_names = []
        print(f"{timestamp()} No DT STM available; skipping DT densification.")

    # AML changes accepted
    # The three products stay independent from this point onward.
    print(
        f"{timestamp()} Step 4 done: global PS="
        f"{ps_stm_densified.sizes['space']}, regional PS="
        f"{0 if regional_ps_stm_densified is None else regional_ps_stm_densified.sizes['space']}, "
        f"DT={0 if dt_stm_densified is None else dt_stm_densified.sizes['space']}."
    )
    # AML changes accepted end

else:
    # AML changes accepted
    # - Reload each independent Step 4 product for checkpoint restarts.
    print(f"{timestamp()} Loading Step 4 product checkpoints...")

    ps_stm_densified = xr.open_zarr(ps_stm_5dens_path)

    if regional_ps_enabled() and os.path.isdir(regional_ps_stm_5dens_path):
        regional_ps_stm_densified = xr.open_zarr(
            regional_ps_stm_5dens_path
        )
    else:
        regional_ps_stm_densified = None

    if dt_csv_exists() and os.path.isdir(dt_stm_5dens_path):
        dt_stm_densified = xr.open_zarr(dt_stm_5dens_path)
        dt_stm_densified = make_dt_stm_zarr_safe(dt_stm_densified)
    else:
        dt_stm_densified = None
    # AML changes accepted end
        
# ========================================= #
#   5. Geocoding, STC, filtering, geometry  #
# ========================================= #

if start_at_checkpoint < 5:
    print("========== Step 5: geocoding, filtering, and viewing geometry ...")

    # - Load Step 4 products.
    print(f"{timestamp()} Loading global PS checkpoint...")
    ps_stm_densified = xr.open_zarr(ps_stm_5dens_path)

    # AML changes accepted
    # Regional PS are loaded as a separate product. Their STC and filtering
    # settings below are independent from the global PS settings.
    if regional_ps_enabled() and os.path.isdir(regional_ps_stm_5dens_path):
        print(f"{timestamp()} Loading regional PS densified checkpoint...")
        regional_ps_stm_densified = xr.open_zarr(
            regional_ps_stm_5dens_path
        )
    else:
        regional_ps_stm_densified = None
    # AML changes accepted end

    if dt_csv_exists() and os.path.isdir(dt_stm_5dens_path):
        print(f"{timestamp()} Loading DT densified checkpoint...")
        dt_stm_densified = xr.open_zarr(dt_stm_5dens_path)
        dt_stm_densified = make_dt_stm_zarr_safe(dt_stm_densified)
    else:
        dt_stm_densified = None

    # <>
    # Support Step 4 checkpoints created before displacement variables were
    # introduced. The helper is idempotent for newly generated checkpoints.
    ps_stm_densified = add_displacement_quality_metrics(
        ps_stm_densified,
        wavelength,
    )
    if regional_ps_stm_densified is not None:
        regional_ps_stm_densified = add_displacement_quality_metrics(
            regional_ps_stm_densified,
            wavelength,
        )
    if dt_stm_densified is not None:
        dt_stm_densified = add_displacement_quality_metrics(
            dt_stm_densified,
            wavelength,
        )
    # <|>

    # - Rebuild model parameter layer names.
    model_parameter_layer_names = []
    for model_name in MODEL_NAMES_PARAMS:
        if model_name in model_types:
            model_parameter_layer_names.extend(MODEL_NAMES_PARAMS[model_name])

    # - Geocode PS points.
    print(f"{timestamp()} Starting PS geocoding...")
    latlonh = radar_to_latlonh(
        azimuth_coords=ps_stm_densified.azimuth.values,
        range_coords=ps_stm_densified.range.values,
        elevation=ps_stm_densified.pnt_height.values,
        metadata=metadata,
    )

    ps_stm_densified["lat"].data = latlonh[0].flatten()
    ps_stm_densified["lon"].data = latlonh[1].flatten()
    ps_stm_densified["h"].data = latlonh[2].flatten()

    # - Calculate STC for global PS only when configured.
    print(f"{timestamp()} Calculating PS unwrapped phase time series...")
    ps_stm_densified["unwrapped_phase"].data = ps_stm_densified["unwrapped_phase"].values

    if ps_compute_stc:
        print(f"{timestamp()} Calculating PS spatiotemporal consistency...")
        ps_stm_densified = compute_spatiotemporal_consistency(
            ps_stm_densified,
            min_dist=stc_min_dist,
            max_dist=stc_max_dist,
            x_crd_layer_name="lon",
            y_crd_layer_name="lat",
            coordinate_type="geographic",
            # <>
            # DePSI STC multiplies its input scatter by 1000 and therefore
            # expects displacement in metres, not unwrapped phase in radians.
            key_unwrapped_phase="displacement_los_m",
            # <|>
        )

    # - Filter PS scatterers only.
    print(f"{timestamp()} Filtering PS scatterers...")
    for layer in ps_filter_dict:
        # <>
        ps_stm_densified = filter_stm_preserve_reference(
            stm=ps_stm_densified,
            layer_to_filter=layer,
            vmin=ps_filter_dict[layer][0],
            vmax=ps_filter_dict[layer][1],
        )
        # <|>

    # <>
    reference_lat, reference_lon = get_reference_coordinates(ps_stm_densified)
    ps_stm_densified.attrs["reference_lat"] = reference_lat
    ps_stm_densified.attrs["reference_lon"] = reference_lon
    # <|>

    # - Add PS viewing geometry.
    print(f"{timestamp()} Adding PS viewing geometry...")
    ps_stm_densified = add_local_viewing_geometry(
        stm=ps_stm_densified,
        orbit_config_file=orbit_file,
        orbit_res=orbit_resolution,
        orbit_mode=orbit_mode,
        orbit=f"{satellite}_{direction}_t{track:0>3d}",
    )

    # - Project PS phase and model parameters onto vertical.
    print(f"{timestamp()} Projecting PS phases and models onto vertical...")
    ps_stm_densified["unwrapped_phase_pov"] = ps_stm_densified["unwrapped_phase"] / np.cos(
        np.radians(ps_stm_densified["local_incidence_angle"])
    )

    for param in model_parameter_layer_names:
        ps_stm_densified[param].data = ps_stm_densified[param].values
        ps_stm_densified[f"{param}_pov"] = ps_stm_densified[param] / np.cos(
            np.radians(ps_stm_densified["local_incidence_angle"])
        )

    # - Save PS result.
    print(f"{timestamp()} Saving PS result checkpoint...")
    ps_stm_densified = ps_stm_densified.chunk({"time": -1, "space": "auto"})
    save_zarr_v2(ps_stm_densified, ps_stm_result_path)

    # AML changes accepted
    # - Process regional PS independently from global PS.
    #   This preserves the focus-region product even when different STC or
    #   scatterer filters are desired for the two PS populations.
    if regional_ps_stm_densified is not None:
        print(f"{timestamp()} Starting regional PS geocoding...")
        latlonh = radar_to_latlonh(
            azimuth_coords=regional_ps_stm_densified.azimuth.values,
            range_coords=regional_ps_stm_densified.range.values,
            elevation=regional_ps_stm_densified.pnt_height.values,
            metadata=metadata,
        )

        regional_ps_stm_densified["lat"].data = latlonh[0].flatten()
        regional_ps_stm_densified["lon"].data = latlonh[1].flatten()
        regional_ps_stm_densified["h"].data = latlonh[2].flatten()
        regional_ps_stm_densified["unwrapped_phase"].data = (
            regional_ps_stm_densified["unwrapped_phase"].values
        )

        if regional_ps_compute_stc:
            print(
                f"{timestamp()} Calculating regional PS spatiotemporal "
                "consistency..."
            )
            regional_ps_stm_densified = compute_spatiotemporal_consistency(
                regional_ps_stm_densified,
                min_dist=stc_min_dist,
                max_dist=stc_max_dist,
                x_crd_layer_name="lon",
                y_crd_layer_name="lat",
                coordinate_type="geographic",
                # <>
                # Keep regional STC in millimetres by supplying LOS
                # displacement in metres to the unchanged DePSI function.
                key_unwrapped_phase="displacement_los_m",
                # <|>
            )

        print(f"{timestamp()} Filtering regional PS scatterers...")
        for layer in regional_ps_filter_dict:
            regional_ps_stm_densified = stm_point_filter(
                stm=regional_ps_stm_densified,
                layer_to_filter=layer,
                vmin=regional_ps_filter_dict[layer][0],
                vmax=regional_ps_filter_dict[layer][1],
                return_removed=False,
            )

        print(f"{timestamp()} Adding regional PS viewing geometry...")
        regional_ps_stm_densified = add_local_viewing_geometry(
            stm=regional_ps_stm_densified,
            orbit_config_file=orbit_file,
            orbit_res=orbit_resolution,
            orbit_mode=orbit_mode,
            orbit=f"{satellite}_{direction}_t{track:0>3d}",
        )

        print(
            f"{timestamp()} Projecting regional PS phases and models onto "
            "vertical..."
        )
        regional_ps_stm_densified["unwrapped_phase_pov"] = (
            regional_ps_stm_densified["unwrapped_phase"]
            / np.cos(
                np.radians(
                    regional_ps_stm_densified["local_incidence_angle"]
                )
            )
        )

        for param in model_parameter_layer_names:
            regional_ps_stm_densified[param].data = (
                regional_ps_stm_densified[param].values
            )
            regional_ps_stm_densified[f"{param}_pov"] = (
                regional_ps_stm_densified[param]
                / np.cos(
                    np.radians(
                        regional_ps_stm_densified["local_incidence_angle"]
                    )
                )
            )

        print(f"{timestamp()} Saving regional PS result checkpoint...")
        regional_ps_stm_densified = regional_ps_stm_densified.chunk({
            "time": -1,
            "space": "auto",
        })
        save_zarr_v2(
            regional_ps_stm_densified,
            regional_ps_stm_result_path,
        )
    # AML changes accepted end

    # - Process DT result if available.
    #   DTs are geocoded and get viewing geometry, but STC and scatter filtering
    #   are skipped.
    if dt_stm_densified is not None:
        print(f"{timestamp()} Starting DT geocoding...")
        latlonh = radar_to_latlonh(
            azimuth_coords=dt_stm_densified.azimuth.values,
            range_coords=dt_stm_densified.range.values,
            elevation=dt_stm_densified.pnt_height.values,
            metadata=metadata,
        )

        dt_stm_densified["lat"].data = latlonh[0].flatten()
        dt_stm_densified["lon"].data = latlonh[1].flatten()
        dt_stm_densified["h"].data = latlonh[2].flatten()

        print(f"{timestamp()} Adding DT viewing geometry...")
        dt_stm_densified = add_local_viewing_geometry(
            stm=dt_stm_densified,
            orbit_config_file=orbit_file,
            orbit_res=orbit_resolution,
            orbit_mode=orbit_mode,
            orbit=f"{satellite}_{direction}_t{track:0>3d}",
        )

        print(f"{timestamp()} Projecting DT phases and models onto vertical...")
        dt_stm_densified["unwrapped_phase_pov"] = dt_stm_densified["unwrapped_phase"] / np.cos(
            np.radians(dt_stm_densified["local_incidence_angle"])
        )

        for param in model_parameter_layer_names:
            dt_stm_densified[param].data = dt_stm_densified[param].values
            dt_stm_densified[f"{param}_pov"] = dt_stm_densified[param] / np.cos(
                np.radians(dt_stm_densified["local_incidence_angle"])
            )

        print(f"{timestamp()} Saving DT result checkpoint...")
        dt_stm_densified = make_dt_stm_zarr_safe(dt_stm_densified)
        dt_stm_densified = dt_stm_densified.chunk({"time": -1, "space": "auto"})
        save_zarr_v2(dt_stm_densified, dt_stm_result_path)

    print(f"{timestamp()} Step 5 done.")

else:
    print(f"{timestamp()} Loading Step 5 result checkpoints...")

    ps_stm_densified = xr.open_zarr(ps_stm_result_path)

    # AML changes accepted
    if regional_ps_enabled() and os.path.isdir(regional_ps_stm_result_path):
        regional_ps_stm_densified = xr.open_zarr(
            regional_ps_stm_result_path
        )
    else:
        regional_ps_stm_densified = None
    # AML changes accepted end

    if dt_csv_exists() and os.path.isdir(dt_stm_result_path):
        dt_stm_densified = xr.open_zarr(dt_stm_result_path)
        dt_stm_densified = make_dt_stm_zarr_safe(dt_stm_densified)
    else:
        dt_stm_densified = None
        
# ========================================= #
#              6. Export results            #
# ========================================= #

if start_at_checkpoint < 6:
    print("========== Step 6: export PS, regional PS, and optional DT results ...")

    # - Load Step 5 result products.
    print(f"{timestamp()} Loading PS result checkpoint...")
    ps_stm_densified = xr.open_zarr(ps_stm_result_path)

    # AML changes accepted
    if regional_ps_enabled() and os.path.isdir(regional_ps_stm_result_path):
        print(f"{timestamp()} Loading regional PS result checkpoint...")
        regional_ps_stm_densified = xr.open_zarr(
            regional_ps_stm_result_path
        )
    else:
        regional_ps_stm_densified = None
    # AML changes accepted end

    if dt_csv_exists() and os.path.isdir(dt_stm_result_path):
        print(f"{timestamp()} Loading DT result checkpoint...")
        dt_stm_densified = xr.open_zarr(dt_stm_result_path)
        dt_stm_densified = make_dt_stm_zarr_safe(dt_stm_densified)
    else:
        dt_stm_densified = None

    # <>
    # Reconstruct derived variables when exporting older Step 5 checkpoints.
    # Scientific calculations remain separate from export schema preparation.
    ps_stm_densified = add_displacement_quality_metrics(
        ps_stm_densified,
        wavelength,
    )
    if regional_ps_stm_densified is not None:
        regional_ps_stm_densified = add_displacement_quality_metrics(
            regional_ps_stm_densified,
            wavelength,
        )
    if dt_stm_densified is not None:
        dt_stm_densified = add_displacement_quality_metrics(
            dt_stm_densified,
            wavelength,
        )
    # <|>

    # - Rebuild model parameter layer names.
    model_parameter_layer_names = []
    for model_name in MODEL_NAMES_PARAMS:
        if model_name in model_types:
            model_parameter_layer_names.extend(MODEL_NAMES_PARAMS[model_name])

    # ========================================= #
    #                6.1 PS export              #
    # ========================================= #

    # <>
    reference_lat, reference_lon = get_reference_coordinates(ps_stm_densified)
    print(
        f"{timestamp()} Portal reference point: "
        f"{reference_lat:.8f}, {reference_lon:.8f}."
    )
    # <|>

    print(f"{timestamp()} Exporting PS results...")
    ps_stm_export = prepare_stm_for_standard_exports(
        ps_stm_densified,
        "PS",
    )

    if "csv" in output_types:
        print(f"{timestamp()} Exporting PS CSV...")
        export_to_csv(
            stm=ps_stm_export,
            save_path=ps_csv_save_path,
            model_parameter_layer_names=model_parameter_layer_names,
            ts_proj=csv_ts_proj,
            point_annotation_label=csv_point_annotation_label,
        )

    if "csv_web_portal" in output_types:
        print(f"{timestamp()} Exporting PS portal CSV...")
        export_to_skygeo_portal_with_reference(
            stm=ps_stm_export,
            reference_lat=reference_lat,
            reference_lon=reference_lon,
            save_path=ps_csv_web_save_path,
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
        print(f"{timestamp()} Exporting PS shapefile...")
        export_to_shapefile(
            stm=ps_stm_export,
            save_path=ps_shape_save_path,
            projection=shape_projection,
            model_parameter_layer_names=model_parameter_layer_names,
            point_annotation_label=shape_point_annotation_label,
        )

    if "convex_hull" in output_types:
        print(f"{timestamp()} Exporting PS convex hull...")
        export_convex_hull_to_shapefile(
            stm=ps_stm_densified,
            save_path=ps_chull_save_path,
            projection=chull_projection,
        )

    # AML changes accepted
    # ========================================= #
    #           6.2 Regional PS export          #
    # ========================================= #

    # Regional PS are exported independently, preserving the focus point and
    # the regional filtering choices made in Step 5.
    if regional_ps_stm_densified is not None:
        print(f"{timestamp()} Exporting regional PS results...")
        regional_ps_stm_export = prepare_stm_for_standard_exports(
            regional_ps_stm_densified,
            "Regional PS",
        )

        if "csv" in output_types:
            print(f"{timestamp()} Exporting regional PS CSV...")
            export_to_csv(
                stm=regional_ps_stm_export,
                save_path=regional_ps_csv_save_path,
                model_parameter_layer_names=model_parameter_layer_names,
                ts_proj=csv_ts_proj,
                point_annotation_label=regional_ps_output_label,
            )

        if "csv_web_portal" in output_types:
            print(f"{timestamp()} Exporting regional PS portal CSV...")
            export_to_skygeo_portal_with_reference(
                stm=regional_ps_stm_export,
                reference_lat=reference_lat,
                reference_lon=reference_lon,
                save_path=regional_ps_csv_web_save_path,
                ts_proj=csv_web_ts_proj,
                point_annotation_label=regional_ps_output_label,
                satellite=satellite,
                asc_dsc=direction,
                azimuth_spacing=metadata["azimuth_pixel_spacing"],
                range_spacing=metadata["range_pixel_spacing"],
                model_names=model_types,
                model_parameter_layer_names=model_parameter_layer_names,
            )

        if "shapefile" in output_types:
            print(f"{timestamp()} Exporting regional PS shapefile...")
            export_to_shapefile(
                stm=regional_ps_stm_export,
                save_path=regional_ps_shape_save_path,
                projection=shape_projection,
                model_parameter_layer_names=model_parameter_layer_names,
                point_annotation_label=regional_ps_output_label,
            )

        if "convex_hull" in output_types:
            print(f"{timestamp()} Exporting regional PS convex hull...")
            export_convex_hull_to_shapefile(
                stm=regional_ps_stm_densified,
                save_path=regional_ps_chull_save_path,
                projection=chull_projection,
            )
    else:
        print(f"{timestamp()} No regional PS result available; skipping it.")
    # AML changes accepted end

    # ========================================= #
    #                6.3 DT export              #
    # ========================================= #

    if dt_stm_densified is not None:
        print(f"{timestamp()} Exporting DT results...")
        dt_stm_export = prepare_stm_for_standard_exports(
            dt_stm_densified,
            "DT",
        )

        if "csv" in output_types:
            print(f"{timestamp()} Exporting DT CSV...")
            export_to_csv(
                stm=dt_stm_export,
                save_path=dt_csv_save_path,
                model_parameter_layer_names=model_parameter_layer_names,
                ts_proj=csv_ts_proj,
                point_annotation_label=csv_point_annotation_label,
            )

        if "csv_web_portal" in output_types:
            print(f"{timestamp()} Exporting DT portal CSV...")
            export_to_skygeo_portal_with_reference(
                stm=dt_stm_export,
                reference_lat=reference_lat,
                reference_lon=reference_lon,
                save_path=dt_csv_web_save_path,
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
            print(f"{timestamp()} Exporting DT shapefile...")
            export_to_shapefile(
                stm=dt_stm_export,
                save_path=dt_shape_save_path,
                projection=shape_projection,
                model_parameter_layer_names=model_parameter_layer_names,
                point_annotation_label=shape_point_annotation_label,
            )

        if "convex_hull" in output_types:
            print(f"{timestamp()} Exporting DT convex hull...")
            export_convex_hull_to_shapefile(
                stm=dt_stm_densified,
                save_path=dt_chull_save_path,
                projection=chull_projection,
            )
    else:
        print(f"{timestamp()} No DT result available; skipping DT export.")

    print(f"{timestamp()} Step 6 export finished.")


# AML changes accepted
# ========================================= #
#   7. DT unwrapped phase + amplitude plots #
# ========================================= #

# Regional PS have already been saved and exported as their own product.
# This final diagnostic remains DT-specific because existing_flag describes
# the individual DT existence periods shown in these plots.
# - Create one subplot for every available DT.
#   The non-interactive backend saves the figure without requiring X11, which
#   keeps this diagnostic compatible with batch and Slurm execution.
if dt_stm_densified is not None and dt_stm_densified.sizes.get("space", 0) > 0:
    print("========== Step 7: DT diagnostic time-series plots ...")
    print(f"{timestamp()} Loading DT plotting variables into memory...")

    if dt_plot_quantity == "phase":
        series_label = "Unwrapped phase"
        series_axis_label = "Unwrapped phase [rad]"
        figure_title = (
            "Designated target densified unwrapped phase and amplitude "
            "time series"
        )
    elif dt_plot_quantity == "displacement":
        series_label = "LOS displacement"
        series_axis_label = "LOS displacement [mm]"
        figure_title = (
            "Designated target densified LOS displacement and amplitude "
            "time series"
        )
    else:
        raise ValueError(
            "dt_plot_quantity must be either 'phase' or 'displacement'."
        )

    dt_plot = dt_stm_densified.compute()
    times = dt_plot["time"].values
    targets = dt_plot["target"].values.astype(str)
    n_targets = dt_plot.sizes["space"]

    # - Arrange DT panels in at most three columns.
    ncols = min(3, n_targets)
    nrows = math.ceil(n_targets / ncols)

    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(16, 3.9 * nrows),
        sharex=True,
        squeeze=False,
    )
    axes = axes.ravel()

    legend_handles = {}

    # - Plot each DT with the requested phase/displacement quantity on the left
    #   axis and amplitude on the right. Orange amplitude remains behind blue.
    for i, target_name in enumerate(targets):
        ax_phase = axes[i]
        ax_amp = ax_phase.twinx()

        unwrapped_phase = dt_plot["unwrapped_phase"].isel(space=i).values
        if dt_plot_quantity == "phase":
            plot_values = unwrapped_phase
        else:
            plot_values = (
                -wavelength / (4.0 * np.pi) * unwrapped_phase * 1000.0
            )
        amplitude = dt_plot["amplitude"].isel(space=i).values
        existing = dt_plot["existing_flag"].isel(space=i).values == 1
        not_existing = ~existing

        ax_amp.plot(
            times,
            amplitude,
            color="tab:orange",
            linewidth=0.9,
            alpha=0.7,
            label="Amplitude",
            zorder=1,
        )
        ax_amp.scatter(
            times[existing],
            amplitude[existing],
            color="tab:orange",
            s=10,
            alpha=0.8,
            label="Amplitude, existing_flag = 1",
            zorder=2,
        )

        ax_phase.plot(
            times,
            plot_values,
            color="tab:blue",
            linewidth=1.0,
            alpha=0.9,
            label=series_label,
            zorder=4,
        )
        ax_phase.scatter(
            times[existing],
            plot_values[existing],
            color="tab:blue",
            s=14,
            alpha=0.9,
            label=f"{series_label}, existing_flag = 1",
            zorder=5,
        )

        # - Mark epochs where the DT does not exist at the bottom of the phase
        #   axis. Densified phase values at these epochs are expected to be NaN.
        if np.any(not_existing):
            ymin, ymax = ax_phase.get_ylim()
            marker_y = ymin + 0.03 * (ymax - ymin)
            ax_phase.scatter(
                times[not_existing],
                np.full(np.count_nonzero(not_existing), marker_y),
                marker="|",
                color="0.55",
                s=18,
                alpha=0.6,
                label="existing_flag = 0",
                zorder=3,
            )

        # - Keep the phase/displacement axis blue and amplitude only orange.
        ax_phase.set_ylabel(series_axis_label, color="tab:blue")
        ax_phase.tick_params(axis="y", colors="tab:blue")
        ax_phase.spines["left"].set_color("tab:blue")

        ax_amp.set_ylabel("Amplitude", color="tab:orange")
        ax_amp.tick_params(axis="y", colors="tab:orange")
        ax_amp.spines["right"].set_color("tab:orange")

        ax_phase.set_zorder(ax_amp.get_zorder() + 1)
        ax_phase.patch.set_visible(False)

        # - Add radar-coordinate, coherence, and existence information.
        azimuth = int(dt_plot["azimuth"].isel(space=i).item())
        range_coordinate = int(dt_plot["range"].isel(space=i).item())
        subpixel_azimuth = float(
            dt_plot["azimuth_subpixel"].isel(space=i).item()
        )
        subpixel_range = float(
            dt_plot["range_subpixel"].isel(space=i).item()
        )
        local_temp_coh = float(
            dt_plot["local_temp_coh"].isel(space=i).item()
        )

        n_existing = int(existing.sum())
        n_total = existing.size

        ax_phase.set_title(
            f"{target_name}\n"
            f"SLC=({azimuth}, {range_coordinate}), "
            f"DT=({subpixel_azimuth:.2f}, {subpixel_range:.2f}), "
            f"coh={local_temp_coh:.2f}\n"
            f"existing={n_existing}/{n_total}"
        )
        ax_phase.grid(alpha=0.25)
        ax_phase.set_xlabel("Time")

        # - Collect each legend entry once for a common figure legend.
        for axis in (ax_phase, ax_amp):
            handles, labels = axis.get_legend_handles_labels()
            for handle, label in zip(handles, labels, strict=True):
                if label not in legend_handles:
                    legend_handles[label] = handle

    # - Remove panels not used by a DT.
    for j in range(n_targets, len(axes)):
        fig.delaxes(axes[j])

    fig.legend(
        legend_handles.values(),
        legend_handles.keys(),
        loc="upper center",
        ncol=3,
        bbox_to_anchor=(0.5, 1.02),
    )
    fig.suptitle(figure_title, fontsize=16, y=1.06)

    fig.tight_layout()
    fig.savefig(
        dt_timeseries_plot_path,
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(fig)

    print(f"{timestamp()} Saved DT diagnostic plot to {dt_timeseries_plot_path}")
else:
    print(f"{timestamp()} No DT result available; skipping DT diagnostic plots.")
# AML changes accepted end
