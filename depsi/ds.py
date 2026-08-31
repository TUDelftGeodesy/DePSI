import itertools
import multiprocessing as mp
import os
from datetime import datetime

import geopandas as gpd
import matplotlib as mpl
import numpy as np
import pandas as pd
import shapely.geometry as sg
import xarray as xr
from scipy import stats
from scipy.spatial import cKDTree

from depsi.densification import _determine_dens_connections

DS_REQUIRED_ATTR_KEYS = ("stack_id", "wavelength", "prf", "az_bw", "r_fs", "r_bw")


def assign_parcel_id(
    slc_stack: xr.Dataset, path_to_shapefile: str, ds_min_cells: int, ps_stm: xr.Dataset
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Assign parcel id to radar pixels.

    This function performs point in polygon test,
    assigns parcel id to each radar coords, and
    computes centroid of each parcel.

    Parameters
    ----------
    slc_stack : xr.Dataset
        SLC stack with three variables: complex, amplitude, phase
        and two coordinates: space (lon, lat, azimuth, range) and time.
    path_to_shapefile : str
        Path to parcel shapefile with attributes id/parcel_id, cropcode, soilcode, knmi_id.
    ds_min_cells : int
        Minimum number of radar pixels inside a parcel polygon.
    ps_stm : xr.Dataset
        STM containing information of selected PS.

    Returns
    -------
    ds_mask : 2D array
        Mask array with the same 2D shape as SLC data containing ID of each parcel.
    ds_stm : xr.Dataset
        STM containing information of DS.
    ds_id : 1D array
        ID of each parcel.
    centroid : 2D array (lon, lat)
        Array containing centroid coordinates of each parcel.
    az_centroid : 1D array
        Array containing azimuth index of centroid of each parcel.
    rg_centroid : 1D array
        Array containing range index of centroid of each parcel.
    """
    ## Load variables
    nlines = slc_stack.sizes["azimuth"]
    npixels = slc_stack.sizes["range"]
    az = slc_stack["azimuth"].values
    rg = slc_stack["range"].values
    lon = slc_stack["lon"].values
    lat = slc_stack["lat"].values
    pts = np.vstack((lon.flatten(), lat.flatten())).T

    ## Allocate empty variable
    ds_mask = np.full(nlines * npixels, np.nan)
    ds_id = []
    centroid = []
    az_centroid = []
    rg_centroid = []

    ## Assign parcel_id to each pixel
    features = gpd.read_file(path_to_shapefile)
    print(f"Original number of parcels: {len(features)}")

    for j, feature in features.iterrows():
        s = sg.shape(feature.geometry)
        if s.is_valid and s.is_simple:
            if s.geom_type == "Polygon":
                coords = list(s.exterior.coords)
            elif s.geom_type == "MultiPolygon":
                item = s.geoms[0]
                coords = list(item.exterior.coords)
            r = sg.LinearRing(coords)
            poly = mpl.path.Path(coords)
            mask = poly.contains_points(pts)
            mask2d = mask.reshape(nlines, npixels)

        if np.count_nonzero(mask) >= ds_min_cells:
            mask_id = np.argwhere(mask)
            if "id" in feature.keys():
                ds_mask[mask_id] = int(feature["id"])
                ds_id.append(int(feature["id"]))
            else:
                raise ValueError("Parcel shapefile must contain 'id' attribute.")
            center = sg.Polygon(r).centroid
            centroid.append([center.coords.xy[0][0], center.coords.xy[1][0]])

            mask2d_id = np.argwhere(mask2d)
            az_centroid.append(np.median(az[np.unique(mask2d_id[:, 0])]).astype(int))
            rg_centroid.append(np.median(rg[np.unique(mask2d_id[:, 1])]).astype(int))

        ## Counter
        if np.mod(j, round(len(features) / (len(features) / 10))) == 0:
            print(f"Done {j} polygons ({round(j / len(features) * 100, 2)} %)")

    ds_mask = ds_mask.reshape(nlines, npixels)
    ds_id = np.array(ds_id)
    centroid = np.array(centroid)

    ## Remove PS pixels inside a parcel in ds_mask
    ps_mask = np.zeros(shape=(nlines, npixels))
    ps_az_idx = np.searchsorted(az, ps_stm["azimuth"].values)
    ps_rg_idx = np.searchsorted(rg, ps_stm["range"].values)
    ps_mask[ps_az_idx, ps_rg_idx] = 1
    ds_mask[ps_mask == 1] = np.nan

    return ds_mask, ds_id, centroid, az_centroid, rg_centroid


def ds_phase_estimation(
    slc_stack: xr.Dataset,
    metadata: dict,
    stack_id: str,
    wavelength: float,
    centroid: np.ndarray,
    az_centroid: np.ndarray,
    rg_centroid: np.ndarray,
    parcel_df: pd.DataFrame,
    ds_mask: np.ndarray,
    ds_ids: np.ndarray,
    ds_id_list: list | None = None,
    shp_test: str | None = None,
    min_seg_len: int = 10,
    coh_threshold: float = 0.2,
    mother_epoch: datetime | None = None,
    ds_filepath: str = None,
) -> xr.DataTree:
    """Estimate equivalent single mother (ESM) phase to each parcel.

    This function estimates ESM phase from a full complex coherence matrix
    after multilooking interferograms based on parcel id.

    Parameters
    ----------
    slc_stack : xr.Dataset
        SLC stack with three variables: complex, amplitude, phase
        and two coordinates: space (lon, lat, azimuth, range) and time.
    metadata : dict
        Metadata of the SLC stack.
    stack_id : str
        Unique ID of the SLC stack.
    wavelength : float
        Radar wavelength in m.
    centroid : 2D array (lon, lat)
        Array containing centroid coordinates of each parcel.
    az_centroid : 1D array
        Array containing azimuth index of centroid of each parcel.
    rg_centroid : 1D array
        Array containing range index of centroid of each parcel.
    parcel_df : gpd.GeoDataFrame or pd.DataFrame
        Parcel dataframe containing attributes id, cropcode, soilcode, meteo_id
    ds_mask : 2D array
        Mask array with the same 2D shape as SLC data containing parcel_id.
    ds_ids : np.ndarray
        Array of parcel ids to perform phase estimation.
    ds_id_list : list or None, optional
        List of parcel ids to recompute phase estimation of particular parcels, by default None.
    shp_test : str or None, optional
        Perform statistical homogeneous pixel (shp) test, by default None.
        Current options: None or "ks_test".
    min_seg_len : int, optional
        Minimum segment length to perform phase linking, by default 10.
    coh_threshold : float, optional
        Coherence threshold to perform phase linking, by default 0.2.
    mother_epoch : datetime or None, optional
        Mother epoch to assign the ESM phase estimates, by default None which will take the first epoch.
    ds_filepath : str, optional
        Filepath to check whether the datatree exists, by default None.

    Returns
    -------
    ds_dtree : xr.DataTree
        DataTree containing the STM and coherence matrix datasets.
    """
    assert min_seg_len >= 2, "Minimum segment length must be at least 2."

    ## Initialize or load the datatree
    if ds_filepath is not None and os.path.isdir(ds_filepath):
        ds_dtree = _open_datatree_compat(ds_filepath)
    else:
        ds_dtree = _ds_dt_init(slc_stack, metadata, stack_id, wavelength, ds_ids, centroid, az_centroid, rg_centroid)

    stm = ds_dtree["ds_stm"].to_dataset()
    cpx_coh = ds_dtree["ds_cpx_coh"].to_dataset()

    for key in DS_REQUIRED_ATTR_KEYS:
        if key not in stm.attrs:
            raise ValueError(f"Missing required attribute '{key}' in stm.")

    wavelength = stm.attrs["wavelength"]
    prf = stm.attrs["prf"]
    az_bw = stm.attrs["az_bw"]
    r_fs = stm.attrs["r_fs"]
    r_bw = stm.attrs["r_bw"]

    if ds_id_list is None:
        ds_ids = stm["ds_id"].values
    else:
        ds_ids = ds_id_list

    for ds_id in ds_ids:
        ds_id = int(ds_id)
        print(f"DS #{ds_id}")

        ## Find parcel index in the stm
        mask = (stm["ds_id"] == ds_id).compute()
        idx = (
            mask.where(mask, drop=True)["space"].values[0]
            if mask.where(mask, drop=True)["space"].values.size > 0
            else None
        )
        if idx is None:
            print("Warning: not found in the STM dataset.")
            continue

        ## Coherence matrix estimation
        ds_npixels, ds_mean_ia, ds_mean_p, ds_mean_amp, ds_cpx_coh = _coherence_matrix(
            slc_stack, ds_mask, ds_id, shp_test
        )
        cpx_coh["ds_cpx_coh"].values[idx] = ds_cpx_coh

        if np.isnan(ds_cpx_coh).all():
            print("Warning: all coherence values are NaN.")
            continue
        else:
            ds_mean_amp = np.sqrt(ds_mean_p)
            ds_enl = _enl(ds_npixels, prf, az_bw, r_fs, r_bw)
            ds_coh = np.abs(ds_cpx_coh)
            ds_coh_dc = np.insert(np.diag(ds_coh, k=1), 0, np.nan)

            ## Full complex coherence-based phase estimation (ESM)
            cpxopt = _phase_linking(ds_cpx_coh, regularization=0, estimator="emi")
            mother_idx = np.where(slc_stack.time.values == mother_epoch)[0][0] if mother_epoch is not None else 0
            ds_phi_esm_full = np.angle(cpxopt[:, mother_idx])  # Gives form S - M (S * conj(M))

            ## h2ph variables
            h2ph_stack = slc_stack["h2ph"].values
            sel = np.where(ds_mask == ds_id)
            h2ph_mean = np.nanmean(h2ph_stack[sel], axis=0)

            ## z2ph constants
            if ds_mean_ia is None or np.isnan(ds_mean_ia):
                print(
                    "Warning: incidence angle is not available in the stack, so z2ph is set to NaN. \
It will be computed later if orbit file config is provided."
                )
                z2ph = np.nan
            else:
                z2ph = (-4 * np.pi) / (wavelength) * np.cos(np.radians(ds_mean_ia))

            ## Segmentation
            nsegments, blocks_idx = _segmentation(
                ds_coh,
                min_seg_len=min_seg_len,
                threshold=coh_threshold,
            )

            ## Coherent-block-based phase estimation (ESM)
            ntime = stm.sizes["time"]
            ds_segments = np.zeros((ntime,), dtype=int)
            ds_phi_esm_block = np.full((ntime,), np.nan)

            for s in range(nsegments):
                rmin, rmax = blocks_idx[s][0], blocks_idx[s][1]
                ds_segments[rmin + 1 : rmax + 1] = 1
                cpxcoh_block = ds_cpx_coh[rmin : rmax + 1, rmin : rmax + 1]
                cpxopt_block = _phase_linking(cpxcoh_block, regularization=0, estimator="emi")
                ds_phi_esm_block[rmin : rmax + 1] = np.angle(cpxopt_block[0])

            ## Contextual data (meteo_id, crop_id, peil_id, soil_id)
            if "id" in parcel_df.columns:
                crop_id = parcel_df[parcel_df["id"] == ds_id]["cropcode"].values[0]
                peil_id = parcel_df[parcel_df["id"] == ds_id]["peilgebied"].values[0]
                soil_id = parcel_df[parcel_df["id"] == ds_id]["soilcode"].values[0]
                meteo_id = parcel_df[parcel_df["id"] == ds_id]["meteo_id"].values[0]
            else:
                raise ValueError("Parcel attributes must contain either 'id' attribute.")

            ## Fill the values to STM
            stm["ds_npixels"].values[idx] = ds_npixels
            stm["ds_enl"].values[idx] = ds_enl
            stm["ds_mean_p"].values[idx] = ds_mean_p
            stm["ds_mean_amp"].values[idx] = ds_mean_amp
            stm["ds_coh_dc"].values[idx] = ds_coh_dc
            stm["ds_phi_esm_full"].values[idx] = ds_phi_esm_full
            stm["ds_phi_esm_block"].values[idx] = ds_phi_esm_block
            stm["ds_nsegments"].values[idx] = nsegments
            stm["ds_segments"].values[idx] = ds_segments
            stm["h2ph"].values[idx] = h2ph_mean
            stm["local_incidence_angle"].values[idx] = ds_mean_ia
            stm["z2ph"].values[idx] = z2ph
            stm["meteo_id"].values[idx] = meteo_id
            stm["crop_id"].values[idx] = crop_id
            stm["peil_id"].values[idx] = peil_id
            stm["soil_id"].values[idx] = soil_id

    ds_dtree["ds_stm"] = stm
    ds_dtree["ds_cpx_coh"] = cpx_coh
    return ds_dtree


def _ds_dt_init(
    slc_stack: xr.Dataset,
    metadata: dict,
    stack_id: str,
    wavelength: float,
    ds_ids: np.ndarray,
    centroid: np.ndarray,
    az_centroid: np.ndarray,
    rg_centroid: np.ndarray,
) -> xr.DataTree:
    """Initialize xr.DataTree for DS targets.

    This function initializes ds_dt dataset with basic information of each parcel
    and allocate empty variables for ds_phase_estimation.

    Parameters
    ----------
    slc_stack : xr.Dataset
        SLC stack with (at least) three variables: complex, amplitude, phase
    metadata : dict
        Metadata of the SLC stack.
    stack_id : str
        Unique ID of the SLC stack.
    wavelength : float
        Radar wavelength in m.
    ds_ids : 1D array
        List of ID of each parcel.
    centroid : 2D array (lon, lat)
        Array containing centroid coordinates of each parcel.
    az_centroid : 1D array
        Array containing azimuth index of centroid of each parcel.
    rg_centroid : 1D array
        Array containing range index of centroid of each parcel.

    Returns
    -------
    ds_dtree : xr.DataTree
        DataTree containing two nodes for STM and coherence matrix.
    """
    stm = xr.Dataset(
        data_vars=dict(
            pnt_class=(["space"], np.repeat(3, ds_ids.size)),
            ds_id=(["space"], ds_ids),
        ),
        coords=dict(
            space=(["space"], np.arange(ds_ids.size)),
            time=(["time"], slc_stack.time.values),
            lat=(["space"], centroid[:, 1]),
            lon=(["space"], centroid[:, 0]),
            azimuth=(["space"], az_centroid),
            range=(["space"], rg_centroid),
        ),
    )

    stm = stm.assign_attrs(
        stack_id=stack_id,
        wavelength=wavelength,
        prf=metadata["pulse_repetition_frequency"],
        az_bw=metadata["total_azimuth_bandwidth"],
        r_fs=metadata["range_sampling_rate"],
        r_bw=metadata["total_range_bandwidth"],
    )

    ## Allocate empty variables
    nspace = stm.sizes["space"]
    ntime = stm.sizes["time"]

    stm = stm.assign(
        ds_npixels=(["space"], np.full((nspace,), np.nan)),
        ds_enl=(["space"], np.full((nspace,), np.nan)),
        ds_mean_p=(["space", "time"], np.full((nspace, ntime), np.nan)),
        ds_mean_amp=(["space", "time"], np.full((nspace, ntime), np.nan)),
        ds_coh_dc=(["space", "time"], np.full((nspace, ntime), np.nan)),
        ds_phi_esm_full=(["space", "time"], np.full((nspace, ntime), np.nan)),
        ds_phi_esm_block=(["space", "time"], np.full((nspace, ntime), np.nan)),
        ds_nsegments=(["space"], np.zeros((nspace,), dtype=int)),
        ds_segments=(["space", "time"], np.zeros((nspace, ntime), dtype=int)),
        h2ph=(["space", "time"], np.full((nspace, ntime), np.nan)),
        local_incidence_angle=(["space"], np.full((nspace,), np.nan)),
        z2ph=(["space"], np.full((nspace,), np.nan)),
        meteo_id=(["space"], np.full((nspace,), np.nan)),
        crop_id=(["space"], np.full((nspace,), np.nan)),
        peil_id=(["space"], np.full((nspace,), np.nan)),
        soil_id=(["space"], np.full((nspace,), None, dtype=object)),
    )

    cpx_coh = xr.Dataset(
        data_vars=dict(
            ds_id=(["space"], ds_ids),
        ),
        coords=dict(
            space=(["space"], np.arange(ds_ids.size)),
            time1=(["time"], slc_stack.time.values),
            time2=(["time"], slc_stack.time.values),
        ),
    )

    cpx_coh = cpx_coh.assign(
        ds_cpx_coh=(
            ["space", "time1", "time2"],
            np.full((nspace, ntime, ntime), np.nan + 1j * np.nan, dtype=np.complex64),
        ),
    )

    ds_dtree = xr.DataTree.from_dict(
        {
            "ds_stm": stm,
            "ds_cpx_coh": cpx_coh,
        },
    )

    return ds_dtree


def _coherence_matrix(slc_stack, ds_mask, ds_id, shp_test=None):
    """Complex coherence matrix estimation.

    Parameters
    ----------
    slc_stack : xr.Dataset
        SLC stack with three variables: complex, amplitude, phase
    ds_mask : 2D array
        Mask array with the same 2D shape as SLC data containing parcel_id.
    ds_id : int
        ID of the parcel for which to estimate coherence
    shp_test : bool, optional
        Whether to perform shape test on the selected pixels, by default True

    Returns
    -------
    npixels : int
        Number of pixels inside the parcel.
    mean_ia : float
        Mean incidence angle of the selected pixels inside the parcel.
    mean_p : 1D array
        Mean intensity of the selected pixels inside the parcel.
    mean_amp : 1D array
        Mean amplitude of the selected pixels inside the parcel.
    cpx_coh : 2D array
        Complex coherence matrix of the selected pixels inside the parcel.
    """
    ## Extract data from the stack
    cpx_stack = slc_stack["complex"].values
    if "incidence_angle" in slc_stack.data_vars:
        ia = slc_stack["incidence_angle"].values
    else:
        ia = np.nan

    ## Select pixels inside parcel
    sel = np.where(ds_mask == ds_id)
    cpx_sel = cpx_stack[sel]
    ia_sel = ia[sel] if not np.isnan(ia).all() else np.nan

    ## Remove nan rows from cpx_sel (indicating unwanted pixel, e.g. water)
    cpx_sel = cpx_sel[~np.isnan(cpx_sel).all(axis=1)]
    print(f"Selected data size = {cpx_sel.shape}")
    if cpx_sel.size == 0:
        print("Warning: no pixels was selected for this parcel! Parcel data values will become 0 or NaNs.")
    else:
        ## Perform brotherhood selection to each parcel
        if shp_test is not None:
            cpx_sel = _shp_test(cpx_sel, shp_test)

    npixels = cpx_sel.shape[0]

    ## Mean intensity, amplitude, and incidence angle
    mean_p = np.nanmean(np.abs(cpx_sel) ** 2, axis=0) if not np.isnan(cpx_sel).all() else np.nan
    mean_amp = np.nanmean(np.abs(cpx_sel), axis=0) if not np.isnan(cpx_sel).all() else np.nan
    mean_ia = np.nanmean(ia_sel) if not np.isnan(ia).all() else np.nan

    ## Calculate coherence matrix and multilooking
    cpx_sum = cpx_sel.T @ np.conj(cpx_sel)
    abs_sum = np.sum(np.power(np.abs(cpx_sel), 2), axis=0)
    abs_sum = abs_sum.reshape(1, abs_sum.shape[0])
    denom = np.sqrt(abs_sum.T @ abs_sum)
    cpx_coh = cpx_sum / denom

    return npixels, mean_ia, mean_p, mean_amp, cpx_coh


def _kstest(x):
    """Kolmogorov-Smirnov test.

    Parameters
    ----------
    x : float
        Amplitude of each SLC pixel within a parcel.

    Returns
    -------
    float
        P-value.
    """
    return stats.kstest(x[0], x[1]).pvalue


def _shp_test(data, method="ks-test"):
    """Statistical homogeneous pixel (SHP) test.

    This function performs brotherhood selection of selected pixels within the specified extent
    by testing whether pixels come from the same distribution.

    Parameters
    ----------
    data : 1D array
        Selected radar pixels inside parcel polygon.
    method : str, optional
        Method to perform the SHP test, by default "ks-test".

    Returns
    -------
    data : 1D array
        Selected radar pixels inside parcel polygon after removing outliers.

    Raises
    ------
    NotImplementedError
        Other methods that are not yet implemented.
    """
    if method == "ks-test" or method == "":
        ## Option 1: K-S Test (Kolmogorov-Smirnov)
        ## One means two samples have the same distribution (null hypothesis)
        ksmat = np.zeros(data.shape[0] * data.shape[0])
        with mp.Pool() as pool:
            iter = itertools.product(np.sort(np.abs(data)), np.sort(np.abs(data)))
            pvals = np.array(pool.map(_kstest, iter))
        idx = np.argwhere(pvals > 0.05)
        if idx.size != 0:
            ksmat[idx] = 1
        ksmat = ksmat.reshape(data.shape[0], data.shape[0])
        ksmat_sum = np.sum(ksmat, axis=1)
        idx = np.argwhere(ksmat_sum == np.max(ksmat_sum))[0, 0]
        data_sorted_idx = np.argsort(np.abs(data))
        data_sorted = data[data_sorted_idx]
        data = data_sorted[ksmat[idx].astype(bool)]

    else:
        raise NotImplementedError("This module is not yet implemented.")

    return data


def _phase_linking(data, regularization=0, estimator="emi"):
    """Phase linking.

    This function provides different phase linking methods for
    equivalent single mother (ESM) phase estimation.

    Parameters
    ----------
    data : 2D array
        Complex coherence.
    regularization : int, optional
        Scaling (spectral regularization), by default 0
    estimator : str, optional
        Method to perform phase linking a.k.a. ESM phase estimation, by default "emi"
        EMI (Eigendecomposition-based Maximum-likelihood-estimator of Interferometric phase)

    Returns
    -------
    cpx_phase : 2D array
        Complex phase estimates.

    Raises
    ------
    NotImplementedError
        Other methods that are not yet implemented.
    """
    ## Data cleaning
    data = np.nan_to_num(data, nan=0.0, posinf=0.0, neginf=0.0)

    ## Spectral regularization
    if regularization == 1:
        beta = 0.5
        data = (1 - beta) * data + beta * np.eye(data.shape[0])

    ## Phase linking
    if estimator == "emi" or estimator == "":
        ## Option 1: EMI (Eigendecomposition-based Maximum-likelihood-estimator of Interferometric phase)
        ## Reference: Ansar et al. (2018) Efficient phase estimation in interferogram stacks.

        ## Implementation using Singular Value Decomposition (svd)
        u, s, vh = np.linalg.svd(np.linalg.pinv(np.abs(data)) * data)
        cpx_opt = u[:, -1:] * s[-1] @ np.conj(u[:, -1:].T)

    else:
        raise NotImplementedError("This module is not yet implemented.")

    return cpx_opt


def _enl(npixels, prf, az_bw, r_fs, r_bw):
    """Compute the equivalent number of looks.

    Parameters
    ----------
    npixels : int
        Number of multilooking pixels.
    prf : float
        Pulse repetition frequency.
    az_bw : float
        azimuth bandwidth.
    r_fs : float
        Range sampling rate.
    r_bw : float
        Range bandwidth.

    Returns
    -------
    int
        Equivalent number of looks.
    """
    osr = prf / az_bw * r_fs / r_bw
    enl = np.floor(npixels / osr)

    return enl


def _segmentation(data, min_seg_len=10, threshold=0.2):
    """Coherent segment block identification.

    Split the full coherence matrix into separate blocks using daisy chain coherence threshold.

    Parameters
    ----------
    data : 2D array
        The full coherence matrix of a parcel.
    min_seg_len : int, optional
        Minimum segment length, by default 10.
    threshold : float, optional
        Threshold for daisy chain coherence, by default 0.2.

    Returns
    -------
    nsegments: int
        Number of segments in the full time series data.
    blocks_idx: list of tuple
        Block indices [start, stop]. The blocks are square.
    """
    # - Mask nan values with 1e-8
    data = np.nan_to_num(data, nan=1e-8)

    # - Identify cut positions based on the threshold
    coh_dc = np.diag(data, k=1)
    cut_positions = np.where(coh_dc <= threshold)[0]
    blocks_idx = []
    if cut_positions.size == 0:
        blocks_idx.append((0, data.shape[0]))
    else:
        start = 0
        for cp in cut_positions:
            stop = cp
            if (stop + 1 - start) >= min_seg_len:
                blocks_idx.append((start, stop))
            start = cp + 1

    nsegments = len(blocks_idx)

    return nsegments, blocks_idx


def _open_datatree_compat(zarr_path: str) -> xr.DataTree:
    """Open DS datatree with a compatibility fallback for known Zarr reader issues."""
    try:
        return xr.open_datatree(zarr_path, engine="zarr")
    except AttributeError as exc:
        if "read_only" not in str(exc):
            raise

        root_ds = xr.open_zarr(zarr_path, consolidated=False)
        dt_dict = {"/": root_ds}
        for group_name in ("ds_stm", "ds_cpx_coh"):
            if os.path.isdir(os.path.join(zarr_path, group_name)):
                dt_dict[group_name] = xr.open_zarr(zarr_path, group=group_name, consolidated=False)
        return xr.DataTree.from_dict(dt_dict)


def select_common_fop_ref(
    ps_stm_list: list, proj_crs: int, dist_ub: float = 20.0, quality_var: str = "full_ts_nad"
) -> tuple[np.ndarray, np.ndarray]:
    """Select a common reference point and common points from the first order points across tracks.

    Parameters
    ----------
    ps_stm_list : list
        List of PS STM datasets from different tracks.
    proj_crs : int
        The EPSG code for the coordinate reference system used for projecting the geographic coordinates.
        The PS STM datasets must contain projected coordinates in this CRS as
        "x_euclidean_proj_epsg{proj_crs}" and "y_euclidean_proj_epsg{proj_crs}".
    dist_ub : float, optional
        The upper bound for the distance between PS points across tracks, by default 20.0.
        Points with distance below this value are considerd common points across tracks.
    quality_var : str, optional
        The name of the quality variable in the PS STM datasets to consider for selecting the reference point,
        by default "full_ts_nad".

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        Tuple containing the selected common reference point indices,
        and the indices of the common points in each PS STM dataset.
    """
    for i in range(len(ps_stm_list)):
        ps_stm = ps_stm_list[i]
        if (
            f"x_euclidean_proj_epsg{proj_crs}" not in ps_stm.coords
            or f"y_euclidean_proj_epsg{proj_crs}" not in ps_stm.coords
        ):
            raise ValueError(
                f"PS STM dataset does not contain projected coordinates in EPSG:{proj_crs}. "
                f"Expected variables: 'x_euclidean_proj_epsg{proj_crs}' and 'y_euclidean_proj_epsg{proj_crs}'."
            )

    # - Store the original point indices before filtering out invalid points (NaN lon/lat).
    # - This is necessary to maintain the correct reference point indices after filtering
    # - since we do not save the filtered PS stms.
    pnt_idxs_ori = []
    for i in range(len(ps_stm_list)):
        ps_stm = ps_stm_list[i]
        mask = np.isfinite(ps_stm.lon) & np.isfinite(ps_stm.lat)
        mask = mask.compute()
        pnt_idxs_ori.append(np.nonzero(mask.values)[0])
        stm_clean = ps_stm.isel(space=mask)
        ps_stm_list[i] = stm_clean

    ps_locs = []
    for i in range(len(ps_stm_list)):
        ps_x = ps_stm_list[i][f"x_euclidean_proj_epsg{proj_crs}"].values
        ps_y = ps_stm_list[i][f"y_euclidean_proj_epsg{proj_crs}"].values
        ps_loc = np.stack((ps_x, ps_y), axis=-1)
        ps_locs.append(ps_loc)

    common_locs = ps_locs[0]
    for i in range(len(ps_stm_list)):
        ps_loc = ps_locs[i]
        tree = cKDTree(common_locs)
        distances, idxs = tree.query(ps_loc, distance_upper_bound=dist_ub)
        valid = np.isfinite(distances) & (idxs < common_locs.shape[0])
        common_locs = common_locs[idxs[valid]]
        if common_locs.shape[0] == 0:
            raise ValueError("No common PS was found! Adjust tolerance.")

    pnt_idx_candidates = np.zeros((len(ps_stm_list), common_locs.shape[0]), dtype=int)
    for i in range(len(ps_stm_list)):
        ps_loc = ps_locs[i]
        tree = cKDTree(ps_loc)
        distances, idxs = tree.query(common_locs, distance_upper_bound=dist_ub)
        valid = np.isfinite(distances)
        pnt_idx_candidates[i] = idxs[valid]

    pnt_q_sum = np.zeros((common_locs.shape[0],))
    for i in range(len(ps_stm_list)):
        pnt_q = ps_stm_list[i][quality_var].values[pnt_idx_candidates[i, :]]
        pnt_q_sum += pnt_q
    min_q_idx = np.argmin(pnt_q_sum)
    ref_idxs = pnt_idx_candidates[:, min_q_idx]

    ref_idxs_ori = []
    pnt_idx_candidates_ori = []
    for pnt_idx_ori, ref_idx, pnt_idx_candidate in zip(pnt_idxs_ori, ref_idxs, pnt_idx_candidates, strict=True):
        ref_idxs_ori.append(pnt_idx_ori[ref_idx])
        pnt_idx_candidates_ori.append(pnt_idx_ori[pnt_idx_candidate])

    return ref_idxs_ori, pnt_idx_candidates_ori


def ps_ds_arc(
    ps_stm: xr.Dataset,
    ds_dtree: xr.DataTree,
    n_connections: int = 1,
    key_xcoord: str = "azimuth",
    key_ycoord: str = "range",
    key_h2ph: str = "h2ph",
    key_sdphase_ps: str = "sd_phase",
    key_sdphase_ds: str = "sd_phase",
) -> xr.DataTree:
    """Arc PS and DS targets based on the nearest neighbor search.

    Parameters
    ----------
    ps_stm : xr.Dataset
        STM dataset of first order points.
    ds_dtree : xr.DataTree
        DataTree containing the DS targets.
    n_connections : int, optional
        Number of connections, by default 1.
    key_xcoord : str, optional
        Key name for the x-coordinate in the STM datasets, by default "azimuth".
    key_ycoord : str, optional
        Key name for the y-coordinate in the STM datasets, by default "range".
    key_h2ph : str, optional
        Key name for the height-to-phase variable in the STM datasets, by default "h2ph".
    key_sdphase_ps : str, optional
        Key name for the single-difference of phase variable in the PS STM datasets, by default "sd_phase".
    key_sdphase_ds : str, optional
        Key name for the single-difference of phase variable in the DS STM datasets, by default "sd_phase".

    Returns
    -------
    xr.DataTree
        DataTree containing the PS and DS targets and their connections.
    """
    assert n_connections == 1, "Currently only n_connections=1 is supported."

    ds_stm = ds_dtree["ds_stm"].to_dataset()

    if key_h2ph not in ps_stm.data_vars or key_h2ph not in ds_stm.data_vars:
        raise ValueError(f"Missing required variable '{key_h2ph}' in the PS or DS STM datasets.")
    if key_sdphase_ps not in ps_stm.data_vars:
        raise ValueError(f"Missing required variable '{key_sdphase_ps}' in the PS STM datasets.")
    if key_sdphase_ds not in ds_stm.data_vars:
        raise ValueError(f"Missing required variable '{key_sdphase_ds}' in the DS STM datasets.")

    # - Query densification connections
    idx_ds_pnts, idx_ps_pnts = _determine_dens_connections(ds_stm, ps_stm, n_connections, key_xcoord, key_ycoord)

    # - Take the mean for arc h2ph
    h2ph = (ds_stm[key_h2ph].isel(space=idx_ds_pnts).data + ps_stm[key_h2ph].isel(space=idx_ps_pnts).data) / 2

    # - Double difference phase (dd_phase) calculation = W{target - source}
    dd_phase = (
        ds_stm[key_sdphase_ds].isel(space=idx_ds_pnts).data
        - ps_stm[key_sdphase_ps].isel(space=idx_ps_pnts).data
        + np.pi
    ) % (2 * np.pi) - np.pi

    arc_stm = xr.Dataset(
        coords={
            "idx_dens": (("space",), idx_ds_pnts),
            "idx_network": (("space",), idx_ps_pnts),
        },
        data_vars={
            "h2ph": (("space", "time"), h2ph),
            "dd_phase": (("space", "time"), dd_phase),
        },
        attrs={"wavelength": ds_stm.attrs["wavelength"]},
    )

    ds_cpx_coh = ds_dtree["ds_cpx_coh"].to_dataset()
    arc_dtree = xr.DataTree.from_dict(
        {
            "arc_stm": arc_stm,
            "ps_stm": ps_stm,
            "ds_stm": ds_stm,
            "ds_cpx_coh": ds_cpx_coh,
        },
    )

    return arc_dtree
