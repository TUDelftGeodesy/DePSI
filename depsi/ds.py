import itertools
import multiprocessing as mp
import os
from datetime import datetime

import fiona
import h5py
import matplotlib as mpl
import numpy as np
import shapely.geometry as sg
import xarray as xr
from scipy import stats


def ds_selection(
    slc_stack,
    stack_id,
    nlines,
    npixels,
    slc_dates,
    mother_date,
    ds_min_cells,
    ds_shp_test,
    path_to_shapefile,
    path_to_stm,
    path_to_pe,
):
    """Function to read an slc stack, identify id pixel, perform multilooking,
    and estimate equivalent single mother phase.

    Parameters
    ----------
    slc_stack : xr.Dataset
        SLC stack with three variables: (complex, amplitude, phase)
        and two coordinates: space (lat, lon, azimuth, range) and time.
    stack_id : str
        Equals to track id (e.g. s1_asc_t088)
    nlines : int
        Number of lines (azimuth).
    npixels : int
        Number of pixels (range).
    slc_dates : list
        List of integer of slc dates in the format yyyyMMdd.
    mother_date : int
        Mother date in the format yyyyMMdd.
    ds_min_cells : int
        Minimum number of radar pixels inside the specified extent.
    ds_shp_test : str, optional
        Statistical homogeneous pixel, by default "yes".
    path_to_shapefile : str
        Path to parcel shapefile with attributes int_id, cropcode, soilcode, knmi_id.
    path_to_stm : str
        Path to STM directory.
    path_to_pe : str
        Path to parcel phase estimation directory.

    Returns
    -------
    ds_stm : xr.Dataset
        Selected DS as virtial PS in the form of STM with two variables: space and time.
    """  # noqa: D401, D205
    ## Assign id pixel
    print(
        "Check whether radar pixel_id has been assigned according to parcels.\n \
        If not, assign pixel_id first before parcel phase estimation."
    )
    filename = os.path.join(path_to_pe, "id_pixel_" + stack_id + ".h5")

    if os.path.isfile(filename):
        print("Radar pixel_id has been assigned, load pixel_id from {} ...".format("id_pixel_" + stack_id + ".h5"))
        with h5py.File(filename, "r") as f:
            pixel_id = f["pixel_id"][()]
        ds_stm = xr.open_zarr(os.path.join(path_to_stm, "ds_stm_" + stack_id + ".zarr"))

    else:
        print("Assigning radar pixel_id to the corresponding parcel ...")
        pixel_id, ds_stm = assign_id_pixel(
            slc_stack,
            nlines,
            npixels,
            path_to_shapefile,
            ds_min_cells,
        )

        fileout = os.path.join(path_to_stm, "ds_stm_" + stack_id + ".zarr")
        ds_stm.to_zarr(fileout)

        print("Saving pixel_id into an HDF file ...")
        export_to_hdf(
            dataset_name=["pixel_id"],
            dataset=[pixel_id],
            out_dir=path_to_pe,
            filename="id_pixel_" + stack_id,
        )

    ## Parcel phase estimation
    print("Check if esm phase has been estimated.")
    filename = os.path.join(path_to_pe, "stack_data_" + stack_id + ".h5")

    if os.path.isfile(filename):
        print("ESM phases have been estimated. Load STM DS ...")
        ds_stm = xr.open_zarr(os.path.join(path_to_stm, "ds_stm_" + stack_id + ".zarr"))

    else:
        print("Multilooking and ESM phase estimation ...")
        ds_stm, parcel_id, cpx_coh = parcel_phase_estimation(
            slc_stack,
            pixel_id,
            ds_stm,
            slc_dates,
            mother_date,
            ds_shp_test,
        )

        fileout = os.path.join(path_to_stm, "ds_stm_" + stack_id + ".zarr")
        ds_stm.to_zarr(fileout, mode="a")

        print("Saving stack_data into an HDF file ...")
        export_to_hdf(
            dataset_name=["parcel_id", "cpx_coh"],
            dataset=[parcel_id, cpx_coh],
            out_dir=path_to_pe,
            filename="stack_data_" + stack_id,
        )

    return ds_stm


def assign_id_pixel(slc_stack, nlines, npixels, path_to_shapefile, ds_min_cells):  # noqa: D417
    """Function that performs point in polygon test,
    assigns parcel id to each radar coords, and
    computes centroid of each parcel.

    Parameters
    ----------
    slc_stack : xr.Dataset
        SLC stack with three variables: complex, amplitude, phase
        and two coordinates: space (lon, lat, azimuth, range) and time.
    nlines : int
        Number of lines (azimuth).
    npixels : int
        Number of pixels (range).
    path_to_shapefile : str
        Path to parcel shapefile with attributes int_id, cropcode, soilcode, knmi_id.
    ds_min_cells : int
        Minimum number of radar pixels inside a parcel polygon.

    Returns
    -------
    pixel_id : 2D array
        Array with the same 2D shape as SLC data containing parcel_id.
    ds_stm : xr.Dataset
        STM containing information of DS.
    """  # noqa: D205, D401
    ## Load variables
    az = slc_stack["azimuth"].values
    rg = slc_stack["range"].values
    lon = slc_stack["lon"].values
    lat = slc_stack["lat"].values
    pts = np.vstack((lon.flatten(), lat.flatten())).T

    ## Allocate empty variable
    pixel_id = np.empty(shape=(nlines * npixels))
    pixel_id[:] = np.nan
    parcel_id = []
    centroid = []
    az_centroid = []
    rg_centroid = []

    ## Assign parcel_id to each pixel
    with fiona.open(path_to_shapefile) as src:
        features = [feature for feature in src]
    print(f"Original number of parcels: {len(features)}")

    for j, feature in enumerate(features):
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
            mask_id = np.argwhere(mask == True)  # noqa: E712
            pixel_id[mask_id] = int(feature["properties"]["id"])
            center = sg.Polygon(r).centroid
            centroid.append([center.coords.xy[0][0], center.coords.xy[1][0]])
            parcel_id.append(int(feature["properties"]["id"]))

            mask2d_id = np.argwhere(mask2d == True)  # noqa: E712
            az_centroid.append(np.median(az[np.unique(mask2d_id[:, 0])]).astype(int))
            rg_centroid.append(np.median(rg[np.unique(mask2d_id[:, 1])]).astype(int))

        ## Counter
        if np.mod(j, round(len(features) / (len(features) / 10))) == 0:
            print(f"Done {j} polygons ({round(j/len(features)*100, 2)} %)")

    pixel_id = pixel_id.reshape(nlines, npixels)
    parcel_id = np.array(parcel_id)
    centroid = np.array(centroid)

    print(f"Current number of parcels: {len(parcel_id)}")
    print("Create ds_stm Xarray and save as a zarr file ...")
    ds_stm = xr.Dataset(
        data_vars=dict(
            pnt_class=(["space"], np.repeat(3, parcel_id.size)),
            pnt_id=(["space"], parcel_id),
        ),
        coords=dict(
            space=(["space"], np.arange(parcel_id.size)),
            lat=(["space"], centroid[:, 1]),
            lon=(["space"], centroid[:, 0]),
            azimuth=(["space"], az_centroid),
            range=(["space"], rg_centroid),
        ),
    )

    return pixel_id, ds_stm


def export_to_hdf(dataset_name, dataset, out_dir, filename):  # noqa: D417
    """Funtion that exports single or multiple dataset into HDF file.

    Parameters
    ----------
    dataset_name : list
        List of strings containing the name of dataset to be stored.
    dataset : list
        List of arrays containing the data to be stored.
    out_dir : str
        Path to directory.
    filename : str
        Output filename without extention.

    Returns
    -------
    None
    """
    out_file = os.path.join(out_dir, filename + ".h5")
    data_dict = {}
    for j, dset_name in enumerate(dataset_name):
        data_dict.update({dset_name: dataset[j]})
    with h5py.File(out_file, "w") as f:
        for dset_name in data_dict:
            f.create_dataset(dset_name, data=data_dict[dset_name])


def parcel_phase_estimation(slc_stack, pixel_id, ds_stm, slc_dates, mother_date, ds_shp_test):  # noqa: D417
    """Function that estimates equivalent single mother phase from a full complex coherence
    using multilooking interferogram based on parcel.

    Parameters
    ----------
    slc_stack : xr.Dataset
        SLC stack with three variables: complex, amplitude, phase
        and two coordinates: space (lon, lat, azimuth, range) and time.
    pixel_id : 2D array
        Array with the same 2D shape as SLC data containing parcel_id.
    ds_stm : xr.Dataset
        STM containing information of DS.
    slc_dates : list
        List of integer of slc dates in the format yyyyMMdd.
    mother_date : int
        Mother date in the format yyyyMMdd.
    ds_shp_test : str, optional
        Statistical homogeneous pixel, by default "yes".

    Returns
    -------
    ds_stm : xr.Dataset
        STM containing information of DS.
    parcel_id : 1D array
        ID of each parcel.
    cpx_coh : 2D array
        Complex coherence of each parcel.
    NOTE: should we store the cpx_coh?
    """  # noqa: D205, D401
    ## Extract data from the stack
    cpx_stack = slc_stack["complex"].values
    parcel_id = ds_stm["pnt_id"].values

    ## Allocate empty variable
    nlooks = [None] * parcel_id.size
    mean_amp = [None] * parcel_id.size
    amp_disp = [None] * parcel_id.size
    ml_ifg = [None] * parcel_id.size
    cpx_coh = [None] * parcel_id.size
    dc_coh = [None] * parcel_id.size
    esm_phase = [None] * parcel_id.size

    ## Iterate the parcels
    for i, pid in enumerate(parcel_id):
        ## Select pixel stack inside the current parcel
        sel = np.where(pixel_id == pid)
        cpx_sel = cpx_stack[sel]

        ## Remove nan rows from cpx_sel (indicating unwanted pixel, e.g. water)
        cpx_sel = cpx_sel[~np.isnan(cpx_sel).all(axis=1)]

        ## Perform brotherhood selection to each parcel
        if ds_shp_test == "yes":
            cpx_sel = shp_test(cpx_sel, "ks-test")
        nlooks[i] = cpx_sel.shape[0]

        ## Mean amplitude and amplitude dispersion
        mean_amp[i] = np.nanmean(np.abs(cpx_sel), axis=0)
        amp_disp[i] = np.nanstd(mean_amp[i]) / np.nanmean(mean_amp[i])

        ## Calculate coherence matrix and multilooking
        cpx_sum = cpx_sel.T @ np.conj(cpx_sel)
        abs_sum = np.sum(np.power(np.abs(cpx_sel), 2), axis=0)
        abs_sum = abs_sum.reshape(1, abs_sum.shape[0])
        ml_ifg[i] = abs_sum / cpx_sel.shape[0]
        denom = np.sqrt(abs_sum.T @ abs_sum)
        cpx_coh[i] = cpx_sum / denom
        dc_coh[i] = np.abs(np.diag(cpx_coh[i], 1))
        dc_coh[i] = np.insert(dc_coh[i], 0, 0)

        ## Equivalent single mother (ESM) phase estimation
        mother_idx = slc_dates.index(mother_date)
        if np.isfinite(cpx_coh[i]).all():
            cpx_phases = phase_linking(cpx_coh[i], regularization=0, estimator="emi")
            esm_phase[i] = np.angle(cpx_phases[mother_idx])
        else:
            esm_phase[i] = np.full((len(slc_dates),), np.nan)

        ## Counter
        if np.mod(i, round(len(parcel_id) / (len(parcel_id) / 10))) == 0:
            print(f"Done {i} polygons ({round(i/len(parcel_id)*100, 2)} %)")

    print("Update ds_stm file ...")
    ds_stm["time"] = [datetime.strptime(str(date_int), "%Y%m%d") for date_int in slc_dates]
    ds_stm = ds_stm.assign(
        pnt_nlooks=(["space"], nlooks),
        mean_amp=(["space", "time"], mean_amp),
        pnt_ampdisp=(["space"], amp_disp),
        dc_coh=(["space", "time"], dc_coh),
        esm_phase=(["space", "time"], esm_phase),
    )

    return ds_stm, parcel_id, cpx_coh


def kstest(x):  # noqa: D103
    return stats.kstest(x[0], x[1]).pvalue


def shp_test(data, method="ks-test"):  # noqa: D417
    """Function that performs brotherhood selection of selected pixels within the specified extent
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
    """  # noqa: D401, D205
    if method == "ks-test" or method == "":
        ## Option 1: K-S Test (Kolmogorov-Smirnov)
        ## One means two samples have the same distribution (null hypothesis)
        ksmat = np.zeros(data.shape[0] * data.shape[0])
        # iter        = itertools.product(np.sort(np.abs(data)), np.sort(np.abs(data)))
        # pvals       = np.array([stats.kstest(x[0], x[1]).pvalue for x in iter])
        with mp.Pool() as pool:
            iter = itertools.product(np.sort(np.abs(data)), np.sort(np.abs(data)))
            pvals = np.array(pool.map(kstest, iter))
        idx = np.argwhere(pvals > 0.05)
        if idx.size != 0:
            ksmat[idx] = 1
        ksmat = ksmat.reshape(data.shape[0], data.shape[0])
        ksmat_sum = np.sum(ksmat, axis=1)
        idx = np.argwhere(ksmat_sum == np.max(ksmat_sum))[0, 0]
        data = data[ksmat[idx].astype(bool)]

    if method == "ad-test":
        ## Option 2: A-D Test (Anderson-Darling)
        raise NotImplementedError("This module is not yet implemented.")

    if method == "ampmean-test":
        ## Option 3: Testing the amplitude mean, assuming the Rayleigh distribution
        ## Reference: Jiang et al. (2014) Fast SHP selection for covariance matrix estimation for multitemporal InSAR.
        raise NotImplementedError("This module is not yet implemented.")

    if method == "ampvar-test":
        ## Option 4: Testing the amplitude variance
        raise NotImplementedError("This module is not yet implemented.")

    return data


def phase_linking(data, regularization=0, estimator="emi"):
    """Function that provides different phase linking methods for
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
    """  # noqa: D205, D401
    ## Spectral regularization
    if regularization == 1:
        beta = 0.5
        data = (1 - beta) * data + beta * np.eye(data.shape[0])

    ## Phase linking
    if estimator == "emi" or estimator == "":
        ## Option 1: EMI (Eigendecomposition-based Maximum-likelihood-estimator of Interferometric phase)
        ## Reference: Ansar et al. (2018) Efficient phase estimation in interferogram stacks.

        ## Implementation using Singular Value Decomposition (svd)
        u, s, vh = np.linalg.svd(np.linalg.inv(np.abs(data)) * data)
        cpx_phases = u[:, -1:] * s[-1] @ np.conj(u[:, -1:].T)

    if estimator == "pta":
        raise NotImplementedError("This module is not yet implemented.")

    return cpx_phases
