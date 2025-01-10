import os
from datetime import datetime
import numpy as np
import h5py
import fiona
import shapely.geometry as sg
import matplotlib as mpl
from scipy import stats
import multiprocessing as mp
import itertools
import xarray as xr
import numba
    

def ds_selection(settings, slc_stack, stack_meta):
    '''
    Function that reads slc stack, identifies id pixel, performs multilooking, 
    and estimate equivalent single mother phase.

    Parameters:
        - settings   : json file containing required variables
        - slc_stack  : Xarray contains slc data
        - stack_meta : stack information of the corresponding track
    
    Returns:
        - ds_stm     : space-time-matrix of distributed scatterers
    '''
    ## Assign id pixel
    print('Check whether radar pixel_id has been assigned according to parcels. If not, assign pixel_id first before parcel phase estimation.')
    filename = os.path.join(settings['phase_est_dir'], 'id_pixel_'+stack_meta['stack_id']+'.h5')

    if os.path.isfile(filename):
        print('Radar pixel_id has been assigned, load pixel_id from {} ...'
              .format('id_pixel_'+stack_meta['stack_id']+'.h5'))
        with h5py.File(filename, 'r') as f:
            pixel_id  = f['pixel_id'][()]
        ds_stm = xr.open_zarr(
            os.path.join(settings['stm_dir'], 'ds_stm_' + stack_meta['stack_id'] + '.zarr')
        )

    else:
        print('Assigning radar pixel_id to the corresponding parcel ...')
        pixel_id, ds_stm = assign_id_pixel(settings, slc_stack, stack_meta)
    
    ## Equivalent single mother phase estimation
    print('Check if esm phase has been estimated.')
    filename = os.path.join(settings['phase_est_dir'], 'stack_data_'+stack_meta['stack_id']+'.h5')

    if os.path.isfile(filename):
        print('Load stack_data from {} ...'
              .format('stack_data_'+stack_meta['stack_id']+'.h5'))
        with h5py.File(filename, 'r') as f:
            parcel_id = f['parcel_id'][()]
            cpx_coh   = f['cpx_coh'][()]
        ds_stm = xr.open_zarr(
            os.path.join(settings['stm_dir'], 'ds_stm_' + stack_meta['stack_id'] + '.zarr')
        )

    else:
        print('Multilooking and ESM phase estimation ...')
        ds_stm = multilooking(settings, slc_stack, stack_meta, pixel_id, ds_stm)

    return ds_stm


# @numba.jit(nogil=True)
def assign_id_pixel(settings, slc_stack, stack_meta):
    '''
    Function that performs point in polygon test,
    assigns parcel id to each radar coords, and 
    computes centroid of each parcel.

    Parameters:
        - settings   : json file containing required variables
        - stack_meta : stack information of the corresponding track
    
    Returns:
        - pixel_id : radar pixel containing parcel_id
        - ds_stm   : contains centroid of each parcel, pnt_class, and pnt_id
    '''
    ## Load variables
    nlines  = stack_meta['nlines']
    npixels = stack_meta['npixels']
    az      = slc_stack['azimuth'].values
    rg      = slc_stack['range'].values
    lon     = slc_stack['lon'].values
    lat     = slc_stack['lat'].values
    pts     = np.vstack((lon.flatten(), lat.flatten())).T

    ## Allocate empty variable
    pixel_id    = np.empty(shape=(nlines*npixels))
    pixel_id[:] = np.nan
    parcel_id   = []
    centroid    = []
    az_centroid = []
    rg_centroid = []

    ## Assign parcel_id to each pixel
    f_aux = settings['parcel_shapefile']
    with fiona.open(f_aux) as src:
        features = [feature for feature in src]
    print('Original number of parcels: {}'.format(len(features)))

    for j, feature in enumerate(features):
        s = sg.shape(feature.geometry)
        if s.is_valid and s.is_simple:
            if s.geom_type == 'Polygon':
                coords = list(s.exterior.coords)
            elif s.geom_type == 'MultiPolygon':
                item   = s.geoms[0]
                coords = list(item.exterior.coords)
            r      = sg.LinearRing(coords)
            poly   = mpl.path.Path(coords)
            mask   = poly.contains_points(pts)
            mask2d = mask.reshape(nlines, npixels)

        if np.count_nonzero(mask) >= settings['ds_min_cells']:
            mask_id           = np.argwhere(mask == True)
            pixel_id[mask_id] = int(feature['properties']['id'])
            center            = sg.Polygon(r).centroid
            centroid.append([center.coords.xy[0][0], center.coords.xy[1][0]])
            parcel_id.append(int(feature['properties']['id']))

            mask2d_id = np.argwhere(mask2d == True)
            az_centroid.append(np.median(az[np.unique(mask2d_id[:,0])]).astype(int))
            rg_centroid.append(np.median(rg[np.unique(mask2d_id[:,1])]).astype(int))
        
        ## Counter
        if np.mod(j, round(len(features)/(len(features)/15))) == 0:
            print('Done {} polygons ({} %)'.format(j, round(j/len(features)*100, 2)))
    
    pixel_id  = pixel_id.reshape(nlines,npixels)
    parcel_id = np.array(parcel_id)
    centroid  = np.array(centroid)

    print('Create ds_stm Xarray and save as a zarr file ...')
    ds_stm = xr.Dataset(
        data_vars = dict(
            pnt_class = (["space"], np.repeat(3, parcel_id.size)),
            pnt_id = (["space"], parcel_id),
        ),
        coords = dict(
            space = (["space"], np.arange(parcel_id.size)),
            lat = (["space"], centroid[:,1]),
            lon = (["space"], centroid[:,0]),
            azimuth = (["space"], az_centroid),
            range = (["space"], rg_centroid),
        )
    )
    fileout = os.path.join(settings['stm_dir'], 'ds_stm_' + stack_meta['stack_id'] + '.zarr')
    ds_stm.to_zarr(fileout)

    print('Saving pixel_id into an HDF file ...')
    export_to_hdf(dataset_name = ['pixel_id'], 
                  dataset      = [pixel_id], 
                  out_dir      = settings['phase_est_dir'], 
                  filename     = 'id_pixel_' + stack_meta['stack_id'])
    
    print('Current number of parcels: {}'.format(len(parcel_id)))

    return pixel_id, ds_stm


def export_to_hdf(dataset_name, dataset, out_dir, filename):
    '''
    Funtion that exports single or multiple dataset into HDF file.

    Parameters:
        - dataset_name  : list of dataset name (string)
        - dataset       : list of dataset to be stored
        - out_dir       : output directory (listed in settings)
        - filename      : output filename without extention
    
    Returns:
        - none
    '''

    out_file = os.path.join(out_dir, filename+'.h5')
    data_dict  = {}
    for j, dset_name in enumerate(dataset_name):
        data_dict.update({dset_name : dataset[j]})
    with h5py.File(out_file, 'w') as f:
        for dset_name in data_dict:
            f.create_dataset(dset_name,  data=data_dict[dset_name])


def multilooking(settings, slc_stack, stack_meta, pixel_id, ds_stm):
    '''
    Function that estimates equivalent single master phase from a full complex coherence
    using multilooking interferogram based on parcel.

    Parameters:
        - settings   : json file containing required variables
        - slc_stack
        - stack_meta : stack information of the corresponding track
        - pixel_id   : radar pixels containing parcel_id
        - ds_stm     : space-time-matrix of distributed scatterers
    
    Returns:
        - ds_stm     : full complex coherence matrix for all possible interferogram stack
    TODO: should we store the cpx_coh? >> require a large storage
    '''
    ## Load variables
    master_date = stack_meta['master_date']
    slc_dates   = stack_meta['slc_dates']
    
    ## Extract data from the stack
    cpx_stack = slc_stack["complex"].values
    parcel_id = ds_stm["pnt_id"].values

    ## Allocate empty variable
    nlooks    = [None] * parcel_id.size
    mean_amp  = [None] * parcel_id.size
    amp_disp  = [None] * parcel_id.size
    ml_ifg    = [None] * parcel_id.size
    cpx_coh   = [None] * parcel_id.size
    dc_coh    = [None] * parcel_id.size
    esm_phase = [None] * parcel_id.size

    ## Iterate the parcels
    for i, pid in enumerate(parcel_id):
        ## Select pixel stack inside the current parcel
        sel     = np.where(pixel_id == pid)
        cpx_sel = cpx_stack[sel]

        ## Remove nan rows from cpx_sel (indicating unwanted pixel, e.g. water)
        cpx_sel  = cpx_sel[~np.isnan(cpx_sel).all(axis=1)]

        ## Perform brotherhood selection to each parcel
        if settings['ds_shp_test'] == 'yes':
            cpx_sel = shp_test(cpx_sel, 'ks-test')
        nlooks[i]   = cpx_sel.shape[0]

        ## Mean amplitude and amplitude dispersion
        mean_amp[i] = np.nanmean(np.abs(cpx_sel), axis=0)
        amp_disp[i] = np.nanstd(mean_amp[i]) / np.nanmean(mean_amp[i])

        ## Calculate coherence matrix and multilooking
        cpx_sum    = cpx_sel.T @ np.conj(cpx_sel)
        abs_sum    = np.sum(np.power(np.abs(cpx_sel), 2), axis=0)
        abs_sum    = abs_sum.reshape(1, abs_sum.shape[0])
        ml_ifg[i]  = abs_sum / cpx_sel.shape[0]
        denom      = np.sqrt(abs_sum.T @ abs_sum)
        cpx_coh[i] = cpx_sum / denom
        dc_coh[i]  = np.abs(np.diag(cpx_coh[i],1))
        dc_coh[i]  = np.insert(dc_coh[i], 0, 0)

        ## Equivalent single master (ESM) phase estimation
        master_idx = slc_dates.index(master_date)
        if np.isfinite(cpx_coh[i]).all():
            cpx_phases = phase_linking(cpx_coh[i], regularization=0, estimator='emi')
            esm_phase[i] = np.angle(cpx_phases[master_idx])
        else:
            esm_phase[i] = np.full((len(slc_dates),), np.nan)

        ## Counter
        if np.mod(i, round(len(parcel_id)/(len(parcel_id)/15))) == 0:
            print('Done {} polygons ({} %)'.format(i, round(i/len(parcel_id)*100, 2)))
    
    print('Update ds_stm file ...')
    ds_stm["time"] = [datetime.strptime(str(date_int), '%Y%m%d') for date_int in slc_dates]
    ds_stm = ds_stm.assign(
        pnt_nlooks = (["space"], nlooks),
        mean_amp = (["space", "time"], mean_amp),
        pnt_ampdisp = (["space"], amp_disp),
        dc_coh = (["space", "time"], dc_coh),
        esm_phase = (["space", "time"], esm_phase),
    )
    fileout = os.path.join(settings['stm_dir'], 'ds_stm_' + stack_meta['stack_id'] + '.zarr')
    ds_stm.to_zarr(fileout, mode='a')

    print('Saving stack_data into an HDF file ...')
    export_to_hdf(dataset_name = ['parcel_id', 'cpx_coh'],
                  dataset      = [parcel_id, cpx_coh],
                  out_dir      = settings['phase_est_dir'],
                  filename     = 'stack_data_' + stack_meta['stack_id'])

    return ds_stm


def task(x):
    return stats.kstest(x[0], x[1]).pvalue


def shp_test(data, method='ks-test'):
    '''
    Function that performs brotherhood selection of selected pixels in a parcel by
    testing whether pixels come from the same distribution.

    Parameters:
        - data  : selected pixels in a parcel (array: number_of_pixels x epochs)
        - method: testing method, either parametric or non-parametric (default: ks-test)
    
    Returns:
        - data  : final pixel selection
    '''

    if method == 'ks-test' or method == '':
        ## Option 1: K-S Test (Kolmogorov-Smirnov)
        ## One means two samples have the same distribution (null hypothesis)
        ksmat       = np.zeros((data.shape[0] * data.shape[0]))
        # iter        = itertools.product(np.sort(np.abs(data)), np.sort(np.abs(data)))
        # pvals       = np.array([stats.kstest(x[0], x[1]).pvalue for x in iter])
        with mp.Pool() as pool:
            iter        = itertools.product(np.sort(np.abs(data)), np.sort(np.abs(data)))
            pvals       = np.array(pool.map(task, iter))
        idx         = np.argwhere(pvals>0.05)
        if idx.size != 0:
            ksmat[idx] = 1
        ksmat       = ksmat.reshape(data.shape[0], data.shape[0])
        ksmat_sum   = np.sum(ksmat, axis=1)
        idx         = np.argwhere(ksmat_sum == np.max(ksmat_sum))[0,0]
        data        = data[ksmat[idx].astype(bool)]
    
    if method == 'ad-test':
        ## Option 2: A-D Test (Anderson-Darling)
        pass

    if method == 'ampmean-test':
        ## Option 3: Testing the amplitude mean, assuming the Rayleigh distribution
        ## Reference: Jiang et al. (2014) Fast SHP selection for covariance matrix estimation for multitemporal InSAR.
        pass

    if method == 'ampvar-test':
        ## Option 4: Testing the amplitude variance
        pass

    return data


def phase_linking(data, regularization=0, estimator='emi'):
    '''
    Function that provides different phase linking methods for esm phase estimation.

    Input:
        - data      : full complex coherence for one parcel
        - estimator : phase linking estimators (default: emi)
    
    Output:
        - esm phase : esm phase estimation for one parcel
    '''

    ## Spectral regularization
    if regularization == 1:
        beta = 0.5
        data = (1-beta)*data + beta*np.eye(data.shape[0])
    
    ## Phase linking
    if estimator == 'emi' or estimator == '':
        ## Option 1: EMI (Eigendecomposition-based Maximum-likelihood-estimator of Interferometric phase)
        ## Reference: Ansar et al. (2018) Efficient phase estimation in interferogram stacks.

        ## Implementation using Singular Value Decomposition (svd)
        u, s, vh = np.linalg.svd(np.linalg.inv(np.abs(data)) * data)
        cpx_phases = u[:,-1:] * s[-1] @ np.conj(u[:,-1:].T)
    
    if estimator == 'pta':
        pass
    
    return cpx_phases