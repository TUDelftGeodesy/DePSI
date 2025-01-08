import os
import numpy as np
import h5py
import fiona
import shapely.geometry as sg
import matplotlib as mpl
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
    
    ## Multilooking and equivalent single mother phase estimation
    ## In progress ...

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
    lon     = slc_stack['lon'].values
    lat     = slc_stack['lat'].values
    pts     = np.vstack((lon.flatten(), lat.flatten())).T

    ## Allocate empty variable
    pixel_id    = np.empty(shape=(nlines*npixels))
    pixel_id[:] = np.nan
    centroid    = []
    parcel_id   = []

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
                r      = sg.LinearRing(coords)
                poly   = mpl.path.Path(coords)
                mask   = poly.contains_points(pts)
            elif s.geom_type == 'MultiPolygon':
                item   = s.geoms[0]
                coords = list(item.exterior.coords)
                r      = sg.LinearRing(coords)
                poly   = mpl.path.Path(coords)
                mask   = poly.contains_points(pts)

        if np.count_nonzero(mask) >= settings['ds_min_cells']:
            mask_id           = np.argwhere(mask == True)
            pixel_id[mask_id] = int(feature['properties']['id'])
            center            = sg.Polygon(r).centroid
            centroid.append([center.coords.xy[0][0], center.coords.xy[1][0]])
            parcel_id.append(int(feature['properties']['id']))
        
        ## Counter
        if np.mod(j, round(len(features)/(len(features)/15))) == 0:
            print('Done {} polygons ({} %)'.format(j, round(j/len(features)*100, 2)))
    
    parcel_id = np.array(parcel_id)
    centroid  = np.array(centroid)

    print('Create ds_stm Xarray and save as a zarr file ...')
    ds_stm = xr.Dataset(
        data_vars = dict(
            pnt_class = (["space"], np.repeat(3, parcel_id.size)),
            pnt_id = (["space"], parcel_id)
        ),
        coords = dict(
            lat = ("space", centroid[:,1]),
            lon = ("space", centroid[:,0]),
        )
    )
    fileout = os.path.join(settings['stm_dir'], 'ds_stm_' + stack_meta['stack_id'] + '.zarr')
    ds_stm.to_zarr(fileout)

    print('Saving pixel_id into an HDF file ...')
    export_to_hdf(dataset_name = ['pixel_id'], 
                  dataset      = [pixel_id.reshape(nlines,npixels)], 
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