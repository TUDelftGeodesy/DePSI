import os, re
from datetime import datetime
import xml.etree.ElementTree as ET
from pathlib import Path
import geopandas as gpd
import json
import numpy as np


def extract_master_date(xml_file):
    """
    Function to extract master date as an integer from doris_input.xml

    Parameters:
        - xml_file  : xml containing the setup for doris input, inc. the master date
    
    Returns:
        - master date as an integer with yyyyMMdd format
    """
    # Parse the XML file
    tree = ET.parse(xml_file)
    root = tree.getroot()
    
    # Find the master_date element
    master_date_element = root.find('.//master_date')
    
    if master_date_element is not None:
        # Extract the date string
        master_date_str = master_date_element.text
        
        # Convert the date string to a datetime object
        master_date_obj = datetime.strptime(master_date_str, '%Y-%m-%d')
        
        # Format the date as an integer in yyyyMMdd format
        master_date_int = int(master_date_obj.strftime('%Y%m%d'))
        
        return master_date_int
    else:
        raise ValueError("master_date element not found in the XML file.")
    

def identify_stacks(settings):
    """
    Function to itendify stack data information and store it to a metadata file

    Parameters:
        - settings : a json file containing initial parameter settings and metadata
    
    Returns:
        - stack_meta : an updated json file containing additional stack information
    """
    ## Detect number of master images / tracks
    dirlist         = os.listdir(settings['stack_root_dir'])
    pattern         = re.compile(rf'^{settings["stack_prefix"]}_?s1_[ad]sc_t\d{{3}}$')
    stack_dirs      = sorted([os.path.join(settings['stack_root_dir'], i) for i in dirlist if pattern.match(i)])
    num_tracks_init = len(stack_dirs)

    print('Found {} DORIS v5 stacks in the specified location'.format(num_tracks_init))

    assert num_tracks_init >= 1, 'Could not find SAR data in the specified directory'

    ## Check if stack intersects the region of interest
    num_tracks = 0
    for i in range(num_tracks_init):
        gdf_stack = gpd.read_file(os.path.join(stack_dirs[i], 'stackburst_coverage.shp'))
        gdf_aoi   = gpd.read_file(settings['aoi_shapefile'])

        ## Ensure aoi shapefiles are in the WGS84 coordinate system
        if gdf_aoi.crs != 'EPSG:4326':
            gdf_aoi = gdf_aoi.to_crs(gdf_stack.crs)
        
        ## Dissolve all features in each shapefile into single geometries
        boundary_stack = gdf_stack.unary_union
        boundary_aoi   = gdf_aoi.unary_union

        ## Check for intersection
        if boundary_stack.intersects(boundary_aoi):
            num_tracks += 1
        else:
            print(f'Discarding stack {stack_id[i]}')
            stack_dirs.remove(stack_dirs[i])
    
    assert num_tracks >= 1, 'Could not find SAR data in the specified directory'

    ## Add parameter to settings and save it
    settings['num_tracks'] = num_tracks
    with open(os.path.join(settings['proj_dir'], settings['settings_filename']), 'w') as file:
        json.dump(settings, file, indent=2)

    ## Extract the metadata of each stack
    ## Read variables
    start_date  = settings['start_date']
    end_date    = settings['end_date']

    stack_meta_list = []
    for i in range(num_tracks):
        master_date = extract_master_date(os.path.join(stack_dirs[i], 'doris_input.xml'))

        assert master_date << start_date or master_date >> end_date, 'Master image is outside specified date range, add it later'

        # for root, dirs, files in os.walk(stack_dirs[i]):
        #     if str(master_date) in dirs:
        #         master_dir = os.path.join(root, str(master_date))
        
        if settings['processor'] == 'flinsar':
            master_dir = os.path.join(stack_dirs[i], str(master_date))
            full_dates = np.loadtxt(os.path.join(stack_dirs[i], 'dates.txt'), dtype=int)
        if settings['processor'] == 'caroline':
            master_dir = os.path.join(stack_dirs[i], 'stack', str(master_date))
            full_dates = np.loadtxt(os.path.join(stack_dirs[i], 'stack', 'dir.txt'), dtype=int)
        sid         = np.argwhere(full_dates >= start_date)[0][0]
        eid         = np.argwhere(full_dates <= end_date)[-1][0]
        slc_dates   = full_dates[sid:eid+1].tolist()
        if master_date < start_date:
            slc_dates.insert(0, master_date)
        if master_date > end_date:
            slc_dates.append(master_date)
        nslcs       = len(slc_dates)

        ## Extract the name of the stack folder
        stack_name = stack_dirs[i].split('/')[-1]
        stack_id = '_'.join(stack_name.split('_')[-3:])

        ## Create list of dates without the mother
        ifg_dates = slc_dates.copy()
        ifg_dates.remove(master_date)
        nifgs     = len(ifg_dates)

        print('Track: {}  Master: {}  nifgs: {}'.format(stack_id, master_date, nifgs))

        ## Create a list of stack paths
        stack_dir  = os.path.dirname(master_dir)
        master_idx = slc_dates.index(master_date)

        if settings['reslc'] == 'yes':
            stack_type = 'cint'
            filelist = sorted([fp for fp in Path(stack_dir).glob('**/cint_srd.raw') 
                           if 'swath' not in str(fp) or 'burst' not in str(fp)])
            filelist = [str(path) for path in filelist]
            slc_paths = [path for path in filelist if extract_date(path) in ifg_dates]
            if settings['processor'] == 'caroline':
                slc_paths.insert(master_idx , os.path.join(master_dir, 'slave_rsmp_reramped.raw'))
            if settings['processor'] == 'flinsar':
                slc_paths.insert(master_idx , os.path.join(master_dir, 'slc_srd.raw'))
            
        if settings['reslc'] == 'no':
            stack_type = 'slc'
            if settings['processor'] == 'caroline':
                filelist = sorted(Path(stack_dir).glob('**/slave_rsmp_reramped.raw'))
            if settings['processor'] == 'flinsar':
                filelist = sorted(Path(stack_dir).glob('**/slc_srd.raw'))
            slc_paths = [path for path in filelist if extract_date(path) in slc_dates]

        ## Read data from the master res file
        if settings['processor'] == 'caroline':
            with open(os.path.join(master_dir, 'master.res'), 'r') as f:
                lines         = f.readlines()
                swath         = lines[37].strip().split()[-1]
                mode          = lines[38].strip().split()[-1]
                r_px_spacing  = float(lines[46].strip().split()[-1])
                az_px_spacing = float(lines[47].strip().split()[-1])
                npixels_res   = int(lines[98].strip().split()[-1])
                nlines_res    = int(lines[99].strip().split()[-1])
        elif settings['processor'] == 'flinsar':
            with open(os.path.join(stack_dirs[i], 'nlines_crp.txt'), 'r') as f:
                lines         = f.readlines()
                nlines_res   = int(lines[0])
            with open(os.path.join(stack_dirs[i], 'npixels_crp.txt'), 'r') as f:
                lines         = f.readlines()
                npixels_res   = int(lines[0])
            swath         = ''
            mode          = ''
            r_px_spacing  = 'n/a'
            az_px_spacing = 'n/a'

        ## Store variables to the metadata file
        stack_meta = {
            'processor'     : settings['processor'],
            'aoi_path'      : settings['aoi_shapefile'],
            'stack_prefix'  : settings['stack_prefix'],
            'meta_dir'      : os.path.join(settings['meta_dir'], stack_id),
            'do_reslc'      : settings['reslc'],
            'master_date'   : master_date,
            'stack_dir'     : stack_dir,
            'master_dir'    : master_dir,
            'slc_paths'     : slc_paths,
            'stack_type'    : stack_type,
            'stack_id'      : stack_id,
            'nslcs'         : nslcs,
            'nifgs'         : nifgs,
            'slc_dates'     : slc_dates,
            'ifg_dates'     : ifg_dates,
            'r_px_spacing'  : r_px_spacing,
            'az_px_spacing' : az_px_spacing,
            'npixels_res'   : npixels_res,
            'nlines_res'    : nlines_res,
        }
        if not os.path.exists(stack_meta['meta_dir']):
            os.makedirs(stack_meta['meta_dir'])
        with open(os.path.join(stack_meta['meta_dir'], 'stack_meta_' + stack_id + '.json'), 'w') as file:
            json.dump(settings, file, indent=2)
        
        stack_meta_list.append(stack_meta)

    return settings, stack_meta_list


def create_processing_folders(settings):
    '''
    Function to create directories if not already there

    Parameters:
        - settings
    
    Returns:
        settings
    '''
    ## Do not modify
    settings['run_dir']        = os.path.join(settings['proj_dir'], settings['run_name'])
    settings['meta_dir']       = os.path.join(settings['run_dir'], 'metadata/')
    settings['phase_est_dir']  = os.path.join(settings['run_dir'], 'phase_estimation/')

    if not os.path.exists(settings['run_dir']):
        os.makedirs(settings['run_dir'])
    
    if not os.path.exists(settings['meta_dir']):
        os.makedirs(settings['meta_dir'])

    if not os.path.exists(settings['phase_est_dir']):
        os.makedirs(settings['phase_est_dir'])
    
    return settings


def extract_date(path):
    match = re.search(r'/(\d{8})/', path)
    return int(match.group(1)) if match else None

