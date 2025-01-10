import os
import sys
import json
import sarxarray
import xarray as xr

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append('/Users/ylumbangaol/Documents/S3/software/caroline_dev/DePSI_group/depsi/')
from ds_utils import identify_stacks, create_processing_folders, load_slc_stack
from classification import ps_selection
from ds import ds_selection

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
PROJ_DIR = '/Users/ylumbangaol/Documents/S3/projects/nieuwolda/insar/'


def main():
    ## Load initial settings from the JSON file
    with open(os.path.join(PROJ_DIR, 'ds_settings.json'), 'r') as file:
        settings = json.load(file)

    ## Create and add directory lists
    settings = create_processing_folders(settings)

    ## Identify stacks and variables
    settings, stack_meta_list = identify_stacks(settings)

    ## Iterate tracks
    for i in range(settings['num_tracks']):
        ## Read stack metadata
        stack_meta = stack_meta_list[i]

        ## Load SLCs
        # slc_stack = sarxarray.from_dataset(slcs)
        slc_stack, stack_meta = load_slc_stack(stack_meta, chunks=(500,500))
        print(f"SLCs stack {stack_meta['stack_id']} was loaded")

        ## PS selection
        ps_stm = ps_selection(slcs=slc_stack,
                              threshold=settings['ps_threshold'],
                              method = settings['ps_method'],
                              output_chunks = 500,
                              )
        fileout = os.path.join(settings['stm_dir'], 'ps_stm_' + stack_meta['stack_id'] + '.zarr')
        ps_stm.to_zarr(fileout, mode='w')

        ## DS selection
        ds_stm = ds_selection(settings, slc_stack, stack_meta)

        ## Merge ps and ds stm
        psds_stm = xr.concat([ps_stm, ds_stm], dim='space')
        print(psds_stm)
        
        print('To be continued ...')


if __name__ == "__main__":
    main()