import os
import sys
import json
import sarxarray

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append('/Users/ylumbangaol/Documents/S3/software/caroline_dev/DePSI_group/depsi/')
from ds_utils import identify_stacks, create_processing_folders
from ds import *

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

    ## ## Interferogram formation
    for i in range(settings['num_tracks']):
        ## Read stack metadata
        stack_meta = stack_meta_list[i]

        ## Load SLCs
        slc_stack = sarxarray.from_dataset(slcs)
        print(f"SLCs stack {stack_meta['stack_id']} was loaded")

        ## PS selection
        # ps_stm = ps_selection(settings, slc_stack, stack_meta)

        ## DS selection
        ds_stm = ds_selection(settings, slc_stack, stack_meta)
        print(f'Stack #{i}')


if __name__ == "__main__":
    main()