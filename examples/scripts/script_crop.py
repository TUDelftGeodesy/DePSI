import sys

import xarray as xr

sys.path.append('/project/caroline/Share/users/caroline-admin/crop_groningen_zarr_stack/DePSI_group/depsi/')

from depsi.utils import crop_slc_spacetime

# The path where the SLCs are saved
slc_path = '/project/caroline/Share/zarr_stacks/nl_groningen_s1_dsc_t037.zarr'

# The path where the cropped SLC is to be saved
slc_save_path = \
    '/project/caroline/Share/users/caroline-admin/crop_groningen_zarr_stack/nl_groningen_s1_dsc_t037_cropped.zarr'

# Crop in space
aoi_file = '/project/caroline/Share/users/caroline-alapadat/projects/exampleGroningenCrop_desc_t037/haren_igrs_CROP.shp'

# Load the SLCs from ZARR
stack = xr.open_zarr(slc_path)

# Crop the SLCs in space
cropped_stack = crop_slc_spacetime(stack,
                                   aoi_filename=aoi_file)


cropped_stack.to_zarr(slc_save_path, mode="w")
