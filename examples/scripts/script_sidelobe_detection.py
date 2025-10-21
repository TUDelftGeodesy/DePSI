"""This script runs through the sidelobe detection algorithm. In this example, the sidelobes are discarded.
One can also choose to add mask_sidelobes as an additional layer in the STM."""

import numpy as np
import xarray as xr

from depsi.classification import detect_side_lobes

# The original STM is output from the script script_ps_point_selection.py
stm_original = (
    "/Users/sanvandiepen/PycharmProjects/workingEnvironment2/test_zarr/nl_amsterdam_s1_dsc_t110_stm_base.zarr"
)
stm_v2 = (
    "/Users/sanvandiepen/PycharmProjects/workingEnvironment2/test_zarr/nl_amsterdam_s1_dsc_t110_base_corr.zarr"
)
stm_save_path = (
    "/Users/sanvandiepen/PycharmProjects/workingEnvironment2/test_zarr/nl_amsterdam_s1_dsc_t110_stm_nosl.zarr"
)

max_pixel_dist = 2
min_correlation = 0.90

# first read the STM
stm = xr.open_zarr(stm_original)
raster = xr.open_zarr(stm_v2)
import pdb; pdb.set_trace()

side_lobes_array, _ = detect_side_lobes(stm, max_pixel_dist, min_correlation)

mask_sidelobes = np.ones(len(stm.space), dtype=bool)
mask_sidelobes[side_lobes_array] = False

stm_res = stm.isel(space=mask_sidelobes)

print(f"Removed {len(stm.space)-np.sum(mask_sidelobes)} sidelobes from the STM")

stm.to_zarr(stm_save_path, mode='w')
