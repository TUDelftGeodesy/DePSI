import os
import json
import sarxarray
import xarray as xr
import h5py

from depsi.classification import ps_selection
from depsi.ds_utils import identify_stacks, create_processing_folders, load_slc_stack
from depsi.ds import assign_id_pixel, parcel_phase_estimation, export_to_hdf

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
                              threshold = settings['ps_threshold'],
                              method = settings['ps_method'],
                              output_chunks = 500,
                              )
        fileout = os.path.join(settings['stm_dir'], 'ps_stm_' + stack_meta['stack_id'] + '.zarr')
        ps_stm.to_zarr(fileout, mode='w')

        ## DS selection
        ## Assign id pixel
        print(
            "Check whether radar pixel_id has been assigned according to parcels.\n \
            If not, assign pixel_id first before parcel phase estimation."
        )
        filename = os.path.join(settings["phase_est_dir"], "id_pixel_" + stack_meta["stack_id"] + ".h5")

        if os.path.isfile(filename):
            print("Radar pixel_id has been assigned, load pixel_id from {} ...".format("id_pixel_" + stack_meta["stack_id"] + ".h5"))
            with h5py.File(filename, "r") as f:
                pixel_id = f["pixel_id"][()]
            ds_stm = xr.open_zarr(os.path.join(settings["stm_dir"], "ds_stm_" + stack_meta["stack_id"] + ".zarr"))

        else:
            print("Assigning radar pixel_id to the corresponding parcel ...")
            pixel_id, ds_stm = assign_id_pixel(
                slc_stack,
                nlines=stack_meta["nlines"],
                npixels=stack_meta["npixels"],
                path_to_shapefile=settings["parcel_shapefile"],
                ds_min_cells=settings["ds_min_cells"],
            )

            fileout = os.path.join(settings["stm_dir"], "ds_stm_" + stack_meta["stack_id"] + ".zarr")
            ds_stm.to_zarr(fileout)

            print("Saving pixel_id into an HDF file ...")
            export_to_hdf(
                dataset_name=["pixel_id"],
                dataset=[pixel_id],
                out_dir=settings["phase_est_dir"],
                filename="id_pixel_" + stack_meta["stack_id"],
            )

        ## Parcel phase estimation
        print("Check if esm phase has been estimated.")
        filename = os.path.join(settings["phase_est_dir"], "stack_data_" + stack_meta["stack_id"] + ".h5")

        if os.path.isfile(filename):
            print("ESM phases have been estimated. Load STM DS ...")
            ds_stm = xr.open_zarr(os.path.join(settings["stm_dir"], "ds_stm_" + stack_meta["stack_id"] + ".zarr"))

        else:
            print("Multilooking and ESM phase estimation ...")
            ds_stm, parcel_id, cpx_coh = parcel_phase_estimation(
                slc_stack,
                pixel_id,
                ds_stm,
                slc_dates=stack_meta["slc_dates"],
                mother_date=stack_meta["mother_date"],
                ds_shp_test=settings["ds_shp_test"],
            )

            fileout = os.path.join(settings["stm_dir"], "ds_stm_" + stack_meta["stack_id"] + ".zarr")
            ds_stm.to_zarr(fileout, mode="a")

            print("Saving stack_data into an HDF file ...")
            export_to_hdf(
                dataset_name=["parcel_id", "cpx_coh"],
                dataset=[parcel_id, cpx_coh],
                out_dir=settings["phase_est_dir"],
                filename="stack_data_" + stack_meta["stack_id"],
            )

        ## Merge ps and ds stm
        psds_stm = xr.concat([ps_stm, ds_stm], dim='space')
        print(psds_stm)
        
        print('To be continued ...')


if __name__ == "__main__":
    main()