# Distributed Scatterers Module

This is a Python implementation of DElft Contextually Aided Distributed scatterer Environment ([DECADE](https://bitbucket.org/grsradartudelft/decade/src/main/)) on MATLAB. To run the module:

```bash
python ds_main.py
```

The example script performs the following things in order:
1. PS and DS selection
2. Atmospheric phase screen
3. Mixed arc formation (under development)
4. Displacement parameter estimation (under development)

The first implementation will output one `xarray.Dataset` for PS and one `xarray.DataTree` for DS. The DS DataTree output will include two `xarray.Dataset` for the space-time matrix and the full complex coherence matrix. The DS space-time matrix has the following variables:
- space                 : idices, lat, lon, azimuth, and range
- time                  : datetime format
- ds_id                 : equivalent to id of each parcel (attribute id inside parcel shapefile)
- pnt_class             : 3 (class for DS)
- ds_npixels            : number of pixels inside the polygon
- ds_enl                : equivalent number of looks
- ds_mean_p             : mean intensity
- ds_mean_amp           : mean amplitude
- ds_coh_dc             : daisy chain coherence (first off diagonal)
- ds_phi_esm_full       : equivalent single mother phase based on the full coherence matrix
- ds_phi_esm_block      : equivalent single mother phase based on coherent block matrix
- ds_nsegments          : number of coherent segments/blocks
- ds_segments           : start and stop indices in time for each segment
- local_incident_angle  : indicence angle based on parcel centroid
- z2ph                  : conversion from vertical height to line-of-sight phase
- meteo_id              : ID of the closest meteorological stations
- crop_id               : ID of croptype
- peil_id               : ID of peilgebied
- soil_id               : ID of soiltype

Also attributes:
- stack_id              : orbit name
- wavelength            : radar wavelength in meter
- prf                   : pulse repetition frequency
- az_bw                 : azimuth bandwidth
- r_fs                  : range sampling rate
- r_bw                  : range bandwidth

In addition to the output, there will be an intermediate output stored in HDF format, such as:
- `id_pixel_{stack_id}`   : h5 file including 2D array masking ds_id to each radar pixel, id of each parcel, centroid coordinates, azimuth centroid, and range centroid

## Dependencies

- Parcel shapefile with attributes id, cropcode, soilcode, peilgebied, and meteo_id
- [KNMI](https://www.knmi.nl/nederland-nu/klimatologie/daggegevens) data files with format etmgeg_<knmi_id>.txt for the Netherlands
- [DePSI](https://github.com/TUDelftGeodesy/DePSI)

## Set up input settings

Users need to specify several input parameters desired for the analysis. The settings should be specified in the `ds_main.py`.
Follow `depsi_main.py` for PS related settings. For DS settings:

- proj_dir              : project directory path
- run_name              : run name
- log_filename          : filename to save the processing log
- aoi_shp_filepath		: path to shapefile of the area of interest
- parcel_shp_filepath   : path to shapefile of parcels, should have attributes id, cropcode, soilcode, and knmi_id
- meteo_dir		        : path to files of KNMI data (for the Netherlands only)
- stack_root_dir        : directory to available stacks
- stack_prefix		    : the name of the stack (e.g. nl_nieuwolda)
- mission               : satellite mission (options: s1)
- wavelength            : radar wavelength in meter
- stack_ids             : list of stack_id
- metadata_paths        : list of metadata path
- mother_epochs         : list of mother epoch in datetime format
- start_date		    : desired start date of time series for analysis in datetime format
- end_date			    : desired end date of time series for analysis in datetime format
- ds_id_list            : None or list of parcel IDs. Specify only to reestimate esm phase for a subset of parcels
- ds_multilooking_window: boundary for multilooking, current implementation is polygon only
- ds_min_cells			: minimum number of pixels inside the boundary
- ds_shp_test			: option to do statistical homogeneous test
- ds_shp_test_method	: method to do the shp test, current implementation is ks-test only
- ds_coh_threshold		: coherence threshold to specify coherent segment blocks
- ds_min_seg_len		: minimum length of coherent periods
- ds_min_group_size     : minimum number of parcel to form a group
- igrs_codes            : list of IGRS station codes if choose igrs as reference point
- igrs_locs             : list of IGRS locs (lat,lon) if choose igrs as reference point
- ps_stm_save_name      : filename to save PS STM (example: ps_stm_{stack_id}_{checkpoint}.zarr)
- ds_dtree_save_name    : filename to save DS DataTree (example: ds_stm_{stack_id}.zarr)

## Project directory structure

### Example
```bash
nieuwolda
├── contextual_data
│   ├── aoi
│   │   ├── aoi.shp (aoi_shp_filepath)
│   │   ├── ...
│   ├── knmi (meteo_dir)
│   │   ├── etmgeg_280.txt
│   │   ├── ...
│   ├── parcels
│   │   ├── nieuwolda_attributes.shp (parcel_shp_filepath)
│   │   ├── ...
├── [proj_dir] e.g.: insar
│   ├── [run_name] e.g.: testpy
│   │   ├── logs
│   │   │   ├── run_test.txt
│   │   │   ├── ...
│   │   ├── phase_estimation
│   │   │   ├── s1_dsc_t037
│   │   │   │   ├── id_pixel_s1_dsc_t037.h5
│   │   │   │   ├── ...
│   │   │   ├── ...
│   │   │   │   ├── ...
│   │   │   │   ├── ...
│   │   ├── stm
│   │   │   ├── ds_dtree_s1_dsc_t037.zarr
│   │   │   ├── ps_stm_s1_dsc_t037_1sel.zarr
│   │   │   ├── ps_stm_s1_dsc_t037_2atmo.zarr
│   │   │   ├── ...
```
