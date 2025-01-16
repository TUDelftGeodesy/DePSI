# Distributed Scatterers Module

This is a Python implementation of DElft Contextually Aided Distributed scatterer Environment ([DECADE](https://bitbucket.org/grsradartudelft/decade/src/main/)) on MATLAB. To run the module:

```bash
python ds_main.py
```

The current implementation performs the following things:
- load the stack
- do reslc
- parcel based multilooking
- equivalent single master phase estimation

The output will be a space-time matrix (STM) Xarray Dataset with the following attributes:
- space         : idices, lat, lon, azimuth, and range
- time          : datetime format
- pnt_id        : equivalent to parcel_id (attribute id inside parcel shapefile)
- pnt_class     : 2 (class for DS)
- pnt_nlooks    : number of pixels inside the polygon
- pnt_ampdisp   : amplitude dispersion
- mean_amp      : mean amplitude
- dc_coh        : daisy chain coherence
- esm_phase     : equivalent single master phase

In addition to the output, there will be intermediate outputs stored in HDF format, such as:
- id_pixel: 2D array contains parcel_id to each radar pixel
- stack_data: 2D array contains complex coherence of each parcel

## Dependencies

- DORIS v5 processed SLCs
- Parcel shapefile with attributes id, cropcode, soilcode, and knmi_id
- [KNMI](https://www.knmi.nl/nederland-nu/klimatologie/daggegevens) data files with format etmgeg_<knmi_id>.txt for the Netherlands
- [DePSI](https://github.com/TUDelftGeodesy/DePSI)
- [sarxarray](https://github.com/TUDelftGeodesy/sarxarray)
- [fiona](https://pypi.org/project/fiona/)
- [shapely](https://pypi.org/project/shapely/)
- [h5py](https://www.h5py.org)

## Set up input settings

Users need to specify several input parameters desired for the analysis. The settings are strored in a json format.

- settings_filename	: file name of this input parameters
- run_name			: unique name for this run, users can have multiple runs with different names
- proj_dir			: the root directory of a project
- stack_root_dir	: path to the root directory of all stacks (the stack list can cover multiple regions), output from DORIS
- aoi_shapefile		: path to shapefile of the area of interest
- parcel_shapefile	: path to shapefile of parcels, should have attributes id, cropcode, soilcode, and knmi_id
- meteo_file		: path to files of KNMI data (for the Netherlands only)
- stack_prefix		: the name of the stack (e.g. nl_nieuwolda or leave empty if no prefix)
- start_date		: desired start date of time series for analysis
- end_date			: desired end date of time series for analysis
- processor			: type of DORIS processor (caroline or flinsar, by default caroline)
- reslc				: option to recompute SLC (yes or no, by default yes)
- ds_ml_window		: extent to be used for multilooking (polygon or grid, by default polygon)
- ds_min_cells		: minimum pixels inside the extent for multilooking (by default 40)
- ds_shp_test		: statistically homogeneous pixels test (yes or no, by default yes)
- pe_method			: phase linking method (by default emi)
- ps_method			: ps selection method (nad or nmad, by default nad)
- ps_threshold		: ps method threshold (by default 0.3)

## Directory structure

### Stack directory (caroline)
```bash
stacks (stack_root_dir)
├── nl_groningen_s1_dsc_t037 (stack_prefix = nl_groningen)
│   ├── doris_input.xml
│   ├── stackarea_of_interest.xxx
│   ├── stackburst_coverage.xxx
│   ├── stack
│   │   ├── yyyyMMdd (master)
│   │   │   ├── cint.raw
│   │   │   ├── dem_radar.raw
│   │   │   ├── h2ph_srd.raw
│   │   │   ├── lam.raw
│   │   │   ├── master.res
│   │   │   ├── phi.raw
│   │   │   ├── slave_rsmp_reramped.raw
│   │   │   ├── ...
│   │   ├── yyyyMMdd (slave)
│   │   │   ├── cint.raw
│   │   │   ├── h2ph_srd.raw
│   │   │   ├── lam.raw
│   │   │   ├── slave.res
│   │   │   ├── phi.raw
│   │   │   ├── slave_rsmp_reramped.raw
│   │   │   ├── ...
│   │   ├── dir.txt
```

### Stack directory (flinsar)
```bash
stacks (stack_root_dir)
├── nl_groningen_s1_dsc_t037 (stack_prefix = nl_groningen)
│   ├── doris_input.xml
│   ├── dates.txt
│   ├── stackarea_of_interest.xxx
│   ├── stackburst_coverage.xxx
│   ├── yyyyMMdd (master)
│   │   ├── cint_srd.raw
│   │   ├── h2ph.raw
│   │   ├── master.res
│   │   ├── slc_srd.raw
│   ├── yyyyMMdd (slave)
│   │   ├── dem_radar.raw
│   │   ├── lam.raw
│   │   ├── master.res
│   │   ├── phi.raw
│   │   ├── slc_srd.raw
```

### Project directory (example)
```bash
nieuwolda
├── contextual_data
│   ├── aoi
│   │   ├── aoi.shp (aoi_shapefile)
│   │   ├── ...
│   ├── knmi (meteo_file)
│   │   ├── etmgeg_280.txt
│   │   ├── ...
│   ├── parcels
│   │   ├── nieuwolda_attributes.shp (parcel_shapefile)
│   │   ├── ...
├── insar (proj_dir); directories below will be created during the run
│   ├── testpy
│   │   ├── metadata
│   │   │   ├── s1_dsc_t037
│   │   │   │   ├── stack_meta_s1_dsc_t037.json
│   │   ├── phase_estimation
│   │   │   ├── id_pixel_s1_dsc_t037.h5
│   │   │   ├── stack_data_s1_dsc_t037.h5
│   │   ├── stm
│   │   │   ├── ds_stm_s1_dsc_t037.zarr
│   │   │   ├── ps_stm_s1_dsc_t037.zarr
```

## Future implementation

- Add meteorological data.
- 
