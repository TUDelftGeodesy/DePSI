# Starting Point of DePSI

The DePSI package is designed to work with coregistered interferogram stacks or reconstructed SLC stacks. It expects this stack data to be in the form of `xarray.Dataset`, with specific attributes and dimensions. Below is a minimum example interferogram stack in `xarray.Dataset` that DePSI can work with:

```python
<xarray.Dataset> Size: 2GB
Dimensions:    (azimuth: 2000, range: 4000, time: 10)
Coordinates:
  * azimuth    (azimuth) int64 16kB 0 1 2 3 4 5 ... 1995 1996 1997 1998 1999
  * range      (range) int64 32kB 0 1 2 3 4 5 ... 3994 3995 3996 3997 3998 3999
  * time       (time) datetime64[us] 80B 2018-09-20 2018-09-26 ... 2018-11-13
Data variables:
    complex    (azimuth, range, time) complex64 640MB dask.array<chunksize=(500, 500, 1), meta=np.ndarray>
    amplitude  (azimuth, range, time) float32 320MB dask.array<chunksize=(500, 500, 1), meta=np.ndarray>
    phase      (azimuth, range, time) float32 320MB dask.array<chunksize=(500, 500, 1), meta=np.ndarray>
```

This is an corregistered Sentinel-1 interferogram stack with 10 interferograms, each having dimensions of 2000 (azimuth) x 4000 (range) pixels. DePSI expects the `xr.Dataset` to match the following criteria:
- It should contain dimensions for `azimuth`, `range`, and `time` (number of interferograms or SLCs).
- It should include data variables for `complex`, `amplitude`, and `phase` corresponding to each interferogram or SLC.

These criteria are the bare minimum requirements for DePSI to function correctly. Additional variables are need for certain functionalities in DePSI, such as latitude (`lat`), longitude (`lon`), will be needed for some distance-related calculations. The Height-to-phase conversion factor (`h2ph`) is required for deformation phase to displacement conversion.

## Loading a binary interferogram stack with `sarxarray`

We recommend using the [`sarxarray`](https://tudelftgeodesy.github.io/sarxarray/) package to read binary interferogram stacks into `xarray.Dataset`. This package provides efficient reading and handling of large SAR datasets using `dask` for lazy loading. For exmaple, the `xr.Dataset` above is created using the following code:   

```python
import sarxarray

# Path to the interferogram dataset
path = Path('nl_amsterdam_s1_asc_t088')

# Make a list of interferograms to read
list_ifgs = [p for p in path.rglob('*_cint_srd.raw')]
list_ifgs.sort()

shape = (2000, 4000) # Assume (azimuth, range) shape known beforehand
reading_chunks = (500, 500) # Reading chunks for Dask lazy loading
stack = sarxarray.from_binary(list_ifgs, shape, dtype=np.complex64, chunks=reading_chunks)
```

The file structure of the dataset is as follows:
```sh
nl_amsterdam_s1_asc_t088/
├── 20180920_cint_srd.raw
├── 20180926_cint_srd.raw
├── 20181002_cint_srd.raw
├── 20181008_cint_srd.raw
├── 20181014_cint_srd.raw
├── 20181020_cint_srd.raw
├── 20181026_cint_srd.raw
├── 20181101_cint_srd.raw
├── 20181107_cint_srd.raw
├── 20181113_cint_srd.raw
```
where each interferogram binary file is named with its acquisition date followed by `_cint_srd.raw`. We assumen the metadata such as the shape of the interferograms and data type are known beforehand.