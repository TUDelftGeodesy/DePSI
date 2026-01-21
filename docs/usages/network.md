# Network Estimation

DePSI selects a sparse set of reliable Persistent Scatterer (PS) points to form a network and estimates the network deformation phase with strict Hypothesis Testing (HT) criteria. This step is crucial for building a robust benchmark for estimating deformation time series of all selected PS points.

## Step 1: Select network points

One can use the `depsi.network_stm_selection` method to select network points from the pre-selected PS points based on their spatial distribution and coherence properties. The method sorts the PS points by a customizable quality metric, such as NMAD or NAD, and selects points with a spatial separation greater than a specified minimum distance. An example of using `depsi.network_stm_selection` to select network points is shown below:

```python
stm_network_pnts = depsi.network_stm_selection(
                  stm_preselention_pnts, # results from depsi.ps_selection, with ~20K PS points
                  min_dist=300, # Minimum distance (in meters) between selected network points
                  azimuth_spacing=10, # Azimuth spacing (in meters)
                  range_spacing=5, # Range spacing (in meters)
                  sortby_var="full_ts_nmad", # Sorting variable
                  )
```
```output
<xarray.Dataset> Size: 492kB
Dimensions:              (space: 1431, time: 10)
Coordinates:
  * space                (space) int64 11kB 0 1 2 9 ... 7405 13727 13426 10879
  * time                 (time) datetime64[us] 80B 2018-09-20 ... 2018-11-13
    lat                  (space) float32 6kB 52.24 52.25 52.26 ... 52.42 52.39
    lon                  (space) float32 6kB 4.878 4.95 5.02 ... 4.819 4.804
    azimuth              (space) int64 11kB 0 0 0 2 2 ... 930 1538 1501 1258
    range                (space) int64 11kB 1197 2456 3672 ... 1471 1166 754
    decimal_years        (time) float64 80B 0.0 0.01643 ... 0.1314 0.1478
    local_x              (space) float64 11kB 3.9e+06 3.899e+06 ... 3.887e+06
    local_y              (space) float64 11kB 3.328e+05 3.377e+05 ... 3.267e+05
Data variables:
    complex              (space, time) complex64 114kB dask.array<chunksize=(1431, 10), meta=np.ndarray>
    amplitude            (space, time) float32 57kB dask.array<chunksize=(1431, 10), meta=np.ndarray>
    phase                (space, time) float32 57kB dask.array<chunksize=(1431, 10), meta=np.ndarray>
    h2ph                 (space, time) float64 114kB dask.array<chunksize=(1431, 10), meta=np.ndarray>
    time_selection_nmad  (space) float32 6kB dask.array<chunksize=(1431,), meta=np.ndarray>
    full_ts_nmad         (space) float32 6kB dask.array<chunksize=(1431,), meta=np.ndarray>
    full_ts_nad          (space) float32 6kB dask.array<chunksize=(1431,), meta=np.ndarray>
    pnt_class            (space) float32 6kB 1.0 1.0 1.0 1.0 ... 1.0 1.0 1.0 1.0
Attributes:
    ps_selection_start_date:  20180920
    ps_selection_end_date:    20181113
```

As shown in the output above, `depsi.network_stm_selection` returns a subset of the input PS points. It uses the approximate radar coordinates with `azimuth_spacing` and `range_spacing` parameters to ensure the sparseness of the selected network points. The points with higher quality metrics (e.g., lower NMAD values) are prioritized during selection.

## Step 2: Form a network

One can use the `depsi.form_network` method to form a network from the selected network points. This method returns an Space-Time Matrix (STM) representing the arcs of the network. An example of using `depsi.form_network` to form a network is shown below:

```python
stm_network_arcs = depsi.form_network(
                stm_network_pnts,
                key_phase='sd_phase', # Variable for single-differenced phase, i.e. phase w.r.t. a reference epoch
                key_h2ph='h2ph', # Hight to phase conversion variable, needed for relative height estimation
                key_Btemporal='decimal_years', # Temporal baseline. Should be a decimal year vector with reference epoch as 0.
                key_xcrds='local_x', # X coordinates of points used for distance calculations
                key_ycrds='local_y', # Y coordinates of points used for distance calculations
                network_method = 'redundant', # Network formation method: 'redundant' or 'delaunay'
                max_length=800,  # Maximum length of network arcs, unit identical to the unit of xcrds and ycrds
            )
```
```output
<xarray.Dataset> Size: 536kB
Dimensions:  (space: 3722, time: 10)
Coordinates:
    source   (space) int64 30kB 0 0 1 1 1 1 2 ... 1368 1372 1376 1379 1386 1405
    target   (space) int64 30kB 12 1057 6 7 18 22 ... 1373 1408 1406 1393 1414
    uid      (space) int64 30kB 10013 11058 20007 ... 13801407 13871394 14061415
Dimensions without coordinates: space, time
Data variables:
    d_phase  (space, time) float32 149kB dask.array<chunksize=(1431, 10), meta=np.ndarray>
    Btemp    (time) float64 80B 0.0 0.01643 0.03285 ... 0.115 0.1314 0.1478
    h2ph     (space, time) float64 298kB dask.array<chunksize=(1431, 10), meta=np.ndarray>
Attributes:
    ps_selection_start_date:  20180920
    ps_selection_end_date:    20181113
    wavelength:               0.055465763
    mother_idx:               0
```

In the example above, we assumed the existence of several variables in the input dataset. In practice, one can load them from the stack data, metadata, or compute them as needed. Below is an example of how these variables can be prepared after the preselection step. In this way, the selected PS points will inherit these variables for subsequent network formation.

```python
# Assign metadata attributes
WAVELENGTH = 0.055465763 # Wavelength in meters (example value for Sentinel-1)
# Index of the mother (reference)acquisition
# Here we assume the first acquisition as the mother
# In practice this depends on the acquisition chosen for the stack corregistration
MOTHER_IDX = 0
stm_preselention_pnts = stm_preselention_pnts.assign_attrs({"wavelength": WAVELENGTH})
stm_preselention_pnts = stm_preselention_pnts.assign_attrs({"mother_idx": MOTHER_IDX})

# Compute single differential phase w.r.t. the mother acquisition
decimal_years = (stm_preselention_pnts["time"] - stm_preselention_pnts["time"].isel(time=stm_preselention_pnts.attrs["mother_idx"])).dt.days / 365.25
stm_preselention_pnts = stm_preselention_pnts.assign_coords({"decimal_years": (("time"), decimal_years.data)})
sd_phase = stm_preselention_pnts["phase"] - stm_preselention_pnts["phase"].isel(time=stm_preselention_pnts.attrs["mother_idx"])
stm_preselention_pnts = stm_preselention_pnts.assign(sd_phase=sd_phase)

# Using pyproj covert WGS84 to UTM coordinates
from depsi.transformations import latlonh_to_xyz
latlonh = np.array([
    stm_preselention_pnts["lat"].values,
    stm_preselention_pnts["lon"].values,
    np.zeros_like(stm_preselention_pnts["lat"].values)
]).T
xyz = latlonh_to_xyz(latlonh).T
stm_preselention_pnts = stm_preselention_pnts.assign_coords({
    "local_x": (("space"), xyz[:, 0]),
    "local_y": (("space"), xyz[:, 1]),
})
```

## Step 3: Estimate network ambiguities

After forming the network arcs, the ambiguity of the arcs can be estimated using the `depsi.periodogram` method. This method performs ambiguity estimation in time for each arc. An example is shown below:

```python
_, ambiguities, _, _, ens_coh = depsi.periodogram(
                stm_network_arcs, # formed arcs STM
                key_dphase="d_phase", # arc phase. The variable "d_phase" is a standard output from depsi.form_network
                key_h2ph="h2ph", #  arc height to phase conversion variable, a standard output from depsi.form_network
                key_Btemporal="Btemp", # temporal baseline variable, a standard output from depsi.form_network
            )
# Assining the results back to the stm_network_arcs dataset
stm_network_arcs["ambiguities"] = ambiguities
stm_network_arcs["temp_coh"] = ens_coh
```
```output
<xarray.Dataset> Size: 864kB
Dimensions:      (space: 3722, time: 10)
Coordinates:
    source       (space) int64 30kB 0 0 1 1 1 1 ... 1372 1376 1379 1386 1405
    target       (space) int64 30kB 12 1057 6 7 18 ... 1373 1408 1406 1393 1414
    uid          (space) int64 30kB 10013 11058 20007 ... 13871394 14061415
Dimensions without coordinates: space, time
Data variables:
    d_phase      (space, time) float32 149kB dask.array<chunksize=(1431, 10), meta=np.ndarray>
    Btemp        (time) float64 80B 0.0 0.01643 0.03285 ... 0.115 0.1314 0.1478
    h2ph         (space, time) float64 298kB dask.array<chunksize=(1431, 10), meta=np.ndarray>
    ambiguities  (space, time) float64 298kB dask.array<chunksize=(1431, 10), meta=np.ndarray>
    temp_coh     (space) float64 30kB dask.array<chunksize=(1431,), meta=np.ndarray>
Attributes:
    ps_selection_start_date:  20180920
    ps_selection_end_date:    20181113
    wavelength:               0.055465763
    mother_idx:               0
```

The new variable `ambiguities` are the estimated ambiguities for each arc in the network. The variable `temp_coh` represents the temporal coherence of each arc, which can be used to assess the quality of the ambiguity estimates.

## Step 4: Spatial integration of ambiguities

The arc ambiguities need to be spatially integrated to obtain the absolute ambiguities for each network point. This is done using the Multiple Hypothesis Testing (MHT) approach as introduced by [van Leijen, 2014](https://doi.org/10.4233/uuid:5dba48d7-ee26-4449-b674-caa8df93e71e). This step is implemented in the `depsi.spatial_integration`. An example is shown below:

```python
stm_network_arcs_solved, stm_network_pnts_solved = depsi.spatial_integration(
    stm_network_pnts,
    stm_network_arcs,
    key_arc_quality = "temp_coh",
    threshold_arc_quality = 0.5,
)
```
```output
<xarray.Dataset> Size: 370kB
Dimensions:              (space: 832, time: 10)
Coordinates:
  * space                (space) int64 7kB 1 19 22 53 ... 9107 6737 13727 13426
  * time                 (time) datetime64[us] 80B 2018-09-20 ... 2018-11-13
    lat                  (space) float32 3kB 52.25 52.25 52.25 ... 52.43 52.42
    lon                  (space) float32 3kB 4.95 4.943 4.956 ... 4.835 4.819
    azimuth              (space) int64 7kB 0 4 4 12 12 ... 1076 874 1538 1501
    range                (space) int64 7kB 2456 2334 2561 ... 1980 1471 1166
    decimal_years        (time) float64 80B 0.0 0.01643 ... 0.1314 0.1478
    local_x              (space) float64 7kB 3.899e+06 3.899e+06 ... 3.884e+06
    local_y              (space) float64 7kB 3.377e+05 3.372e+05 ... 3.274e+05
Data variables:
    complex              (space, time) complex64 67kB 0j ... (-1615.9874+1488...
    amplitude            (space, time) float32 33kB 0.0 0.0 ... 2.197e+03
    phase                (space, time) float32 33kB 0.0 0.0 0.0 ... 2.212 2.397
    h2ph                 (space, time) float64 67kB 0.0005848 ... 0.0002541
    time_selection_nmad  (space) float32 3kB 0.0 0.0 0.0 ... 0.02991 0.02992
    full_ts_nmad         (space) float32 3kB 0.0 0.0 0.0 ... 0.02991 0.02992
    full_ts_nad          (space) float32 3kB 0.0 0.0 0.0 ... 0.157 0.1395 0.3579
    pnt_class            (space) float32 3kB 1.0 1.0 1.0 1.0 ... 1.0 1.0 1.0 1.0
    sd_phase             (space, time) float32 33kB 0.0 0.0 0.0 ... 4.918 5.103
    ambiguities          (space, time) int16 17kB 0 0 0 0 0 0 0 ... 0 0 0 0 0 0
    unwrapped_phase      (space, time) float64 67kB 0.0 0.0 0.0 ... 4.918 5.103
Attributes:
    ps_selection_start_date:  20180920
    ps_selection_end_date:    20181113
    wavelength:               0.055465763
    mother_idx:               0
    idx_refpnt:               0
```

This adds `ambiguities` and `unwrapped_phase` variables to the network points dataset. One more metadata attribute `idx_refpnt` is also added to indicate the index of the reference point in the `stm_network_pnts`.

## References

Van Leijen, Frederik Johannes. "Persistent scatterer interferometry based on geodetic estimation theory." (2014).