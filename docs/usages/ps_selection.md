# Select Persistent Scatterers (PS)

The first step in the DePSI workflow is to identify candidate Persistent Scatterer (PS) points from the interferogram stack. PS points are characterized by their coherent amplitude properties over time.

One can use `depsi.ps_selection` method to perform PS selection based on various criteria. Currently, the following methods are implemented:

- **Normalized Median Absolute Deviation (NMAD) (Default)**, calculated as the ratio of the median absolute deviation to the median of the amplitude time series. Lower NMAD values indicate more stable scatterers with less temporal variability. Compared to NAD, NMAD uses the median instead of the mean, making it more robust to outliers and better suited for identifying high-quality PS points in the presence of noise or anomalous measurements. 

- **Normalized Amplitude Dispersion (NAD)**, calculated as the ratio of the standard deviation to the mean of the amplitude time series. Lower NAD values indicate more stable scatterers.

For more details on these metrics, please refer to [W. S. Brouwer, 2025](https://doi.org/10.1109/TGRS.2025.3600893).

An example of using `depsi.ps_selection` to select PS points based on NMAD is shown below:

```python
from depsi.ps_selection import ps_selection

stm_preselention_pnts = depsi.ps_selection(stack, # xr.Dataset of interferogram stack
                   method="nmad", # Selection method: "nmad" or "nad"
                   threshold=0.05, # Points with NMAD below this threshold are selected
                   output_chunks=5000, # Chunk size for output Dask array
                   )
```
```output
<xarray.Dataset> Size: 6MB
Dimensions:              (time: 10, space: 16510)
Coordinates:
  * time                 (time) datetime64[us] 80B 2018-09-20 ... 2018-11-13
  * space                (space) int64 132kB 0 1 2 3 ... 16506 16507 16508 16509
    lat                  (space) float32 66kB dask.array<chunksize=(5000,), meta=np.ndarray>
    lon                  (space) float32 66kB dask.array<chunksize=(5000,), meta=np.ndarray>
    azimuth              (space) int64 132kB dask.array<chunksize=(5000,), meta=np.ndarray>
    range                (space) int64 132kB dask.array<chunksize=(5000,), meta=np.ndarray>
Data variables:
    complex              (space, time) complex64 1MB dask.array<chunksize=(5000, 10), meta=np.ndarray>
    amplitude            (space, time) float32 660kB dask.array<chunksize=(5000, 10), meta=np.ndarray>
    phase                (space, time) float32 660kB dask.array<chunksize=(5000, 10), meta=np.ndarray>
    h2ph                 (space, time) float64 1MB dask.array<chunksize=(5000, 10), meta=np.ndarray>
    time_selection_nmad  (space) float32 66kB dask.array<chunksize=(5000,), meta=np.ndarray>
    full_ts_nmad         (space) float32 66kB dask.array<chunksize=(5000,), meta=np.ndarray>
    full_ts_nad          (space) float32 66kB dask.array<chunksize=(5000,), meta=np.ndarray>
    pnt_class            (space) int8 17kB 1 1 1 1 1 1 1 1 1 ... 1 1 1 1 1 1 1 1
Attributes:
    ps_selection_start_date:  20180920
    ps_selection_end_date:    20181113
```

As shown in the output above, `depsi.ps_selection` returns the selected PS as a new `xarray.Dataset` containing two main dimensions: `space` and `time`, with multiple coordinates and data variables. It attaches the following data variables to the output dataset:

- **full_ts_nmad**: NMAD values calculated over the full time series for each point.
- **full_ts_nad**: NAD values calculated over the full time series for each point.
- **time_selection_nmad**: NMAD values calculated over the selected time series for each point.
- **pnt_class**: Point classification label to distinguish different types of points. Since this is a PS selection step, all selected points are labeled as 1 (PS).

Note that in this case, `time_selection_nmad` is equivalent to `full_ts_nmad` since the selection is performed over the full time series. However, if a subset of the time series is specified for selection, `time_selection_nmad` will reflect the NMAD values calculated over that subset. One can perform PS selection on a specific time subset by providing `ps_selection_start_date` and `ps_selection_end_date` parameters to the `depsi.ps_selection` method:

```python
stm_preselention_pnts = depsi.ps_selection(stack,
                   method="nmad",
                   threshold=0.05,
                   output_chunks=5000,
                   ps_selection_start_date="20181002", # Start date for selection
                   ps_selection_end_date="20181101",   # End date for selection
                   )
```

In this case, although the `time` dimension of the output is still identical to the input stack, i.e. no epochs are removed, the selection NMAD values (`time_selection_nmad`) are calculated only over the specified date range (2018-10-02 to 2018-11-01).

Another note is that by default, `depsi.ps_selection` does not perform in-memory persistence of the output dataset. Instead, it builds a Dask computation graph for the output, as shown in the output above, `time_selection_nmad`, `full_ts_nmad`, and `full_ts_nad` are all Dask arrays. This allows for efficient handling of large datasets without loading everything into memory at once. However, if one wants to persist the output dataset in memory for faster access in subsequent steps, one can set the `mem_persist` parameter to `True`:

```python
stm_preselention_pnts = depsi.ps_selection(stack,
                   method="nmad",
                   threshold=0.05, 
                   output_chunks=5000,
                   mem_persist=True, # Persist the output in memory for faster access later
                   )
```

This will trigger Dask to compute and store the output dataset in memory immediately after selection.

## References
W. S. Brouwer and R. F. Hanssen, "On the Definition of an Independent Stochastic Model for InSAR Time Series," in IEEE Transactions on Geoscience and Remote Sensing, vol. 63, pp. 1-11, 2025, Art no. 4508311, doi: 10.1109/TGRS.2025.3600893.

