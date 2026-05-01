# Densification

After the network PS points have been estimated, and the atmospheric phase screen estimated from the network residual phases, the next step is to estimate the deformation time series at all PS points from the network PS points. This step is called densification, as it densifies the deformation estimates from the sparse network PS points to all PS points.

This is done using the `depsi.densification` method. It constructs arcs between each non-network PS point and N nearest network PS points, and uses the same arc estimation approach as in the network estimation step to estimate the deformation time series at each non-network PS point.

An example of using `depsi.densification` is shown below:

```python
stm_densified = depsi.densification(stm_preselention_pnts, # xr.Dataset of pre-selected PS points
                                    stm_network_pnts_solved, # xr.Dataset of network PS points with estimated ambiguities
                                    key_h2ph='h2ph', # Variable name for height-to-phase conversion
                                    key_Btemporal='decimal_years' # Temporal baseline in decimal years, with reference epoch at 0
                                    )
```
```output
<xarray.Dataset> Size: 9MB
Dimensions:              (space: 16510, time: 10)
Coordinates:
  * space                (space) int64 132kB 1 19 22 53 ... 16506 16508 16509
  * time                 (time) datetime64[us] 80B 2018-09-20 ... 2018-11-13
    lat                  (space) float32 66kB 52.25 52.25 52.25 ... 52.5 52.5
    lon                  (space) float32 66kB 4.95 4.943 4.956 ... 4.943 4.962
    azimuth              (space) int64 132kB 0 4 4 12 12 ... 1999 1999 1999 1999
    range                (space) int64 132kB 2456 2334 2561 ... 3095 3648 3990
    decimal_years        (time) float64 80B 0.0 0.01643 ... 0.1314 0.1478
    local_x              (space) float64 132kB 3.899e+06 3.899e+06 ... 3.876e+06
    local_y              (space) float64 132kB 3.377e+05 3.372e+05 ... 3.365e+05
    pnt_uid              (space) uint64 132kB 17210692856312994489 ... 144259...
Data variables:
    complex              (space, time) complex64 1MB 0j ... (-125328.234+6352...
    amplitude            (space, time) float32 660kB 0.0 0.0 ... 1.405e+05
    phase                (space, time) float32 660kB 0.0 0.0 ... -1.854 2.672
    h2ph                 (space, time) float64 1MB 0.0005848 ... 0.0004645
    time_selection_nmad  (space) float32 66kB 0.0 0.0 0.0 ... 0.04473 0.03029
    full_ts_nmad         (space) float32 66kB 0.0 0.0 0.0 ... 0.04473 0.03029
    full_ts_nad          (space) float32 66kB 0.0 0.0 0.0 ... 0.1557 0.1225
    pnt_class            (space) float32 66kB 1.0 1.0 1.0 1.0 ... 1.0 1.0 1.0
    sd_phase             (space, time) float32 660kB 0.0 0.0 ... -0.5955 3.931
    ambiguities          (space, time) float64 1MB 0.0 0.0 0.0 ... -2.0 -3.0
    unwrapped_phase      (space, time) float64 1MB 0.0 0.0 0.0 ... -13.16 -14.92
    local_temp_coh       (space) float64 132kB nan nan nan ... 0.5748 0.6928
Attributes:
    ps_selection_start_date:  20180920
    ps_selection_end_date:    20181113
    wavelength:               0.055465763
    mother_idx:               0
    idx_refpnt:               0
```

This step adds `ambiguities` and `unwrapped_phase` data variables to the output dataset, which contain the estimated integer ambiguities and unwrapped phases at each PS point, respectively. Besides, the `local_temp_coh` variable is also added, which indicates the local temporal coherence, representing the best temporal coherence of all arcs connecting to the given PS point.