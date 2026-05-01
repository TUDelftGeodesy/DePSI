# Atmospheric Phase Screen estimation

In DePSI, Atmospheric Phase Screen (APS) is estimated using Kriging methods from the residual phase of network points after removing modeled deformation, assuming the APS is a spatially correlated and temporally uncorrelated signal. In background, we use the `pykrige` library to perform Kriging interpolation.

## Step 1: Estimate model parameters for network points

Before estimating APS, the model parameters should be estimated for the later removal of modeled deformation from the observed phase. Below is an example of estimating model parameters after performing the spatial intergration of the network, i.e. the ambiguities and unwrapped phases are available per network point in `stm_network_pnts_solved`:

```python
stm_network_pnts_solved, model_parameter_layer_names = estimate_model_params(
    stm=stm_network_pnts_solved,
    models=["offset","linear", "height"],
    key_observations="unwrapped_phase",
    key_h2ph="sd_h2ph",
    key_time="temporal_baseline"
)
```

```python
# Temporarily rename the model parameter layers for later identification
stm_pnts_output = stm_pnts_output.rename_vars({
    "pnt_offset": "mother_atmosphere",
})
# Temporally rename the space dimension since it will otherwise clash with the space dimension of the predicted coordinates
stm_pnts_output = stm_pnts_output.rename_dims({"space": "space_fon"})
```


## Step 2: Prepare prediction coordinates

Kriging system is based on distance. Therefore we need Euclidean coordinates (in meter units) in both the network points dataset and the prediction points dataset. Below is an example of converting geographic coordinates (longitude and latitude in degrees) to Euclidean coordinates using `depsi.utils.convert_geographic_coords_to_euclidean` function:


```python
from depsi.utils import convert_geographic_coords_to_euclidean

# RD New, suitable for Dutch region.
# You may change it to other EPSG codes as needed.
target_epsg = 28992

# Add to network points dataset
crd_x, crd_y = convert_geographic_coords_to_euclidean(
    stm["lon"].values, stm["lat"].values, target_crs=f"EPSG:{target_epsg}"
)
stm_pnts_output = stm_pnts_output.assign_coords(
    {
        "x_euclidean": ("space_fon", crd_x),
        "y_euclidean": ("space_fon", crd_y),
    }
)

# Prepare prediction coordinates for predicting points
# For example, we predict APS for all preselected points
crd_x_preselection, crd_y_preselection = convert_geographic_coords_to_euclidean(
    stm_preselection["lon"].values, stm_preselection["lat"].values, target_crs=f"EPSG:{target_epsg}"
)
atmosphere_prediction_coords = xarray.Dataset(
    coords={
        f"x_euclidean": ("space", crd_x_preselection),
        f"y_euclidean": ("space", crd_y_preselection),
    }
)
```

## Step 3: Estimate atmospheric phase screens

Then we can estimate the atmospheric phase screens for the prediction coordinates using the `estimate_atmosphere_phase` function, and add the estimated APS into the preseletion STM:

```python
atmospheric_phase_screens = estimate_atmosphere_phase(
    stm=stm_pnts_output,
    prediction_coords=atmosphere_prediction_coords,
    psc_phase_residuals="phase_residuals", # Data variables of phase residuals, added by `estimate_model_params` step
    atmosphere_mother="mother_atmosphere", # Mannually assigned atmosphere for mother (reference) acquisition
    key_Btemporal="temporal_baseline", # Temporal baseline coordinates in decimal years
)

# create the zero layer for the mother atmosphere to be appended to the interferometric atmospheric phase screens
mother_atmo_stm = xarray.Dataset(
    data_vars={
        "atmosphere_predicted": ("space", np.zeros((len(stm.space), ))),
        "atmosphere_sigmasq": ("space", np.zeros((len(stm.space), ))),},
    coords=stm["phase"].sel(time=stm.ps_sd_mother).coords
)

# Concatenate the mother atmosphere to the estimated atmospheric phase screens
atmosphere = xarray.concat([atmospheric_phase_screens, mother_atmo_stm], dim="time").sortby("time")

# Add the atmosphere into the original STM
stm["atmosphere_predicted"] = atmosphere["atmosphere_predicted"].transpose("space", "time")
stm["atmosphere_sigmasq"] = atmosphere["atmosphere_sigmasq"].transpose("space", "time")
```

Note that for the mother acquisition, the atmospheric phase screen is set to zero.

