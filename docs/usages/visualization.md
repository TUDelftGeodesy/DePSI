# Visualization

The `depsi.visualization` provides a set of functions for visualizing DePSI (intermediate) results.

## Plot a Mean Reflectivity Map

`plot_mrm` accepts a single `xarray.DataArray` containing the MRM. The array is expected to be arranged in radar coordinates, typically with dimensions or coordinates named `azimuth` and `range`.

```python
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from depsi.visualization import plot_mrm

# Creating a dummy 4x5 MRM for demonstration purposes
# In practice, mrm can be computed from a stack of radar images
mrm = xr.DataArray(
    np.arange(20, dtype=float).reshape(4, 5),
    dims=("azimuth", "range"),
    coords={"azimuth": np.arange(4), "range": np.arange(5)},
)

ax = plot_mrm(mrm)
ax.set_title("Mean Reflectivity Map")
plt.show()
```

This produces a grayscale image in which the color scale is defined by the 2nd and 95th percentiles of the input data by default. The function returns the Matplotlib axes it used, which makes it easy to overlay scatter or arc plots afterward.

## Plot point STMs

`plot_points` is used to visualize a point-based STM as a radar-coordinate scatter plot. It expects an `xarray.Dataset` with radar coordinates (`range` and `azimuth`), and colors points according to a data variable, such as `velocity`.

```python
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from depsi.visualization import plot_mrm, plot_points

mrm = xr.DataArray(
    np.arange(20, dtype=float).reshape(4, 5),
    dims=("azimuth", "range"),
    coords={"azimuth": np.arange(4), "range": np.arange(5)},
)

stm_points = xr.Dataset(
    data_vars={
        "velocity": ("space", np.array([1.0, 2.0, 3.0, 4.0])),
    },
    coords={
        "range": ("space", np.array([0, 1, 2, 3])),
        "azimuth": ("space", np.array([0, 1, 2, 3])),
    },
)

fig, ax = plt.subplots()
plot_points(stm_points, color_by="velocity", mrm=mrm, ax=ax)
```

## Plot arc STMs

`plot_arcs` is used to visualize arc-based STMs as line segments between source and target points, colored by a scalar arc variable such as temporal coherence. Input data are expected to include integer `source` and `target` indices into a point STM, along with `range` and `azimuth` coordinates on the point dataset.

```python
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from depsi.visualization import plot_arcs, plot_mrm

mrm = xr.DataArray(
    np.arange(20, dtype=float).reshape(4, 5),
    dims=("azimuth", "range"),
    coords={"azimuth": np.arange(4), "range": np.arange(5)},
)

stm_points = xr.Dataset(
    data_vars={},
    coords={
        "range": ("space", np.array([0, 1, 2, 3])),
        "azimuth": ("space", np.array([0, 1, 2, 3])),
    },
)

stm_arcs = xr.Dataset(
    data_vars={
        "source": ("space", np.array([0, 1, 2])),
        "target": ("space", np.array([1, 2, 3])),
        "temp_coh": ("space", np.array([0.9, 0.1, 0.5])),
    },
)

fig, ax = plt.subplots()
plot_arcs(stm_points, stm_arcs, color_by="temp_coh", mrm=mrm, ax=ax)
plt.show()
```

## Use MRM as a Basemap

The MRM as an xarray.DataArray can serve as a basemap for overlaying points or arcs. This allows you to maintain the radar coordinate system while visualizing additional data on top of the reflectivity map.

```python
from depsi.visualization import plot_arcs, plot_mrm, plot_points

fig, ax = plt.subplots()
plot_points(stm_points, color_by="velocity", mrm=mrm, ax=ax)
plot_arcs(stm_points, stm_arcs, color_by="temp_coh", mrm=mrm, ax=ax)
```
