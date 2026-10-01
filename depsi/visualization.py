"""Visualization utilities for point and arc STMs.

These functions replace the plotting boilerplate repeated across the example notebooks: drawing a Mean
Reflectivity Map (MRM) as a radar-coordinate basemap, then overlaying point STMs (as a scatter, optionally
coloured by a data variable) or arc STMs (as line segments between their source/target points, optionally
coloured by a data variable) on top of it.

All three functions plot in radar coordinates (`range`/`azimuth` by default), matching the coordinate system
of the MRM raster - not `lat`/`lon`, which does not align with the MRM's pixel grid.
"""

import matplotlib.colors as plc
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from matplotlib.axes import Axes


def plot_mrm(
    mrm: xr.DataArray,
    pctl_min: int = 2,
    pctl_max: int = 95,
    ax: Axes | None = None,
) -> Axes:
    """Plot a Mean Reflectivity Map (MRM) as a radar-coordinate basemap.

    Parameters
    ----------
    mrm: xarray.DataArray
        Mean Reflectivity Map, e.g. as returned by `stack.slcstack.mrm()`. If it is a lazy (dask-backed)
        array, it is computed before plotting.
    pctl_min: int
        Minimum percentile for the MRM image. Default is 2.
    pctl_max: int
        Maximum percentile for the MRM image. Default is 95.
    ax: matplotlib.axes.Axes, optional
        Axes to plot on. A new figure and axes are created if not given.

    Returns
    -------
    matplotlib.axes.Axes
        The axes the MRM was plotted on.
    """
    if ax is None:
        _, ax = plt.subplots()
    mrm_values = mrm.data.flatten()
    pctl_min = np.percentile(mrm_values, pctl_min)
    pctl_max = np.percentile(mrm_values, pctl_max)

    mrm.plot(ax=ax, vmin=pctl_min, vmax=pctl_max, cmap="gray", add_colorbar=False)

    ax.axis("off")

    return ax


def plot_points(
    stm_points: xr.Dataset,
    color_by: str,
    mrm: xr.DataArray | None = None,
    ax: Axes | None = None,
    colorbar: bool = True,
    mrm_kwargs: dict | None = None,
    **scatter_kwargs,
) -> Axes:
    """Plot a point STM as a scatter, optionally over an MRM basemap.

    Parameters
    ----------
    stm_points: xarray.Dataset
        Point STM with `range` and `azimuth` coordinates/data variables (used as x/y), plus
        `color_by` as a data variable used for colouring.
    color_by: str
        Name of a data variable on `stm_points` used as `hue` in
        `stm_points.plot.scatter(..., hue=color_by)`.
    mrm: xarray.DataArray, optional
        Mean Reflectivity Map to draw as a basemap before the scatter, via `plot_mrm`.
    ax: matplotlib.axes.Axes, optional
        Axes to plot on. A new figure and axes are created if not given.
    colorbar: bool
        Whether to add a colorbar when `color_by` is given. Default is True.
    mrm_kwargs: dict, optional
        Extra keyword arguments forwarded to `plot_mrm` when `mrm` is given.
    **scatter_kwargs
        Extra keyword arguments forwarded to `stm_points.plot.scatter` (e.g. `s`, `marker`, `vmin`,
        `vmax`, `norm`). Defaults set by this function are `cmap="jet_r"`, `s=1`,
        `marker="s"`, and `edgecolor="none"` unless provided explicitly.

    Returns
    -------
    matplotlib.axes.Axes
        The axes the points were plotted on.
    """
    if ax is None:
        _, ax = plt.subplots()

    # Plot the MRM basemap first (if any)
    if mrm is not None:
        plot_mrm(mrm, ax=ax, **(mrm_kwargs or {}))

    # Set default scatter kwargs
    scatter_kwargs.setdefault("cmap", "jet_r")
    scatter_kwargs.setdefault("s", 1)
    scatter_kwargs.setdefault("marker", "s")
    scatter_kwargs.setdefault("edgecolor", "none")

    # Plot the scatter points
    stm_points.plot.scatter(
        x="range",
        y="azimuth",
        hue=color_by,
        ax=ax,
        add_colorbar=colorbar,
        **scatter_kwargs,
    )

    return ax


def plot_arcs(
    stm_points: xr.Dataset,
    stm_arcs: xr.Dataset,
    color_by: str,
    mrm: xr.DataArray | None = None,
    cmap: str = "jet_r",
    vmin: float = 0.0,
    vmax: float = 1.0,
    colorbar: bool = True,
    linewidth: float = 0.5,
    ax: Axes | None = None,
    mrm_kwargs: dict | None = None,
) -> Axes:
    """Plot arc STMs as line segments between their source and target points, optionally over an MRM basemap.

    Parameters
    ----------
    stm_points: xarray.Dataset
        Point STM the arcs' `source`/`target` indices refer to, with `x_coord`/`y_coord` (default
        `range`/`azimuth`) coordinates or data variables on its `space` dimension.
    stm_arcs: xarray.Dataset
        Arc STM with a `space` dimension, and `source`/`target` data variables holding integer positional
        indices into `stm_points`'s `space` dimension.
    color_by: str
        Name of a data variable on `stm_arcs` used to color the arcs. Its absolute value is used
        (e.g. for complex-valued temporal coherence).
    mrm: xarray.DataArray, optional
        Mean Reflectivity Map to draw as a basemap before the arcs, via `plot_mrm`.
    cmap: str
        Colormap used for arc colors. Default is "jet_r".
    vmin, vmax: float, optional
        Color scale limits for `color_by`. If not provided, the data min/max are used.
    colorbar: bool
        Whether to add a colorbar for `color_by`. Default is True.
    linewidth: float
        Line width for the arcs. Default is 0.5.
    ax: matplotlib.axes.Axes, optional
        Axes to plot on. A new figure and axes are created if not given.
    mrm_kwargs: dict, optional
        Extra keyword arguments forwarded to `plot_mrm` when `mrm` is given.

    Returns
    -------
    matplotlib.axes.Axes
        The axes the arcs were plotted on.
    """
    if ax is None:
        _, ax = plt.subplots()

    if mrm is not None:
        plot_mrm(mrm, ax=ax, **(mrm_kwargs or {}))

    source = stm_arcs["source"].values
    target = stm_arcs["target"].values
    n_arcs = stm_arcs.sizes["space"]

    xx = np.stack([stm_points["range"].values[source], stm_points["range"].values[target]]).T
    yy = np.stack([stm_points["azimuth"].values[source], stm_points["azimuth"].values[target]]).T

    values = np.abs(stm_arcs[color_by].values)
    norm = plc.Normalize(
        vmin=vmin if vmin is not None else values.min(),
        vmax=vmax if vmax is not None else values.max(),
    )
    cmap_obj = plt.get_cmap(cmap)
    colors = cmap_obj(norm(values))

    for i in range(n_arcs):
        ax.plot(xx[i], yy[i], color=colors[i], linewidth=linewidth)

    if color_by is not None and colorbar:
        sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])
        plt.colorbar(sm, ax=ax, label=color_by)

    return ax
