"""Visualization utilities for point and arc STMs.

These functions replace the plotting boilerplate repeated across the example notebooks: drawing a Mean
Reflectivity Map (MRM) as a radar-coordinate basemap, then overlaying point STMs (as a scatter, optionally
coloured by a data variable) or arc STMs (as line segments between their source/target points, optionally
coloured by a data variable) on top of it.

All three functions plot in radar coordinates (`range`/`azimuth` by default), matching the coordinate system
of the MRM raster - not `lat`/`lon`, which does not align with the MRM's pixel grid.
"""

from typing import Literal

import matplotlib.colors as plc
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from matplotlib.axes import Axes


def plot_mrm(
    mrm: xr.DataArray,
    ax: Axes | None = None,
    clim: tuple[float, float] = (0, 40000),
    cmap: str = "gray",
    aspect: float | Literal["auto", "equal"] | None = None,
) -> Axes:
    """Plot a Mean Reflectivity Map (MRM) as a radar-coordinate basemap.

    Parameters
    ----------
    mrm: xarray.DataArray
        Mean Reflectivity Map, e.g. as returned by `stack.slcstack.mrm()`. If it is a lazy (dask-backed)
        array, it is computed before plotting.
    ax: matplotlib.axes.Axes, optional
        Axes to plot on. A new figure and axes are created if not given.
    clim: tuple of float
        Colour limits passed to the MRM image, as (vmin, vmax). Default is (0, 40000).
    cmap: str
        Colormap for the MRM image. Default is "gray".
    aspect: float, "auto", "equal", or None
        Axes aspect ratio, forwarded to `ax.set_aspect`. Left as the matplotlib default (None) unless given:
        the range/azimuth pixel spacing ratio is product-specific, so there is no universally correct default.

    Returns
    -------
    matplotlib.axes.Axes
        The axes the MRM was plotted on.
    """
    if ax is None:
        _, ax = plt.subplots()

    if getattr(mrm, "chunks", None) is not None:
        mrm = mrm.compute()

    im = mrm.plot(ax=ax, cmap=cmap, add_colorbar=False)
    im.set_clim(clim)

    if aspect is not None:
        ax.set_aspect(aspect)
    ax.axis("off")

    return ax


def plot_points(
    stm_points: xr.Dataset,
    mrm: xr.DataArray | None = None,
    color_by: str | None = None,
    ax: Axes | None = None,
    cmap: str = "jet_r",
    colorbar: bool = True,
    x_coord: str = "range",
    y_coord: str = "azimuth",
    mrm_kwargs: dict | None = None,
    **scatter_kwargs,
) -> Axes:
    """Plot a point STM as a scatter, optionally over an MRM basemap.

    Parameters
    ----------
    stm_points: xarray.Dataset
        Point STM with a `space` dimension, and `x_coord`/`y_coord` (default `range`/`azimuth`) coordinates
        or data variables on that dimension.
    mrm: xarray.DataArray, optional
        Mean Reflectivity Map to draw as a basemap before the scatter, via `plot_mrm`.
    color_by: str, optional
        Name of a data variable on `stm_points` to colour the points by. If not given, `scatter_kwargs["c"]`
        is used if present, otherwise all points are drawn in a single default colour.
    ax: matplotlib.axes.Axes, optional
        Axes to plot on. A new figure and axes are created if not given.
    cmap: str
        Colormap used when `color_by` is given. Default is "jet_r".
    colorbar: bool
        Whether to add a colorbar when `color_by` is given. Default is True.
    x_coord, y_coord: str
        Names of the coordinates/data variables on `stm_points` to use for the point positions.
        Default is "range" and "azimuth".
    mrm_kwargs: dict, optional
        Extra keyword arguments forwarded to `plot_mrm` when `mrm` is given.
    **scatter_kwargs
        Extra keyword arguments forwarded to `ax.scatter` (e.g. `s`, `marker`, `vmin`, `vmax`, `norm`).

    Returns
    -------
    matplotlib.axes.Axes
        The axes the points were plotted on.
    """
    if ax is None:
        _, ax = plt.subplots()

    if mrm is not None:
        plot_mrm(mrm, ax=ax, **(mrm_kwargs or {}))

    if color_by is not None:
        scatter_kwargs["c"] = stm_points[color_by].values
        scatter_kwargs.setdefault("cmap", cmap)

    scatter = ax.scatter(stm_points[x_coord].values, stm_points[y_coord].values, **scatter_kwargs)

    if color_by is not None and colorbar:
        plt.colorbar(scatter, ax=ax, label=color_by)

    return ax


def plot_arcs(
    stm_points: xr.Dataset,
    stm_arcs: xr.Dataset,
    mrm: xr.DataArray | None = None,
    color_by: str | None = None,
    ax: Axes | None = None,
    cmap: str = "jet_r",
    vmin: float | None = None,
    vmax: float | None = None,
    colorbar: bool = True,
    linewidth: float = 0.5,
    default_color: str = "tab:blue",
    x_coord: str = "range",
    y_coord: str = "azimuth",
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
    mrm: xarray.DataArray, optional
        Mean Reflectivity Map to draw as a basemap before the arcs, via `plot_mrm`.
    color_by: str, optional
        Name of a data variable on `stm_arcs` to colour the arcs by (its absolute value is used, e.g. for a
        complex-valued temporal coherence). If not given, all arcs are drawn in `default_color`.
    ax: matplotlib.axes.Axes, optional
        Axes to plot on. A new figure and axes are created if not given.
    cmap: str
        Colormap used when `color_by` is given. Default is "jet_r".
    vmin, vmax: float, optional
        Colour scale limits used when `color_by` is given. Default to the data's own min/max.
    colorbar: bool
        Whether to add a colorbar when `color_by` is given. Default is True.
    linewidth: float
        Line width for the arcs. Default is 0.5.
    default_color: str
        Colour used for all arcs when `color_by` is not given. Default is "tab:blue".
    x_coord, y_coord: str
        Names of the coordinates/data variables on `stm_points` to use for the arc endpoint positions.
        Default is "range" and "azimuth".
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

    xx = np.stack([stm_points[x_coord].values[source], stm_points[x_coord].values[target]]).T
    yy = np.stack([stm_points[y_coord].values[source], stm_points[y_coord].values[target]]).T

    norm = None
    if color_by is not None:
        values = np.abs(stm_arcs[color_by].values)
        norm = plc.Normalize(
            vmin=vmin if vmin is not None else values.min(),
            vmax=vmax if vmax is not None else values.max(),
        )
        cmap_obj = plt.get_cmap(cmap)
        colors = cmap_obj(norm(values))
    else:
        colors = [default_color] * n_arcs

    for i in range(n_arcs):
        ax.plot(xx[i], yy[i], color=colors[i], linewidth=linewidth)

    if color_by is not None and colorbar:
        sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])
        plt.colorbar(sm, ax=ax, label=color_by)

    return ax
