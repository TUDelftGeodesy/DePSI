import dask.array as da
import matplotlib.pyplot as plt
import numpy as np
import pytest
import xarray as xr
from matplotlib.collections import PathCollection, QuadMesh

from depsi.visualization import plot_arcs, plot_mrm, plot_points


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


@pytest.fixture
def mrm():
    # Real MRM DataArrays carry evenly-spaced range/azimuth coordinates, which is what makes
    # xarray's `.plot()` render an AxesImage (imshow) rather than falling back to a QuadMesh
    # (pcolormesh) for coordinate-less 2D data.
    data = np.arange(20, dtype=float).reshape(4, 5)
    return xr.DataArray(data, dims=("azimuth", "range"), coords={"azimuth": np.arange(4), "range": np.arange(5)})


@pytest.fixture
def stm_points():
    return xr.Dataset(
        data_vars={
            "velocity": ("space", np.array([1.0, 2.0, 3.0, 4.0])),
        },
        coords={
            "range": ("space", np.array([0, 1, 2, 3])),
            "azimuth": ("space", np.array([0, 1, 2, 3])),
        },
    )


@pytest.fixture
def stm_arcs():
    return xr.Dataset(
        data_vars={
            "source": ("space", np.array([0, 1, 2])),
            "target": ("space", np.array([1, 2, 3])),
            "temp_coh": ("space", np.array([0.9, 0.1, 0.5])),
        },
    )


def test_plot_mrm_returns_axes_with_image(mrm):
    # xarray's DataArray.plot() renders 2D data as a QuadMesh (pcolormesh), not an AxesImage.
    ax = plot_mrm(mrm)
    assert ax is not None
    assert len(ax.collections) == 1
    assert ax.collections[0].get_clim() == (0, 40000)


def test_plot_mrm_custom_clim_and_cmap(mrm):
    ax = plot_mrm(mrm, clim=(0, 10), cmap="viridis")
    quadmesh = ax.collections[0]
    assert quadmesh.get_clim() == (0, 10)
    assert quadmesh.get_cmap().name == "viridis"


def test_plot_mrm_computes_dask_backed_array(mrm):
    lazy_mrm = mrm.copy()
    lazy_mrm.data = da.from_array(mrm.data, chunks=(2, 5))

    ax = plot_mrm(lazy_mrm)
    assert len(ax.collections) == 1


def test_plot_mrm_no_colorbar(mrm):
    ax = plot_mrm(mrm)
    assert len(ax.figure.axes) == 1  # no colorbar axes added


def test_plot_points_positions(stm_points):
    ax = plot_points(stm_points)
    offsets = ax.collections[0].get_offsets()
    expected = np.column_stack([stm_points["range"].values, stm_points["azimuth"].values])
    assert np.allclose(offsets, expected)


def test_plot_points_color_by_adds_colorbar(stm_points):
    ax = plot_points(stm_points, color_by="velocity")
    scatter = ax.collections[0]
    assert np.allclose(scatter.get_array(), stm_points["velocity"].values)
    assert len(ax.figure.axes) == 2  # scatter axes + colorbar axes


def test_plot_points_no_color_by_no_colorbar(stm_points):
    ax = plot_points(stm_points)
    assert len(ax.figure.axes) == 1


def test_plot_points_over_mrm(stm_points, mrm):
    ax = plot_points(stm_points, mrm=mrm)
    assert sum(isinstance(c, QuadMesh) for c in ax.collections) == 1
    assert sum(isinstance(c, PathCollection) for c in ax.collections) == 1


def test_plot_arcs_line_positions(stm_points, stm_arcs):
    ax = plot_arcs(stm_points, stm_arcs)
    assert len(ax.lines) == 3

    source = stm_arcs["source"].values
    target = stm_arcs["target"].values
    for i, line in enumerate(ax.lines):
        expected_x = [stm_points["range"].values[source[i]], stm_points["range"].values[target[i]]]
        expected_y = [stm_points["azimuth"].values[source[i]], stm_points["azimuth"].values[target[i]]]
        assert np.allclose(line.get_xdata(), expected_x)
        assert np.allclose(line.get_ydata(), expected_y)


def test_plot_arcs_default_color(stm_points, stm_arcs):
    ax = plot_arcs(stm_points, stm_arcs)
    colors = {line.get_color() for line in ax.lines}
    assert colors == {"tab:blue"}
    assert len(ax.figure.axes) == 1  # no colorbar without color_by


def test_plot_arcs_color_by_adds_colorbar(stm_points, stm_arcs):
    ax = plot_arcs(stm_points, stm_arcs, color_by="temp_coh")
    colors = {tuple(np.atleast_1d(line.get_color())) for line in ax.lines}
    assert len(colors) == 3  # three distinct temp_coh values -> three distinct colors
    assert len(ax.figure.axes) == 2


def test_plot_arcs_color_by_uses_absolute_value(stm_points, stm_arcs):
    stm_arcs["temp_coh"] = ("space", np.array([-0.9, 0.1, 0.5]))
    ax_negative = plot_arcs(stm_points, stm_arcs, color_by="temp_coh", vmin=0, vmax=1)

    stm_arcs["temp_coh"] = ("space", np.array([0.9, 0.1, 0.5]))
    ax_positive = plot_arcs(stm_points, stm_arcs, color_by="temp_coh", vmin=0, vmax=1)

    assert np.array_equal(ax_negative.lines[0].get_color(), ax_positive.lines[0].get_color())


def test_plot_arcs_over_mrm(stm_points, stm_arcs, mrm):
    ax = plot_arcs(stm_points, stm_arcs, mrm=mrm)
    assert sum(isinstance(c, QuadMesh) for c in ax.collections) == 1
    assert len(ax.lines) == 3
