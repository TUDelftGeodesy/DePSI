import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
import xarray as xr
from numpy.testing import assert_allclose
from shapely.geometry import MultiPolygon, Polygon
from xarray import DataTree

from depsi.ds import (
    _coherence_matrix,
    _phase_linking,
    _segmentation,
    _shp_test,
    assign_parcel_id,
    ds_phase_estimation,
    ps_ds_arc,
    select_common_fop_ref,
)


@pytest.fixture
def slc_stack_for_assign_id():
    az = np.array([0.0, 1.0, 2.0])
    rg = np.array([0.0, 1.0, 2.0])
    lon = np.array(
        [
            [0.5, 1.5, 2.5],
            [0.5, 1.5, 2.5],
            [0.5, 1.5, 2.5],
        ]
    )
    lat = np.array(
        [
            [0.5, 0.5, 0.5],
            [1.5, 1.5, 1.5],
            [2.5, 2.5, 2.5],
        ]
    )
    return xr.Dataset(
        coords={
            "azimuth": ("azimuth", az),
            "range": ("range", rg),
            "lon": (("azimuth", "range"), lon),
            "lat": (("azimuth", "range"), lat),
        }
    )


@pytest.fixture
def parcel_shapefile_with_id(tmp_path):
    shp_path = tmp_path / "parcel.shp"
    gdf = gpd.GeoDataFrame(
        {"id": [7]},
        geometry=[
            Polygon(
                [
                    (0.0, 0.0),
                    (3.0, 0.0),
                    (3.0, 3.0),
                    (0.0, 3.0),
                    (0.0, 0.0),
                ]
            )
        ],
        crs="EPSG:4326",
    )
    gdf.to_file(
        shp_path,
        driver="ESRI Shapefile",
        index=False,
    )
    return str(shp_path)


@pytest.fixture
def parcel_shapefile_without_id(tmp_path):
    shp_path = tmp_path / "parcel_no_id.shp"
    gdf = gpd.GeoDataFrame(
        {
            "label": ["parcel"],
        },
        geometry=[
            Polygon(
                [
                    (0.0, 0.0),
                    (3.0, 0.0),
                    (3.0, 3.0),
                    (0.0, 3.0),
                    (0.0, 0.0),
                ]
            )
        ],
        crs="EPSG:4326",
    )
    gdf.to_file(
        shp_path,
        driver="ESRI Shapefile",
        index=False,
    )
    return str(shp_path)


def test_assign_parcel_id_assigns_parcels_and_centroids(slc_stack_for_assign_id, parcel_shapefile_with_id):
    """The function should assign the same ID to all pixels inside one polygon and return a centroid."""
    ps_stm = xr.Dataset(
        coords={"azimuth": ("space", np.array([], dtype=float)), "range": ("space", np.array([], dtype=float))}
    )

    ds_mask, ds_id, centroid, az_centroid, rg_centroid = assign_parcel_id(
        slc_stack_for_assign_id,
        parcel_shapefile_with_id,
        ds_min_cells=1,
        ps_stm=ps_stm,
    )

    assert ds_mask.shape == (3, 3)
    np.testing.assert_array_equal(ds_id, np.array([7]))
    np.testing.assert_array_equal(ds_mask, np.full((3, 3), 7.0))
    np.testing.assert_allclose(centroid, np.array([[1.5, 1.5]]))
    assert az_centroid == [1]
    assert rg_centroid == [1]


def test_assign_parcel_id_removes_ps_pixels_from_mask(slc_stack_for_assign_id, parcel_shapefile_with_id):
    """PS pixels already flagged in the STM should be excluded from parcel assignment."""
    ps_stm = xr.Dataset(coords={"azimuth": ("space", np.array([1.0])), "range": ("space", np.array([1.0]))})

    ds_mask, _, _, _, _ = assign_parcel_id(
        slc_stack_for_assign_id,
        parcel_shapefile_with_id,
        ds_min_cells=1,
        ps_stm=ps_stm,
    )

    expected = np.full((3, 3), 7.0)
    expected[1, 1] = np.nan
    np.testing.assert_array_equal(ds_mask, expected)


def test_assign_parcel_id_requires_id_attribute(slc_stack_for_assign_id, parcel_shapefile_without_id):
    """Parsing a parcel shapefile without an id field should fail clearly."""
    ps_stm = xr.Dataset(
        coords={"azimuth": ("space", np.array([], dtype=float)), "range": ("space", np.array([], dtype=float))}
    )

    with pytest.raises(ValueError, match="Parcel shapefile must contain 'id' attribute"):
        assign_parcel_id(
            slc_stack_for_assign_id,
            parcel_shapefile_without_id,
            ds_min_cells=1,
            ps_stm=ps_stm,
        )


@pytest.fixture
def parcel_shapefile_multipolygon_with_id(tmp_path):
    shp_path = tmp_path / "parcel_multipolygon.shp"
    gdf = gpd.GeoDataFrame(
        {"id": [11]},
        geometry=[
            MultiPolygon(
                [
                    Polygon(
                        [
                            (0.0, 0.0),
                            (1.2, 0.0),
                            (1.2, 1.2),
                            (0.0, 1.2),
                            (0.0, 0.0),
                        ]
                    ),
                    Polygon(
                        [
                            (2.0, 2.0),
                            (3.0, 2.0),
                            (3.0, 3.0),
                            (2.0, 3.0),
                            (2.0, 2.0),
                        ]
                    ),
                ]
            )
        ],
        crs="EPSG:4326",
    )
    gdf.to_file(
        shp_path,
        driver="ESRI Shapefile",
        index=False,
    )
    return str(shp_path)


def test_assign_parcel_id_filters_polygons_below_min_cells(slc_stack_for_assign_id, parcel_shapefile_with_id):
    """Parcels below ds_min_cells should be skipped and return empty ids/centroids."""
    ps_stm = xr.Dataset(
        coords={"azimuth": ("space", np.array([], dtype=float)), "range": ("space", np.array([], dtype=float))}
    )

    ds_mask, ds_id, centroid, az_centroid, rg_centroid = assign_parcel_id(
        slc_stack_for_assign_id,
        parcel_shapefile_with_id,
        ds_min_cells=100,
        ps_stm=ps_stm,
    )

    assert ds_mask.shape == (3, 3)
    assert np.isnan(ds_mask).all()
    assert ds_id.size == 0
    assert centroid.size == 0
    assert az_centroid == []
    assert rg_centroid == []


def test_assign_parcel_id_supports_multipolygon_geometry(
    slc_stack_for_assign_id,
    parcel_shapefile_multipolygon_with_id,
):
    """A multipolygon parcel should be parsed and assigned from its first polygon geometry."""
    ps_stm = xr.Dataset(
        coords={"azimuth": ("space", np.array([], dtype=float)), "range": ("space", np.array([], dtype=float))}
    )

    ds_mask, ds_id, centroid, az_centroid, rg_centroid = assign_parcel_id(
        slc_stack_for_assign_id,
        parcel_shapefile_multipolygon_with_id,
        ds_min_cells=1,
        ps_stm=ps_stm,
    )

    np.testing.assert_array_equal(ds_id, np.array([11]))
    assert np.count_nonzero(~np.isnan(ds_mask)) == 1
    np.testing.assert_allclose(centroid, np.array([[0.6, 0.6]]))
    assert az_centroid == [0]
    assert rg_centroid == [0]


@pytest.fixture
def ds_phase_inputs():
    slc_stack = xr.Dataset(
        data_vars={
            "complex": (
                ("azimuth", "range", "time"),
                np.array(
                    [
                        [[1 + 0j, 1 + 0j, 1 + 0j], [1 + 0j, 1 + 0j, 1 + 0j]],
                        [[1 + 0j, 1 + 0j, 1 + 0j], [1 + 0j, 1 + 0j, 1 + 0j]],
                    ],
                    dtype=np.complex64,
                ),
            ),
            "h2ph": (
                ("azimuth", "range", "time"),
                np.array(
                    [
                        [[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],
                        [[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],
                    ],
                    dtype=np.float32,
                ),
            ),
        },
        coords={
            "azimuth": np.array([0, 1]),
            "range": np.array([0, 1]),
            "time": np.array([0, 1, 2]),
        },
    )
    metadata = {
        "pulse_repetition_frequency": 10.0,
        "total_azimuth_bandwidth": 5.0,
        "range_sampling_rate": 20.0,
        "total_range_bandwidth": 8.0,
    }
    ds_ids = np.array([10, 20])
    centroid = np.array([[5.0, 52.0], [5.1, 52.1]])
    az_centroid = np.array([0, 1])
    rg_centroid = np.array([0, 1])
    ds_mask = np.array([[10, 20], [10, 20]], dtype=float)
    parcel_df = pd.DataFrame(
        {
            "id": [10, 20],
            "cropcode": [7, 8],
            "peilgebied": [3, 4],
            "soilcode": [4, 5],
            "meteo_id": [9, 10],
        }
    )
    return {
        "slc_stack": slc_stack,
        "metadata": metadata,
        "stack_id": "stack_001",
        "wavelength": 55.0,
        "centroid": centroid,
        "az_centroid": az_centroid,
        "rg_centroid": rg_centroid,
        "parcel_df": parcel_df,
        "ds_mask": ds_mask,
        "ds_ids": ds_ids,
    }


def test_ds_phase_estimation_returns_datatree_and_fills_outputs(monkeypatch, ds_phase_inputs):
    """Valid parcel inputs should fill ds_stm and ds_cpx_coh nodes in the returned datatree."""

    def fake_coherence_matrix(_slc_stack, _mask, ds_id_value, shp_test=None):
        scale = 1.0 if ds_id_value == 10 else 2.0
        return (
            2,
            30.0,
            np.array([1.0, 4.0, 9.0]) * scale,
            np.array([1.0, 2.0, 3.0]) * scale,
            np.eye(3, dtype=complex),
        )

    monkeypatch.setattr("depsi.ds._coherence_matrix", fake_coherence_matrix)

    result = ds_phase_estimation(**ds_phase_inputs)
    assert isinstance(result, xr.DataTree)

    stm = result["ds_stm"].to_dataset()
    cpx = result["ds_cpx_coh"].to_dataset()

    np.testing.assert_allclose(stm["ds_npixels"].values, np.array([2.0, 2.0]))
    np.testing.assert_allclose(stm["ds_mean_p"].values[0], np.array([1.0, 4.0, 9.0]))
    np.testing.assert_allclose(stm["ds_mean_p"].values[1], np.array([2.0, 8.0, 18.0]))
    np.testing.assert_allclose(cpx["ds_cpx_coh"].values[0], np.eye(3, dtype=complex))
    np.testing.assert_allclose(cpx["ds_cpx_coh"].values[1], np.eye(3, dtype=complex))
    np.testing.assert_array_equal(stm["crop_id"].values, np.array([7.0, 8.0]))
    np.testing.assert_array_equal(stm["meteo_id"].values, np.array([9.0, 10.0]))


def test_ds_phase_estimation_respects_ds_id_list_subset(monkeypatch, ds_phase_inputs):
    """When ds_id_list is provided, only the selected parcel ids should be recomputed."""
    calls = []

    def fake_coherence_matrix(_slc_stack, _mask, ds_id_value, shp_test=None):
        calls.append(int(ds_id_value))
        return 2, 30.0, np.array([1.0, 4.0, 9.0]), np.array([1.0, 2.0, 3.0]), np.eye(3, dtype=complex)

    monkeypatch.setattr("depsi.ds._coherence_matrix", fake_coherence_matrix)

    result = ds_phase_estimation(**ds_phase_inputs, ds_id_list=[20])
    stm = result["ds_stm"].to_dataset()

    assert calls == [20]
    assert np.isnan(stm["ds_npixels"].values[0])
    assert stm["ds_npixels"].values[1] == 2


def test_ds_phase_estimation_requires_id_field_in_parcel_dataframe(monkeypatch, ds_phase_inputs):
    """A parcel dataframe without an id column should fail clearly during STM fill."""

    def fake_coherence_matrix(_slc_stack, _mask, _ds_id_value, shp_test=None):
        return 2, 30.0, np.array([1.0, 4.0, 9.0]), np.array([1.0, 2.0, 3.0]), np.eye(3, dtype=complex)

    monkeypatch.setattr("depsi.ds._coherence_matrix", fake_coherence_matrix)
    parcel_df_no_id = ds_phase_inputs["parcel_df"].drop(columns=["id"])
    params = dict(ds_phase_inputs)
    params["parcel_df"] = parcel_df_no_id

    with pytest.raises(ValueError, match="Parcel attributes must contain either 'id' attribute"):
        ds_phase_estimation(**params)


def test_ds_phase_estimation_raises_when_required_attrs_missing(tmp_path):
    """Loading a preexisting datatree without required STM attrs should fail fast."""
    ds_stm_missing_attrs = xr.Dataset(
        data_vars={"ds_id": ("space", np.array([10]))},
        coords={"space": np.array([0]), "time": np.array([0, 1, 2])},
        attrs={
            "stack_id": "stack_001",
            "wavelength": 55.0,
            "prf": 10.0,
            "az_bw": 5.0,
            "r_fs": 20.0,
            # r_bw intentionally missing
        },
    )
    cpx = xr.Dataset(
        data_vars={"ds_cpx_coh": (("space", "time1", "time2"), np.zeros((1, 3, 3), dtype=np.complex64))},
        coords={"space": [0], "time1": [0, 1, 2], "time2": [0, 1, 2]},
    )
    ds_path = tmp_path / "ds_dtree_stack_001.zarr"
    xr.DataTree.from_dict({"ds_stm": ds_stm_missing_attrs, "ds_cpx_coh": cpx}).to_zarr(ds_path, mode="w")

    slc_stack = xr.Dataset(
        data_vars={"complex": (("azimuth", "range", "time"), np.ones((1, 1, 3), dtype=np.complex64))},
        coords={"azimuth": [0], "range": [0], "time": [0, 1, 2]},
    )
    metadata = {
        "pulse_repetition_frequency": 10.0,
        "total_azimuth_bandwidth": 5.0,
        "range_sampling_rate": 20.0,
        "total_range_bandwidth": 8.0,
    }

    with pytest.raises(ValueError, match="Missing required attribute 'r_bw' in stm"):
        ds_phase_estimation(
            slc_stack=slc_stack,
            metadata=metadata,
            stack_id="stack_001",
            wavelength=55.0,
            centroid=np.array([[5.0, 52.0]]),
            az_centroid=np.array([0]),
            rg_centroid=np.array([0]),
            parcel_df=pd.DataFrame({"id": [10], "cropcode": [1], "peilgebied": [1], "soilcode": [1], "meteo_id": [1]}),
            ds_mask=np.array([[10]], dtype=float),
            ds_ids=np.array([10]),
            ds_filepath=str(ds_path),
        )


@pytest.fixture
def simple_slc_stack():
    """
    Stack dimensions:
      azimuth: 2
      range:   3
      time:    2

    The `complex` variable has shape:
      (azimuth, range, time)

    Selected parcel pixels are:
      (0, 0): [1+0j, 1+0j]
      (0, 1): [1+0j, 1+0j]
      (1, 0): [nan+nanj, nan+nanj]  -> should be removed
    """
    cpx_data = np.array(
        [
            [
                [1.0 + 0.0j, 1.0 + 0.0j],
                [1.0 + 0.0j, 1.0 + 0.0j],
                [2.0 + 0.0j, 2.0 + 0.0j],
            ],
            [
                [np.nan + 1j * np.nan, np.nan + 1j * np.nan],
                [3.0 + 0.0j, 3.0 + 0.0j],
                [4.0 + 0.0j, 4.0 + 0.0j],
            ],
        ],
        dtype=np.complex128,
    )

    incidence_angle = np.array(
        [
            [30.0, 31.0, 32.0],
            [33.0, 34.0, 35.0],
        ]
    )

    return xr.Dataset(
        data_vars={
            "complex": (
                ("azimuth", "range", "time"),
                cpx_data,
            ),
            "incidence_angle": (
                ("azimuth", "range"),
                incidence_angle,
            ),
        },
        coords={
            "azimuth": [0, 1],
            "range": [0, 1, 2],
            "time": [0, 1],
        },
    )


def test_coherence_matrix_selected_pixels_and_statistics(simple_slc_stack):
    ds_id = 10

    ds_mask = np.array(
        [
            [10, 10, 0],
            [10, 0, 0],
        ]
    )

    npixels, mean_ia, mean_p, mean_amp, cpx_coh = _coherence_matrix(
        simple_slc_stack,
        ds_mask=ds_mask,
        ds_id=ds_id,
        shp_test=None,
    )

    # Three mask pixels were selected, but one is all-NaN and must be removed.
    assert npixels == 2

    # Note: mean_ia is currently calculated from all selected mask pixels,
    # including the all-NaN complex pixel at location (1, 0).
    assert_allclose(mean_ia, np.mean([30.0, 31.0, 33.0]))

    # Two retained complex rows are both [1+0j, 1+0j].
    expected_mean_p = np.array([1.0, 1.0])
    expected_mean_amp = np.array([1.0, 1.0])

    assert_allclose(mean_p, expected_mean_p)
    assert_allclose(mean_amp, expected_mean_amp)

    # Both acquisitions are perfectly coherent.
    expected_coh = np.array(
        [
            [1.0 + 0.0j, 1.0 + 0.0j],
            [1.0 + 0.0j, 1.0 + 0.0j],
        ]
    )

    assert cpx_coh.shape == (2, 2)
    assert_allclose(cpx_coh, expected_coh, atol=1e-12)


def test_coherence_matrix_without_incidence_angle():
    cpx_data = np.array(
        [
            [
                [1.0 + 0j, 1.0 + 0j],
                [1.0 + 0j, 1.0 + 0j],
            ]
        ]
    )

    slc_stack = xr.Dataset(
        data_vars={
            "complex": (("azimuth", "range", "time"), cpx_data),
        },
        coords={
            "azimuth": [0],
            "range": [0, 1],
            "time": [0, 1],
        },
    )

    ds_mask = np.array([[7, 7]])

    npixels, mean_ia, mean_p, mean_amp, cpx_coh = _coherence_matrix(
        slc_stack,
        ds_mask=ds_mask,
        ds_id=7,
        shp_test=None,
    )

    assert npixels == 2
    assert np.isnan(mean_ia)
    assert_allclose(mean_p, [1.0, 1.0])
    assert_allclose(mean_amp, [1.0, 1.0])
    assert_allclose(cpx_coh, np.ones((2, 2), dtype=complex))


def test_coherence_matrix_no_matching_parcel_returns_empty_statistics():
    cpx_data = np.ones((2, 2, 2), dtype=np.complex128)

    slc_stack = xr.Dataset(
        data_vars={
            "complex": (("azimuth", "range", "time"), cpx_data),
        }
    )

    ds_mask = np.array(
        [
            [1, 1],
            [1, 1],
        ]
    )

    npixels, mean_ia, mean_p, mean_amp, cpx_coh = _coherence_matrix(
        slc_stack,
        ds_mask=ds_mask,
        ds_id=999,
        shp_test=None,
    )

    assert npixels == 0
    assert np.isnan(mean_ia)

    # np.nanmean over an empty selection produces NaNs.
    assert np.all(np.isnan(mean_p))
    assert np.all(np.isnan(mean_amp))

    # No acquisitions are specified in the test input interpretation;
    # coherence retains the temporal dimensions: 2 x 2.
    assert cpx_coh.shape == (2, 2)
    assert np.all(np.isnan(cpx_coh))


def expected_emi(data, regularization=0):
    """Independent reference implementation of the current EMI branch."""
    cleaned = np.nan_to_num(data, nan=0.0, posinf=0.0, neginf=0.0)

    if regularization == 1:
        beta = 0.5
        cleaned = (1 - beta) * cleaned + beta * np.eye(cleaned.shape[0])

    u, s, _ = np.linalg.svd(np.linalg.pinv(np.abs(cleaned)) * cleaned)
    return u[:, -1:] * s[-1] @ np.conj(u[:, -1:].T)


@pytest.fixture
def coherence_matrix():
    """Small deterministic complex coherence-like matrix."""
    return np.array(
        [
            [1.0 + 0.0j, 0.70 + 0.20j, 0.40 - 0.10j],
            [0.70 - 0.20j, 1.0 + 0.0j, 0.60 + 0.30j],
            [0.40 + 0.10j, 0.60 - 0.30j, 1.0 + 0.0j],
        ],
        dtype=np.complex128,
    )


def test_phase_linking_emi_matches_expected_result(coherence_matrix):
    """The default estimator should produce the EMI result."""
    result = _phase_linking(coherence_matrix)
    expected = expected_emi(coherence_matrix)

    assert result.shape == coherence_matrix.shape
    assert np.iscomplexobj(result)
    assert_allclose(result, expected, rtol=1e-12, atol=1e-12)


def test_phase_linking_empty_estimator_uses_emi(coherence_matrix):
    """An empty estimator string should use the EMI branch."""
    result = _phase_linking(coherence_matrix, estimator="")
    expected = _phase_linking(coherence_matrix, estimator="emi")

    assert_allclose(result, expected, rtol=1e-12, atol=1e-12)


def test_phase_linking_regularization_matches_expected_result(coherence_matrix):
    """regularization=1 should apply beta=0.5 spectral regularization."""
    result = _phase_linking(coherence_matrix, regularization=1)
    expected = expected_emi(coherence_matrix, regularization=1)

    assert_allclose(result, expected, rtol=1e-12, atol=1e-12)

    unregularized = _phase_linking(coherence_matrix, regularization=0)
    assert not np.allclose(result, unregularized)


def test_phase_linking_replaces_nan_and_inf_with_zero(coherence_matrix):
    """NaN and infinite values should be handled exactly as zero-valued inputs."""
    dirty_data = coherence_matrix.copy()
    dirty_data[0, 1] = np.nan + 0.0j
    dirty_data[1, 0] = np.inf + 0.0j
    dirty_data[2, 1] = -np.inf + 0.0j

    cleaned_data = np.nan_to_num(
        dirty_data,
        nan=0.0,
        posinf=0.0,
        neginf=0.0,
    )

    result_dirty = _phase_linking(dirty_data)
    result_clean = _phase_linking(cleaned_data)

    assert np.all(np.isfinite(result_dirty))
    assert_allclose(result_dirty, result_clean, rtol=1e-12, atol=1e-12)


@pytest.mark.parametrize("estimator", ["sequential", "eig", "invalid", "EMI"])
def test_phase_linking_rejects_unsupported_estimators(
    coherence_matrix,
    estimator,
):
    """Only lowercase 'emi' and an empty string are currently supported."""
    with pytest.raises(
        NotImplementedError,
        match="This module is not yet implemented.",
    ):
        _phase_linking(coherence_matrix, estimator=estimator)


def coherence_matrix_from_daisy_chain(coh_dc):
    """
    Create a square test matrix whose first upper diagonal contains coh_dc.

    Parameters
    ----------
    coh_dc : sequence of float
        Coherence values between adjacent acquisitions. For N acquisitions,
        provide N - 1 values.

    Returns
    -------
    np.ndarray
        Square coherence matrix of shape (N, N).
    """
    n = len(coh_dc) + 1
    data = np.eye(n, dtype=float)

    coh_dc = np.asarray(coh_dc, dtype=float)
    row_idx = np.arange(n - 1)

    # Upper and lower first off-diagonals
    data[row_idx, row_idx + 1] = coh_dc
    data[row_idx + 1, row_idx] = coh_dc

    return data


def test_segmentation_no_cuts_returns_full_block():
    """All daisy-chain coherence values exceed threshold."""
    data = coherence_matrix_from_daisy_chain([0.8] * 14)

    nsegments, blocks_idx = _segmentation(
        data,
        min_seg_len=10,
        threshold=0.2,
    )

    assert nsegments == 1
    assert blocks_idx == [(0, 15)]


def test_segmentation_finds_two_valid_blocks():
    """
    A coherence <= threshold between positions 9 and 10 produces a cut.

    With 20 acquisitions:
      - first block: indices 0..9, returned as (0, 9)
      - second block: expected to start at index 10

    Note: this test reflects the CURRENT implementation, which does not
    append the final segment after the last cut.
    """
    coh_dc = np.full(19, 0.8)
    coh_dc[9] = 0.2  # threshold equality must create a cut

    data = coherence_matrix_from_daisy_chain(coh_dc)

    nsegments, blocks_idx = _segmentation(
        data,
        min_seg_len=10,
        threshold=0.2,
    )

    # Current function appends the first valid block only.
    assert nsegments == 1
    assert blocks_idx == [(0, 9)]


def test_segmentation_nan_is_treated_as_low_coherence():
    """
    NaN is replaced by 1e-8, which is below the default threshold,
    and thus causes a cut.
    """
    coh_dc = np.full(14, 0.8)
    coh_dc[9] = np.nan

    data = coherence_matrix_from_daisy_chain(coh_dc)

    nsegments, blocks_idx = _segmentation(
        data,
        min_seg_len=10,
        threshold=0.2,
    )

    assert nsegments == 1
    assert blocks_idx == [(0, 9)]


@pytest.mark.parametrize(
    ("coherence", "threshold", "expected_nsegments", "expected_blocks"),
    [
        # Value equal to threshold: cut because comparison is <=.
        (0.2, 0.2, 1, [(0, 9)]),
        # Value above threshold: no cut, one whole block.
        (0.200001, 0.2, 1, [(0, 20)]),
    ],
)
def test_segmentation_threshold_boundary(
    coherence,
    threshold,
    expected_nsegments,
    expected_blocks,
):
    coh_dc = np.full(19, 0.8)
    coh_dc[9] = coherence
    data = coherence_matrix_from_daisy_chain(coh_dc)

    nsegments, blocks_idx = _segmentation(
        data,
        min_seg_len=10,
        threshold=threshold,
    )

    assert nsegments == expected_nsegments
    assert blocks_idx == expected_blocks


def test_segmentation_discards_block_shorter_than_minimum():
    """
    A cut after index 4 makes the first candidate segment too short.
    The present implementation also does not append the final segment.
    """
    coh_dc = np.full(19, 0.8)
    coh_dc[4] = 0.1

    data = coherence_matrix_from_daisy_chain(coh_dc)

    nsegments, blocks_idx = _segmentation(
        data,
        min_seg_len=10,
        threshold=0.2,
    )

    assert nsegments == 0
    assert blocks_idx == []


def make_ps_stm(
    x,
    y,
    *,
    lon=None,
    lat=None,
    quality=None,
    h2ph=None,
    sd_phase=None,
):
    """Create a minimal PS STM-like Dataset for unit testing."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n_space = len(x)

    if lon is None:
        lon = x + 4.0
    if lat is None:
        lat = y + 52.0
    if quality is None:
        quality = np.ones(n_space)
    if h2ph is None:
        h2ph = np.ones((n_space, 2))
    if sd_phase is None:
        sd_phase = np.zeros((n_space, 2))

    return xr.Dataset(
        data_vars={
            "full_ts_nad": ("space", np.asarray(quality, dtype=float)),
            "h2ph": (("space", "time"), np.asarray(h2ph, dtype=float)),
            "sd_phase": (("space", "time"), np.asarray(sd_phase, dtype=float)),
        },
        coords={
            "space": np.arange(n_space),
            "time": np.array(["2024-01-01", "2024-01-02"], dtype="datetime64[ns]"),
            "lon": ("space", np.asarray(lon, dtype=float)),
            "lat": ("space", np.asarray(lat, dtype=float)),
            "x_euclidean_proj_epsg28992": ("space", x),
            "y_euclidean_proj_epsg28992": ("space", y),
        },
        attrs={"wavelength": 0.056},
    )


def make_ds_dtree(ds_stm, ds_cpx_coh=None):
    """Wrap a DS STM Dataset in the DataTree structure expected by ps_ds_arc."""
    if ds_cpx_coh is None:
        ds_cpx_coh = xr.Dataset({"coherence": (("space", "time"), np.ones_like(ds_stm["sd_phase"].values))})

    return DataTree.from_dict(
        {
            "ds_stm": ds_stm,
            "ds_cpx_coh": ds_cpx_coh,
        }
    )


def test_select_common_fop_ref_returns_common_points_and_best_reference():
    # All three tracks have points within 20 m of each other.
    # The second candidate has the smallest summed quality:
    # 5 + 1 + 2 = 8, versus 10 + 8 + 9 = 27 for the first candidate.
    track_1 = make_ps_stm(
        x=[0.0, 100.0],
        y=[0.0, 100.0],
        quality=[10.0, 5.0],
    )
    track_2 = make_ps_stm(
        x=[2.0, 102.0],
        y=[1.0, 98.0],
        quality=[8.0, 1.0],
    )
    track_3 = make_ps_stm(
        x=[-2.0, 99.0],
        y=[1.0, 101.0],
        quality=[9.0, 2.0],
    )

    ref_idxs, pnt_idx_candidates = select_common_fop_ref(
        ps_stm_list=[track_1, track_2, track_3],
        proj_crs=28992,
        dist_ub=20.0,
    )

    expected_candidates = np.array(
        [
            [0, 1],
            [0, 1],
            [0, 1],
        ]
    )

    np.testing.assert_array_equal(pnt_idx_candidates, expected_candidates)
    np.testing.assert_array_equal(ref_idxs, np.array([1, 1, 1]))


def test_select_common_fop_ref_restores_original_indices_after_nan_filtering():
    # The middle point of track_1 has invalid geographic coordinates.
    # It is removed before matching, but returned indices must refer to
    # the original unfiltered Dataset.
    track_1 = make_ps_stm(
        x=[0.0, 50.0, 100.0],
        y=[0.0, 50.0, 100.0],
        lon=[4.0, np.nan, 4.1],
        lat=[52.0, 52.1, 52.2],
        quality=[10.0, 999.0, 1.0],
    )
    track_2 = make_ps_stm(
        x=[1.0, 101.0],
        y=[1.0, 99.0],
        quality=[5.0, 2.0],
    )

    ref_idxs, pnt_idx_candidates = select_common_fop_ref(
        ps_stm_list=[track_1, track_2],
        proj_crs=28992,
        dist_ub=10.0,
    )

    expected_candidates = np.array(
        [
            [0, 2],
            [0, 1],
        ]
    )

    np.testing.assert_array_equal(pnt_idx_candidates, expected_candidates)
    np.testing.assert_array_equal(ref_idxs, np.array([2, 1]))


def test_select_common_fop_ref_raises_when_no_points_are_common():
    track_1 = make_ps_stm(x=[0.0], y=[0.0])
    track_2 = make_ps_stm(x=[1000.0], y=[1000.0])

    with pytest.raises(ValueError, match="No common PS was found"):
        select_common_fop_ref(
            ps_stm_list=[track_1, track_2],
            proj_crs=28992,
            dist_ub=20.0,
        )


def test_ps_ds_arc_builds_expected_arc_dataset(monkeypatch):
    ps_stm = make_ps_stm(
        x=[0.0, 10.0, 20.0],
        y=[0.0, 10.0, 20.0],
        h2ph=[
            [2.0, 4.0],
            [10.0, 14.0],
            [100.0, 200.0],
        ],
        sd_phase=[
            [0.0, 0.0],
            [0.2, -0.2],
            [0.0, 0.0],
        ],
    )

    ds_stm = make_ps_stm(
        x=[1.0, 11.0],
        y=[1.0, 11.0],
        h2ph=[
            [6.0, 8.0],
            [14.0, 18.0],
        ],
        sd_phase=[
            [0.1, 0.3],
            [-0.3, 0.4],
        ],
    )

    ds_cpx_coh = xr.Dataset(
        {
            "coh": (
                ("space", "time"),
                np.array([[0.9, 0.8], [0.7, 0.6]]),
            )
        }
    )
    ds_dtree = make_ds_dtree(ds_stm, ds_cpx_coh)

    def fake_determine_connections(*args, **kwargs):
        return np.array([0, 1]), np.array([1, 0])

    monkeypatch.setattr(
        "depsi.ds._determine_dens_connections",
        fake_determine_connections,
    )

    result = ps_ds_arc(
        ps_stm=ps_stm,
        ds_dtree=ds_dtree,
        n_connections=1,
    )

    arc_stm = result["arc_stm"].to_dataset()

    assert set(result.children) == {"arc_stm", "ps_stm", "ds_stm", "ds_cpx_coh"}
    assert arc_stm.attrs["wavelength"] == pytest.approx(0.056)

    np.testing.assert_array_equal(arc_stm["idx_dens"].values, [0, 1])
    np.testing.assert_array_equal(arc_stm["idx_network"].values, [1, 0])

    # Arc h2ph is the average of connected DS and PS h2ph values.
    expected_h2ph = np.array(
        [
            [(6.0 + 10.0) / 2, (8.0 + 14.0) / 2],
            [(14.0 + 2.0) / 2, (18.0 + 4.0) / 2],
        ]
    )
    np.testing.assert_allclose(arc_stm["h2ph"].values, expected_h2ph)

    # dd_phase = wrapped(ds_phase - ps_phase), constrained to [-pi, pi).
    expected_dd_phase = np.array(
        [
            [0.1 - 0.2, 0.3 - (-0.2)],
            [-0.3 - 0.0, 0.4 - 0.0],
        ]
    )
    np.testing.assert_allclose(arc_stm["dd_phase"].values, expected_dd_phase)


def test_ps_ds_arc_wraps_double_difference_phase(monkeypatch):
    ps_stm = make_ps_stm(
        x=[0.0],
        y=[0.0],
        h2ph=[[1.0, 1.0]],
        sd_phase=[[3.0, -3.0]],
    )
    ds_stm = make_ps_stm(
        x=[1.0],
        y=[1.0],
        h2ph=[[3.0, 3.0]],
        sd_phase=[[-3.0, 3.0]],
    )

    def fake_determine_connections(*args, **kwargs):
        return np.array([0]), np.array([0])

    monkeypatch.setattr(
        "depsi.ds._determine_dens_connections",
        fake_determine_connections,
    )

    result = ps_ds_arc(ps_stm, make_ds_dtree(ds_stm))
    dd_phase = result["arc_stm"].to_dataset()["dd_phase"].values

    raw_difference = np.array([[-3.0 - 3.0, 3.0 - (-3.0)]])
    expected = (raw_difference + np.pi) % (2 * np.pi) - np.pi

    np.testing.assert_allclose(dd_phase, expected)
    assert np.all(dd_phase >= -np.pi)
    assert np.all(dd_phase < np.pi)


def test_ps_ds_arc_rejects_multiple_connections():
    ps_stm = make_ps_stm(x=[0.0], y=[0.0])
    ds_stm = make_ps_stm(x=[1.0], y=[1.0])

    with pytest.raises(AssertionError, match="n_connections=1"):
        ps_ds_arc(
            ps_stm=ps_stm,
            ds_dtree=make_ds_dtree(ds_stm),
            n_connections=2,
        )


def test_ps_ds_arc_raises_for_missing_ps_sd_phase():
    ps_stm = make_ps_stm(x=[0.0], y=[0.0]).drop_vars("sd_phase")
    ds_stm = make_ps_stm(x=[1.0], y=[1.0])

    with pytest.raises(ValueError, match="Missing required variable 'sd_phase'"):
        ps_ds_arc(
            ps_stm=ps_stm,
            ds_dtree=make_ds_dtree(ds_stm),
        )


def test_ps_ds_arc_raises_for_missing_ds_sd_phase():
    ps_stm = make_ps_stm(x=[0.0], y=[0.0])
    ds_stm = make_ps_stm(x=[1.0], y=[1.0])
    ds_dtree = make_ds_dtree(ds_stm)
    ds_stm = ds_stm.drop_vars("sd_phase")
    ds_dtree["ds_stm"] = ds_stm

    with pytest.raises(ValueError, match="Missing required variable 'sd_phase'"):
        ps_ds_arc(
            ps_stm=ps_stm,
            ds_dtree=ds_dtree,
        )


class DummyPool:
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        return False

    def map(self, function, iterable):
        return [function(item) for item in iterable]


def test_shp_test_keeps_pixels_matching_selected_reference(monkeypatch):
    """Test deterministic SHP selection.

    We prescribe a 3x3 p-value matrix:

        [[0.90, 0.80, 0.01],
         [0.80, 0.90, 0.01],
         [0.01, 0.01, 0.90]]

    With threshold p > 0.05:
      - pixel 0 matches pixels 0 and 1: total = 2
      - pixel 1 matches pixels 0 and 1: total = 2
      - pixel 2 matches only itself   : total = 1

    np.argmax selects the first maximum: row 0.
    Thus, pixels at original positions 0 and 1 are retained.
    """
    data = np.array(
        [
            3.0 + 4.0j,  # magnitude = 5
            1.0 + 0.0j,  # magnitude = 1
            2.0 + 0.0j,  # magnitude = 2
        ]
    )

    pvalue_lookup = {
        (1.0, 1.0): 0.90,
        (1.0, 2.0): 0.80,
        (1.0, 5.0): 0.01,
        (2.0, 1.0): 0.80,
        (2.0, 2.0): 0.90,
        (2.0, 5.0): 0.01,
        (5.0, 1.0): 0.01,
        (5.0, 2.0): 0.01,
        (5.0, 5.0): 0.90,
    }

    def fake_kstest(pair):
        sample_1, sample_2 = pair
        key = (float(sample_1), float(sample_2))
        return pvalue_lookup[key]

    monkeypatch.setattr("depsi.ds.mp.Pool", DummyPool)
    monkeypatch.setattr("depsi.ds._kstest", fake_kstest)

    result = _shp_test(data, method="ks-test")

    expected = np.array(
        [
            1.0 + 0.0j,  # magnitude = 1
            2.0 + 0.0j,  # magnitude = 2
        ]
    )

    np.testing.assert_array_equal(result, expected)


def test_shp_test_empty_method_uses_ks_test_branch(monkeypatch):
    """method='' is documented/implemented as an alias for 'ks-test'.

    All p-values are > 0.05, so all pixels remain selected.
    """
    data = np.array(
        [
            3.0 + 4.0j,  # magnitude = 5
            1.0 + 0.0j,  # magnitude = 1
            2.0 + 0.0j,  # magnitude = 2
        ]
    )

    monkeypatch.setattr("depsi.ds.mp.Pool", DummyPool)
    monkeypatch.setattr("depsi.ds._kstest", lambda pair: 1.0)

    result = _shp_test(data, method="")
    expected = data[np.argsort(np.abs(data))]

    np.testing.assert_array_equal(result, expected)


def test_shp_test_rejects_unknown_method():
    """An unimplemented SHP method should raise the documented error."""
    data = np.array(
        [
            3.0 + 4.0j,  # magnitude = 5
            1.0 + 0.0j,  # magnitude = 1
            2.0 + 0.0j,  # magnitude = 2
        ]
    )

    with pytest.raises(
        NotImplementedError,
        match="This module is not yet implemented",
    ):
        _shp_test(data, method="anderson-darling")
