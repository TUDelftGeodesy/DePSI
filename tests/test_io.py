import csv
import os
from pathlib import Path

import dask.array as da
import numpy as np
import pandas as pd
import pytest
import xarray as xr

from depsi import io


def test_read_metadata_lines_pixels():
    with pytest.warns(DeprecationWarning):
        metadata = io.read_metadata("tests/data/example.res")
        assert metadata["n_lines"] == 14065
        assert metadata["n_pixels"] == 49843


def test_get_targets_from_slc():
    # Grid dimensions like the real dataset
    az = np.arange(65)
    rg = np.arange(222)
    times = np.array([np.datetime64(f"2025-01-{i + 1:02d}") for i in range(5)])

    # Use dask arrays and correct types
    slc_stack = xr.Dataset(
        data_vars={
            "h2ph": (("azimuth", "range", "time"), da.ones((65, 222, 5), dtype=np.float32)),
            "lat": (("azimuth", "range"), da.linspace(53.0, 53.1, 65 * 222, dtype=np.float32).reshape(65, 222)),
            "lon": (("azimuth", "range"), da.linspace(6.0, 6.1, 65 * 222, dtype=np.float32).reshape(65, 222)),
            "complex": (("azimuth", "range", "time"), da.ones((65, 222, 5), dtype=np.complex64)),
            "amplitude": (("azimuth", "range", "time"), da.ones((65, 222, 5), dtype=np.float32)),
            "phase": (("azimuth", "range", "time"), da.zeros((65, 222, 5), dtype=np.float32)),
        },
        coords={"azimuth": az, "range": rg, "time": times},
    )

    # 3 targets (arbitrary positions within grid)
    n_targets = 3
    target_az = np.array([1, 20, 50])
    target_rg = np.array([3, 100, 200])

    targets = xr.Dataset(
        data_vars={
            "detection_flag": (("space", "time"), np.random.randint(0, 2, size=(n_targets, len(times)))),
        },
        coords={
            "space": np.arange(n_targets),
            "time": times,
            "azimuth": ("space", target_az),
            "range": ("space", target_rg),
            "target": ("space", ["A", "B", "C"]),
            "lat": ("space", [53.01, 53.05, 53.07]),
            "lon": ("space", [6.01, 6.05, 6.07]),
            "height": ("space", [40.0, 42.0, 41.0]),
        },
    )

    # Call function
    res = io.get_targets_from_slc(slc_stack, targets)

    # Assertions
    assert res.sizes["space"] == n_targets
    assert res.sizes["time"] == len(times)
    assert "detection_flag" in res
    for var in ["h2ph", "lat", "lon", "complex", "amplitude", "phase"]:
        assert var in res


def test_read_rcs_csv():
    # Define the CSV path inside the test
    data_path = os.path.join(os.path.dirname(__file__), "data", "s1_dsc037_RC.csv")

    # Read CSV using the function under test
    dataset = io.read_rcs_csv(data_path)

    # Check that the '*****' row existed in the original file
    with open(data_path) as f:
        content = f.read()
    assert "*****" in content, "Metadata markers '*****' not found in CSV file"

    assert isinstance(dataset, xr.Dataset)
    assert set(dataset.dims) == {"space", "time"}
    assert dataset.sizes["space"] == 24
    assert dataset.sizes["time"] == 422

    expected_coords = {"space", "time", "target", "azimuth_subpixel", "range_subpixel"}
    assert expected_coords.issubset(dataset.coords)

    expected_data_vars = {"lat", "lon", "height", "validation", "existing_flag"}
    assert expected_data_vars.issubset(dataset.data_vars)

    for coord in ["azimuth_subpixel", "range_subpixel"]:
        assert dataset[coord].dims == ("space",)
        assert np.issubdtype(dataset[coord].dtype, np.floating)

    for var in ["lat", "lon", "height"]:
        assert dataset[var].dims == ("space",)
        assert np.issubdtype(dataset[var].dtype, np.floating)

    assert dataset["validation"].dims == ("space",)
    assert np.isin(dataset["validation"].values, [-1, 0, 1]).all()

    assert dataset["existing_flag"].dims == ("space", "time")
    assert np.isin(dataset["existing_flag"].values, [0, 1]).all()
    assert np.issubdtype(dataset["time"].dtype, np.datetime64)


def test_read_knmi_txt():
    filepath = os.path.join(os.path.dirname(__file__), "data", "knmi", "etmgeg_280.txt")

    df = io.read_knmi_txt(filepath)

    required_cols = ["meteo_id", "datum", "pr", "pet"]
    for col in required_cols:
        assert col in df.columns, f"Required column '{col}' missing from DataFrame"

    assert pd.api.types.is_numeric_dtype(df["meteo_id"]), "Column 'meteo_id' is not numeric"
    assert pd.api.types.is_datetime64_any_dtype(df["datum"]), "Column 'datum' is not datetime"
    assert pd.api.types.is_numeric_dtype(df["pr"]), "Column 'pr' is not numeric"
    assert pd.api.types.is_numeric_dtype(df["pet"]), "Column 'pet' is not numeric"


def test_read_slc_stack_zarr_engine_uses_open_zarr_and_from_dataset(monkeypatch):
    opened = xr.Dataset({"h2ph": (("azimuth", "range", "time"), np.ones((2, 3, 1), dtype=np.float32))})
    converted = xr.Dataset({"complex": (("azimuth", "range", "time"), np.ones((2, 3, 1), dtype=np.complex64))})

    calls = {}

    def fake_open_zarr(filename):
        calls["open_zarr_filename"] = filename
        return opened

    def fake_from_dataset(dataset):
        calls["from_dataset_input"] = dataset
        return converted

    monkeypatch.setattr(io.xr, "open_zarr", fake_open_zarr)
    monkeypatch.setattr(io.sarxarray, "from_dataset", fake_from_dataset)

    result = io.read_slc_stack("dummy_stack.zarr", engine="zarr")

    assert result is converted
    assert calls["open_zarr_filename"] == "dummy_stack.zarr"
    assert calls["from_dataset_input"] is opened


def test_read_slc_stack_doris_requires_nlines_and_npixels():
    with pytest.raises(ValueError, match="'nlines' and 'npixels' must be provided"):
        io.read_slc_stack("dummy_stack_dir", engine="doris")


def test_read_slc_stack_doris_raises_when_no_slc_files(monkeypatch):
    monkeypatch.setattr(io, "glob", lambda pattern: [])

    with pytest.raises(FileNotFoundError, match=r"No SLC files found in .*matching pattern .*/slc_srd.raw"):
        io.read_slc_stack("dummy_stack_dir", engine="doris", nlines=10, npixels=20)


def test_read_slc_stack_doris_calls_from_binary_with_expected_arguments(monkeypatch):
    stack_files = ["/tmp/acq1/slc_srd.raw", "/tmp/acq2/slc_srd.raw"]
    expected = xr.Dataset({"complex": (("azimuth", "range", "time"), np.ones((2, 3, 2), dtype=np.complex64))})
    calls = {}

    monkeypatch.setattr(io, "glob", lambda pattern: stack_files)

    def fake_from_binary(file_list, shape, dtype, chunks):
        calls["file_list"] = file_list
        calls["shape"] = shape
        calls["dtype"] = dtype
        calls["chunks"] = chunks
        return expected

    monkeypatch.setattr(io.sarxarray, "from_binary", fake_from_binary)

    result = io.read_slc_stack(Path("dummy_stack_dir"), engine="doris", nlines=65, npixels=222, chunks=(64, 64))

    assert result is expected
    assert calls["file_list"] == stack_files
    assert calls["shape"] == (65, 222)
    assert calls["dtype"] == np.complex64
    assert calls["chunks"] == (64, 64)


def test_read_slc_stack_doris_wraps_from_binary_errors(monkeypatch):
    monkeypatch.setattr(io, "glob", lambda pattern: ["/tmp/acq1/slc_srd.raw"])

    def fake_from_binary(*args, **kwargs):
        raise RuntimeError("boom")

    monkeypatch.setattr(io.sarxarray, "from_binary", fake_from_binary)

    with pytest.raises(RuntimeError, match="Failed to load the SLC stack"):
        io.read_slc_stack("dummy_stack_dir", engine="doris", nlines=10, npixels=10)


def test_read_slc_stack_rejects_unsupported_engine():
    with pytest.raises(ValueError, match="Unsupported engine 'netcdf'. Use 'zarr' or 'doris'."):
        io.read_slc_stack("dummy_stack", engine="netcdf")


def _make_minimal_stm_for_csv_export():
    space = np.arange(2)
    time = np.array([np.datetime64("2025-01-01"), np.datetime64("2025-01-02")])

    return xr.Dataset(
        data_vars={
            "azimuth": ("space", np.array([10, 20], dtype=np.int32)),
            "range": ("space", np.array([100, 200], dtype=np.int32)),
            "lat": ("space", np.array([53.123456789, 53.223456789])),
            "lon": ("space", np.array([6.123456789, 6.223456789])),
            "h": ("space", np.array([40.1234, 41.5678])),
            "stc": ("space", np.array([1.23456, 2.34567])),
            "amplitude": (("space", "time"), np.array([[10.1234, 11.9876], [12.1111, 13.2222]])),
            "unwrapped_phase": (("space", "time"), np.array([[0.123456, 0.654321], [1.111111, 1.999999]])),
            "linear": ("space", np.array([3.141592, 2.718281])),
            "unwrapped_phase_pov": (("space", "time"), np.array([[9.100001, 9.200002], [8.300003, 8.400004]])),
            "linear_pov": ("space", np.array([7.555551, 6.444449])),
            "x_euclidean_proj_epsg28992": ("space", np.array([120000.123, 120010.987])),
            "y_euclidean_proj_epsg28992": ("space", np.array([480000.321, 480010.654])),
        },
        coords={"space": space, "time": time},
    )


def test_export_to_csv_writes_expected_headers_and_los_values(tmp_path):
    stm = _make_minimal_stm_for_csv_export()
    output_path = tmp_path / "stm_los.csv"

    io.export_to_csv(
        stm=stm,
        save_path=str(output_path),
        model_parameter_layer_names=["linear"],
        ts_proj="los",
        point_annotation_label="TEST",
    )

    with output_path.open(newline="") as f:
        rows = list(csv.reader(f))

    header = rows[0]
    first_row = rows[1]

    expected_tail = ["20250101", "20250102", "a_20250101", "a_20250102"]
    assert header[-4:] == expected_tail
    assert "linear" in header

    idx_linear = header.index("linear")
    idx_20250101 = header.index("20250101")
    idx_a_20250102 = header.index("a_20250102")

    assert first_row[0] == "TEST_az00000010r00000100"
    assert float(first_row[idx_linear]) == pytest.approx(3.14159, abs=1e-5)
    assert float(first_row[idx_20250101]) == pytest.approx(0.12346, abs=1e-5)
    assert float(first_row[idx_a_20250102]) == pytest.approx(11.988, abs=1e-3)


def test_export_to_csv_vertical_uses_pov_layers(tmp_path):
    stm = _make_minimal_stm_for_csv_export()
    output_path = tmp_path / "stm_vertical.csv"

    io.export_to_csv(
        stm=stm,
        save_path=str(output_path),
        model_parameter_layer_names=["linear"],
        ts_proj="vertical",
        point_annotation_label="TEST",
    )

    with output_path.open(newline="") as f:
        rows = list(csv.reader(f))

    header = rows[0]
    first_row = rows[1]

    idx_linear = header.index("linear")
    idx_20250101 = header.index("20250101")

    assert float(first_row[idx_linear]) == pytest.approx(7.55555, abs=1e-5)
    assert float(first_row[idx_20250101]) == pytest.approx(9.1, abs=1e-5)


def test_export_to_csv_rejects_non_csv_path(tmp_path):
    stm = _make_minimal_stm_for_csv_export()

    with pytest.raises(AssertionError, match="is not a csv"):
        io.export_to_csv(
            stm=stm,
            save_path=str(tmp_path / "not_csv.txt"),
            model_parameter_layer_names=["linear"],
            ts_proj="los",
            point_annotation_label="TEST",
        )


def test_export_to_csv_raises_for_unknown_insert_marker(tmp_path, monkeypatch):
    stm = _make_minimal_stm_for_csv_export()
    output_path = tmp_path / "stm_invalid.csv"

    monkeypatch.setattr(io, "CSV_FIELD_NAMES", ["ID", "FUNC_INSERTS_UNKNOWN_HERE"])

    with pytest.raises(ValueError, match="Function insert FUNC_INSERTS_UNKNOWN_HERE requested but not defined"):
        io.export_to_csv(
            stm=stm,
            save_path=str(output_path),
            model_parameter_layer_names=["linear"],
            ts_proj="los",
            point_annotation_label="TEST",
        )


def _make_minimal_stm_for_shapefile_export():
    space = np.arange(2)
    return xr.Dataset(
        data_vars={
            "azimuth": ("space", np.array([10, 20], dtype=np.int32)),
            "range": ("space", np.array([100, 200], dtype=np.int32)),
            "lat": ("space", np.array([53.123456789, 53.223456789])),
            "lon": ("space", np.array([6.123456789, 6.223456789])),
            "h": ("space", np.array([40.1234, 41.5678])),
            "stc": ("space", np.array([1.23456, 2.34567])),
            "x_euclidean_proj_epsg28992": ("space", np.array([120000.123, 120010.987])),
            "y_euclidean_proj_epsg28992": ("space", np.array([480000.321, 480010.654])),
            "long_parameter_name": ("space", np.array([3.141592, 2.718281])),
        },
        coords={"space": space},
    )


def test_export_to_shapefile_builds_expected_properties_and_calls_to_file(monkeypatch, tmp_path):
    stm = _make_minimal_stm_for_shapefile_export()

    captured = {}

    class FakeGeoDataFrame:
        def __init__(self, properties, crs):
            captured["properties"] = properties
            captured["crs"] = crs

        def to_file(self, save_path):
            captured["save_path"] = save_path

    monkeypatch.setattr(io.gpd, "GeoDataFrame", FakeGeoDataFrame)

    save_path = str(tmp_path / "stm_points.shp")
    io.export_to_shapefile(
        stm=stm,
        save_path=save_path,
        projection="RD",
        model_parameter_layer_names=["long_parameter_name"],
        point_annotation_label="TEST",
    )

    props = captured["properties"]
    assert captured["crs"] == "EPSG:28992"
    assert captured["save_path"] == save_path

    assert props["ID"] == ["TEST_az00000010r00000100", "TEST_az00000020r00000200"]
    assert props["Azimuth"] == [10, 20]
    assert props["Range"] == [100, 200]
    assert props["STC [mm]"] == [1.235, 2.346]
    assert props["long_param"] == [3.14159, 2.71828]
    assert len(props["geometry"]) == 2


def test_export_to_shapefile_rejects_unknown_projection(tmp_path):
    stm = _make_minimal_stm_for_shapefile_export()

    with pytest.raises(AssertionError, match="Unknown requested projection"):
        io.export_to_shapefile(
            stm=stm,
            save_path=str(tmp_path / "stm_points.shp"),
            projection="UTM",
            model_parameter_layer_names=["long_parameter_name"],
            point_annotation_label="TEST",
        )


def test_export_to_shapefile_rejects_non_shp_path(tmp_path):
    stm = _make_minimal_stm_for_shapefile_export()

    with pytest.raises(AssertionError, match="is not a shapefile"):
        io.export_to_shapefile(
            stm=stm,
            save_path=str(tmp_path / "stm_points.geojson"),
            projection="RD",
            model_parameter_layer_names=["long_parameter_name"],
            point_annotation_label="TEST",
        )


def test_export_to_shapefile_raises_for_unknown_insert_marker(monkeypatch, tmp_path):
    stm = _make_minimal_stm_for_shapefile_export()

    monkeypatch.setattr(io, "SHAPEFILE_FIELD_NAMES", {"geometry": "Point", "properties": ["FUNC_INSERTS_UNKNOWN"]})

    with pytest.raises(ValueError, match="Function insert FUNC_INSERTS_UNKNOWN requested but not defined"):
        io.export_to_shapefile(
            stm=stm,
            save_path=str(tmp_path / "stm_points.shp"),
            projection="RD",
            model_parameter_layer_names=["long_parameter_name"],
            point_annotation_label="TEST",
        )


def _make_minimal_stm_for_convex_hull_export():
    space = np.arange(5)
    return xr.Dataset(
        data_vars={
            # Rectangle with one interior point.
            "x_euclidean_proj_epsg28992": ("space", np.array([0.0, 2.0, 2.0, 0.0, 1.0])),
            "y_euclidean_proj_epsg28992": ("space", np.array([0.0, 0.0, 1.0, 1.0, 0.5])),
            "lon": ("space", np.array([6.0, 6.2, 6.2, 6.0, 6.1])),
            "lat": ("space", np.array([53.0, 53.0, 53.1, 53.1, 53.05])),
        },
        coords={"space": space},
    )


def test_export_convex_hull_to_shapefile_builds_polygon_and_calls_to_file(monkeypatch, tmp_path):
    stm = _make_minimal_stm_for_convex_hull_export()

    captured = {}

    class FakeGeoDataFrame:
        def __init__(self, properties, crs):
            captured["properties"] = properties
            captured["crs"] = crs

        def to_file(self, save_path):
            captured["save_path"] = save_path

    monkeypatch.setattr(io.gpd, "GeoDataFrame", FakeGeoDataFrame)

    save_path = str(tmp_path / "hull.shp")
    io.export_convex_hull_to_shapefile(stm=stm, save_path=save_path, projection="RD")

    assert captured["crs"] == "EPSG:28992"
    assert captured["save_path"] == save_path
    assert "geometry" in captured["properties"]
    assert len(captured["properties"]["geometry"]) == 1

    polygon = captured["properties"]["geometry"][0]
    coords = list(polygon.exterior.coords)
    assert len(coords) >= 5
    assert coords[0] == coords[-1]


def test_export_convex_hull_to_shapefile_rejects_unknown_projection(tmp_path):
    stm = _make_minimal_stm_for_convex_hull_export()

    with pytest.raises(AssertionError, match="Unknown requested projection"):
        io.export_convex_hull_to_shapefile(stm=stm, save_path=str(tmp_path / "hull.shp"), projection="UTM")


def test_export_convex_hull_to_shapefile_rejects_non_shp_path(tmp_path):
    stm = _make_minimal_stm_for_convex_hull_export()

    with pytest.raises(AssertionError, match="is not a shapefile"):
        io.export_convex_hull_to_shapefile(stm=stm, save_path=str(tmp_path / "hull.geojson"), projection="RD")
