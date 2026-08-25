import os

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
    DATA_PATH = os.path.join(os.path.dirname(__file__), "data", "s1_dsc037_RC.csv")

    # Read CSV using the function under test
    df, dates = io.read_rcs_csv(DATA_PATH)

    # 1️⃣ Check that the '*****' row existed in the original file
    with open(DATA_PATH) as f:
        content = f.read()
    assert "*****" in content, "Metadata markers '*****' not found in CSV file"

    # 2️⃣ Check required columns exist
    required_cols = ["Range", "Azimuth", "Lat", "Lon", "Height"]
    for col in required_cols:
        assert col in df.columns, f"Required column '{col}' missing from DataFrame"

    # 3️⃣ Check required columns are numeric (float)
    for col in required_cols:
        assert pd.api.types.is_numeric_dtype(df[col]), f"Column '{col}' is not numeric"

    # 4️⃣ Check date columns are present and numeric
    for date_col in dates:
        assert date_col in df.columns, f"Date column '{date_col}' missing from DataFrame"
        assert pd.api.types.is_numeric_dtype(df[date_col]), f"Date column '{date_col}' is not numeric"


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
