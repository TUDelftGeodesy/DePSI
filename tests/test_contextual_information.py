import geopandas as gpd
import numpy as np
import pytest
import xarray as xr
from shapely.geometry import box

from depsi.contextual_information import (
    add_contextual_information_from_grid,
    add_contextual_information_from_point,
    add_contextual_information_from_polygon,
)

rng = np.random.default_rng(seed=42)


class TestAddContextualInformationFromPoint:
    def test_broadcasts_one_station_to_every_point(self):
        """A single station's data should be added to every STM point, and the STM should not be mutated in place."""
        n_space = 5
        lat = rng.uniform(52.35, 52.38, n_space)
        lon = rng.uniform(4.89, 4.92, n_space)
        time = np.array(["2026-07-01", "2026-07-13", "2026-07-25"], dtype="datetime64[ns]")

        stm = xr.Dataset(
            data_vars={"sd_h2ph": (["space", "time"], rng.uniform(0.5, 1.5, (n_space, len(time))))},
            coords={"lat": ("space", lat), "lon": ("space", lon), "time": time},
        )
        original_stm = stm.copy(deep=True)

        stm_out = add_contextual_information_from_point(
            stm,
            selected_attributes=["temperature"],
            attribute_sources={
                "temperature": (
                    52.3676,
                    4.9041,
                    {
                        "2026-07-01": {"temperature": 18.0},
                        "2026-07-13": {"temperature": 20.0},
                        "2026-07-25": {"temperature": 19.0},
                    },
                ),
            },
        )

        assert "temperature" in stm_out.data_vars
        expected = np.broadcast_to([18.0, 20.0, 19.0], (n_space, 3))
        np.testing.assert_array_equal(stm_out["temperature"].values, expected)

        # the input STM must not be mutated, even though this attribute was added via stmtools
        assert stm.identical(original_stm)

    def test_raises_value_error_when_selected_attribute_has_no_source(self):
        """A selected attribute without an entry in attribute_sources should be rejected."""
        time = np.array(["2026-07-01"], dtype="datetime64[ns]")
        stm = xr.Dataset(
            data_vars={"sd_h2ph": (["space", "time"], np.ones((1, 1)))},
            coords={"lat": ("space", [52.365]), "lon": ("space", [4.895]), "time": time},
        )

        with pytest.raises(ValueError):
            add_contextual_information_from_point(
                stm,
                selected_attributes=["temperature"],
                attribute_sources={
                    "temperature": (
                        52.3676,
                        4.9041,
                        {
                            "2026-07-01": {"wrong_name": 18.0},
                            "2026-07-13": {"wrong_name": 20.0},
                            "2026-07-25": {"temperature": 19.0},
                        },
                    ),
                },
            )

    def test_raises_value_error_when_selected_attribute_has_wrong_attribute_name(self):
        """A selected attribute with a wrong attribute name in attribute_sources should be rejected."""
        time = np.array(["2026-07-01"], dtype="datetime64[ns]")
        stm = xr.Dataset(
            data_vars={"sd_h2ph": (["space", "time"], np.ones((1, 1)))},
            coords={"lat": ("space", [52.365]), "lon": ("space", [4.895]), "time": time},
        )

        with pytest.raises(ValueError, match="No data source given for 'temperature'"):
            add_contextual_information_from_point(
                stm,
                selected_attributes=["temperature"],
                attribute_sources={},
            )


class TestAddContextualInformationFromGrid:
    def test_matches_each_point_to_nearest_station(self):
        """Points should be matched to their nearest station, and missing dates should be filled with NaN."""
        time = np.array(["2026-07-01", "2026-07-13", "2026-07-25"], dtype="datetime64[ns]")

        # two STM points: one right on top of the north station, one right on top of the south station
        stm = xr.Dataset(
            data_vars={"sd_h2ph": (["space", "time"], rng.uniform(0.5, 1.5, (2, len(time))))},
            coords={"lat": ("space", [52.40, 52.30]), "lon": ("space", [4.90, 4.90]), "time": time},
        )
        original_stm = stm.copy(deep=True)

        stm_out = add_contextual_information_from_grid(
            stm,
            selected_attributes=["temperature"],
            attribute_sources={
                "temperature": [
                    (
                        52.40,
                        4.90,
                        {"2026-07-01": {"temperature": 14.0}, "2026-07-13": {"temperature": 15.0}},
                        # 2026-07-25 deliberately missing for the north station
                    ),
                    (
                        52.30,
                        4.90,
                        {
                            "2026-07-01": {"temperature": 22.0},
                            "2026-07-13": {"temperature": 23.0},
                            "2026-07-25": {"temperature": 24.0},
                        },
                    ),
                ],
            },
        )

        assert "temperature" in stm_out.data_vars
        north_point, south_point = stm_out["temperature"].values
        np.testing.assert_array_equal(north_point, [14.0, 15.0, np.nan])
        np.testing.assert_array_equal(south_point, [22.0, 23.0, 24.0])

        assert stm.identical(original_stm)

    def test_raises_value_error_when_selected_attribute_has_no_source(self):
        """A selected attribute with a wrong attribute name."""
        time = np.array(["2026-07-01"], dtype="datetime64[ns]")
        stm = xr.Dataset(
            data_vars={"sd_h2ph": (["space", "time"], np.ones((1, 1)))},
            coords={"lat": ("space", [52.365]), "lon": ("space", [4.895]), "time": time},
        )

        with pytest.raises(ValueError):
            add_contextual_information_from_grid(
                stm,
                selected_attributes=["temperature"],
                attribute_sources={
                    "temperature": [],
                },
            )

    def test_raises_value_error_corrupted_location_coords(self):
        """A selected attribute with corrupted location coordinates."""
        time = np.array(["2026-07-01"], dtype="datetime64[ns]")
        stm = xr.Dataset(
            data_vars={"sd_h2ph": (["space", "time"], np.ones((1, 1)))},
            coords={"lat": ("space", [52.365]), "lon": ("space", [4.895]), "time": time},
        )

        with pytest.raises(ValueError):
            add_contextual_information_from_grid(
                stm,
                selected_attributes=["temperature"],
                attribute_sources={
                    "temperature": [
                        (
                            52.40,  # 1D coords which is corrupted
                            {"2026-07-01": {"temperature": 14.0}, "2026-07-13": {"temperature": 15.0}},
                            # 2026-07-25 deliberately missing for the north station
                        ),
                    ],
                },
            )


class TestAddContextualInformationFromPolygon:
    def test_matches_points_by_containment_and_leaves_outside_points_none(self):
        """Points inside a polygon should get its field value; points outside every polygon should get None."""
        stm = xr.Dataset(
            coords={
                "lat": ("space", [52.365, 52.365, 52.500]),
                "lon": ("space", [4.895, 4.905, 4.500]),
                "time": np.array(["2026-07-01"], dtype="datetime64[ns]"),
            }
        ).chunk({"space": -1, "time": -1})
        original_stm = stm.copy(deep=True)

        polygons = gpd.GeoDataFrame(
            {"land_use": ["urban"]},
            geometry=[box(4.890, 52.360, 4.900, 52.370)],
            crs="EPSG:4326",
        )

        stm_out = add_contextual_information_from_polygon(
            stm,
            selected_attributes=["land_use"],
            attribute_sources={"land_use": (polygons, "land_use")},
        )

        assert "land_use" in stm_out.data_vars
        assert "time" not in stm_out["land_use"].dims
        assert list(stm_out["land_use"].values) == ["urban", None, None]

        assert stm.identical(original_stm)

    def test_renames_field_to_attribute_name_when_different(self):
        """When the polygon's field name differs from the attribute name, the result should be renamed."""
        stm = xr.Dataset(
            coords={
                "lat": ("space", [52.365]),
                "lon": ("space", [4.895]),
                "time": np.array(["2026-07-01"], dtype="datetime64[ns]"),
            }
        ).chunk({"space": -1, "time": -1})

        polygons = gpd.GeoDataFrame(
            {"LU_CODE": ["urban"]},
            geometry=[box(4.890, 52.360, 4.900, 52.370)],
            crs="EPSG:4326",
        )

        stm_out = add_contextual_information_from_polygon(
            stm,
            selected_attributes=["land_use"],
            attribute_sources={"land_use": (polygons, "LU_CODE")},
        )

        assert "land_use" in stm_out.data_vars
        assert "LU_CODE" not in stm_out.data_vars

    def test_raises_value_error_when_selected_attribute_has_no_source(self):
        """A selected attribute without an entry in attribute_sources should be rejected."""
        stm = xr.Dataset(
            coords={
                "lat": ("space", [52.365]),
                "lon": ("space", [4.895]),
                "time": np.array(["2026-07-01"], dtype="datetime64[ns]"),
            }
        ).chunk({"space": -1, "time": -1})

        with pytest.raises(ValueError, match="No data source given for 'land_use'"):
            add_contextual_information_from_polygon(
                stm,
                selected_attributes=["land_use"],
                attribute_sources={},
            )


class TestLoadStationSeriesFromFile:
    def test_reads_temperature_from_a_csv_file(self, tmp_path):
        """add_contextual_information_from_point should read a station's data from a CSV path, not just a dict."""
        temperature_csv = tmp_path / "temperature.csv"
        temperature_csv.write_text("Date,temperature\n2026-07-01,18.0\n2026-07-13,20.0\n2026-07-25,19.0\n")

        time = np.array(["2026-07-01", "2026-07-13", "2026-07-25"], dtype="datetime64[ns]")
        stm = xr.Dataset(
            data_vars={"sd_h2ph": (["space", "time"], rng.uniform(0.5, 1.5, (3, len(time))))},
            coords={"lat": ("space", [52.365, 52.366, 52.367]), "lon": ("space", [4.895, 4.896, 4.897]), "time": time},
        )

        stm_out = add_contextual_information_from_point(
            stm,
            selected_attributes=["temperature"],
            attribute_sources={"temperature": (52.3676, 4.9041, str(temperature_csv))},
        )

        assert "temperature" in stm_out.data_vars
        expected = np.broadcast_to([18.0, 20.0, 19.0], (3, 3))
        np.testing.assert_array_equal(stm_out["temperature"].values, expected)
