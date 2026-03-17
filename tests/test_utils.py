import numpy as np
import pytest
import xarray as xr

from depsi.utils import (
    compute_phase_difference,
    concatenate_stms,
    convert_geographic_coords_to_euclidean,
    generate_pnt_uids,
)

rng = np.random.default_rng(seed=42)


class TestGeneratePntUids:
    @pytest.mark.parametrize(
        "num_points, max_coordinates",
        [
            (411, 1e3),
            (1213, 1e4),
            (19, 1e3),
        ],
    )
    def test_generate_pnt_uids_robust(self, num_points, max_coordinates):
        """same radar coordinates should give same pnt_uid"""
        rng = np.random.default_rng(seed=12)
        # Create a sample STM dataset with azimuth and range coordinates
        azimuth_coords = rng.integers(0, max_coordinates, size=num_points)
        range_coords = rng.integers(0, max_coordinates, size=num_points)
        stm = xr.Dataset(coords={"azimuth": ("space", azimuth_coords), "range": ("space", range_coords)})

        # Generate unique point identifiers
        stm_with_uids_1 = generate_pnt_uids(stm)
        stm_with_uids_2 = generate_pnt_uids(stm)

        # Assert that the unique identifiers are the same for both generations
        np.testing.assert_array_equal(stm_with_uids_1["pnt_uid"].values, stm_with_uids_2["pnt_uid"].values)

    @pytest.mark.parametrize(
        "num_points, max_coordinates",
        [
            (97, 1e3 + 369),
            (1021, 1e4 - 1341),
            (5129, 1e5 + 35231),
        ],
    )
    def test_generate_pnt_uids_random_input(self, num_points, max_coordinates):
        """Test the generation of unique point identifiers."""
        rng = np.random.default_rng(seed=42)
        # Create a sample STM dataset with azimuth and range coordinates
        azimuth_coords = rng.integers(0, max_coordinates, size=num_points)
        range_coords = rng.integers(0, max_coordinates, size=num_points)
        stm = xr.Dataset(coords={"azimuth": ("space", azimuth_coords), "range": ("space", range_coords)})

        # Generate unique point identifiers
        stm_with_uids = generate_pnt_uids(stm)

        # results should be 1) unique 2) same shape as input and 3) integer
        pnt_uids = stm_with_uids["pnt_uid"].values
        assert len(np.unique(pnt_uids)) == num_points
        assert pnt_uids.shape == (num_points,)
        assert np.issubdtype(pnt_uids.dtype, np.integer)

    @pytest.mark.parametrize("overwrite_flag", [True, False])
    def test_generate_pnt_uids_overwrite(self, overwrite_flag, caplog):
        """Test the overwrite functionality of unique point identifier generation."""
        # Create a sample STM dataset with azimuth and range coordinates
        azimuth_coords = np.array([0, 1, 2, 3, 4])
        range_coords = np.array([10, 11, 12, 13, 14])
        stm = xr.Dataset(coords={"azimuth": ("space", azimuth_coords), "range": ("space", range_coords)})

        # First generation of unique point identifiers
        stm_with_uids = generate_pnt_uids(stm)

        # Store the first set of unique identifiers
        first_uids = stm_with_uids["pnt_uid"].values.copy()

        if overwrite_flag:
            # Generate unique point identifiers again
            stm_with_uids_overwritten = generate_pnt_uids(stm_with_uids, overwrite=overwrite_flag)
            # Assert that the unique identifiers remain the same after overwriting
            np.testing.assert_array_equal(stm_with_uids_overwritten["pnt_uid"].values, first_uids)
        else:
            # Attempt to generate unique point identifiers again without overwrite
            with caplog.at_level("WARNING"):
                _ = generate_pnt_uids(stm_with_uids, overwrite=overwrite_flag)

    def test_generate_pnt_uids_duplicate_detection(self):
        """Test that duplicate unique point identifiers are detected."""
        # Create a sample STM dataset with azimuth and range coordinates that will produce duplicates
        azimuth_coords = np.array([0, 0, 1, 1, 2])
        range_coords = np.array([10, 10, 11, 11, 12])
        stm = xr.Dataset(coords={"azimuth": ("space", azimuth_coords), "range": ("space", range_coords)})

        # Attempt to generate unique point identifiers and expect a ValueError
        with pytest.raises(ValueError):
            _ = generate_pnt_uids(stm)


class TestComputePhaseDifference:
    def test_compute_phase_difference(self):
        rng = np.random.default_rng(seed=42)
        complex = rng.uniform(-1, 1, (17, 3)) + 1j * rng.uniform(-1, 1, (17, 3))
        phase = np.angle(complex)

        # Phase difference between same phases should be zero
        assert np.allclose(compute_phase_difference(phase, phase, method="subtract"), 0.0)

        # Phase difference between complex and itself should be zero
        assert np.allclose(compute_phase_difference(complex, complex, method="conjmult"), 0.0)

        # Phase difference between complex and its negative should be pi or -pi
        assert np.allclose(np.abs(compute_phase_difference(complex, -complex, method="conjmult")), np.pi)

    def test_compute_phase_difference_errors(self):
        rng = np.random.default_rng(seed=42)
        complex = rng.uniform(-1, 1, (17, 3)) + 1j * rng.uniform(-1, 1, (17, 3))
        phase = np.angle(complex)

        with pytest.raises(NotImplementedError):
            _ = compute_phase_difference(phase, phase, method="invalid_method")

        with pytest.raises(ValueError):
            _ = compute_phase_difference(phase, phase, method="conjmult")  # conjmult requires complex inputs


def test_convert_geographic_coords_to_euclidean():
    """Test the conversion of geographic coordinates to Euclidean coordinates."""
    # Create a sample STM dataset with latitude and longitude coordinates
    latitudes = np.array([34.0, 34.1, 34.2])
    longitudes = np.array([-118.0, -118.1, -118.2])
    stm = xr.Dataset(coords={"latitude": ("space", latitudes), "longitude": ("space", longitudes)})

    # Convert geographic coordinates to Euclidean coordinates
    x, y = convert_geographic_coords_to_euclidean(stm["longitude"], stm["latitude"])

    # Assert that the output coordinates have the correct shape
    assert x.shape == (3,)
    assert y.shape == (3,)


def test_concatenate_multiple_stms():
    """Test concatenation of multiple STM datasets."""
    # Make three sample STM datasets
    space_dim_size = [3, 7, 4]
    time_var_common = rng.uniform(10, 20, 5)
    list_stms = []
    sp_coords_start = 0
    for idx, sp_size in enumerate(space_dim_size):
        stm_i = xr.Dataset(
            data_vars={
                "phase": (("space", "time"), rng.uniform(-np.pi, np.pi, (sp_size, 5))),
                "h2ph": (("space", "time"), rng.uniform(0, 1, (sp_size, 5))),
                f"time_var_{idx}": (("time"), rng.uniform(0, 10, 5)),  # time-only variable stm specific
                "time_var_common": (("time"), time_var_common),  # time-only variable common
            },
            coords={
                "space": ("space", rng.integers(sp_coords_start, sp_coords_start + 100, sp_size)),
                "time": ("time", np.arange(5)),
            },
        )
        sp_coords_start += 100
        list_stms.append(stm_i)

    # Concatenate the STM datasets
    concatenated_stm = concatenate_stms(list_stms)

    # Assert that the concatenated dataset has the correct shape
    assert concatenated_stm["phase"].shape == (14, 5)
    assert concatenated_stm["h2ph"].shape == (14, 5)
    # Assert that time-only variables are preserved correctly
    for idx in range(3):
        assert f"time_var_{idx}" in concatenated_stm
        np.testing.assert_array_equal(
            concatenated_stm[f"time_var_{idx}"].values,
            list_stms[idx][f"time_var_{idx}"].values,
        )
    np.testing.assert_array_equal(
        concatenated_stm["time_var_common"].values,
        list_stms[0]["time_var_common"].values,
    )
