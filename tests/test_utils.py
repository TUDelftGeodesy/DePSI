import math

import numpy as np
import pytest
import xarray as xr

from depsi.utils import generate_pnt_uids


class TestGeneratePntUids:
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
        scale = 10 ** (math.floor(math.log10(range_coords.max())) + 1)
        stm = xr.Dataset(coords={"azimuth": ("space", azimuth_coords), "range": ("space", range_coords)})

        # Generate unique point identifiers
        stm_with_uids = generate_pnt_uids(stm)

        # Expected unique identifiers
        # 10**5 is determined based on the maximum range coordinate scale which is 1e4
        expected_uids = (azimuth_coords.astype(np.uint64) + 1) * scale + (
            range_coords.astype(np.uint64) + 1
        )

        # Assert that the generated unique identifiers match the expected values
        np.testing.assert_array_equal(stm_with_uids["pnt_uid"].values, expected_uids)

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
