import numpy as np
import pytest
import xarray as xr

from depsi.viewing_geometry import add_cross_range

rng = np.random.default_rng(seed=42)


class TestAddCrossRange:
    def test_computes_crossrange_from_existing_incidence_angle_and_sd_h2ph(self):
        """add_cross_range should compute sd_cr2ph = sd_h2ph * sin(local_incidence_angle) from existing STM data."""
        n_space, n_time = 3, 2
        sd_h2ph = rng.uniform(0.5, 1.5, (n_space, n_time))
        local_incidence_angle = np.array([30.0, 35.0, 40.0])

        stm = xr.Dataset(
            data_vars={
                "sd_h2ph": (["space", "time"], sd_h2ph),
                "local_incidence_angle": (["space"], local_incidence_angle),
            },
            coords={"time": np.array(["2026-07-01", "2026-07-13"], dtype="datetime64[ns]")},
        )
        original_stm = stm.copy(deep=True)

        stm_out = add_cross_range(stm)

        assert "sd_cr2ph" in stm_out.data_vars
        expected = sd_h2ph * np.sin(np.radians(local_incidence_angle[:, np.newaxis]))
        np.testing.assert_allclose(stm_out["sd_cr2ph"].values, expected)

        # the input STM must not be mutated
        assert stm.identical(original_stm)

    def test_rejects_stm_missing_sd_h2ph(self):
        """Missing sd_h2ph should raise a clear ValueError naming it."""
        stm = xr.Dataset(
            data_vars={"local_incidence_angle": (["space"], [30.0, 35.0])},
            coords={"time": np.array(["2026-07-01"], dtype="datetime64[ns]")},
        )

        with pytest.raises(ValueError, match="sd_h2ph"):
            add_cross_range(stm)

    def test_rejects_stm_missing_local_incidence_angle(self):
        """Missing local_incidence_angle should raise a ValueError hinting to call add_local_viewing_geometry."""
        stm = xr.Dataset(
            data_vars={"sd_h2ph": (["space", "time"], rng.uniform(0.5, 1.5, (2, 1)))},
            coords={"time": np.array(["2026-07-01"], dtype="datetime64[ns]")},
        )

        with pytest.raises(ValueError, match="add_local_viewing_geometry"):
            add_cross_range(stm)
