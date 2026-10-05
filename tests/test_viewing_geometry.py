import numpy as np
import pytest
import xarray as xr

from depsi.viewing_geometry import (
    add_cross_range,
    add_local_viewing_geometry,
    estimate_plane_viewing_geometry,
    fit_plane_viewing_geometry,
)

rng = np.random.default_rng(seed=42)


@pytest.fixture
def sample_stm():
    return xr.Dataset(
        data_vars={
            "lon": (["space"], np.array([0.0, 1.0, 2.0])),
            "lat": (["space"], np.array([10.0, 10.5, 11.0])),
        }
    )


@pytest.fixture
def mock_orbit_geometry(monkeypatch):
    def _mock(orbits, footprints=None, asc_inc=None, desc_inc=None, asc_alpha=None, desc_alpha=None):
        monkeypatch.setattr(
            "depsi.viewing_geometry.identify_s1_orbits_in_aoi", lambda lon, lat: (orbits, footprints or {})
        )
        monkeypatch.setattr("depsi.viewing_geometry.cfg.ConfigFile", lambda file_path: file_path)
        monkeypatch.setattr("depsi.viewing_geometry.SARModeFromCfg", lambda cfg_obj, orbit_mode: orbit_mode)
        monkeypatch.setattr(
            "depsi.viewing_geometry.nd.viewing_geometry",
            lambda *args, **kwargs: (None, None, asc_inc, desc_inc, asc_alpha, desc_alpha),
        )

    return _mock


class TestAddLocalViewingGeometry:
    def test_adds_local_incidence_and_azimuth_for_matching_orbit(self, sample_stm, mock_orbit_geometry):
        """The function should assign the orbit-specific incidence and azimuth arrays to the STM."""
        original_stm = sample_stm.copy(deep=True)
        desc_inc = np.array([[0.2, 0.4], [0.6, 0.7], [0.8, 0.9]], dtype=float)
        desc_alpha = np.array([[1.1, 1.4], [1.6, 1.7], [1.8, 1.9]], dtype=float)
        mock_orbit_geometry(("s1_dsc_t110",), asc_inc=None, desc_inc=desc_inc, asc_alpha=None, desc_alpha=desc_alpha)

        stm_out = add_local_viewing_geometry(sample_stm, "dummy.cfg", 0.01, "IWS", "s1_dsc_t110")

        assert set(stm_out.data_vars) == {"lon", "lat", "local_incidence_angle", "local_azimuth_angle"}
        np.testing.assert_allclose(stm_out["local_incidence_angle"].values, np.degrees(desc_inc[:, 0]))
        np.testing.assert_allclose(stm_out["local_azimuth_angle"].values, np.degrees(desc_alpha[:, 0]))
        assert sample_stm.identical(original_stm)

    @pytest.mark.parametrize(
        (
            "orbit",
            "orbits",
            "footprints",
            "asc_inc",
            "desc_inc",
            "asc_alpha",
            "desc_alpha",
            "expected_inc",
            "expected_alpha",
        ),
        [
            (
                "s1_asc_t110",
                ("s1_asc_t110", "s1_asc_t111"),
                {
                    "s1_asc_t110": [np.array([[1.0, 2.0], [2.0, 3.0], [3.0, 4.0]])],
                    "s1_asc_t111": [np.array([[2.5, 3.5], [3.5, 4.5], [4.5, 5.5]])],
                },
                np.array([[0.7, 0.3], [0.8, 0.5], [0.9, 0.6]], dtype=float),
                None,
                np.array([[1.1, 1.2], [1.3, 1.4], [1.5, 1.6]], dtype=float),
                None,
                np.array([[0.7, 0.3], [0.8, 0.5], [0.9, 0.6]], dtype=float)[:, 0],
                np.array([[1.1, 1.2], [1.3, 1.4], [1.5, 1.6]], dtype=float)[:, 0],
            ),
            (
                "s1_asc_t110",
                ("s1_asc_t110", "s1_asc_t111"),
                {
                    "s1_asc_t110": [np.array([[4.0, 9.0], [5.0, 10.0], [6.0, 11.0]])],
                    "s1_asc_t111": [np.array([[3.0, 8.0], [4.0, 9.0], [5.0, 10.0]])],
                },
                np.array([[0.2, 0.3], [0.4, 0.5], [0.6, 0.7]], dtype=float),
                None,
                np.array([[1.1, 1.2], [1.3, 1.4], [1.5, 1.6]], dtype=float),
                None,
                np.array([[0.2, 0.3], [0.4, 0.5], [0.6, 0.7]], dtype=float)[:, 0],
                np.array([[1.1, 1.2], [1.3, 1.4], [1.5, 1.6]], dtype=float)[:, 0],
            ),
            (
                "s1_dsc_t110",
                ("s1_dsc_t110", "s1_dsc_t111"),
                {
                    "s1_dsc_t110": [np.array([[5.0, 10.0], [6.0, 11.0], [7.0, 12.0]])],
                    "s1_dsc_t111": [np.array([[2.0, 7.0], [3.0, 8.0], [4.0, 9.0]])],
                },
                None,
                np.array([[0.9, 0.1], [0.8, 0.2], [0.7, 0.3]], dtype=float),
                None,
                np.array([[1.9, 1.0], [1.8, 1.1], [1.7, 1.2]], dtype=float),
                np.array([[0.9, 0.1], [0.8, 0.2], [0.7, 0.3]], dtype=float)[:, 0],
                np.array([[1.9, 1.0], [1.8, 1.1], [1.7, 1.2]], dtype=float)[:, 0],
            ),
        ],
    )
    def test_selects_correct_overlap_branch(
        self,
        sample_stm,
        mock_orbit_geometry,
        orbit,
        orbits,
        footprints,
        asc_inc,
        desc_inc,
        asc_alpha,
        desc_alpha,
        expected_inc,
        expected_alpha,
    ):
        """The overlap branch should pick the correct same-direction orbit column for each geometry case."""
        mock_orbit_geometry(orbits, footprints, asc_inc, desc_inc, asc_alpha, desc_alpha)

        stm_out = add_local_viewing_geometry(sample_stm, "dummy.cfg", 0.01, "IWS", orbit)

        np.testing.assert_allclose(stm_out["local_incidence_angle"].values, np.degrees(expected_inc))
        np.testing.assert_allclose(stm_out["local_azimuth_angle"].values, np.degrees(expected_alpha))

    def test_rejects_orbit_not_found_in_aoi(self, sample_stm, monkeypatch):
        """The selected orbit must be among the orbit footprints identified in the AOI."""
        monkeypatch.setattr("depsi.viewing_geometry.identify_s1_orbits_in_aoi", lambda lon, lat: (("s1_dsc_t111",), {}))

        with pytest.raises(AssertionError, match="Orbit provided is"):
            add_local_viewing_geometry(sample_stm, "dummy.cfg", 0.01, "IWS", "s1_dsc_t110")


class TestPlaneViewingGeometry:
    def test_fit_plane_viewing_geometry_recovers_known_plane(self):
        """The least-squares fit should recover the exact plane coefficients used to generate the angle field."""
        x = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
        y = np.array([0.0, 2.0, 1.0, 3.0, 2.0])
        expected_coeffs = np.array([2.5, -1.25, 10.0], dtype=float)
        angle = expected_coeffs[0] * x + expected_coeffs[1] * y + expected_coeffs[2]

        coeffs = fit_plane_viewing_geometry(x, y, angle)

        np.testing.assert_allclose(coeffs, expected_coeffs, rtol=0.0, atol=1e-12)
        np.testing.assert_allclose(estimate_plane_viewing_geometry(x, y, coeffs), angle, rtol=0.0, atol=1e-12)

    def test_fit_plane_viewing_geometry_uses_least_squares_for_noisy_data(self):
        """With noise, the fit should stay close to the underlying plane instead of returning arbitrary values."""
        x = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
        y = np.array([0.0, 2.0, 1.0, 3.0, 5.0, 4.0])
        true_coeffs = np.array([1.5, -0.75, 2.5], dtype=float)
        noise = np.array([0.15, -0.10, 0.05, -0.20, 0.10, -0.05], dtype=float)
        angle = true_coeffs[0] * x + true_coeffs[1] * y + true_coeffs[2] + noise

        coeffs = fit_plane_viewing_geometry(x, y, angle)

        np.testing.assert_allclose(coeffs, true_coeffs, rtol=2e-2, atol=5e-2)
        np.testing.assert_allclose(estimate_plane_viewing_geometry(x, y, coeffs), angle, rtol=5e-2, atol=2.5e-1)


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
