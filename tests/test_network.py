"""test_network.py"""

import numpy as np
import pytest
import xarray as xr

from depsi.network import _compute_phase_difference, generate_arcs


@pytest.fixture
def stm_random():
    """Fixture to create a random STM dataset."""
    rng = np.random.default_rng(42)
    Npoints = 12  # Number of points
    Ntimes = 31  # Number of epochs
    # Coordinates and time
    lat = np.linspace(51.14, 51.15, Npoints)
    lon = rng.uniform(6.9, 7.0, Npoints)
    time = np.arange(Ntimes)
    # Data
    complex = rng.uniform(-1, 1, (Npoints, Ntimes)) + 1j * rng.uniform(-1, 1, (Npoints, Ntimes))
    phase = np.angle(complex)
    h2ph = rng.uniform(1e3, 1e4, (Npoints, Ntimes))
    # Create the xarray Dataset
    stm = xr.Dataset(
        data_vars={
            "phase": (("space", "time"), phase),
            "h2ph": (("space", "time"), h2ph),
            "complex": (("space", "time"), complex),
        },
        coords={
            "space": ("space", np.arange(Npoints)),
            "time": ("time", time),
            "lat": ("space", lat),
            "lon": ("space", lon),
        },
    )

    return stm


class TestNetwork:
    @pytest.mark.parametrize("method", ["subtract", "conjmult"])
    def test_compute_phase_difference(self, stm_random, method):
        arcs = generate_arcs(stm_random, key_phase="phase", key_h2ph="h2ph", key_Btemp="time")
        d_phase_subtract_0_0 = _compute_phase_difference(
            stm_random, arcs["source"], arcs["source"], "phase", "complex", method=method
        )
        d_phase_subtract_0_1 = _compute_phase_difference(
            stm_random, arcs["source"], arcs["target"], "phase", "complex", method=method
        )
        # Phase differences should be zero for the same source.
        assert d_phase_subtract_0_0 == pytest.approx(np.zeros(d_phase_subtract_0_0.shape), abs=1e-7)
        # Phase difference should be within the range of -2*pi to 2*pi for different sources.
        assert d_phase_subtract_0_1 == pytest.approx(np.zeros(d_phase_subtract_0_1.shape), abs=2 * np.pi + 1e-7)

    def test_stm_to_arcs_subtract(self, stm_random):
        # Generate arcs of a Delaunay network with subtracted phase differences.
        stm_arcs = generate_arcs(
            stm_random,
            key_phase="phase",
            key_h2ph="h2ph",
            key_Btemp="time",
            network_method="delaunay",
            max_length=0.05,
            dphase_method="subtract",
        )

        assert all(
            [all([-2 * np.pi <= phase <= 2 * np.pi for phase in phases]) for phases in stm_arcs["d_phase"].values]
        )

    def test_stm_to_arcs_conjmult(self, stm_random):
        # Generate arcs of a Delaunay network with conjugate multiplication phase differences.
        stm_arcs = generate_arcs(
            stm_random,
            key_phase="phase",
            key_h2ph="h2ph",
            key_Btemp="time",
            network_method="delaunay",
            max_length=0.05,
            dphase_method="conjmult",
        )

        assert all([all([-np.pi <= phase <= np.pi for phase in phases]) for phases in stm_arcs["d_phase"].values])

    def test_stm_to_arcs_fail(self, stm_random):
        # Test incorrect method fail.
        with pytest.raises(NotImplementedError):
            generate_arcs(
                stm_random,
                key_phase="phase",
                key_h2ph="h2ph",
                key_Btemp="time",
                network_method="unknown",
                dphase_method="subtract",
            )
        with pytest.raises(NotImplementedError):
            generate_arcs(
                stm_random,
                key_phase="phase",
                key_h2ph="h2ph",
                key_Btemp="time",
                network_method="delaunay",
                dphase_method="unknown",
            )
