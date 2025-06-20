"""test_network.py"""

import numpy as np
import pytest
import xarray as xr

from depsi.network import _compute_phase_difference, arc_selection, form_network, remove_isolated_stm


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


@pytest.fixture
def arcs_random(stm_random):
    """Fixture to create a random STM dataset."""
    # Fully connected arcs
    arcs = form_network(stm_random, key_phase="phase", key_h2ph="h2ph", key_Btemp="time")

    # All arcs has quality 0.9, exept the last two are 0.0
    real_ens_coh = np.zeros((arcs.sizes["space"],))  # Put all values in real, all imaginary are 0
    real_ens_coh[:-2] = 0.9
    real_ens_coh[:5] = 0.99
    arcs["ens_coh"] = (("space"), real_ens_coh + 1j * np.zeros((arcs.sizes["space"],)))

    return arcs


class TestNetworkFormation:
    @pytest.mark.parametrize("method", ["subtract", "conjmult"])
    def test_compute_phase_difference(self, stm_random, method):
        arcs = form_network(stm_random, key_phase="phase", key_h2ph="h2ph", key_Btemp="time")
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
        stm_arcs = form_network(
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
        stm_arcs = form_network(
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
            form_network(
                stm_random,
                key_phase="phase",
                key_h2ph="h2ph",
                key_Btemp="time",
                network_method="unknown",
                dphase_method="subtract",
            )
        with pytest.raises(NotImplementedError):
            form_network(
                stm_random,
                key_phase="phase",
                key_h2ph="h2ph",
                key_Btemp="time",
                network_method="delaunay",
                dphase_method="unknown",
            )


class TestArcSelection:
    @pytest.mark.parametrize("thres, min_n_connection", [(0.99, 0), (0.5, 999)])
    def test_select_arcs_return_zero(self, arcs_random, thres, min_n_connection):
        """Should return zero arcs, two high threshold or too high min_n_connection."""
        # Select arcs based on ens_coh threshold.
        selected_arcs = arc_selection(
            arcs_random,
            threshold=thres,
            selection_method="ens_coh",
            min_n_connection=min_n_connection,
        )

        assert selected_arcs.sizes["space"] == 0

    @pytest.mark.parametrize("thres, min_n_connection", [(0.5, 2), (0.5, 1)])
    def test_select_arcs_discard_two(self, arcs_random, thres, min_n_connection):
        """Should only discard two arcs, with ens_coh < 0.5."""
        # Select arcs based on ens_coh threshold.
        selected_arcs = arc_selection(
            arcs_random,
            threshold=thres,
            selection_method="ens_coh",
            min_n_connection=min_n_connection,
        )

        assert selected_arcs.sizes["space"] == arcs_random.sizes["space"] - 2

    def test_select_arcs_non_connected(self, arcs_random, caplog):
        """Should only discard two arcs, with ens_coh < 0.5."""
        # this should raise a logger warning
        with caplog.at_level("WARNING"):
            _ = arc_selection(
                arcs_random,
                threshold=0.99,
                selection_method="ens_coh",
                min_n_connection=0,
            )

    def test_remove_isolated_stm_keep_all_pnts(self, stm_random, arcs_random):
        """Should remove isolated STM points."""
        stm_updated, arcs_updated = remove_isolated_stm(stm_random, arcs_random)

        assert stm_updated.sizes["space"] == stm_random.sizes["space"]
        assert arcs_updated.sizes["space"] == arcs_random.sizes["space"]

    def test_remove_isolated_stm_discard_one(self, stm_random, arcs_random):
        """Should remove isolated STM points."""
        # remove arcs with source or target == 1
        arcs = arcs_random.copy(deep=True)
        arcs = arcs.where((arcs["source"] != 1) & (arcs["target"] != 1), drop=True)

        stm_updated, arcs_updated = remove_isolated_stm(stm_random, arcs)

        assert stm_updated.sizes["space"] == stm_random.sizes["space"] - 1
