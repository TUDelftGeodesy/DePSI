"""test_network.py"""

import numpy as np
import pytest
import xarray as xr

from depsi.network import (
    _compute_phase_difference,
    _network_relation_matrix,
    arc_selection,
    form_network,
    remove_network_points_min_connections,
)


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
def stm_random_grid():
    """STM points forming a 10x10 grid"""
    N_time = 50
    grid_shape = 10
    x_grid, y_grid = np.meshgrid(np.arange(0, 100, grid_shape), np.arange(0, 100, grid_shape))
    N_points = x_grid.flatten().shape[0]

    stm = xr.Dataset(
        coords={
            "space": (["space"], np.arange(N_points)),
            "time": (["time"], np.arange(N_time)),
            "x": (["space"], x_grid.flatten()),
            "y": (["space"], y_grid.flatten()),
        },
        data_vars={
            "phase": (["space", "time"], np.random.uniform(0, 1, (N_points, N_time))),
            "h2ph": (["space", "time"], np.random.uniform(0, 1, (N_points, N_time))),
            "ambiguity": (["space", "time"], np.random.choice([-1, 0, 1], (N_points, N_time), p=[0.02, 0.96, 0.02])),
        },
    )

    return stm


@pytest.fixture
def arcs_random(stm_random):
    """Fixture of fully connected arcs from stm_random."""
    # Fully connected arcs
    # Defaul method is redundant
    # No max_length, so all points are connected
    arcs = form_network(stm_random, key_phase="phase", key_h2ph="h2ph", key_Btemp="time")

    # Most arcs has quality 0.9
    # Except the last two are 0.0
    # The first five are 0.99
    real_ens_coh = np.zeros((arcs.sizes["space"],))  # Put all values in real, all imaginary are 0
    real_ens_coh[:-2] = 0.9
    real_ens_coh[:5] = 0.99
    arcs["ens_coh"] = (("space"), real_ens_coh + 1j * np.zeros((arcs.sizes["space"],)))

    return arcs


class TestNetworkFormation:
    def test_form_network_simulated_grid(self, stm_random_grid):
        arcs = form_network(
            stm_random_grid,
            key_phase="phase",
            key_h2ph="h2ph",
            key_Btemp="time",
            key_xlabel="x",
            key_ylabel="y",
            max_length=25,
            n_links=8,
            num_partitions=8,
        )

        source = arcs["source"].values
        target = arcs["target"].values

        assert arcs.sizes["space"] == 428  # nr arcs should be 428 with a 10x10 grid setting
        assert np.all(np.diff(source) >= 0)  # check if source is mono-increasing
        assert np.all(source < target)  # check if all sources < targets
        assert (
            np.unique(np.column_stack((source, target)), axis=0).shape[0] == source.shape[0]
        )  # check if all (source, target) pairs are unique

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

        # Threshold is 0.5, so only the last two arcs are discarded
        # The min_n_connection should not affect the selection
        assert selected_arcs.sizes["space"] == arcs_random.sizes["space"] - 2

    def test_select_arcs_non_connected(self, arcs_random, caplog):
        """Should keep the first five arcs which are disconnected."""
        # this should raise a logger warning of disconnected arcs
        with caplog.at_level("WARNING"):
            _ = arc_selection(
                arcs_random,
                threshold=0.99,
                selection_method="ens_coh",
                min_n_connection=0,
            )

    def test_remove_network_points_min_connections_nconnection_zero(self, stm_random, arcs_random):
        """Raise error when min_connections <1."""
        with pytest.raises(ValueError):
            remove_network_points_min_connections(stm_random, arcs_random, min_connections=-1)
        with pytest.raises(ValueError):
            remove_network_points_min_connections(stm_random, arcs_random, min_connections=0)

    def test_remove_network_points_min_connections_keep_all_pnts(self, stm_random, arcs_random):
        """No STM points removed since no arc is discarded."""
        stm_updated, arcs_updated = remove_network_points_min_connections(stm_random, arcs_random, min_connections=1)

        assert stm_updated.sizes["space"] == stm_random.sizes["space"]
        assert arcs_updated.sizes["space"] == arcs_random.sizes["space"]

    def test_remove_network_points_min_connections_discard_one(self, stm_random, arcs_random):
        """Remove one STM point."""
        # remove arcs with source or target == 1
        arcs = arcs_random.copy(deep=True)
        arcs = arcs.where((arcs["source"] != 1) & (arcs["target"] != 1), drop=True)

        stm_updated, arcs_updated = remove_network_points_min_connections(stm_random, arcs, min_connections=1)

        # Should remove the point with index 1
        assert stm_updated.sizes["space"] == stm_random.sizes["space"] - 1


class TestNetworkUnwrap:
    @pytest.mark.parametrize(
        ["idx_source", "idx_target", "n_points"],
        [
            (np.array([0, 1, 2]), np.array([1, 2, 3]), 4),  # 4 points, 3 arcs
            (np.array([0, 1, 2]), np.array([1, 2, 3]), 7),  # 7 points, 3 arcs
            (np.array([1, 1, 2, 2]), np.array([0, 2, 1, 3]), 4),  # 4 points, 4 arcs, unsorted
            (np.array([0, 0, 0, 1, 1, 2, 2]), np.array([1, 2, 3, 3, 4, 3, 4]), 5),  # 5 points, 6 arcs
        ],
    )
    def test_init_network_relation_matrix(
        self,
        idx_source,
        idx_target,
        n_points,
    ):
        A = _network_relation_matrix(idx_source, idx_target, n_points)

        # Create expected matrix in a for loop
        A_exp = np.zeros((idx_source.shape[0], n_points), dtype=int)
        for i, (src, tgt) in enumerate(zip(idx_source, idx_target, strict=False)):
            A_exp[i, src] = -1
            A_exp[i, tgt] = 1

        assert A.shape == A_exp.shape
        assert np.all(A.todense() == A_exp)
