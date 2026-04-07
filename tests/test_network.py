"""test_network.py"""

import numpy as np
import pytest
import xarray as xr

from depsi.network import (
    _ensure_network_min_connections,
    _ensure_single_network,
    _network_relation_matrix,
    _remove_network_points_min_connections,
    form_network,
    spatial_integration,
)


@pytest.fixture
def stm_random():
    """Fixture to create a random STM dataset."""
    rng = np.random.default_rng(42)
    Npoints = 12  # Number of points
    Ntimes = 31  # Number of epochs
    # Coordinates and time
    lat = rng.uniform(51.14, 51.15, Npoints)
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
    arcs = form_network(stm_random, key_phase="phase", key_h2ph="h2ph", key_Btemporal="time")

    # Most arcs have quality 0.9
    # Except the last two have quality 0.0
    # And the first five have quality 0.99
    temp_coh = np.zeros((arcs.sizes["space"],))
    temp_coh[:-2] = 0.9
    temp_coh[:5] = 0.99
    arcs["temp_coh"] = (("space"), temp_coh)

    return arcs


def _build_network_components(component_sizes: list[int]) -> tuple[xr.Dataset, xr.Dataset]:
    """Build coordinate only point/arcs STMs from connected-component sizes."""
    n_points = int(sum(component_sizes))
    stm_pnts = xr.Dataset(coords={"space": ("space", np.arange(n_points))})

    source = []
    target = []
    offset = 0
    for size in component_sizes:
        # Build each component as a simple chain graph.
        for idx in range(offset, offset + size - 1):
            source.append(idx)
            target.append(idx + 1)
        offset += size

    n_arcs = len(source)
    stm_arcs = xr.Dataset(
        coords={
            "space": ("space", np.arange(n_arcs)),
            "source": ("space", np.array(source, dtype=int)),
            "target": ("space", np.array(target, dtype=int)),
        }
    )

    return stm_arcs, stm_pnts


class TestNetworkFormation:
    def test_form_network_simulated_grid(self, stm_random_grid):
        arcs = form_network(
            stm_random_grid,
            key_phase="phase",
            key_h2ph="h2ph",
            key_Btemporal="time",
            key_xcrds="x",
            key_ycrds="y",
            max_length=25,
            min_links=8,
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
        assert np.unique(arcs["uid"].values).shape[0] == arcs.sizes["space"]

    def test_stm_to_arcs_subtract(self, stm_random):
        # Generate arcs of a Delaunay network with subtracted phase differences.
        stm_arcs = form_network(
            stm_random,
            key_phase="phase",
            key_h2ph="h2ph",
            key_Btemporal="time",
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
            key_Btemporal="time",
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
                key_Btemporal="time",
                network_method="unknown",
                dphase_method="subtract",
            )
        with pytest.raises(NotImplementedError):
            form_network(
                stm_random,
                key_phase="phase",
                key_h2ph="h2ph",
                key_Btemporal="time",
                network_method="delaunay",
                dphase_method="unknown",
            )


class TestNetworkEnsure:
    @pytest.mark.parametrize("thres, min_n_connections", [(0.5, 2), (0.5, 1)])
    def test_select_arcs_discard_two(self, arcs_random, stm_random, thres, min_n_connections):
        """Should only discard two arcs, with temp_coh < 0.5."""
        # Select arcs based on temp_coh threshold.
        mask = np.abs(arcs_random["temp_coh"]) > thres  # mask as DataArray
        arcs_selected = arcs_random.where(mask, drop=True)
        arcs_results, _ = _ensure_network_min_connections(arcs_selected, stm_random, min_connections=min_n_connections)

        # Threshold is 0.5, so only the last two arcs are discarded
        # The min_n_connections should not affect the selection
        assert arcs_results.sizes["space"] == arcs_random.sizes["space"] - 2

    def test__remove_network_points_min_connections_nconnection_zero(self, stm_random, arcs_random):
        """Raise error when min_connections <1."""
        with pytest.raises(ValueError):
            _remove_network_points_min_connections(stm_random, arcs_random, min_connections=-1)
        with pytest.raises(ValueError):
            _remove_network_points_min_connections(stm_random, arcs_random, min_connections=0)

    def test__remove_network_points_min_connections_keep_all_pnts(self, stm_random, arcs_random):
        """No STM points removed since no arc is discarded."""
        stm_updated, arcs_updated = _remove_network_points_min_connections(stm_random, arcs_random, min_connections=1)

        assert stm_updated.sizes["space"] == stm_random.sizes["space"]
        assert arcs_updated.sizes["space"] == arcs_random.sizes["space"]

    def test__remove_network_points_min_connections_discard_one(self, stm_random, arcs_random):
        """Remove one STM point."""
        # remove arcs with source or target == 1
        arcs = arcs_random.copy(deep=True)
        arcs = arcs.where((arcs["source"] != 1) & (arcs["target"] != 1), drop=True)

        stm_updated, arcs_updated = _remove_network_points_min_connections(stm_random, arcs, min_connections=1)

        # Should remove the point with index 1
        assert stm_updated.sizes["space"] == stm_random.sizes["space"] - 1

    @pytest.mark.parametrize("component_sizes", [[12]])
    def test_ensure_single_network_no_separated_part(self, component_sizes):
        """Keep the network unchanged when there is only one connected component."""
        stm_arcs, stm_pnts = _build_network_components(component_sizes)

        stm_arcs_out, stm_pnts_out = _ensure_single_network(stm_arcs, stm_pnts)

        assert stm_pnts_out.sizes["space"] == stm_pnts.sizes["space"]
        assert stm_arcs_out.sizes["space"] == stm_arcs.sizes["space"]
        assert np.array_equal(stm_arcs_out["source"].values, stm_arcs["source"].values)
        assert np.array_equal(stm_arcs_out["target"].values, stm_arcs["target"].values)

    @pytest.mark.parametrize(
        "component_sizes",
        [
            [9, 1],
            [17, 2, 1],
            [41, 3, 2, 2, 2],
        ],
    )
    def test_ensure_single_network_keep_largest_significant(self, component_sizes):
        """Keep only the largest component when it is significant enough."""
        stm_arcs, stm_pnts = _build_network_components(component_sizes)
        largest_size = max(component_sizes)

        stm_arcs_out, stm_pnts_out = _ensure_single_network(stm_arcs, stm_pnts)

        assert stm_pnts_out.sizes["space"] == largest_size
        assert stm_arcs_out.sizes["space"] == largest_size - 1
        assert np.all(stm_arcs_out["source"].values >= 0)
        assert np.all(stm_arcs_out["target"].values >= 0)
        assert np.all(stm_arcs_out["source"].values < largest_size)
        assert np.all(stm_arcs_out["target"].values < largest_size)

    @pytest.mark.parametrize(
        "component_sizes",
        [
            [6, 4],
            [8, 7, 1],
            [10, 9, 1, 1],
        ],
    )
    def test_ensure_single_network_raise_when_largest_not_significant(self, component_sizes):
        """Raise when the largest component is not clearly dominant."""
        stm_arcs, stm_pnts = _build_network_components(component_sizes)

        with pytest.raises(RuntimeError):
            _ensure_single_network(stm_arcs, stm_pnts)


class TestNetworkUnwrap:
    @pytest.mark.parametrize(
        ["id_ref", "idx_err_space", "idx_err_time", "error_values", "skip_network_adjustment"],
        [
            (3, [], [], [], False),  # No error
            (3, [2, 11], [7, 13], [-1, 1], False),  # Two errors in arc ambiguities
            (3, [2, 11], [7, 13], [-1, 1], True),  # Two errors, skip network adjustment, should still be corrected
            (9, [0, 4, 8], [5, 10, 15], [1, -100, 1], False),  # Three errors, one large, but should be corrected
        ],
    )
    def test_spatial_integration(self, id_ref, idx_err_space, idx_err_time, error_values, skip_network_adjustment):
        """Test spatial unwrapping based on arc ambiguities.

        Build points with true value of ambiguities.
        Construct arcs with arc ambiguities derived from true ambiguities.
        Add tiny errors to arc ambiguities at certain space/time indices.

        Then perform spatial unwrapping with a specified reference point.

        The spatial unwrapping should be able to solve the point ambiguities correctly.
        The solved ambiguities should w.r.t. the reference point.
        """
        # Set up test parameters
        rng = np.random.default_rng(42)
        Npoints = 17  # Number of points
        Ntimes = 29  # Number of epochs
        time = np.arange(Ntimes)
        complex = rng.uniform(-1, 1, (Npoints, Ntimes)) + 1j * rng.uniform(-1, 1, (Npoints, Ntimes))
        phase = np.angle(complex)
        h2ph = rng.uniform(1e3, 1e4, (Npoints, Ntimes))

        # Create the points
        stm_pnts = xr.Dataset(
            data_vars={
                "phase": (("space", "time"), phase),
                "h2ph": (("space", "time"), h2ph),
                "complex": (("space", "time"), complex),
                "ambiguities_true": (
                    ("space", "time"),
                    np.round(rng.normal(0, 0.5, (Npoints, Ntimes))).astype(int).clip(-1, 1),
                ),
            },
            coords={
                "space": ("space", np.arange(Npoints)),
                "time": ("time", time),
                "azimuth": ("space", np.round(rng.normal(0, 10, (Npoints))).astype(int)),
                "range": ("space", np.round(rng.normal(0, 10, (Npoints))).astype(int)),
            },
            attrs={"wavelength": 0.056},  # Wavelength in meters
        )

        # Construct arcs based on true ambiguities
        # All arcs by default have 0.99 temp_coh
        stm_arcs = form_network(
            stm_pnts,
            key_xcrds="azimuth",
            key_ycrds="range",
            key_phase="phase",
            key_h2ph="h2ph",
            key_Btemporal="time",
            network_method="redundant",
            max_length=30,
        )
        temp_coh = np.zeros((stm_arcs.sizes["space"],)) + 0.99
        stm_arcs["temp_coh"] = (("space"), temp_coh)

        # Compute arc ambiguities from true point ambiguities
        ambigs = (
            stm_pnts["ambiguities_true"].values[stm_arcs["target"].values, :]
            - stm_pnts["ambiguities_true"].values[stm_arcs["source"].values, :]
        )
        # Introduce some errors in ambiguities
        ambigs_errors = np.zeros_like(ambigs)
        for idx_s, idx_t, err in zip(idx_err_space, idx_err_time, error_values, strict=False):
            ambigs_errors[idx_s, idx_t] += err
        stm_arcs["ambiguities"] = (("space", "time"), ambigs + ambigs_errors)

        stm_arcs_output, stm_pnts_output = spatial_integration(
            stm_pnts, stm_arcs, idx_refpnt=id_ref, key_sdphase="phase", skip_network_adjustment=skip_network_adjustment
        )

        # Verify output dimensions, no points should be rejected
        assert stm_pnts_output.sizes["space"] == stm_pnts.sizes["space"]

        # Check that the solved ambiguities match the true ambiguities w.r.t. the reference point
        assert np.allclose(
            stm_pnts_output["ambiguities"].values
            - stm_pnts["ambiguities_true"].values
            + np.tile(stm_pnts["ambiguities_true"].isel(space=id_ref).values, (stm_pnts.sizes["space"], 1)),
            0,
        )

        # Check that the unwrapped phase is correct
        ref_phase = stm_pnts["phase"].isel(space=id_ref).values
        relative_ambiguities = stm_pnts["ambiguities_true"].values - np.tile(
            stm_pnts["ambiguities_true"].isel(space=id_ref).values, (stm_pnts.sizes["space"], 1)
        )  # true ambiguities relative to reference point
        unwrapped_phase_expected = (
            stm_pnts["phase"].values
            + relative_ambiguities * 2 * np.pi
            - np.tile(ref_phase, (stm_pnts.sizes["space"], 1))
        )
        assert np.allclose(stm_pnts_output["unwrapped_phase"].values, unwrapped_phase_expected)

        # Check that the reference point index remains the same
        assert stm_pnts_output.attrs["idx_refpnt"] == id_ref

    @pytest.mark.parametrize("idx_ref", [0, 5, 10, 16])
    def test_spatial_integration_ref_pnt_removed(self, idx_ref):
        """Raise error when reference point is removed"""
        # Set up test parameters
        rng = np.random.default_rng(42)
        Npoints = 17  # Number of points
        Ntimes = 29  # Number of epochs
        time = np.arange(Ntimes)
        complex = rng.uniform(-1, 1, (Npoints, Ntimes)) + 1j * rng.uniform(-1, 1, (Npoints, Ntimes))
        phase = np.angle(complex)
        h2ph = rng.uniform(1e3, 1e4, (Npoints, Ntimes))

        # Create the points
        stm_pnts = xr.Dataset(
            data_vars={
                "phase": (("space", "time"), phase),
                "h2ph": (("space", "time"), h2ph),
                "complex": (("space", "time"), complex),
                "ambiguities_true": (
                    ("space", "time"),
                    np.round(rng.normal(0, 0.5, (Npoints, Ntimes))).astype(int).clip(-1, 1),
                ),
            },
            coords={
                "space": ("space", np.arange(Npoints)),
                "time": ("time", time),
                "azimuth": ("space", np.round(rng.normal(0, 10, (Npoints))).astype(int)),
                "range": ("space", np.round(rng.normal(0, 10, (Npoints))).astype(int)),
            },
            attrs={"wavelength": 0.056},
        )

        # Construct arcs based on true ambiguities
        # All arcs by default have 0.99 temp_coh
        stm_arcs = form_network(
            stm_pnts,
            key_xcrds="azimuth",
            key_ycrds="range",
            key_phase="phase",
            key_h2ph="h2ph",
            key_Btemporal="time",
            network_method="redundant",
            max_length=30,
        )
        temp_coh = np.zeros((stm_arcs.sizes["space"],)) + 0.99
        stm_arcs["temp_coh"] = (("space"), temp_coh)

        # Set temp_coh of all arcs connected to reference point to 0.01
        mask_ref_arcs = (stm_arcs["source"] == idx_ref) | (stm_arcs["target"] == idx_ref)
        stm_arcs["temp_coh"] = stm_arcs["temp_coh"].where(~mask_ref_arcs, other=0.01)

        # Compute arc ambiguities from true point ambiguities
        ambigs = (
            stm_pnts["ambiguities_true"].values[stm_arcs["target"].values, :]
            - stm_pnts["ambiguities_true"].values[stm_arcs["source"].values, :]
        )
        stm_arcs["ambiguities"] = (("space", "time"), ambigs)

        with pytest.raises(ValueError):
            spatial_integration(stm_pnts, stm_arcs, idx_refpnt=idx_ref, key_sdphase="phase")

    @pytest.mark.parametrize(
        ["idx_source", "idx_target", "n_points", "idx_refpnt"],
        [
            (np.array([0, 1, 2]), np.array([1, 2, 3]), 4, 0),  # 4 points, 3 arcs
            (np.array([0, 1, 2]), np.array([1, 2, 3]), 7, 0),  # 7 points, 3 arcs
            (np.array([1, 1, 2, 2]), np.array([0, 2, 1, 3]), 4, 2),  # 4 points, 4 arcs, unsorted
            (np.array([0, 0, 0, 1, 1, 2, 2]), np.array([1, 2, 3, 3, 4, 3, 4]), 5, 3),  # 5 points, 6 arcs
        ],
    )
    def test_init_network_relation_matrix(
        self,
        idx_source,
        idx_target,
        n_points,
        idx_refpnt,
    ):
        A = _network_relation_matrix(idx_source, idx_target, n_points, idx_refpnt)

        # Create expected matrix in a for loop
        A_exp = np.zeros((idx_source.shape[0], n_points), dtype=int)
        for i, (src, tgt) in enumerate(zip(idx_source, idx_target, strict=False)):
            A_exp[i, src] = -1
            A_exp[i, tgt] = 1
        A_exp = np.delete(A_exp, idx_refpnt, axis=1)  # Remove reference point column

        assert A.shape == A_exp.shape
        assert np.all(A == A_exp)

    @pytest.mark.parametrize(
        ["idx_source", "idx_target", "n_points", "idx_refpnt"],
        [
            (np.array([0, 1, 2]), np.array([1, 2, 3]), 4, 0),  # 4 points, 3 arcs
            (np.array([0, 1, 2]), np.array([1, 2, 3]), 7, 0),  # 7 points, 3 arcs
            (np.array([1, 1, 2, 2]), np.array([0, 2, 1, 3]), 4, 2),  # 4 points, 4 arcs, unsorted
            (np.array([0, 0, 0, 1, 1, 2, 2]), np.array([1, 2, 3, 3, 4, 3, 4]), 5, 3),  # 5 points, 6 arcs
        ],
    )
    def test_init_network_relation_matrix_sparse(
        self,
        idx_source,
        idx_target,
        n_points,
        idx_refpnt,
    ):
        A = _network_relation_matrix(idx_source, idx_target, n_points, idx_refpnt, sparse_mode=True)

        # Create expected matrix in a for loop
        A_exp = np.zeros((idx_source.shape[0], n_points), dtype=int)
        for i, (src, tgt) in enumerate(zip(idx_source, idx_target, strict=False)):
            A_exp[i, src] = -1
            A_exp[i, tgt] = 1
        A_exp = np.delete(A_exp, idx_refpnt, axis=1)  # Remove reference point column

        assert A.shape == A_exp.shape
        assert np.all(A.todense() == A_exp)
