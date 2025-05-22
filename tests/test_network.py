"""test_network.py"""

import numpy as np
import pytest
import xarray as xr

from depsi.network import _compute_phase_difference, generate_arcs, stm_to_arcs


@pytest.fixture
def stm_sparse():
    # A sparse STM.
    return xr.open_zarr("tests/data/stm_sparse.zarr")


@pytest.fixture
def stm_sparse_arcs_unzipped(stm_sparse):
    # A sparse network.
    _, arcs = generate_arcs(stm_sparse, method="delaunay", max_length=0.05)

    # Compute the phase differences between points and themselves or their neighbors using different methods.
    arcs_unzipped = list(zip(*arcs, strict=False))
    arcs_unzipped = [list(arcs_unzipped[0]), list(arcs_unzipped[1])]

    return arcs_unzipped


class TestNetwork:
    def test_generate_arcs_fail(self, stm_sparse):
        # Test min_links must be strictly positive.
        result = generate_arcs(stm_sparse, method="redundant", min_links=0)
        assert result is None

        # Test num_partitions must be strictly positive.
        result = generate_arcs(stm_sparse, method="redundant", num_partitions=0)
        assert result is None

        # Test incorrect method fail.
        with pytest.raises(NotImplementedError):
            result = generate_arcs(stm_sparse, method="unknown")

    def test_generate_arcs_delaunay(self, stm_sparse):
        # Generate a Delaunay network with long edges removed.
        coordinates, arcs = generate_arcs(stm_sparse, method="delaunay", max_length=0.05)

        assert len(coordinates) == 156
        assert len(arcs) == 442

    def test_generate_arcs_redundant(self, stm_sparse):
        # Generate a 'redundant' network with long edges removed.
        coordinates, arcs = generate_arcs(
            stm_sparse, method="redundant", max_length=0.05, min_links=8, num_partitions=4
        )

        assert len(coordinates) == 156
        assert len(arcs) == 797

    def test_compute_phase_difference_subtract(self, stm_sparse, stm_sparse_arcs_unzipped):
        arcs = stm_sparse_arcs_unzipped
        d_phase_subtract_0_0 = _compute_phase_difference(stm_sparse, arcs[0], arcs[0], method="subtract")
        d_phase_subtract_0_1 = _compute_phase_difference(stm_sparse, arcs[0], arcs[1], method="subtract")

        # Test these phase differences.
        assert d_phase_subtract_0_0.values == pytest.approx(np.zeros(d_phase_subtract_0_0.shape), abs=1e-7)
        assert d_phase_subtract_0_1.values == pytest.approx(np.zeros(d_phase_subtract_0_1.shape), abs=2 * np.pi + 1e-7)

    def test_compute_phase_difference_conjmult(self, stm_sparse, stm_sparse_arcs_unzipped):
        arcs = stm_sparse_arcs_unzipped
        d_phase_conjmult_0_0 = _compute_phase_difference(stm_sparse, arcs[0], arcs[0], method="conjmult")
        d_phase_conjmult_0_1 = _compute_phase_difference(stm_sparse, arcs[0], arcs[1], method="conjmult")

        # Test these phase differences.
        assert d_phase_conjmult_0_0 == pytest.approx(np.zeros(d_phase_conjmult_0_0.shape), abs=1e-7)
        assert d_phase_conjmult_0_1 == pytest.approx(np.zeros(d_phase_conjmult_0_1.shape), abs=np.pi + 1e-7)

    def test_stm_to_srcs_subtract(self, stm_sparse):
        # Generate arcs of a Delaunay network with subtracted phase differences.
        stm_arcs = stm_to_arcs(stm_sparse, network="delaunay", max_length=0.05, difference="subtract")

        assert len(stm_arcs["space"]) == 442
        assert all(
            [all([-2 * np.pi <= phase <= 2 * np.pi for phase in phases]) for phases in stm_arcs["d_phase"].values]
        )

    def test_stm_to_srcs_conjmult(self, stm_sparse):
        # Generate arcs of a Delaunay network with conjugate multiplication phase differences.
        stm_arcs = stm_to_arcs(stm_sparse, network="delaunay", max_length=0.05, difference="conjmult")

        assert len(stm_arcs["space"]) == 442
        assert all([all([-np.pi <= phase <= np.pi for phase in phases]) for phases in stm_arcs["d_phase"].values])

    def test_stm_to_srcs_fail(self, stm_sparse):
        # Test incorrect method fail.
        with pytest.raises(NotImplementedError):
            stm_to_arcs(stm_sparse, network="unknown", difference="subtract")
        with pytest.raises(NotImplementedError):
            stm_to_arcs(stm_sparse, network="delaunay", difference="unknown")
