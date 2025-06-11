"""test_network.py"""

import numpy as np
import pytest
import xarray as xr

from depsi.network import _compute_phase_difference, generate_arcs


@pytest.fixture
def stm_sparse():
    # A sparse STM.
    return xr.open_zarr("tests/data/stm_sparse.zarr")


class TestNetwork:
    def test_compute_phase_difference_subtract(self, stm_sparse):
        arcs = generate_arcs(stm_sparse)
        d_phase_subtract_0_0 = _compute_phase_difference(stm_sparse, arcs["source"], arcs["source"], method="subtract")
        d_phase_subtract_0_1 = _compute_phase_difference(stm_sparse, arcs["source"], arcs["target"], method="subtract")

        # Test these phase differences.
        assert d_phase_subtract_0_0.values == pytest.approx(np.zeros(d_phase_subtract_0_0.shape), abs=1e-7)
        assert d_phase_subtract_0_1.values == pytest.approx(np.zeros(d_phase_subtract_0_1.shape), abs=2 * np.pi + 1e-7)

    def test_compute_phase_difference_conjmult(self, stm_sparse):
        arcs = generate_arcs(stm_sparse)
        d_phase_conjmult_0_0 = _compute_phase_difference(stm_sparse, arcs["source"], arcs["source"], method="conjmult")
        d_phase_conjmult_0_1 = _compute_phase_difference(stm_sparse, arcs["source"], arcs["target"], method="conjmult")

        # Test these phase differences.
        assert d_phase_conjmult_0_0 == pytest.approx(np.zeros(d_phase_conjmult_0_0.shape), abs=1e-7)
        assert d_phase_conjmult_0_1 == pytest.approx(np.zeros(d_phase_conjmult_0_1.shape), abs=np.pi + 1e-7)

    def test_stm_to_srcs_subtract(self, stm_sparse):
        # Generate arcs of a Delaunay network with subtracted phase differences.
        stm_arcs = generate_arcs(stm_sparse, network_method="delaunay", max_length=0.05, difference="subtract")

        assert len(stm_arcs["space"]) == 442
        assert all(
            [all([-2 * np.pi <= phase <= 2 * np.pi for phase in phases]) for phases in stm_arcs["d_phase"].values]
        )

    def test_stm_to_srcs_conjmult(self, stm_sparse):
        # Generate arcs of a Delaunay network with conjugate multiplication phase differences.
        stm_arcs = generate_arcs(stm_sparse, network_method="delaunay", max_length=0.05, difference="conjmult")

        assert len(stm_arcs["space"]) == 442
        assert all([all([-np.pi <= phase <= np.pi for phase in phases]) for phases in stm_arcs["d_phase"].values])

    def test_stm_to_srcs_fail(self, stm_sparse):
        # Test incorrect method fail.
        with pytest.raises(NotImplementedError):
            generate_arcs(stm_sparse, network_method="unknown", difference="subtract")
        with pytest.raises(NotImplementedError):
            generate_arcs(stm_sparse, network_method="delaunay", difference="unknown")
