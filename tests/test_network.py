"""test_network.py"""

import pytest
import xarray as xr

from depsi.network import generate_arcs, stm_to_arcs


@pytest.fixture
def stm_sparse():
    # A sparse STM.
    return xr.open_zarr("tests/data/stm_sparse.zarr")


class TestNetwork:
    def test_generate_arcs_fail(self, stm_sparse):
        # Test min_links must be strictly positive.
        result = generate_arcs(stm_sparse, method="redundant", min_links=0)
        assert result is None

        # Test num_partitions must be strictly positive.
        result = generate_arcs(stm_sparse, method="redundant", num_partitions=0)
        assert result is None

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

    def test_stm_to_srcs_subtract(self, stm_sparse):
        # Generate arcs of a Delaunay network with subtracted phase differences.
        stm_arcs = stm_to_arcs(stm_sparse, network="delaunay", max_length=0.05, difference="subtract")
        assert len(stm_arcs["space"]) == 442

        # TODO(tvl) these stm_to_arcs tests are pretty sparse and only check the size of the network.
        # They should be extended to also test the phase difference computation.
        # However, I don't know how to determine that any of these could be validated as 'correct'.

    def test_stm_to_srcs_conjmult(self, stm_sparse):
        # Generate arcs of a Delaunay network with conjugate multiplication phase differences.
        stm_arcs = stm_to_arcs(stm_sparse, network="delaunay", max_length=0.05, difference="conjmult")
        assert len(stm_arcs["space"]) == 442
