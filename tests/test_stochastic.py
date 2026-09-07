import numpy as np
import pytest
import xarray as xr

from depsi.stochastic import _q_no_mother_atmo, _q_with_mother_atmo, vce_temporal


def simulated_stm(n_ifg, n_points):
    """Function to simulate arcs stm for testing."""
    rng = np.random.default_rng(31)

    # Simulate coordinates
    lat = np.linspace(51.14, 51.15, n_points)
    lon = rng.uniform(6.9, 7.0, n_points)

    # Simulate Btemp
    TIME_STEP = 15  # Time step in days
    years = np.linspace(0, (n_ifg - 1) * TIME_STEP / 365.25, n_ifg)  # years from 0 to n_ifg-1

    # Simulated phase values, linear+ noise
    phase = np.tile(np.linspace(-np.pi, np.pi, n_ifg), (n_points, 1)) + rng.normal(0, 0.1, (n_points, n_ifg))

    # Simulated h2ph values
    h2ph = rng.random((n_points, n_ifg)) * 1e-3  # fixed h2ph values for each arc

    arcs = xr.Dataset(
        data_vars={
            "lon": (("space",), lon),
            "lat": (("space",), lat),
            "phase": (("space", "time"), phase),
            "h2ph": (("space", "time"), h2ph),
            "Btemporal": (("time",), years),
        }
    )

    arcs.attrs["wavelength"] = 0.056  # example wavelength in meters

    return arcs


@pytest.mark.parametrize(
    ["n_ifg", "n_points", "include_mother_atmo"],
    [
        (12, 41, False),
        (19, 107, False),
        (42, 127, True),
    ],
)
def test_vce_temporal(n_ifg, n_points, include_mother_atmo):
    arcs = simulated_stm(n_ifg, n_points)
    sigma2 = vce_temporal(
        arcs, key_phase="phase", key_Btemporal="Btemporal", key_h2ph="h2ph", include_mother_atmo=include_mother_atmo
    )

    assert sigma2.shape[0] == n_ifg + 1
    if include_mother_atmo:
        assert sigma2[0] != 0
    else:
        assert sigma2[0] == 0


@pytest.mark.parametrize("Nifgs", [4, 11, 23])
def test_q_no_mother_atmo(Nifgs):
    """Test the _q_with_mother_atmo function."""
    Qy1, Qy = _q_no_mother_atmo(Nifgs)

    assert Qy1.shape == (Nifgs, Nifgs, Nifgs)
    assert Qy.shape == (Nifgs, Nifgs)
    assert np.all(np.isin(np.unique(Qy1.flatten()), [0, 2]))
    assert np.all(Qy == np.diag(np.diag(Qy)))  # check Qy is a diagonal matrix


@pytest.mark.parametrize("Nifgs", [4, 11, 23])
def test_q_with_mother_atmo(Nifgs):
    """Test the _q_with_mother_atmo function."""
    Qy1, Qy = _q_with_mother_atmo(Nifgs)

    assert Qy1.shape == (Nifgs, Nifgs, Nifgs + 1)
    assert Qy.shape == (Nifgs, Nifgs)
    assert np.all(np.isin(np.unique(Qy1.flatten()), [0, 2]))
    assert np.all(np.isin(np.unique(Qy1[:, :, 0].flatten()), [2]))  # check first slice is all 2
