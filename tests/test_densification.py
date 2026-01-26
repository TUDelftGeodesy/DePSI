import numpy as np
import pytest
import xarray as xr

from depsi.densification import densification
from depsi.network import (
    form_network,
    spatial_integration,
)

# Random number generator
rng = np.random.default_rng(42)

# Constants in this test
wavelength = 0.055465763  # Sentinel-1, in meters
m2ph = -4 * np.pi / wavelength
az_range = (0, 100)  # azimuth coordinate range
rg_range = (0, 100)  # range coordinate range
h2ph_range = (1e-5, 1e-4)  # height to phase conversion factor range
vel_range = (-0.01, 0.01)  # velocity range in m/yr
height_range = (1e-1, 1e0)  # height range in meters
time_range = 3.0  # time span in years


def make_stm_network_pnts(npoints_net, ntime):
    """function for preparing solved network points."""
    # Network points true parameters
    time = np.linspace(0, time_range, ntime)  # Time in years
    pnt_vel_true = rng.uniform(vel_range[0], vel_range[1], (npoints_net,))  # Random velocities in m/yr
    pnt_height_true = rng.uniform(height_range[0], height_range[1], (npoints_net,))  # Random heights in m
    h2ph = rng.uniform(h2ph_range[0], h2ph_range[1], (npoints_net, ntime))  # Random height to phase conversion factors
    # Spatial unwrapping parameters
    id_ref = 3  # Reference point index in the network
    arc_length_max = 70  # unitless based on az and range coords

    # Simulate network points phase and ambiguities
    phase_height = h2ph * m2ph * np.tile(pnt_height_true, (ntime, 1)).T
    phase_vel = m2ph * np.tile(time, (npoints_net, 1)) * np.tile(pnt_vel_true, (ntime, 1)).T
    phase_true = phase_height + phase_vel
    phase = (np.mod(phase_true + np.pi, 2 * np.pi)) - np.pi  # wrap the phase to [-pi, pi]
    sd_phase = phase - np.tile(phase[:, 0], (ntime, 1)).T  # single difference phase reference to first epoch
    ambiguities_true = np.round((phase_true - phase) / (2 * np.pi)).astype(int)

    # Create the network points
    stm_pnts_network = xr.Dataset(
        data_vars={
            "sd_phase": (("space", "time"), sd_phase),
            "phase_true": (("space", "time"), phase_true),
            "h2ph": (("space", "time"), h2ph),
            "ambiguities_true": (
                ("space", "time"),
                ambiguities_true,
            ),
        },
        coords={
            "space": ("space", np.arange(npoints_net)),
            "time": ("time", time),
            "azimuth": ("space", rng.choice(np.arange(az_range[0], az_range[1]), size=npoints_net, replace=False)),
            "range": ("space", rng.choice(np.arange(rg_range[0], rg_range[1]), size=npoints_net, replace=False)),
        },
    )

    stm_pnts_network = stm_pnts_network.assign_attrs({"wavelength": wavelength})

    # Construct arcs based on true ambiguities
    stm_arcs_network = form_network(
        stm_pnts_network,
        key_xcrds="azimuth",
        key_ycrds="range",
        key_phase="sd_phase",
        key_h2ph="h2ph",
        key_Btemporal="time",
        network_method="redundant",
        max_length=arc_length_max,
    )
    temp_coh = np.zeros((stm_arcs_network.sizes["space"],)) + 0.99  # All arcs by default have 0.99 temp_coh
    stm_arcs_network["temp_coh"] = (("space"), temp_coh)

    # Compute arc ambiguities from true point ambiguities
    # This assumes no mistakes in network arc periodogram estimation
    ambigs = (
        stm_pnts_network["ambiguities_true"].values[stm_arcs_network["target"].values, :]
        - stm_pnts_network["ambiguities_true"].values[stm_arcs_network["source"].values, :]
    )
    stm_arcs_network["ambiguities"] = (("space", "time"), ambigs)

    # Network integration based on true arc ambiguities
    stm_arcs_network_output, stm_pnts_network_output = spatial_integration(
        stm_pnts_network, stm_arcs_network, idx_refpnt=id_ref
    )

    return stm_pnts_network_output


def make_stm_pnt_densification(npoints_dens, ntime):
    # Make densification datasets
    time = np.linspace(0, time_range, ntime)  # Time in years
    pnt_vel_true = rng.uniform(vel_range[0], vel_range[1], (npoints_dens,))  # Random velocities in m/yr
    pnt_height_true = rng.uniform(height_range[0], height_range[1], (npoints_dens,))  # Random heights in m
    h2ph = rng.uniform(h2ph_range[0], h2ph_range[1], (npoints_dens, ntime))

    # Simulate phases and ambiguities
    phase_height = h2ph * m2ph * np.tile(pnt_height_true, (ntime, 1)).T
    phase_vel = m2ph * np.tile(time, (npoints_dens, 1)) * np.tile(pnt_vel_true, (ntime, 1)).T
    phase_true = phase_height + phase_vel
    phase = (np.mod(phase_true + np.pi, 2 * np.pi)) - np.pi  # wrap the phase to [-pi, pi]
    sd_phase = phase - np.tile(phase[:, 0], (ntime, 1)).T  # reference to first epoch
    ambiguities_true = np.round((phase_true - phase) / (2 * np.pi)).astype(int)

    # Create densification points dataset
    stm_pnt_densification = xr.Dataset(
        data_vars={
            "sd_phase": (("space", "time"), sd_phase),
            "phase_true": (("space", "time"), phase_true),
            "h2ph": (("space", "time"), h2ph),
            "ambiguities_true": (
                ("space", "time"),
                ambiguities_true,
            ),
        },
        coords={
            "space": ("space", np.arange(npoints_dens)),
            "time": ("time", time),
            "azimuth": ("space", rng.choice(np.arange(az_range[0], az_range[1]), size=npoints_dens, replace=False)),
            "range": ("space", rng.choice(np.arange(rg_range[0], rg_range[1]), size=npoints_dens, replace=False)),
        },
    )

    return stm_pnt_densification


@pytest.mark.parametrize("npoints_net, npoints_dens, ntime", [(17, 7, 24), (21, 13, 16)])
def test_densification(npoints_net, npoints_dens, ntime):
    """Test densification function by comparing unwrapped phase to true double-difference phase."""
    stm_network_pnts = make_stm_network_pnts(npoints_net, ntime)

    stm_pnt_densification = make_stm_pnt_densification(npoints_dens, ntime)

    # Densification parameters
    stm_densified = densification(
        stm_pnt_densification,
        stm_network_pnts,
        n_connections=1,
        key_xcoord="azimuth",
        key_ycoord="range",
        key_Btemporal="time",
        key_h2ph="h2ph",
        key_sdphase="sd_phase",
    )

    # Densification should inherit reference point index attribute from network points
    assert "idx_refpnt" in stm_densified.attrs

    # Following variables should be in the densification output
    assert "ambiguities" in stm_densified.data_vars
    assert "unwrapped_phase" in stm_densified.data_vars
    assert "local_temp_coh" in stm_densified.data_vars

    # Validate on unwrapped phase
    # The densification output unwrapped phase should be close to true double-difference phase
    unw_phase = stm_densified["unwrapped_phase"].data
    true_phase = stm_densified["phase_true"].data
    # Reference phase from network points
    # Reference phase are also true phase therefore it needs to be single-differenced
    reference_phase = stm_network_pnts["phase_true"].data[stm_densified.attrs["idx_refpnt"], :]
    sd_reference_phase = reference_phase - reference_phase[0]
    # Get true double-difference phase
    true_phase_dd = (
        true_phase
        - np.tile(true_phase[:, 0], (unw_phase.shape[1], 1)).T
        - np.tile(sd_reference_phase, (unw_phase.shape[0], 1))
    )
    # Verify the unwrapped phase is close to true double-difference phase
    assert np.allclose(unw_phase - true_phase_dd, 0)


def test_densification_inconsistent_time():
    """Value error raised when time coordinates are inconsistent."""
    ntime_net = 5
    ntime_dens = 7  # Different number of time points
    npoints_net = 11
    npoints_dens = 5

    stm_network_pnts = make_stm_network_pnts(npoints_net, ntime_net)
    stm_pnt_densification = make_stm_pnt_densification(npoints_dens, ntime_dens)

    with pytest.raises(ValueError):
        _ = densification(
            stm_pnt_densification,
            stm_network_pnts,
        )


def test_densification_missing_wavelength():
    """Value error raised when wavelength not provided as either argument or attribute."""
    ntime = 5

    npoints_net = 11
    npoints_dens = 5

    stm_network_pnts = make_stm_network_pnts(npoints_net, ntime)
    stm_pnt_densification = make_stm_pnt_densification(npoints_dens, ntime)

    # Remove wavelength attribute
    stm_network_pnts.attrs.pop("wavelength")

    with pytest.raises(ValueError):
        _ = densification(
            stm_pnt_densification,
            stm_network_pnts,
        )


def test_densification_missing_idx_refpnt():
    """Value error raised when index of reference point idx_refpnt not provided as either argument or attribute."""
    ntime = 5
    npoints_net = 11
    npoints_dens = 5

    stm_network_pnts = make_stm_network_pnts(npoints_net, ntime)
    stm_pnt_densification = make_stm_pnt_densification(npoints_dens, ntime)

    # Remove idx_refpnt attribute
    stm_network_pnts.attrs.pop("idx_refpnt")

    with pytest.raises(ValueError):
        _ = densification(
            stm_pnt_densification,
            stm_network_pnts,
        )


@pytest.mark.parametrize("n_connections", [-1, 4])
def test_densification_wrong_n_connections(n_connections):
    """Value error raised when n_connections is not a positive odd number."""
    ntime = 5
    npoints_net = 11
    npoints_dens = 5

    stm_network_pnts = make_stm_network_pnts(npoints_net, ntime)
    stm_pnt_densification = make_stm_pnt_densification(npoints_dens, ntime)

    with pytest.raises(ValueError):
        _ = densification(
            stm_pnt_densification,
            stm_network_pnts,
            n_connections=n_connections,
        )


def test_densification_time_variables():
    """Time-only variables should be preserved with correct dimensions after densification."""
    stm_network_pnts = make_stm_network_pnts(6, 7)
    stm_pnt_densification = make_stm_pnt_densification(13, 7)

    # Assign some time-only variables
    time_var_common = ("time", rng.uniform(0, 1, size=stm_network_pnts.sizes["time"]))
    stm_network_pnts["time_var_common"] = time_var_common
    stm_pnt_densification["time_var_common"] = time_var_common
    stm_network_pnts["time_var_net"] = ("time", rng.uniform(0, 1, size=stm_network_pnts.sizes["time"]))
    stm_pnt_densification["time_var_dens"] = ("time", rng.uniform(0, 1, size=stm_pnt_densification.sizes["time"]))

    stm_densified = densification(
        stm_pnt_densification,
        stm_network_pnts,
        n_connections=1,
        key_xcoord="azimuth",
        key_ycoord="range",
        key_Btemporal="time",
        key_h2ph="h2ph",
        key_sdphase="sd_phase",
    )

    # Check that time-only variables are present in the output
    for var in ["time_var_common", "time_var_net", "time_var_dens"]:
        assert var in stm_densified.data_vars
        assert stm_densified[var].dims == ("time",)

    assert np.allclose(
        stm_densified["time_var_net"].data,
        stm_network_pnts["time_var_net"].data,
    )
    assert np.allclose(
        stm_densified["time_var_dens"].data,
        stm_pnt_densification["time_var_dens"].data,
    )
    assert np.allclose(
        stm_densified["time_var_common"].data,
        stm_network_pnts["time_var_common"].data,
    )
    assert np.allclose(
        stm_densified["time_var_common"].data,
        stm_pnt_densification["time_var_common"].data,
    )
