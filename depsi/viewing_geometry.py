from typing import Literal

import drama.mission.timeline.no_drama as nd
import numpy as np
import xarray as xr
from drama.io import cfg
from drama.performance.sar import SARModeFromCfg

from depsi.utils import identify_s1_orbits_in_aoi


def add_local_viewing_geometry(
    stm: xr.Dataset, orbit_config_file: str, orbit_res: float, orbit_mode: Literal["IWS"], orbit: str
) -> xr.Dataset:
    """Add the local incidence angle and alpha to each point in the STM.

    Parameters
    ----------
    stm: xr.Dataset
        Space-time matrix with at least the following properties:
        variables: lon, lat
    orbit_config_file: str
        Full path to S1_XTI.cfg (provided with the DePSI package in config/drama)
    orbit_res: float
        Orbital resolution in m (typically 0.01)
    orbit_mode: Literal["IWS"]
        Orbit mode of the satellite. For S1 this is generally IWS. Others not implemented
    orbit: str
        String indicating which orbit we are looking at, of the form `s1_dsc_t110`

    Returns
    -------
    xr.Dataset
        Space-time matrix with the local incidence angle and azimuth angle per point added
    """
    assert orbit_mode == "IWS", f"orbit mode {orbit_mode} is not implemented!"

    # necessary for DRaMA deprecation error avoidance
    np.NaN = np.nan

    # get the footprints of the overlapping AoIs
    orbs, footprints = identify_s1_orbits_in_aoi(stm["lon"].values, stm["lat"].values)
    assert orbit in orbs, f"Orbit provided is {orbit} but only found orbits {orbs}..."

    # get the orbit and the local incidence angles
    mode = SARModeFromCfg(cfg.ConfigFile(orbit_config_file), orbit_mode)
    nr_asc_obs, nr_desc_obs, asc_inc, desc_inc, asc_alpha, desc_alpha = nd.viewing_geometry(
        orbit_config_file, mode, orbit_res, stm["lat"].values, stm["lon"].values
    )

    # figure out which orbit is the correct one
    asc_dsc = orbit.split("_")[1]
    correct_direction = [o for o in orbs if asc_dsc in o]
    other_same_dir_orbit = [c for c in correct_direction if c != orbit]
    if len(other_same_dir_orbit) == 0:
        if asc_dsc == "dsc":
            incs = desc_inc[:, 0]
            alphas = desc_alpha[:, 0]
        else:
            incs = asc_inc[:, 0]
            alphas = asc_alpha[:, 0]
    elif len(other_same_dir_orbit) == 1:
        other_same_dir_orbit = other_same_dir_orbit[0]
        footprint_corr = footprints[orbit][0]
        footprint_other = footprints[other_same_dir_orbit][0]
        leftmost_corr = min([i[0] for i in footprint_corr])
        leftmost_other = min([i[0] for i in footprint_other])
        if asc_dsc == "dsc":
            comp_row = desc_inc[0]
            if leftmost_corr > leftmost_other:  # the orbit of interest is further east. For DSC, this is larger inc
                if comp_row[0] > comp_row[1]:
                    index = 0
                else:
                    index = 1
            else:  # the orbit of interest is further west. For DSC, this is smaller inc angle
                if comp_row[0] > comp_row[1]:
                    index = 1
                else:
                    index = 0
            incs = desc_inc[:, index]
            alphas = desc_alpha[:, index]
        else:  # asc
            comp_row = asc_inc[0]
            if leftmost_corr > leftmost_other:  # the orbit of interest is further east. For ASC, this is smaller inc
                if comp_row[0] > comp_row[1]:
                    index = 1
                else:
                    index = 0
            else:  # the orbit of interest is further west. For ASC, this is larger inc angle
                if comp_row[0] > comp_row[1]:
                    index = 0
                else:
                    index = 1
            incs = asc_inc[:, index]
            alphas = asc_alpha[:, index]
    else:
        raise ValueError(f"Cannot handle more than two orbits in the same direction, found {orbs}!")

    stm = stm.assign({"local_incidence_angle": (["space"], incs * 180 / np.pi)})
    stm = stm.assign({"local_azimuth_angle": (["space"], alphas * 180 / np.pi)})
    return stm


def fit_plane_viewing_geometry(x, y, angle):
    """Fit a plane based on x and y coordinates and the corresponding incidence angle or azimuth of the ZDP.

    Args:
    ----
        x (array-like): x coordinates
        y (array-like): y coordinates
        angle (array-like): incidence angles or azimuth Zero-Doppler plane

    Returns:
    -------
        float: coefficients of the plane
    """
    A = np.column_stack((x, y, np.ones_like(x)))

    # Solve linear system to estimate coefficients a, b, c
    coeffs, _, _, _ = np.linalg.lstsq(A, angle, rcond=None)
    return coeffs


def estimate_plane_viewing_geometry(x, y, coeffs):
    """Estimate value for the incidence angle or azimuth of the ZDP based on x and y coordinates.

    Args:
    ----
        x (array-like): x coordinates
        y (array-like): y coordinates
        coeffs (floats): a,b,c coefficients describing the plane equations

    Returns:
    -------
        array-like: predicted value for the angle given x and y coordinates
    """
    a, b, c = coeffs
    return a * x + b * y + c
