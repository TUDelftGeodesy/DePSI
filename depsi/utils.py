import math
from typing import Literal

import numpy as np

EARTH_RADIUS = 6378136  # m


def _orbit_fit(orbit, verbose=0, der=True):
    """Return a orbit_fit dict.

    Modified from the "orbitFit" function:
    https://github.com/Pbaz98/Caroline-Radar-Coding-Toolbox/blob/main/gecoris/geoUtils.py#L325

    Satellite state vector interpolation using Chebyshev polynomials of
    7th order (according to DLR recommendations). Function returns Chebyshev
    polynomial coefficients. Use these to evaluate orbit state at given time
    via function 'orbitVal'.

    input: snappy 'orbit' object (as read by 'read_metadata' function)

    CHANGE LOG
    - 30/6/2023: Modified to adapt the input to a np.array Nx4 (N number of timesamples)
    - 22/09/23: add the flag for derivative or not
    """
    # parse masterorb:
    t = orbit[:, 0]
    x = orbit[:, 1]
    y = orbit[:, 2]
    z = orbit[:, 3]

    # interpolate orbits using Chebyshev polynomials of 7th order:
    t0 = (min(t) + max(t)) / 2
    px = t - t0  # time argument px (centered around mid interval)
    cx = np.polynomial.chebyshev.chebfit(px, x, 7)  # position
    cy = np.polynomial.chebyshev.chebfit(px, y, 7)
    cz = np.polynomial.chebyshev.chebfit(px, z, 7)

    if der:
        cvx = np.polynomial.chebyshev.chebder(cx)  # velocity
        cvy = np.polynomial.chebyshev.chebder(cy)
        cvz = np.polynomial.chebyshev.chebder(cz)
    else:
        x_vel = orbit[:, 4]
        y_vel = orbit[:, 5]
        z_vel = orbit[:, 6]

        cvx = np.polynomial.chebyshev.chebfit(px, x_vel, 7)  # velocity
        cvy = np.polynomial.chebyshev.chebfit(px, y_vel, 7)
        cvz = np.polynomial.chebyshev.chebfit(px, z_vel, 7)

    cax = np.polynomial.chebyshev.chebder(cvx)  # acceleration
    cay = np.polynomial.chebyshev.chebder(cvy)
    caz = np.polynomial.chebyshev.chebder(cvz)

    if verbose:
        # position fit residuals:
        x_res = np.polynomial.chebyshev.chebval(px, cx) - x
        y_res = np.polynomial.chebyshev.chebval(px, cy) - y
        z_res = np.polynomial.chebyshev.chebval(px, cz) - z
        x_std = np.std(x_res)
        y_std = np.std(y_res)
        z_std = np.std(z_res)
        print(f"Orbit fit position residuals: X {x_std:.4f} m, Y {y_std:.4f} m, Z {z_std:.4f} m. ")
        # velocity residuals:
        vx_res = np.polynomial.chebyshev.chebval(px, np.polynomial.chebyshev.chebder(cx)) - x_vel
        vy_res = np.polynomial.chebyshev.chebval(px, np.polynomial.chebyshev.chebder(cy)) - y_vel
        vz_res = np.polynomial.chebyshev.chebval(px, np.polynomial.chebyshev.chebder(cz)) - z_vel
        vx_std = np.std(vx_res)
        vy_std = np.std(vy_res)
        vz_std = np.std(vz_res)
        print(f"Orbit fit velocity residuals: vX {vx_std:.4f} m/s, vY {vy_std:.4f} m/s, vZ {vz_std:.4f} m/s. ")

    orbit_fit = dict()
    orbit_fit["t0"] = t0
    orbit_fit["cx"] = cx
    orbit_fit["cy"] = cy
    orbit_fit["cz"] = cz
    orbit_fit["cvx"] = cvx
    orbit_fit["cvy"] = cvy
    orbit_fit["cvz"] = cvz
    orbit_fit["cax"] = cax
    orbit_fit["cay"] = cay
    orbit_fit["caz"] = caz

    return orbit_fit


def get_distance(
    s: list | tuple | np.ndarray, t: list | tuple | np.ndarray, mode: Literal["euclidean", "geographic"] = "euclidean"
):
    """Calculate the distance between two points.

    Parameters
    ----------
    s: list | tuple | np.ndarray
        The source point, formatted as (x, y) / (lon, lat)
    t: list | tuple | np.ndarray
        The target point, formatted as (x, y) / (lon, lat)
    mode: Literal["euclidean", "geographic"], default "euclidean"
        Whether the source and target points are given in (x, y) (units meters) or (lon, lat) (units degrees)

    Returns
    -------
        The distance between the two points in meters.
    """
    # TODO(tvl) More complex distance functions can be implemented here.
    #
    # For example, network generation gives better results for square Euclidean distances.
    # Non-square coordinate systems (like non-square image coordinates) may distort the circle properties of Delaunay
    # networks into ellipses.
    # Non-Euclidean coordinate systems (like angular lat-lon systems) may distort these same properties depending on
    # the distance to a pole.
    #
    # Coordinate system transformations may be done on the STM before generating the network.
    # However, there may be cases where it is impossible or undesirable to transform the point coordinates in the STM.
    # In such a case, coordinates may be transformed inside this function.
    if mode == "euclidean":
        return math.dist(s, t)
    elif mode == "geographic":
        # this is the Haversine formula
        lat1 = s[1]
        lat2 = t[1]
        dphi = np.radians(lat1 - lat2)
        dlambda = np.radians(s[0] - t[0])
        dist = (
            2
            * EARTH_RADIUS
            * np.arcsin(
                np.sqrt(
                    (1 - np.cos(dphi) + np.cos(np.radians(lat1)) * np.cos(np.radians(lat2)) * (1 - np.cos(dlambda))) / 2
                )
            )
        )
        return dist
    raise ValueError(f"Unknown mode {mode}, only know euclidean and geographic!")
