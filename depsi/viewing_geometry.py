import numpy as np


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
        aray-like: predicted value for the angle given x and y coordinates
    """
    a, b, c = coeffs
    return a * x + b * y + c
