"""Geocoding module.

Functions for transformations between WGS84 coordinates (lat, lon, ellipsoid height) and
Radar Coordinates (azimuth, range, ellipsoid height), as well as back and forth between ECEF coordinates (X, Y, Z)
and time coordinates (time).

When transferring from Radar to WGS84, the following steps are performed (in radar_to_latlonh):
1. radar_to_xyz (calling radar_to_time internally)
2. xyz_to_latlonh

When transfering from WGS84 to Radar, the following steps are performed (in latlonh_to_radar or latlonh_to_radar_vec):
1. latlonh_to_xyz
2. xyz_to_time
3. time_to_radar
"""

import collections
import datetime
import logging
from typing import Literal

import numpy
import pyproj
from numpy.polynomial.chebyshev import chebval

from depsi.utils import npdatetime64_to_datetime

logger = logging.getLogger(__name__)

SPEED_OF_LIGHT = 299792458.0
MJD_EPOCH = datetime.datetime(2000, 1, 1)
WGS84 = pyproj.CRS.from_epsg(4326).ellipsoid

OrbitFit = collections.namedtuple("OrbitFit", "time0, cx, cy, cz, cvx, cvy, cvz, cax, cay, caz")

VALIDATION_NONORBIT_KEYS = {
    "radar_to_latlonh": [
        "scene_centre_longitude",
        "scene_centre_latitude",
        "first_azimuth_time",
        "pulse_repetition_frequency",
        "first_range_time",
        "range_sampling_rate",
    ],
    "latlonh_to_radar": ["first_azimuth_time", "number_of_rows", "pulse_repetition_frequency"],
}


def seconds_of_day(t: datetime.datetime) -> float:
    """Calculate the number of seconds elapsed since the start of the day for a given epoch.

    Parameters
    ----------
    t: datetime.datetime
        Epoch to be converted

    Returns
    -------
    float
        Number of seconds elapsed since the start of that day
    """
    total_seconds = (t - t.replace(hour=0, minute=0, second=0, microsecond=0)).total_seconds()
    return total_seconds


def latlonh_to_xyz(latlonh: numpy.ndarray) -> numpy.ndarray:
    """Transform WGS84 coordinates to ECEF coordinates.

    Transformation of ellipsoidal to Cartesian geocentric (ECEF) coordinates using WGS84/GRS80 ellipsoid.

    Parameters
    ----------
    latlonh: numpy.ndarray
        latitude, longitude, ellipsoidal height coordinates as a numpy array in [degrees/m] of shape (N, 3).

    Returns
    -------
    xyz: numpy.ndarray
        Cartesian geocentric coordinates in [m].
    """
    transformer = pyproj.Transformer.from_crs("EPSG:4979", "EPSG:4978", always_xy=True)
    x, y, z = transformer.transform(latlonh[:, 1], latlonh[:, 0], latlonh[:, 2])

    return numpy.array([x, y, z]).squeeze()


def xyz_to_latlonh(xyz: numpy.ndarray) -> numpy.ndarray:
    """Transform ECEF Coordinates to WGS84 coordinates.

    Transformation of geocentric Cartesian (ECEF) coordinates to WGS84 ellipsoidal coordinates
    using iteration method.

    Parameters
    ----------
    xyz: numpy.ndarray
        x, y, z Cartesian geocentric coordinates in metres of shape (N, 3).

    Returns
    -------
    latlonh: numpy.ndarray
        Latitude, longitude, ellipsoidal height coordinates as degrees and metres.
    """
    transformer = pyproj.Transformer.from_crs("EPSG:4978", "EPSG:4979", always_xy=True)
    lon, lat, h = transformer.transform(xyz[:, 0], xyz[:, 1], xyz[:, 2])

    return numpy.array([lat, lon, h]).squeeze()


def orbit_fit(orbit_time: numpy.ndarray, xyz_pos: numpy.ndarray, xyz_vel: numpy.ndarray | None) -> OrbitFit:
    """Fit orbit based on state vector using Chebyshev polynomials.

    Satellite state vector interpolation using Chebyshev polynomials of
    7th order. Function returns Chebyshev polynomial coefficients. Use
    these to evaluate orbit state at given time via function 'orbit_eval'.

    Parameters
    ----------
    orbit_time: numpy.ndarray
        Time vector of satellite state data extracted from image metadata
    xyz_pos: numpy.ndarray
        Position vector of satellite state data
    xyz_vel: numpy.ndarray | None
        Velocity vector of satellite state data. If `None`, it is estimated from `xyz_pos`

    Returns
    -------
    collections.namedtuple
        Returns OrbitFit instance with interpolated data.

    """
    # time argument px (centered around mid-interval)
    time0 = (min(orbit_time) + max(orbit_time)) / 2
    px = orbit_time - time0

    # fit position
    cx = numpy.polynomial.chebyshev.chebfit(px, xyz_pos[:, 0], 7)
    cy = numpy.polynomial.chebyshev.chebfit(px, xyz_pos[:, 1], 7)
    cz = numpy.polynomial.chebyshev.chebfit(px, xyz_pos[:, 2], 7)

    if xyz_vel is None:  # Doris does not return velocity state data --> we calculate it instead
        cvx = numpy.polynomial.chebyshev.chebder(cx)
        cvy = numpy.polynomial.chebyshev.chebder(cy)
        cvz = numpy.polynomial.chebyshev.chebder(cz)
    else:
        # fit velocity
        cvx = numpy.polynomial.chebyshev.chebfit(px, xyz_vel[:, 0], 7)
        cvy = numpy.polynomial.chebyshev.chebfit(px, xyz_vel[:, 1], 7)
        cvz = numpy.polynomial.chebyshev.chebfit(px, xyz_vel[:, 2], 7)

    # fit acceleration
    cax = numpy.polynomial.chebyshev.chebder(cvx)
    cay = numpy.polynomial.chebyshev.chebder(cvy)
    caz = numpy.polynomial.chebyshev.chebder(cvz)

    return OrbitFit(time0, cx, cy, cz, cvx, cvy, cvz, cax, cay, caz)


def orbit_eval(orbit: OrbitFit, azimuth_time: numpy.ndarray) -> tuple[numpy.ndarray, numpy.ndarray]:
    """Evaluate the orbit state vector of `orbit_fit`.

    Evaluate orbit state vector (XYZ) at given azimuth time using fitted Chebyshev
    polynomials by function 'orbit_fit'.

    Parameters
    ----------
    orbit: collections.namedtuple
        Orbit fit instance generated by 'orbit_fit' function.
    azimuth_time: numpy.ndarray
        Azimuth time in seconds from start of the day.

    Returns
    -------
    tuple
        Satellite position and velocity arrays in geocentric Cartesian system (XYZ).

    """
    # evaluate position at requested azimuth epoch:
    x = chebval(azimuth_time - orbit.time0, orbit.cx)
    y = chebval(azimuth_time - orbit.time0, orbit.cy)
    z = chebval(azimuth_time - orbit.time0, orbit.cz)

    # velocity at epoch t
    vel_x = chebval(azimuth_time - orbit.time0, orbit.cvx)
    vel_y = chebval(azimuth_time - orbit.time0, orbit.cvy)
    vel_z = chebval(azimuth_time - orbit.time0, orbit.cvz)

    return numpy.array([x, y, z]).squeeze(), numpy.array([vel_x, vel_y, vel_z]).squeeze()


def xyz_to_topocentric(xyz: numpy.ndarray, sat_xyz: numpy.ndarray) -> numpy.ndarray:
    """Transform XYZ coordinates to topocentric coordinates (incidence, azimuth, distance).

    Calculate the topocentric coordinates (incidence, azimuth angle, distance) from
    the Cartesian coordinates (XYZ) of topocenter (target) to the satellite

    Parameters
    ----------
    xyz: numpy.ndarray
        x, y, z Cartesian geocentric coordinates of the topocenter in metres.
    sat_xyz: numpy.ndarray
        x, y, z Cartesian geocentric coordinates of the satellite in metres.

    Returns
    -------
    numpy.ndarray
        Array of incidence angle, azimuth angle (measured from north), in radians,
        and distance, measured in meteres

    """
    dx = sat_xyz - xyz
    latlonh = xyz_to_latlonh(xyz)

    if xyz.ndim < 2:
        dx = dx[:, None]
        latlonh = latlonh[:, None]

    normal = numpy.array(
        [
            numpy.cos(latlonh[0, :]) * numpy.cos(latlonh[1, :]),
            numpy.cos(latlonh[0, :]) * numpy.sin(latlonh[1, :]),
            numpy.sin(latlonh[0, :]),
        ]
    )

    ip = normal[0, :] * dx[0, :] + normal[1, :] * dx[1, :] + normal[2, :] * dx[2, :]

    distance = numpy.sqrt(numpy.sum(numpy.power(dx, 2), axis=0))
    incidence = numpy.arccos(ip / distance)
    azimuth_angle = numpy.arctan2(-normal[1, :] * dx[0, :] + normal[0, :] * dx[1, :], ip * -normal[2, :] + dx[2, :])
    return numpy.array([incidence, azimuth_angle, distance]).squeeze()


def latlonh_to_radar(latlonh: numpy.ndarray, metadata: dict) -> tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]:
    """Transform WGS84 coordinates to radar coordinates (azimuth, range).

    Transformation of ellipsoidal geodetic coordinates to radar coordinates
    (azimuth,range), without auxuliary corrections.

    Parameters
    ----------
    latlonh: numpy.ndarray
        Ellipsoidal geodetic coordinates (latitude, longitude, ellipsoidal height)
        in radians and meters.
    metadata: dict
        Image metadata, at least `orbit_time`, `orbit_position`, `orbit_velocity`, `first_azimuth_time`,
        `number_of_rows` and `pulse_repetition_frequency`

    Returns
    -------
    tuple
        Radar coordinates (azimuth, range).

    """
    metadata_validated = validate_geocoding_metadata(metadata, mode="latlonh_to_radar", orbit_required=True)

    xyz = latlonh_to_xyz(latlonh)

    # get time coords
    azimuth_time, range_time, satellite_vector = xyz_to_time(xyz, metadata_validated)

    # convert to pixels
    return time_to_radar(azimuth_time, range_time, metadata_validated)


def latlonh_to_radar_vec(latlonh: numpy.ndarray, metadata: dict) -> tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]:
    """Transform WGS84 coordinates to radar coordinates (azimuth, range) in vectorized format.

    Transformation of ellipsoidal geodetic coordinates to radar coordinates
    (azimuth,range), without auxuliary corrections.

    Parameters
    ----------
    latlonh: numpy.ndarray
        Ellipsoidal geodetic coordinates (latitude, longitude, ellipsoidal height)
        in radians and meters.
    metadata: dict
        Image metadata, at least `orbit_time`, `orbit_position`, `orbit_velocity`, `first_azimuth_time`,
        `number_of_rows` and `pulse_repetition_frequency`

    Returns
    -------
    tuple
        Radar coordinates (azimuth, range).

    """
    metadata_validated = validate_geocoding_metadata(metadata, mode="latlonh_to_radar", orbit_required=True)

    xyz = latlonh_to_xyz(latlonh)

    # get time coords
    azimuth_time, range_time, satellite_vector = xyz_to_time_vec(xyz, metadata_validated)

    # convert to pixels
    return time_to_radar(azimuth_time, range_time, metadata_validated)


def xyz_to_time(
    xyz: numpy.ndarray, metadata: dict, maxiter: int = 10, criter: float = 1e-10
) -> tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]:
    """Transform ECEF coordinates to radar time coordinates (azimuth time, range time).

    Return azimuth time, range time and satellite vector for a target given in geocentric
    Cartesian coordinates (XYZ) by inverse range-Doppler solution.

    Parameters
    ----------
    xyz: numpy.ndarray
        ECEF coordinates (x, y, z) in meters.
    metadata: dict
        Image metadata, at least `orbit_time`, `orbit_position`, `orbit_velocity`, `first_azimuth_time`,
        `number_of_rows` and `pulse_repetition_frequency`
    maxiter: int, default 10
        How many iterations at most to perform in the optimization
    criter: float, default 1e-10
        The limit below which the solution is accepted as converged

    Returns
    -------
    tuple
        (azimuth time, range time (two-way), satellite vector)

    """
    orbit = orbit_fit(
        metadata["orbit_time"],
        metadata["orbit_position"],
        metadata["orbit_velocity"],
    )

    # initial value for azimuth time:
    t0 = seconds_of_day(metadata["first_azimuth_time"])
    t0 += (numpy.round(metadata["number_of_rows"] / 2) - 1) / metadata["pulse_repetition_frequency"]

    if xyz.ndim < 2:
        xyz = xyz[:, None]

    # loop through points
    idx = 0
    t_azimuth = numpy.empty((xyz.shape[1]), dtype=datetime.datetime)
    t_range = numpy.empty((xyz.shape[1],))
    sat_xyz = numpy.empty_like(xyz)

    for coord in xyz.T:
        i = 1
        while i < maxiter:
            dx = coord[0] - chebval(t0 - orbit.time0, orbit.cx)
            dy = coord[1] - chebval(t0 - orbit.time0, orbit.cy)
            dz = coord[2] - chebval(t0 - orbit.time0, orbit.cz)

            vx = chebval(t0 - orbit.time0, orbit.cvx)
            vy = chebval(t0 - orbit.time0, orbit.cvy)
            vz = chebval(t0 - orbit.time0, orbit.cvz)

            ax = chebval(t0 - orbit.time0, orbit.cax)
            ay = chebval(t0 - orbit.time0, orbit.cay)
            az = chebval(t0 - orbit.time0, orbit.caz)

            # inverse range-Doppler solution:
            dt = -(vx * dx + vy * dy + vz * dz) / (ax * dx + ay * dy + az * dz - vx**2 - vy**2 - vz**2)

            t0 = t0 + dt

            if numpy.abs(dt) < criter:
                break
            if i >= maxiter:
                logger.warning("Warning, range-Doppler solution didn't converge!")

            i += 1

        # compute corresponding range time and satellite position
        x_sat = chebval(t0 - orbit.time0, orbit.cx)
        y_sat = chebval(t0 - orbit.time0, orbit.cy)
        z_sat = chebval(t0 - orbit.time0, orbit.cz)
        sat_xyz[:, idx] = numpy.array([x_sat, y_sat, z_sat])

        t_range[idx] = (numpy.sqrt(numpy.sum(numpy.power(coord - sat_xyz[:, idx], 2))) / SPEED_OF_LIGHT) * 2  # 2-way!
        t_azimuth[idx] = metadata["first_azimuth_time"].replace(
            hour=0, minute=0, second=0, microsecond=0
        ) + datetime.timedelta(seconds=t0)

        idx += 1

    return t_azimuth, t_range, sat_xyz.squeeze()


def xyz_to_time_vec(
    xyz: numpy.ndarray, metadata: dict, maxiter: int = 10, criter: float = 1e-10
) -> tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]:
    """Transform ECEF coordinates to radar time coordinates (azimuth time, range time) in a vectorized manner.

    Return azimuth time, range time and satellite vector for a target given in geocentric
    Cartesian coordinates (XYZ) by inverse range-Doppler solution.

    Parameters
    ----------
    xyz: numpy.ndarray
        ECEF coordinates (x, y, z) in meters. Shape (3,) or (3, N).
    metadata: dict
        Image metadata, at least `orbit_time`, `orbit_position`, `orbit_velocity`, `first_azimuth_time`,
        `number_of_rows` and `pulse_repetition_frequency`
    maxiter: int, default 10
        How many iterations at most to perform in the optimization
    criter: float, default 1e-10
        The limit below which the solution is accepted as converged

    Returns
    -------
    tuple
        (azimuth time, range time (two-way), satellite vector)

        - t_azimuth: 1D array of Python datetimes, shape (N,)
        - t_range:   1D float array, shape (N,)
        - sat_xyz:   (3,) for N=1, or (3, N) for N>1 (like original .squeeze() behavior)
    """
    orbit = orbit_fit(
        metadata["orbit_time"],
        metadata["orbit_position"],
        metadata["orbit_velocity"],
    )

    t0 = seconds_of_day(metadata["first_azimuth_time"])
    t0 += (numpy.round(metadata["number_of_rows"] / 2) - 1) / metadata["pulse_repetition_frequency"]

    xyz = numpy.asarray(xyz)
    if xyz.ndim == 1:
        xyz = xyz[:, None]
    elif xyz.ndim != 2 or xyz.shape[0] != 3:
        raise ValueError("xyz must have shape (3,) or (3, N)")

    N = xyz.shape[1]

    t = numpy.full(N, t0, dtype=float)

    # Newton iterations (vectorized over all points)
    for _ in range(maxiter):
        tau = t - orbit.time0

        # Satellite position, velocity, acceleration at each time
        x_sat = chebval(tau, orbit.cx)
        y_sat = chebval(tau, orbit.cy)
        z_sat = chebval(tau, orbit.cz)

        vx = chebval(tau, orbit.cvx)
        vy = chebval(tau, orbit.cvy)
        vz = chebval(tau, orbit.cvz)

        ax = chebval(tau, orbit.cax)
        ay = chebval(tau, orbit.cay)
        az = chebval(tau, orbit.caz)

        # Offsets target - satellite (shape (N,))
        dx = xyz[0, :] - x_sat
        dy = xyz[1, :] - y_sat
        dz = xyz[2, :] - z_sat

        # Inverse range-Doppler solution (vectorized)
        num = -(vx * dx + vy * dy + vz * dz)
        den = ax * dx + ay * dy + az * dz - (vx**2 + vy**2 + vz**2)

        dt = num / den

        t = t + dt

        if numpy.max(numpy.abs(dt)) < criter:
            break
    else:
        logger.warning("Warning, range-Doppler solution didn't fully converge for all points!")

    # Final satellite positions for converged times
    tau = t - orbit.time0
    x_sat = chebval(tau, orbit.cx)
    y_sat = chebval(tau, orbit.cy)
    z_sat = chebval(tau, orbit.cz)

    sat_xyz = numpy.vstack((x_sat, y_sat, z_sat))

    # Two-way range time
    diff = xyz - sat_xyz
    ranges = numpy.sqrt(numpy.sum(diff**2, axis=0))
    t_range = 2.0 * ranges / SPEED_OF_LIGHT

    # Azimuth times as a 1D array of Python datetime objects
    day_start = metadata["first_azimuth_time"].replace(hour=0, minute=0, second=0, microsecond=0)
    t_azimuth = numpy.array(
        [day_start + datetime.timedelta(seconds=float(sec)) for sec in t],
        dtype=object,
    )

    sat_xyz_out = sat_xyz.squeeze()

    return t_azimuth, t_range, sat_xyz_out


def time_to_radar(
    azimuth_time: datetime.datetime, range_time: float, metadata: dict
) -> tuple[numpy.ndarray, numpy.ndarray]:
    """Convert Radar time representation to radar coordinates (azimuth, range).

    Return radar coordinates (azimuth,range) in pixels from time representation.

    Parameters
    ----------
    azimuth_time: datetime.datetime
        Azimuth time in UTC.
    range_time: float
        Range time in seconds, 2-way.
    metadata: dict
        Image metadata, at least `pulse_repetition_frequency`, `first_azimuth_time`, `range_sampling_rate`,
        `first_range_time`

    Returns
    -------
    tuple
        Radar coordinates (azimuth, range) in decimal pixels.

    """
    azimuth_coord = metadata["pulse_repetition_frequency"] * numpy.array(
        [(az - metadata["first_azimuth_time"]).total_seconds() for az in azimuth_time]
    )
    range_coord = metadata["range_sampling_rate"] * (range_time - metadata["first_range_time"])
    return azimuth_coord, range_coord


def radar_to_time(
    azimuth_coords: numpy.ndarray,
    range_coords: numpy.ndarray,
    metadata: dict,
) -> tuple[numpy.ndarray, numpy.ndarray]:
    """Convert radar coordinates (azimuth, range) to Radar time representation.

    Return time representation of radar coordinates from pixels (azimuth, range).

    Parameters
    ----------
    azimuth_coords: numpy.ndarray
        Azimuth pixel coordinates.
    range_coords: numpy.ndarray
        Range pixel coordinates.
    metadata: dict
        Image metadata, at least `pulse_repetition_frequency`, `first_azimuth_time`, `range_sampling_rate`,
        `first_range_time`

    Returns
    -------
    tuple
        Time radar coordinates, azimuth as datetime.datetime and range in 2-way seconds.

    """
    azimuth_time = npdatetime64_to_datetime(metadata["first_azimuth_time"], tz_aware=False) + numpy.array(
        [datetime.timedelta(seconds=az / metadata["pulse_repetition_frequency"]) for az in azimuth_coords]
    )
    range_time = range_coords / metadata["range_sampling_rate"] + metadata["first_range_time"]
    return azimuth_time, range_time


def radar_to_xyz(
    azimuth_coords: numpy.ndarray,
    range_coords: numpy.ndarray,
    elevation: numpy.ndarray,
    metadata: dict,
    return_satellite_vector: bool = False,
) -> numpy.ndarray | tuple[numpy.ndarray, numpy.ndarray]:
    """Convert radar coordinates (azimuth, range) to ECEF XYZ Cartesian coordinates.

    Transformation from radar coordinates to geocentric Cartesian coordinates (XYZ).

    Parameters
    ----------
    azimuth_coords: numpy.ndarray
        Azimuth pixel radar coordinates.
    range_coords: numpy.ndarray
        Range pixel radar coordinates.
    elevation: numpy.ndarray
        Ellipsoidal elevation.
    metadata: dict
        Image metadata, at least `scene_centre_longitude`, `scene_centre_latitude`, `orbit_time`, `orbit_position`,
        `orbit_velocity`, `pulse_repetition_frequency`, `first_azimuth_time`, `range_sampling_rate`, `first_range_time`
    return_satellite_vector: bool, Optional
        Switch to return satellite state vector.

    Returns
    -------
    numpy.ndarray, tuple
        Return geocentric Cartesian coordinates (XYZ) as a single array or as a tuple
        with satellite state vector (Optional).

    """
    e2 = (WGS84.semi_major_metre**2 - WGS84.semi_minor_metre**2) / WGS84.semi_major_metre**2
    # iteration criterions:
    maxiter = 10
    criter = 1e-6

    # initialize iteration for scene center:
    center_lat = metadata["scene_centre_latitude"] / 180 * numpy.pi
    center_lon = metadata["scene_centre_longitude"] / 180 * numpy.pi
    center_N = WGS84.semi_major_metre / numpy.sqrt(1 - e2 * (numpy.sin(center_lat) ** 2))
    center_X = (center_N + numpy.mean(elevation)) * numpy.cos(center_lat) * numpy.cos(center_lon)
    center_Y = (center_N + numpy.mean(elevation)) * numpy.cos(center_lat) * numpy.sin(center_lon)
    center_Z = (center_N + numpy.mean(elevation) - e2 * center_N) * numpy.sin(center_lat)

    t_azimuth, t_range = radar_to_time(azimuth_coords, range_coords, metadata)
    orbit = orbit_fit(
        metadata["orbit_time"],
        metadata["orbit_position"],
        metadata["orbit_velocity"],
    )
    t_azimuth = numpy.array([seconds_of_day(t) for t in t_azimuth])
    sat_xyz, sat_vel = orbit_eval(orbit, t_azimuth)

    if sat_xyz.ndim < 2:
        sat_xyz = sat_xyz[:, None]
        sat_vel = sat_vel[:, None]

    xyz = numpy.empty(sat_xyz.shape)

    # loop through points:
    idx = 0
    for sat, vel, h, t in zip(sat_xyz.T, sat_vel.T, elevation, t_range, strict=False):
        x = center_X.copy()
        y = center_Y.copy()
        z = center_Z.copy()
        i = 1

        while i < maxiter:
            dx = x - sat[0]
            dy = y - sat[1]
            dz = z - sat[2]
            eq = -1 * numpy.array(
                [
                    vel[0] * dx + vel[1] * dy + vel[2] * dz,
                    dx**2 + dy**2 + dz**2 - (SPEED_OF_LIGHT * t / 2) ** 2,
                    (x**2 + y**2) / ((WGS84.semi_major_metre + h) ** 2) + (z / (WGS84.semi_minor_metre + h)) ** 2 - 1,
                ]
            )
            design = numpy.array(
                [
                    vel,
                    2 * numpy.array([dx, dy, dz]),
                    numpy.array(
                        [
                            2 * x / ((WGS84.semi_major_metre + h) ** 2),
                            2 * y / ((WGS84.semi_major_metre + h) ** 2),
                            2 * z / ((WGS84.semi_minor_metre + h) ** 2),
                        ]
                    ),
                ]
            )

            sol = numpy.linalg.solve(design, eq)
            x += sol[0]
            y += sol[1]
            z += sol[2]

            if numpy.any(numpy.abs(sol) > criter):
                if i < maxiter:
                    i += 1
                    continue
                else:
                    logger.warning(f"radar_to_xyz did not converge after {maxiter} iterations!")
            else:
                break
        xyz[:, idx] = numpy.array([x, y, z])
        idx += 1

    if return_satellite_vector:
        return xyz.squeeze(), sat_xyz.squeeze()

    return xyz.squeeze()


def radar_to_latlonh(
    azimuth_coords: numpy.ndarray,
    range_coords: numpy.ndarray,
    elevation: numpy.ndarray,
    metadata: dict,
) -> numpy.ndarray:
    """Convert radar coordinates (azimuth, range) to WGS84 Lat/lon/h coordinates.

    Parameters
    ----------
    azimuth_coords: numpy.ndarray
        Azimuth pixel radar coordinates.
    range_coords: numpy.ndarray
        Range pixel radar coordinates.
    elevation: numpy.ndarray
        Ellipsoidal elevation.
    metadata: dict
        Image metadata, at least `scene_centre_longitude`, `scene_centre_latitude`,
        `pulse_repetition_frequency`, `first_azimuth_time`, `range_sampling_rate`, `first_range_time`, and
        (`orbit_txyz` OR `orbit_time`, `orbit_position`), optionally `orbit_velocity`

    Returns
    -------
    numpy.ndarray, tuple
        Return latitude/longitude/height coordinates

    """
    metadata_validated = validate_geocoding_metadata(metadata, mode="radar_to_latlonh", orbit_required=True)
    xyz = radar_to_xyz(azimuth_coords, range_coords, elevation, metadata_validated, return_satellite_vector=False)
    latlonh = xyz_to_latlonh(xyz.T)
    return latlonh


def validate_geocoding_metadata(
    metadata: dict, mode: Literal["latlonh_to_radar", "radar_to_latlonh"], orbit_required: bool
) -> dict:
    """Validate that all fields necessary for the geocoding are present, and regulate the orbit metadata.

    Since DORIS v5 outputs `orbit_txyz` instead of orbit_time, orbit_position, and orbit_velocity, this function will
    assign the correct columns of `orbit_txyz` to the correct expected fields.

    Parameters
    ----------
    metadata: dict
        Metadata readout as per `sarxarray.read_metadata` (>=v1.2.2)
    mode: Literal["latlonh_to_radar", "radar_to_latlonh"]
        Which keys to check that exist
    orbit_required: bool
        Whether or not to construct `orbit_time`, `orbit_position` and `orbit_velocity` if they don't exist

    Returns
    -------
    dict
        The validated and corrected metadata dictionary

    """
    assert mode in VALIDATION_NONORBIT_KEYS.keys(), f"Unknown mode {mode}, known are {VALIDATION_NONORBIT_KEYS.keys()}"
    for key in VALIDATION_NONORBIT_KEYS[mode]:
        assert key in metadata.keys(), f"Key {key} is missing from metadata, cannot proceed!"

    if orbit_required:
        # Orbit info from DORIS5 is read in in orbit_txyz instead of orbit_time and orbit_position -->
        # needs to be split up
        if "orbit_time" not in metadata.keys():
            if "orbit_txyz" in metadata.keys():
                metadata["orbit_time"] = metadata["orbit_txyz"][:, 0]
            else:
                raise ValueError("Key orbit_time is missing and cannot be reconstructed from orbit_txyz!")

        if "orbit_position" not in metadata.keys():
            if "orbit_txyz" in metadata.keys():
                metadata["orbit_position"] = metadata["orbit_txyz"][:, 1:]
            else:
                raise ValueError("Key orbit_position is missing and cannot be reconstructed from orbit_txyz!")

        if "orbit_velocity" not in metadata.keys():
            metadata["orbit_velocity"] = None  # this one can be reconstructed when necessary

    return metadata
