import collections
import datetime
import logging
from typing import NamedTuple

import numpy
from numpy.polynomial.chebyshev import chebval

logger = logging.getLogger(__name__)


class Ellipsoid(NamedTuple):
    """Base class for ellipsoidal constants."""

    a: float
    b: float
    f: float


SPEED_OF_LIGHT = 299792458.0
MJD_EPOCH = datetime.datetime(2000, 1, 1)
WGS84 = Ellipsoid(6378137.0, 6356752.3141, 0.003352810681182)

OrbitFit = collections.namedtuple("OrbitFit", "time0, cx, cy, cz, cvx, cvy, cvz, cax, cay, caz")


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
        latitude, longitude, ellipsoidal height coordinates as a numpy array in [radians/m].

    Returns
    -------
    xyz: numpy.ndarray
        Cartesian geocentric coordinates in [m].
    """
    e2 = 2 * WGS84.f - WGS84.f**2

    if latlonh.ndim < 2:
        latlonh = latlonh[:, None]

    n = numpy.divide(WGS84.a, numpy.sqrt(1 - e2 * numpy.power(numpy.sin(latlonh[0, :]), 2)))
    x = (n + latlonh[2, :]) * numpy.cos(latlonh[0, :]) * numpy.cos(latlonh[1, :])
    y = (n + latlonh[2, :]) * numpy.cos(latlonh[0, :]) * numpy.sin(latlonh[1, :])
    z = (n - n * e2 + latlonh[2, :]) * numpy.sin(latlonh[0, :])

    return numpy.array([x, y, z]).squeeze()


def xyz_to_latlonh(xyz: numpy.ndarray) -> numpy.ndarray:
    """Transform ECEF Coordinates to GRS80 coordinates.

    Transformation of geocentric Cartesian (ECEF) coordinates to GRS80 ellipsoidal coordinates
    using iteration method.

    Parameters
    ----------
    xyz: numpy.ndarray
        x, y, z Cartesian geocentric coordinates in metres.

    Returns
    -------
    latlonh: numpy.ndarray
        Latitude, longitude, ellipsoidal height coordinates as radians and metres.
    """
    e2 = 2 * WGS84.f - WGS84.f**2

    if xyz.ndim < 2:
        xyz = xyz[:, None]

    r = numpy.sqrt(numpy.power(xyz[0, :], 2) + numpy.power(xyz[1, :], 2))

    i = 1
    maxiter = 5
    n_p = xyz[2, :]
    while i < maxiter:
        phi = numpy.arctan((xyz[2, :] + e2 * n_p) / r)
        n = WGS84.a / numpy.sqrt(1 - e2 * numpy.power(numpy.sin(phi), 2))
        n_p = n * numpy.sin(phi)
        i += 1

    return numpy.array([phi, numpy.arctan2(xyz[1, :], xyz[0, :]), r / numpy.cos(phi) - n]).squeeze()


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
        xyz_vel = numpy.zeros(xyz_pos.shape)
        xyz_vel[:, 0] = numpy.polynomial.chebyshev.chebder(cx)
        xyz_vel[:, 1] = numpy.polynomial.chebyshev.chebder(cy)
        xyz_vel[:, 2] = numpy.polynomial.chebyshev.chebder(cz)

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
        Image metadata

    Returns
    -------
    tuple
        Radar coordinates (azimuth, range).

    """
    xyz = latlonh_to_xyz(latlonh)

    # get time coords
    azimuth_time, range_time, satellite_vector = xyz_to_time(xyz, metadata)

    # convert to pixels
    return time_to_radar(azimuth_time, range_time, metadata)


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
        Image metadata

    Returns
    -------
    tuple
        Radar coordinates (azimuth, range).

    """
    xyz = latlonh_to_xyz(latlonh)

    # get time coords
    azimuth_time, range_time, satellite_vector = xyz_to_time_vec(xyz, metadata)

    # convert to pixels
    return time_to_radar(azimuth_time, range_time, metadata)


def xyz_to_time(xyz: numpy.ndarray, metadata: dict) -> tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]:
    """Transform ECEF coordinates to radar time coordinates (azimuth time, range time).

    Return azimuth time, range time and satellite vector for a target given in geocentric
    Cartesian coordinates (XYZ) by inverse range-Doppler solution.

    Parameters
    ----------
    xyz: numpy.ndarray
        ECEF coordinates (x, y, z) in meters.
    metadata: dict
        Image metadata

    Returns
    -------
    tuple
        (azimuth time, range time (two-way), satellite vector)

    """
    maxiter = 10
    criter = 1e-10

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


def xyz_to_time_vec(xyz: numpy.ndarray, metadata: dict) -> tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]:
    """Transform ECEF coordinates to radar time coordinates (azimuth time, range time) in a vectorized manner.

    Return azimuth time, range time and satellite vector for a target given in geocentric
    Cartesian coordinates (XYZ) by inverse range-Doppler solution.

    Parameters
    ----------
    xyz: numpy.ndarray
        ECEF coordinates (x, y, z) in meters. Shape (3,) or (3, N).
    metadata: dict
        Image metadata

    Returns
    -------
    tuple
        (azimuth time, range time (two-way), satellite vector)

        - t_azimuth: 1D array of Python datetimes, shape (N,)
        - t_range:   1D float array, shape (N,)
        - sat_xyz:   (3,) for N=1, or (3, N) for N>1 (like original .squeeze() behavior)
    """
    maxiter = 10
    criter = 1e-10

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
        Image metadata

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
        Image metadata

    Returns
    -------
    tuple
        Time radar coordinates, azimuth as datetime.datetime and range in 2-way seconds.

    """
    azimuth_time = metadata["first_azimuth_time"] + numpy.array(
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
        Image metadata
    return_satellite_vector: bool, Optional
        Switch to return satellite state vector.

    Returns
    -------
    numpy.ndarray, tuple
        Return geocentric Cartesian coordinates (XYZ) as a single array or as a tuple
        with satellite state vector (Optional).

    """
    e2 = (WGS84.a**2 - WGS84.b**2) / WGS84.a**2
    # iteration criterions:
    maxiter = 10
    criter = 1e-6

    # initialize iteration for scene center:
    center_lat = metadata["scene_centre_latitude"] / 180 * numpy.pi
    center_lon = metadata["scene_centre_longitude"] / 180 * numpy.pi
    center_N = WGS84.a / numpy.sqrt(1 - e2 * (numpy.sin(center_lat) ** 2))
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
                    (x**2 + y**2) / ((WGS84.a + h) ** 2) + (z / (WGS84.b + h)) ** 2 - 1,
                ]
            )
            design = numpy.array(
                [
                    vel,
                    2 * numpy.array([dx, dy, dz]),
                    numpy.array(
                        [
                            2 * x / ((WGS84.a + h) ** 2),
                            2 * y / ((WGS84.a + h) ** 2),
                            2 * z / ((WGS84.b + h) ** 2),
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
                    logger.warning("radar_to_xyz did not converge after 10 iterations!")
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
        Image metadata

    Returns
    -------
    numpy.ndarray, tuple
        Return latitude/longitude/height coordinates

    """
    xyz = radar_to_xyz(azimuth_coords, range_coords, elevation, metadata, return_satellite_vector=False)
    latlonh = xyz_to_latlonh(xyz)
    return latlonh
