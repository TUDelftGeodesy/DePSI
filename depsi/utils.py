import math
import os
from typing import Literal

import asf_search as asf
import dask.array as da
import pyproj

try:
    from datetime import UTC, datetime
except ImportError:  # UTC can only be imported from Python 3.11 onwards
    import warnings
    from datetime import datetime, timezone

    UTC = timezone.utc
    warnings.warn(
        """
    DePSI uses datetime.UTC, which is only supported from Python 3.11 onwards.
    For older Python versions, datetime.timezone.utc is used.
    This might be deprecated in newer DePSI versions.
    """,
        DeprecationWarning,
        stacklevel=1,  # necessary to start the call stack here.
    )

import geopandas
import numpy as np
import pytz
import xarray as xr


def wrap_phase(phs_abs):
    """Wrap the absolute phase to the range [-pi, pi).

    Parameters
    ----------
    phs_abs : array_like or float
        The absolute phase.

    Returns
    -------
    ndarray or float
        The wrapped phase in the range [-pi, pi).
    """
    phs_wrapped = np.remainder(phs_abs + np.pi, 2 * np.pi) - np.pi

    return phs_wrapped


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
    source: list | tuple | np.ndarray,
    target: list | tuple | np.ndarray,
    mode: Literal["euclidean", "geographic"] = "euclidean",
):
    """Calculate the distance between two points.

    Parameters
    ----------
    source: list | tuple | np.ndarray
        The source point, formatted as (x, y) / (lon, lat)
    target: list | tuple | np.ndarray
        The target point, formatted as (x, y) / (lon, lat)
    mode: Literal["euclidean", "geographic"], default "euclidean"
        Whether the source and target points are given in (x, y) (units meters) or (lon, lat) (units degrees)

    Returns
    -------
        The distance between the two points in meters.
    """
    if mode == "euclidean":
        return math.dist(source, target)
    elif mode == "geographic":
        # this is the Haversine formula
        lat1 = source[1]
        lat2 = target[1]
        dphi = np.radians(lat1 - lat2)
        dlambda = np.radians(source[0] - target[0])
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


def npdatetime64_to_datetime(date: np.datetime64, tz_aware: bool = True) -> datetime:
    """Convert a numpy datetime64 object to a python datetime object.

    Parses the np.datetime64 object into a datetime object.

    Parameters
    ----------
    date : np.datetime64
      the date to be converted
    tz_aware: bool, default True
      whether the returned datetime object should be timezone-aware or not

    Returns
    -------
    datetime.datetime
      The same date converted to a datetime object
    """
    timestamp = (date - np.datetime64("1970-01-01T00:00:00")) / np.timedelta64(1, "s")
    dt_obj = datetime.fromtimestamp(timestamp, UTC)
    if not tz_aware:
        dt_obj = datetime.strptime(dt_obj.strftime("%Y%m%d:%H%M%S"), "%Y%m%d:%H%M%S")
    return dt_obj


def _get_aoi_shapefile_bounding_box(aoi_filename: str) -> tuple:
    """Read a .shp shapefile and return the bounding box.

    The shapefile is read and the area of interest is retrieved. The bounding box is then computed, and a tuple of the
    coordinates of the bounding box is returned.

    Parameters
    ----------
    aoi_filename: str | None
      full path to the AoI shapefile, expects .shp format.

    Returns
    -------
    tuple
      tuple of two lists, the first containing the longitude extent, the second the latitude extent of the
      bounding box.

    Raises
    ------
    AssertionError
        Raised when:
        - the aoi_filename does not exist
        - the aoi_filename does not end in .shp

    ValueError
        Raised when:
        - the provided shapefile contains zero polygons, or more than one polygon
        - the provided shapefile contains an invalid polygon
    """
    assert os.path.exists(aoi_filename), f"The file {aoi_filename} does not exist!"
    assert aoi_filename.split(".")[-1] == "shp", f"The provided file {aoi_filename} is not of .shp type!"

    # open the file, and iterate through the geometry
    shape = geopandas.read_file(aoi_filename)
    # calculate the coordinates of the bounding box of the provided AoI
    bounding_box = shape.total_bounds
    # format as longitude extent (in x), latitude extent (in y)
    bounding_box_formatted = ([bounding_box[0], bounding_box[2]], [bounding_box[1], bounding_box[3]])

    return bounding_box_formatted


def crop_slc_spacetime(
    slcs: xr.Dataset,
    aoi_filename: str | None = None,
    start_date: datetime | str | None = None,
    end_date: datetime | str | int | None = None,
) -> xr.Dataset:
    """Crop an SLC stack in both space and time.

    To crop in space an AoI shapefile is processed, and the SLC stack is cropped to the bounding box.
    To crop in time, a start date is provided, and two options for the end date are available:
    - datetime | str - this directly provides the end date
    - int - this provides a number of SLCs intended to be in the crop. The end date is set automatically.

    If the aoi_filename is not provided and left to None, only a crop in time is performed.
    If the start_date and end_date are not provided and left to None, only a crop in space is performed.

    Parameters
    ----------
    slcs : xr.Dataset
      the SLC stack to be cropped. Requires at least the following coordinates or variables:
      In case of a crop in time:
      - time -> the dates of the images
      In case of a crop in space:
      - lat -> the latitude of the pixels
      - lon -> the longitude of the pixels
    aoi_filename: str | None
      full path to the AoI shapefile, expects .shp format. Set to None if no crop in space is requested.
    start_date : datetime | str | None
      the start date of the crop, in one of three formats:
      - datetime object
      - str object, formatted as YYYYMMDD
      - None, no cropping in time requested
    end_date : datetime | str | int | None
      the end date of the crop, in one of four formats:
      - datetime object
      - str object, formatted as YYYYMMDD
      - int object, which is interpreted as the number of images intended in the crop (including the start date). If
        more images are requested than exist since the start date, all images from start_date until the last image
        are provided.
      - None, no cropping in time requested

    Returns
    -------
    xr.Dataset
      The cropped dataset

    Raises
    ------
    AssertionError
        Raised when:
        - a start_date or end_date is provided in string format, but not in YYYYMMDD format
        - a start_date is provided, but the end_date is set to None
        - an end_date is provided, but the start_date is set to None
        - the aoi_filename does not exist
        - the aoi_filename does not end in .shp
        - one of the required coordinates or variables is not available in slcs

    ValueError
        Raised when:
        - start_date is not of type datetime | str | None
        - end_date is not of type datetime | str | int | None
    """
    # Check the input

    if aoi_filename is not None:
        assert os.path.exists(aoi_filename), f"The file {aoi_filename} does not exist!"
        assert aoi_filename.split(".")[-1] == "shp", f"The provided file {aoi_filename} is not of .shp type!"
        for axis in ["lat", "lon"]:
            assert axis in slcs.keys(), f"Expected axis {axis} in SLCs but it is not present!"

    # convert the input to a timezone-aware datetime object
    if isinstance(start_date, str):
        assert len(start_date) == 8, f"Unknown start_date format {start_date}, expected YYYYMMDD!"
        format_start_date = datetime(
            eval(start_date[:4]), eval(start_date[4:6].lstrip("0")), eval(start_date[6:].lstrip("0")), tzinfo=pytz.UTC
        )
    elif isinstance(start_date, datetime):
        format_start_date = datetime(start_date.year, start_date.month, start_date.day, tzinfo=pytz.UTC)
    elif start_date is None:
        assert end_date is None, f"Start date is None while end date is {end_date} (not None!)"
        format_start_date = None
    else:
        raise ValueError(f'Expected start_date of type "str" | "datetime" | None, got {type(start_date)}!')

    if isinstance(end_date, str):
        assert len(end_date) == 8, f"Unknown end_date format {end_date}, expected YYYYMMDD!"
        format_end_date = datetime(
            eval(end_date[:4]), eval(end_date[4:6].lstrip("0")), eval(end_date[6:].lstrip("0")), tzinfo=pytz.UTC
        )
    elif isinstance(end_date, datetime):
        format_end_date = datetime(end_date.year, end_date.month, end_date.day, tzinfo=pytz.UTC)
    elif isinstance(end_date, int):
        fmt_dates = [npdatetime64_to_datetime(date) for date in slcs["time"].values]
        valid_dates = [date for date in fmt_dates if date >= format_start_date]
        end_idx = fmt_dates.index(valid_dates[0]) + end_date - 1
        end_idx = min(end_idx, len(fmt_dates) - 1)
        format_end_date = fmt_dates[end_idx]
    elif end_date is None:
        assert start_date is None, f"Start date is None while end date is {end_date} (not None!)"
        format_end_date = None
    else:
        raise ValueError(f'Expected end_date of type "str" | "datetime" | "int" | None, got {type(end_date)}!')

    # TIME CROP
    if format_start_date is not None and format_end_date is not None:
        # first the last assertion
        assert "time" in slcs.keys(), "Expected axis 'time' in SLCs but it is not present!"

        fmt_dates = np.array([npdatetime64_to_datetime(date) for date in slcs["time"].values])
        time_mask = (format_start_date <= fmt_dates) & (fmt_dates <= format_end_date)
        slcs = slcs.sel(time=slcs["time"].values[time_mask])

    # SPACE CROP
    if aoi_filename is not None:
        bounding_box = _get_aoi_shapefile_bounding_box(aoi_filename)
        space_mask = (
            (slcs["lat"] >= min(bounding_box[1]))
            & (slcs["lat"] <= max(bounding_box[1]))
            & (slcs["lon"] >= min(bounding_box[0]))
            & (slcs["lon"] <= max(bounding_box[0]))
        )

        comp_space_mask = space_mask.compute()
        az_sum = comp_space_mask.sum(dim="azimuth")
        rg_sum = comp_space_mask.sum(dim="range")

        # first and last non zero
        min_range, max_range = (
            az_sum.where(az_sum > 0, drop=True)["range"].min().values,
            az_sum.where(az_sum > 0, drop=True)["range"].max().values,
        )
        min_azimuth, max_azimuth = (
            rg_sum.where(rg_sum > 0, drop=True)["azimuth"].min().values,
            rg_sum.where(rg_sum > 0, drop=True)["azimuth"].max().values,
        )
        # data at original locations not nan
        slcs = slcs.sel(azimuth=range(min_azimuth, max_azimuth), range=range(min_range, max_range))

    return slcs


def project_stm_coordinates(stm: xr.Dataset, projection: str = "RD") -> xr.Dataset:
    """Project the latitude and longitude of a space-time matrix to another reference frame.

    The latitude and longitude layers are transformed into the desired projection, default Rijksdriehoek or RD.

    Parameters
    ----------
    stm: xr.Dataset
      Space-time matrix with the layers `lat` (latitude) and `lon` (longitude in WGS84 (EPSG:4326), and coordinate
      `space`
    projection: str, optional
      Projection to which the latitude and longitude coordinates should be transformed. "RD" defaults to "EPSG:28992".
      Default "RD"

    Returns
    -------
    xr.Dataset
      Space-time matrix with the added layers `projection_x` and `projection_y`, where projection is the requested
      parameter `projection` in lower case.

    Raises
    ------
    AssertionError
      When layers "lon" or "lat" do not exist in `stm`.
    """
    assert "lon" in stm.keys(), "Expected a space-time matrix with longitude layer named lon but it is not there!"
    assert "lat" in stm.keys(), "Expected a space-time matrix with latitude layer named lat but it is not there!"
    if projection == "RD":
        projection_formatted = "EPSG:28992"
    elif projection[:5] == "EPSG:":
        projection_formatted = projection
    else:
        raise ValueError(f"Invalid projection provided! Expected 'RD' or 'EPSG:###' but got {projection}!")

    wgs84 = pyproj.Transformer.from_crs("EPSG:4326", projection_formatted, always_xy=True).transform
    # Convert Lat and Lon to coordinates
    proj_x, proj_y = wgs84(stm["lon"], stm["lat"])

    # Add coordinates to the dataset
    stm = stm.assign({f"{projection.lower()}_x": (["space"], proj_x)})
    stm = stm.assign({f"{projection.lower()}_y": (["space"], proj_y)})

    return stm


def add_stm_time_deltas(stm: xr.Dataset) -> xr.Dataset:
    """Add the time differences since the first image to a space-time matrix.

    Parameters
    ----------
    stm: xr.Dataset
      the space-time matrix with an axis "time"

    Returns
    -------
    xr.Dataset
      the space-time matrix with two new variables:
      - `days_since_first_img`, the number of days since the first epoch in the STM
      - `years_since_first_img`, the number of years since the first epoch in the STM, assuming 365.2425 days per year

    """
    assert "time" in stm.keys(), "Expected STM to have time axis but it's not there!"
    # Add extra time coordinate variables for time intervals since first image
    days = np.array(
        [
            (npdatetime64_to_datetime(date) - npdatetime64_to_datetime(stm["time"].values[0])).days
            for date in stm["time"].values
        ]
    )
    stm = stm.assign({"days_since_first_img": (["time"], days)})
    stm = stm.assign({"years_since_first_img": (["time"], days / 365.2425)})

    return stm


def stm_compute_single_time_differences(
    stm: xr.Dataset, single_difference_mother: str | datetime = "auto"
) -> xr.Dataset:
    """Compute the single differences of an STM in time with respect to a given mother image.

    This computes the single difference complex value, phase, unnormalized amplitude, and h2ph values with respect
    to the provided single difference mother. The mother image is the first image acquired on or after the provided
    date (if a datetime object or str object is provided), or the mother image of the input dataset (if 'auto' mode
    is selected).

    Parameters
    ----------
    stm: xr.Dataset
      the space-time matrix with an axis "time" and "space", and variables `h2ph` and `complex`
    single_difference_mother: datetime | str
      the date to be used as the mother image for the single difference computations, in one of three formats:
      - 'auto' : will detect the mother image in the input SLC dataset, and use that epoch.
      - datetime object
      - str object, formatted as YYYYMMDD

    Returns
    -------
    xr.Dataset
      the space-time matrix with four new variables:
        - sd_h2ph (space, time): single difference height to phase conversion with respect to single_difference_mother
        - sd_complex (space, time): single difference complex phasor with respect to single_difference_mother
        - sd_amplitude_unnormalized (space, time): single difference complex phasor amplitude to
            single_difference_mother, not normalized
        - sd_phase (space, time): single difference phase with respect to single_difference_mother

    Raises
    ------
    ValueError
      Raised when:
        - single_difference_mother is of an unsupported format
        - the date provided to single_difference_mother is not in the input stack date range
    """
    # Identify the mother image
    if isinstance(single_difference_mother, datetime):
        format_mother_date = datetime(
            single_difference_mother.year, single_difference_mother.month, single_difference_mother.day, tzinfo=pytz.UTC
        )
        mother_index = [
            idx for idx, date in enumerate(stm["time"].values) if format_mother_date <= npdatetime64_to_datetime(date)
        ]  # select all images beyond the mother date
    elif isinstance(single_difference_mother, str):
        if single_difference_mother == "auto":
            mother_index = np.where(abs(stm["h2ph"]).sum(axis=0).values == 0)[0]
        elif len(single_difference_mother) == 8:
            format_mother_date = datetime(
                eval(single_difference_mother[:4]),
                eval(single_difference_mother[4:6].lstrip("0")),
                eval(single_difference_mother[6:].lstrip("0")),
                tzinfo=pytz.UTC,
            )
            mother_index = [
                idx
                for idx, date in enumerate(stm["time"].values)
                if format_mother_date <= npdatetime64_to_datetime(date)
            ]  # select all images beyond the mother date
        else:
            raise ValueError(f'Cannot parse {single_difference_mother}, not of type "auto" or "YYYYMMDD"!')
    else:
        raise ValueError(f"Unknown format {type(single_difference_mother)} for single_difference_mother!")
    if len(mother_index) == 0:
        raise ValueError(
            f"Cannot find provided mother date {single_difference_mother}, "
            "please provide a date that is within the range of the stack! Possible dates: "
            f"{stm.time.values[0]}--{stm.time.values[-1]}"
        )
    sd_mother_index = mother_index[0]  # 0 in case more than 1 image is detected
    # In that case we take the first image that was detected, as this is expected
    sd_mother = npdatetime64_to_datetime(stm["time"].values[sd_mother_index])

    # Format the single difference mother, and save it to the STM
    stm.attrs["ps_sd_mother"] = sd_mother.strftime("%Y%m%d")

    # calculate the h2ph single difference (= daughter - mother)
    sd_h2ph = stm["h2ph"] - stm["h2ph"][:, sd_mother_index]
    stm = stm.assign({"sd_h2ph": (["space", "time"], sd_h2ph.data)})

    # calculate the complex single difference, the amplitude, and the phase
    mother_comp = stm["complex"][:, sd_mother_index].conj()
    sd_complex_transposed = stm["complex"].transpose() * mother_comp
    sd_complex = sd_complex_transposed.transpose()
    sd_phase = da.angle(sd_complex)
    sd_amplitude_unnormalized = da.abs(sd_complex)
    stm = stm.assign({"sd_complex": (["space", "time"], sd_complex.data)})
    stm = stm.assign({"sd_amplitude_unnormalized": (["space", "time"], sd_amplitude_unnormalized.data)})
    stm = stm.assign({"sd_phase": (["space", "time"], sd_phase.data)})

    return stm


def identify_s1_orbits_in_aoi(lon: list | np.ndarray, lat: list | np.ndarray) -> tuple[list[str], dict]:
    """Identify the Sentinel-1 orbit numbers and directions crossing a AoI.

    Parameters
    ----------
    lon: list | np.ndarray
        List of all the longitudes of all the points of interest in the AoI
    lat: list | np.ndarray
        List of all the latitudes of all the points of interest in the AoI

    Returns
    -------
    list
        The orbits overlapping with the AoI
    dict
        The footprints of the overlapping SLCs per track
    """
    bbox = [[np.min(lon), np.max(lon)], [np.min(lat), np.max(lat)]]
    wkt = (
        f"POLYGON(("
        f"{bbox[0][0]} {bbox[1][0]}, "
        f"{bbox[0][1]} {bbox[1][0]}, "
        f"{bbox[0][1]} {bbox[1][1]}, "
        f"{bbox[0][0]} {bbox[1][1]}, "
        f"{bbox[0][0]} {bbox[1][0]}))"
    )
    slcs = None
    counter = 0
    while slcs is None:
        try:
            slcs = asf.geo_search(
                intersectsWith=wkt,
                platform=asf.PLATFORM.SENTINEL1,
                beamMode="IW",
                processingLevel="SLC",
                start="one month ago",
                end="now",
            )
        except (asf.exceptions.ASFSearch5xxError, asf.exceptions.ASFSearchError, TimeoutError):
            counter += 1
            print(f"ASF encountered an internal error. Retrying... (#{counter})")

    orbits = [
        f"s1_{slc.properties['flightDirection'].lower().replace('e', '')[:3]}_t{slc.properties['pathNumber']:0>3d}"
        for slc in slcs
    ]
    filtered_orbits = list(sorted(list(set(orbits))))

    extents = [slc.geojson()["geometry"]["coordinates"][0] for slc in slcs]
    footprints = {}
    for orbit in filtered_orbits:
        footprints[orbit] = []

    for extent in range(len(extents)):
        footprints[orbits[extent]].append(extents[extent])

    return filtered_orbits, footprints
