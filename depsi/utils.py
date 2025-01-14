import os
from datetime import datetime, timezone

import geopandas
import numpy as np
import pytz
import xarray as xr

UTC = timezone.utc


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


def _npdatetime64_to_datetime(date: np.datetime64) -> datetime:
    """Convert a numpy datetime64 object to a python datetime object.

    Parses the np.datetime64 object into a datetime object.

    Parameters
    ----------
    date : np.datetime64
      the date to be converted

    Returns
    -------
    datetime.datetime
      The same date converted to a datetime object
    """
    timestamp = (date - np.datetime64("1970-01-01T00:00:00")) / np.timedelta64(1, "s")
    return datetime.fromtimestamp(timestamp, UTC)


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
      tuple of two np.ndarrays, the first containing the latitudes, the second the longitudes of the bounding box.

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
    bounding_box = shape.envelope.boundary[0].xy

    return bounding_box


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
        fmt_dates = [_npdatetime64_to_datetime(date) for date in slcs["time"].values]
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

        fmt_dates = np.array([_npdatetime64_to_datetime(date) for date in slcs["time"].values])
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
        slcs = slcs.where(space_mask.compute(), drop=True)

    return slcs
