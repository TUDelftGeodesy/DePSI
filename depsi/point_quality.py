from typing import Literal

import dask.array as da
import numpy as np
import ruptures as rpt
import xarray as xr

# The partitioning requires a jump size when using pelt mode. This should always be 5.
# TODO: add documentation as to why this should be 5
PELT_JUMP = 5
# The conversion from NAD and NMAD to standard deviation is done using an empirical approximation. Four approximations
# have been modeled: NAD - mean (the mean of the simulations for NAD), NAD - 2sigma (the mean of the simulations
# plus 2 sigma for NAD), and the same two for NMAD. The function is of the form value = A + B*x + C*x^2 + D*x^3 .
# In this dictionary the values are stored as [A, B, C, D]
NAD_NMAD_TO_SIGMA_CONVERSION = {
    "nad": {
        "mean": [-7.65752941e-03, 1.33360757e00, -3.18428074e00, 9.35392564e00],
        "2sigma": [-0.03222335, 2.02221987, -5.76342934, 14.47118093],
    },
    "nmad": {
        "mean": [-1.44869469e-02, 2.00028682e00, -5.23271341e00, 2.11111801e01],
        "2sigma": [0.01907808, 1.2852969, 1.90052824, 11.60677721],
    },
}


def _estimate_breakpoints(
    amplitude_array: xr.DataArray,
    db_segmentation: bool = False,
    search_method: Literal["pelt", "binseg"] = "pelt",
    cost_model: str = "l2",
    size: int = 27,
) -> [xr.DataArray, xr.DataArray]:
    """Identify breakpoints in the amplitude timeseries of a DataArray of amplitude information.

    Identifies breakpoints in the amplitude timeseries of a set of points, based on a given search method, cost model
    and minimum partition size between the breakpoints. The search can be performed on the normal amplitude data or on
    a dB scale.

    Parameters
    ----------
    amplitude_array: xr.DataArray
      The data array with the amplitude time series of the points.
    db_segmentation: bool, optional
      Toggle to turn partitioning on the dB scale on (True) or off (False). Defaults to False.
    search_method: Literal["pelt", "binseg"], optional
      Which search method to use. Default (and recommended) "pelt"
    cost_model: str, optional
      Which cost model to use. Default (and recommended) "l2"
    size: int, optional
      minimum size of the partitions. Default 27

    Returns
    -------
    xr.DataArray
      boolean data array where True indicates a breakpoint at that epoch
    xr.DataArray
      integer data array where each partition has been assigned a unique identifier.
    """
    match db_segmentation:
        case False:
            amplitude_ts = amplitude_array
        case True:
            amplitude_ts = 10 * da.log10(amplitude_array)
        case _:
            raise ValueError(f"db_segmentation should be False or True but is {db_segmentation}!")

    match search_method:
        case "pelt":
            breakpoints = xr.map_blocks(
                _pelt_block,
                amplitude_ts,
                args=(cost_model, size),
                template=amplitude_ts,
            )
        case "binseg":
            breakpoints = xr.map_blocks(
                _binseg_block,
                amplitude_ts,
                args=(cost_model, size),
                template=amplitude_ts,
            )
        case _:
            raise ValueError(f"search_method should be 'pelt' or 'binseg' but is {search_method}!")

    # Compute the breakpoints (necessary for the identifiers since the chunking breaks)
    computed_breakpoints = breakpoints.values
    # to compute statistics per partition we need to assign each partition a unique identifier.
    # This consists of two parts:
    # 1. a cumulative sum over the entire breakpoint True/False array. This will increase the identifier by 1 for each
    # breakpoint encountered.
    # 2. a point index to increase the identifier by 1 at the start of each new point. Otherwise the last partition of
    # point 1 and the first partition of point 2 will have the same identifier.
    breakpoints.data = computed_breakpoints
    breakpoints_idx_p1 = np.cumsum(breakpoints.data).reshape(breakpoints.shape)
    point_idx = np.arange(0, breakpoints.shape[0]).reshape((breakpoints.shape[0], 1))
    point_idx = np.hstack([point_idx for _ in range(breakpoints.shape[1])])
    breakpoints_idx = breakpoints_idx_p1 + point_idx

    return breakpoints, breakpoints_idx


def _pelt_block(amplitude_ts: xr.Dataset, cost_model: str, size: int) -> xr.Dataset:
    """Compute the pelt cost function for breakpoints in chunks.

    Parameters
    ----------
    amplitude_ts: xr.Dataset
      dataset containing the amplitude values of a chunk of points.
    cost_model: str
      Which cost model to use. Default (and recommended) "l2"
    size: int
      minimum size of the partitions. Default 27

    Returns
    -------
    xr.Dataset
      Boolean dataset of the same size of the amplitude array, where True indicates a detected breakpoint
    """
    groups = amplitude_ts.groupby("space")
    stmat_out = groups.map(_pelt_single_point, cost_model=cost_model, size=size)
    return stmat_out


def _pelt_single_point(amplitude_ts: xr.Dataset, cost_model: str, size: int) -> xr.Dataset:
    """Compute the pelt cost function for breakpoints for a single point.

    Parameters
    ----------
    amplitude_ts: xr.Dataset
      dataset containing the amplitude values of one point.
    cost_model: str
      Which cost model to use. Default (and recommended) "l2"
    size: int
      minimum size of the partitions. Default 27

    Returns
    -------
    xr.Dataset
      Boolean dataset of the same size of the amplitude array, where True indicates a detected breakpoint
    """
    amplitude_computed = amplitude_ts.compute().data
    penalty = np.var(amplitude_computed) * np.log(amplitude_ts.time.shape[0])
    algo = rpt.Pelt(model=cost_model, min_size=size, jump=PELT_JUMP).fit(amplitude_computed.flatten())
    breakpoints = algo.predict(pen=penalty)

    # Convert to the boolean array
    breakpoints_xarray = amplitude_ts.copy()
    breakpoints_xarray.name = "breakpoints"
    breakpoints_list = [False for _ in range(amplitude_ts.time.shape[0])]
    for breakpoint_ in breakpoints:
        if breakpoint_ < amplitude_ts.time.shape[0]:
            breakpoints_list[breakpoint_] = True
    breakpoints_xarray.data = da.array(breakpoints_list).reshape(amplitude_computed.shape)

    return breakpoints_xarray


def _binseg_block(amplitude_ts: xr.Dataset, cost_model: str, size: int) -> xr.Dataset:
    """Compute the BinSeg cost function for breakpoints in chunks.

    Parameters
    ----------
    amplitude_ts: xr.Dataset
      dataset containing the amplitude values of a chunk of points.
    cost_model: str
      Which cost model to use. Default (and recommended) "l2"
    size: int
      minimum size of the partitions. Default 27

    Returns
    -------
    xr.Dataset
      Boolean dataset of the same size of the amplitude array, where True indicates a detected breakpoint
    """
    groups = amplitude_ts.groupby("space")
    stmat_out = groups.map(_binseg_single_point, cost_model=cost_model, size=size)
    return stmat_out


def _binseg_single_point(amplitude_ts: xr.Dataset, cost_model: str, size: int) -> xr.Dataset:
    """Compute the BinSeg cost function for breakpoints for a single point.

    Parameters
    ----------
    amplitude_ts: xr.Dataset
      dataset containing the amplitude values of one point.
    cost_model: str
      Which cost model to use. Default (and recommended) "l2"
    size: int
      minimum size of the partitions. Default 27

    Returns
    -------
    xr.Dataset
      Boolean dataset of the same size of the amplitude array, where True indicates a detected breakpoint
    """
    amplitude_computed = amplitude_ts.compute().data
    penalty = np.var(amplitude_computed) * np.log(amplitude_ts.time.shape[0])
    algo = rpt.Binseg(model=cost_model).fit(amplitude_computed)
    all_breakpoints = algo.predict(pen=penalty)

    # Binseg does not check on partition size, so we do that manually by removing breakpoints too close to the previous
    # detected breakpoint
    if len(all_breakpoints) > 0:  # we detected something
        breakpoints = [all_breakpoints[0]]
        for bkp in all_breakpoints[1:]:
            if bkp - breakpoints[-1] >= size:
                breakpoints.append(bkp)
    else:
        breakpoints = []

    # Convert to the boolean array
    breakpoints_xarray = amplitude_ts.copy()
    breakpoints_xarray.name = "breakpoints"
    breakpoints_list = [False for _ in range(amplitude_ts.time.shape[0])]
    for breakpoint_ in breakpoints:
        if breakpoint_ < amplitude_ts.time.shape[0]:
            breakpoints_list[breakpoint_] = True
    breakpoints_xarray.data = da.array(breakpoints_list).reshape(amplitude_computed.shape)

    return breakpoints_xarray


def _detect_outliers(
    amplitude_array: xr.DataArray, db_outlier_detection: bool = True, window_size: int = 15, n_sigma: int = 3
) -> xr.Dataset:
    """Detect outliers based on a hampel filter.

    Parameters
    ----------
    amplitude_array: xr.DataArray
      DataArray containing the amplitude values
    db_outlier_detection: bool, optional
      whether or not to do outlier detection in dB. Default and advised True
    window_size: int, optional
      window size of the hampel filter used for detection. Default 15
    n_sigma: int, optional
      number of standard deviations difference required before outlier is detected. Default 3


    Returns
    -------
    xr.Dataset
      A boolean dataset with True indicating an outlier detected for that point at that epoch

    """
    match db_outlier_detection:
        case True:
            amplitude_ts = 10 * np.log10(amplitude_array)
        case False:
            amplitude_ts = amplitude_array
        case _:
            raise ValueError(f"db_outlier_detection should be boolean but is {db_outlier_detection}!")

    # Set up the filter value using a hampel filter with the window size
    filter_value = np.zeros((amplitude_ts.shape[0], amplitude_ts.shape[1], window_size))
    for count, shift in enumerate(range(-(window_size - 1) // 2, (window_size - 1) // 2 + 1)):
        filter_value[:, :, count] = np.roll(amplitude_ts, shift, axis=1)
    # Calculate the critical value
    x0 = np.median(filter_value, axis=2)
    critical_value = (
        1.4826 * n_sigma * np.median(np.abs(filter_value - x0.reshape((x0.shape[0], x0.shape[1], 1))), axis=2)
    )
    # Detect the outliers
    outliers = np.abs(amplitude_ts - x0) > critical_value
    # Since np.roll rolls over the end to the start, the first part and last part of the outliers array is based on
    # disconnected data. We turn these to False.
    outliers[:, : (window_size - 1) // 2] = False
    outliers[:, -(window_size - 1) // 2 :] = False
    outliers_xarray = amplitude_ts.copy()
    outliers_xarray.name = "breakpoints"
    outliers_xarray.data = outliers
    return outliers_xarray


def _nad_nmad_quality_metrics(
    nad_nmad: xr.DataArray, input_mode: Literal["nad", "nmad"], output_mode: Literal["mean", "2sigma"] = "2sigma"
) -> xr.DataArray:
    """Estimate the quality metrics for NAD and NMAD.

    These are empirical relations for a cubic function to relate the NMAD and NAD to mean cloud (50%, mean)
    and quality (95%, 2sigma).

    Parameters
    ----------
    nad_nmad: xr.DataArray
        Input array with NAD or NMAD values
    input_mode: Literal["nad", "nmad"]
        Type of data in the input array
    output_mode: Literal["mean", "2sigma"], optional
        Type of quality metric requested. Default 2sigma.

    Returns
    -------
    xr.DataArray
        DataArray with the quality metrics

    Raises
    ------
    AssertionError
      when `input_mode` or `output_mode` do not match the given options
    """
    assert input_mode in ["nad", "nmad"], f"Expected input_mode in ['nad', 'nmad'] but got {input_mode}!"
    assert output_mode in ["mean", "2sigma"], f"Expected output_mode in ['mean', '2sigma'] but got {output_mode}!"

    output = (
        NAD_NMAD_TO_SIGMA_CONVERSION[input_mode][output_mode][0]
        + NAD_NMAD_TO_SIGMA_CONVERSION[input_mode][output_mode][1] * nad_nmad
        + NAD_NMAD_TO_SIGMA_CONVERSION[input_mode][output_mode][2] * nad_nmad**2
        + NAD_NMAD_TO_SIGMA_CONVERSION[input_mode][output_mode][3] * nad_nmad**3
    )

    return output
