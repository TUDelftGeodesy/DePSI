from typing import Literal

import dask.array as da
import numpy as np
import ruptures as rpt
import xarray as xr

PELT_JUMP = 5


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
