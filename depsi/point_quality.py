from typing import Literal

import dask.array as da
import numpy as np
import ruptures as rpt
import xarray as xr

from depsi.classification import _nad_block, _nmad_block
from depsi.utils import crop_slc_spacetime, npdatetime64_to_datetime

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
# The following list is used to check which output variables are valid in the STM partitioning.
POSSIBLE_OUTPUT_VARIABLES_PARTITIONING = [
    "nad",
    "nmad",
    "mad",
    "quality_nmad_2sigma",
    "quality_nmad_mean",
    "quality_nad_2sigma",
    "quality_nad_mean",
    "amplitude_mean",
    "amplitude_sigma",
    "amplitude_median",
]


def stm_partitioning(
    stm: xr.Dataset,
    db_partitioning: bool = False,
    search_method: Literal["pelt", "binseg"] = "pelt",
    cost_model: str = "l2",
    min_partition_size: int = 27,
    amplitude_variable_name: str = "amplitude",
    output_variable_prefix: str = "partition",
    output_variables: tuple = ("nad", "nmad", "quality_nmad_2sigma"),
) -> xr.Dataset:
    """Perform partitioning of a space-time matrix based on amplitude data.

    Identifies breakpoints in the amplitude timeseries of a set of points, based on a given search method, cost model
    and minimum partition size between the breakpoints. The search can be performed on the normal amplitude data or on
    a dB scale. Different output variables (quality metrics) can be defined to be calculated per partition.

    If the layers `breakpoints` and `partition_id` already exist, these are used as the partitions instead.

    Parameters
    ----------
    stm: xr.Dataset
      The dataset with the amplitude of the points and coordinates `space` and `time`
    db_partitioning: bool, optional
      Toggle to turn partitioning on the dB scale on (True) or off (False). Defaults to False.
    search_method: Literal["pelt", "binseg"], optional
      Which search method to use. Default (and recommended) "pelt"
    cost_model: str, optional
      Which cost model to use. Default (and recommended) "l2"
    min_partition_size: int, optional
      minimum size of the partitions. Default 27
    amplitude_variable_name: str, optional
      the name of the variable in which the amplitude information is stored in the STM. Default `amplitude`
    output_variable_prefix: str, optional
      the prefix of the output variables, named `output_variable_prefix`_`output_variable`. Default `partition`
    output_variables: tuple, optional
      the name of the requested output variables per partition. Default ("nad", "nmad", "quality_nmad_2sigma"). Options:
      - nad - the NAD of the partition
      - nmad - the NMAD of the partition
      - mad - the MAD of the partition
      - quality_nmad_2sigma - the quality of the partition based on an empirical relation, see
        `_nad_nmad_quality_metrics`, based on NMAD and 2-sigma
      - quality_nmad_mean - the quality of the partition based on an empirical relation, see
        `_nad_nmad_quality_metrics`, based on NMAD and mean cloud
      - quality_nad_2sigma - the quality of the partition based on an empirical relation, see
        `_nad_nmad_quality_metrics`, based on NAD and 2-sigma
      - quality_nad_mean - the quality of the partition based on an empirical relation, see
        `_nad_nmad_quality_metrics`, based on NAD and mean cloud
      - amplitude_mean - the mean of the amplitude
      - amplitude_sigma - the standard deviation of the amplitude
      - amplitude_median - the median of the amplitude

    Returns
    -------
    xr.Dataset
      the same xr.Dataset with new variables:
      - breakpoints (space, time): boolean array of the breakpoint locations
      - partition_id (space, time): unique identifier for each partition
      - `output_variable_prefix`_`output_variable` (space, time): the requested output variables from
        `output_variables`, where the provided values within a partition are at all epochs equal to the requested
        output variable.

    Raises
    ------
    AssertionError
      - when the amplitude layer name `amplitude_variable_name` is not in the STM
    NotImplementedError
      - when there are variables in `output_variables` that do not exist in the documentation.
    """
    assert amplitude_variable_name in stm.keys(), f"Expected key {amplitude_variable_name} but it is not there!"

    for variable in output_variables:
        if variable not in POSSIBLE_OUTPUT_VARIABLES_PARTITIONING:
            raise NotImplementedError(
                f"Variable {variable} requested but not implemented. Currently implemented are "
                f"{POSSIBLE_OUTPUT_VARIABLES_PARTITIONING}."
            )

    if "breakpoints" not in stm.keys() or "partition_id" not in stm.keys():
        breakpoints, partition_identifiers = _estimate_breakpoints(
            stm[amplitude_variable_name],
            db_partitioning,
            search_method,
            cost_model,
            min_partition_size,
        )
        stm = stm.assign({"breakpoints": (["space", "time"], breakpoints)})
        stm = stm.assign({"partition_id": (["space", "time"], partition_identifiers)})

    if len(output_variables) > 0:
        # persist the amplitude values to memory to facilitate IO during the groups
        groups = stm[amplitude_variable_name]
        groups.data = groups.values
        groups = groups.groupby(stm["partition_id"])

        partition_stats = groups.map(_compute_partition_nad_nmad_amp_stats)
        for output_variable in output_variables:
            if "quality" in output_variable:
                _, var, metric = output_variable.split("_")
                output = _nad_nmad_quality_metrics(partition_stats[f"partition_{var}"].data, var, metric)
            else:
                output = partition_stats[f"partition_{output_variable}"].data
            stm = stm.assign({f"{output_variable_prefix}_{output_variable}": (["space", "time"], output)})
    return stm


def stm_add_incremental_recal_nad_nmad(
    stm: xr.Dataset,
    method: Literal["nmad", "nad"] = "nmad",
    mode: Literal["incremental", "recalibration"] = "recalibration",
    recalibration_jump_size: int = 10,
) -> xr.Dataset:
    """Calculate the incremental or recalibration NAD or NMAD of an STM.

    Incremental NAD or NMAD yields the NAD or NMAD of all images up to and including that epoch.
    Recalibration NAD or NMAD yields the NAD or NMAD of all images up to and including the last recalibration epoch,
    dictated by `recalibration_jump_size`. In case the STM was generated using an initialization period, this is
    automatically detected and alters the output as follows:
    - Any epoch before the initialization period is set to 0
    - Any epoch during the initialization period is set to the NAD or NMAD of the full initialization period
    - Any epoch after the initialization period will use the first epoch of the initialization period as first epoch
      for the NAD and NMAD estimation. The recalibration epochs will start counting after the end of the initialization
      period.

    Parameters
    ----------
    stm: xr.Dataset
      The dataset with the amplitude of the points and coordinates `space` and `time`, attributes
      `ps_selection_start_date` and `ps_selection_end_date`, and variable `amplitude`
    method: Literal["nmad", "nad"], optional
      Which method to use, either NMAD or NAD. Default "nmad"
    mode: Literal["incremental", "recalibration"], optional
      Which output mode to use, either incremental or recalibration. Default "recalibration"
    recalibration_jump_size: int, optional
      The number of epochs the NMAD or NAD should remain constant in recalibration mode. Default 10. This is ignored
      when mode "incremental" is selected.

    Returns
    -------
    xr.Dataset
      The same dataset with one new variable, e.g. `recalibration_nmad`. Always formatted as `mode`_`method`.

    Raises
    ------
    AssertionError
      Raised when:
      - method is not nad or nmad
      - mode is not incremental or recalibration
      - "ps_selection_start_date", "ps_selection_end_date", "time", "space" or "amplitude" is missing from the STM
    """
    assert method in ["nad", "nmad"], f"Method is {method}, expected nad or nmad!"
    assert mode in ["incremental", "recalibration"], f"Mode is {mode}, expected incremental or recalibration"
    for key in ["time", "amplitude", "space"]:
        assert key in stm.keys(), f"Expected STM with key {key} but it is not there!"
    for key in ["ps_selection_start_date", "ps_selection_end_date"]:
        assert key in stm.attrs, f"Expected STM with attribute {key} but it is not there!"

    imgs = []
    recalibration_idx = -1  # start at -1 so that the first addition will trigger a reset of the data layer
    current_image = None

    # detect if an initialization epoch was used in the PS selection
    initialization_mode = True
    first_epoch = npdatetime64_to_datetime(stm["time"].values[0])
    last_epoch = npdatetime64_to_datetime(stm["time"].values[-1])
    if (stm.ps_selection_start_date == first_epoch.strftime("%Y%m%d")) and (
        stm.ps_selection_end_date == last_epoch.strftime("%Y%m%d")
    ):
        initialization_mode = False

    # loop over the time values
    for date in stm["time"].values:
        # determine the time crop and recalibration index for the current image
        if initialization_mode:
            current_epoch = npdatetime64_to_datetime(date)
            current_epoch_formatted = f"{current_epoch.year}{current_epoch.month:0>2d}{current_epoch.day:0>2d}"

            if current_epoch_formatted < stm.ps_selection_start_date:
                start_date = stm.ps_selection_start_date
                end_date = stm.ps_selection_start_date
                if current_epoch_formatted == f"{first_epoch.year}{first_epoch.month:0>2d}{first_epoch.day:0>2d}":
                    recalibration_idx = 0  # compute the first one
                else:
                    recalibration_idx = -1  # just reuse the previous one, they are set to zero
            elif current_epoch_formatted <= stm.ps_selection_end_date:
                start_date = stm.ps_selection_start_date
                end_date = stm.ps_selection_end_date
                if current_epoch_formatted in [stm.ps_selection_start_date, stm.ps_selection_end_date]:
                    recalibration_idx = 0  # only recompute for the first and last image
                else:
                    recalibration_idx = -1  # just reuse the previous one, they are all the same
            else:
                start_date = stm.ps_selection_start_date
                end_date = current_epoch_formatted
                recalibration_idx += 1  # add, and do modulo the jump size, so that it will be 0 every time a new image
                # should be loaded
                recalibration_idx %= recalibration_jump_size
        else:
            start_date = npdatetime64_to_datetime(stm["time"].values[0])
            end_date = npdatetime64_to_datetime(date)
            recalibration_idx += 1  # add, and do modulo the jump size, so that it will be 0 every time a new image
            # should be loaded
            recalibration_idx %= recalibration_jump_size

        # Check if we need to recompute or can use the previous one
        match mode:
            case "incremental":  # always requires an update
                update = True
            case "recalibration":  # only requires an update when the index is 0
                update = recalibration_idx == 0
            case _:
                raise NotImplementedError(f"Unknown mode {mode}, known are recalibration and incremental!")

        if update:  # only recompute if `update` is triggered, otherwise use the previous one
            current_crop = crop_slc_spacetime(stm, start_date=start_date, end_date=end_date)
            match method:
                case "nad":
                    block_func = _nad_block
                case "nmad":
                    block_func = _nmad_block
            current_image = xr.map_blocks(
                block_func,
                current_crop["amplitude"],
                template=current_crop["amplitude"].isel(time=0).drop_vars("time"),
            )

        imgs.append(current_image.copy())

    layer = da.vstack(imgs).T
    stm = stm.assign({f"{mode}_{method}": (["space", "time"], layer)})

    # Rechunk is needed because after calculating incremental NAD/NMAD, the chunk size in time will be inconsistent
    stm = stm.chunk({"time": -1})

    return stm


def _estimate_breakpoints(
    amplitude_array: xr.DataArray,
    db_partitioning: bool = False,
    search_method: Literal["pelt", "binseg"] = "pelt",
    cost_model: str = "l2",
    min_partition_size: int = 27,
) -> tuple[np.ndarray, np.ndarray]:
    """Identify breakpoints in the amplitude timeseries of a DataArray of amplitude information.

    Identifies breakpoints in the amplitude timeseries of a set of points, based on a given search method, cost model
    and minimum partition size between the breakpoints. The search can be performed on the normal amplitude data or on
    a dB scale.

    Parameters
    ----------
    amplitude_array: xr.DataArray
      The data array with the amplitude time series of the points.
    db_partitioning: bool, optional
      Toggle to turn partitioning on the dB scale on (True) or off (False). Defaults to False.
    search_method: Literal["pelt", "binseg"], optional
      Which search method to use. Default (and recommended) "pelt"
    cost_model: str, optional
      Which cost model to use. Default (and recommended) "l2"
    min_partition_size: int, optional
      minimum size of the partitions. Default 27

    Returns
    -------
    np.ndarray
      boolean data array where True indicates a breakpoint at that epoch
    np.ndarray
      integer data array where each partition has been assigned a unique identifier.
    """
    if db_partitioning:
        amplitude_ts = 10 * da.log10(amplitude_array)
    else:
        amplitude_ts = amplitude_array

    match search_method:
        case "pelt":
            breakpoints = xr.map_blocks(
                _pelt_block,
                amplitude_ts,
                args=(cost_model, min_partition_size),
                template=amplitude_ts,
            )
        case "binseg":
            breakpoints = xr.map_blocks(
                _binseg_block,
                amplitude_ts,
                args=(cost_model, min_partition_size),
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
    breakpoints_idx_p1 = np.cumsum(computed_breakpoints).reshape(computed_breakpoints.shape)
    point_idx = np.arange(0, computed_breakpoints.shape[0]).reshape((computed_breakpoints.shape[0], 1))
    point_idx = np.hstack([point_idx for _ in range(computed_breakpoints.shape[1])])
    breakpoints_idx = breakpoints_idx_p1 + point_idx

    return computed_breakpoints, breakpoints_idx


def _pelt_block(amplitude_ts: xr.Dataset, cost_model: str, min_partition_size: int) -> xr.Dataset:
    """Compute the pelt cost function for breakpoints in chunks.

    Parameters
    ----------
    amplitude_ts: xr.Dataset
      dataset containing the amplitude values of a chunk of points.
    cost_model: str
      Which cost model to use
    min_partition_size: int
      minimum size of the partitions

    Returns
    -------
    xr.Dataset
      Boolean dataset of the same size of the amplitude array, where True indicates a detected breakpoint
    """
    groups = amplitude_ts.groupby("space")
    stmat_out = groups.map(_pelt_single_point, cost_model=cost_model, min_partition_size=min_partition_size)
    return stmat_out


def _pelt_single_point(amplitude_ts: xr.Dataset, cost_model: str, min_partition_size: int) -> xr.Dataset:
    """Compute the pelt cost function for breakpoints for a single point.

    Parameters
    ----------
    amplitude_ts: xr.Dataset
      dataset containing the amplitude values of one point.
    cost_model: str
      Which cost model to use
    min_partition_size: int
      minimum size of the partitions

    Returns
    -------
    xr.Dataset
      Boolean dataset of the same size of the amplitude array, where True indicates a detected breakpoint
    """
    amplitude_computed = amplitude_ts.compute().data
    penalty = np.var(amplitude_computed) * np.log(amplitude_ts.time.shape[0])
    algo = rpt.Pelt(model=cost_model, min_size=min_partition_size, jump=PELT_JUMP).fit(amplitude_computed.flatten())
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


def _binseg_block(amplitude_ts: xr.Dataset, cost_model: str, min_partition_size: int) -> xr.Dataset:
    """Compute the BinSeg cost function for breakpoints in chunks.

    Parameters
    ----------
    amplitude_ts: xr.Dataset
      dataset containing the amplitude values of a chunk of points.
    cost_model: str
      Which cost model to use
    min_partition_size: int
      minimum size of the partitions

    Returns
    -------
    xr.Dataset
      Boolean dataset of the same size of the amplitude array, where True indicates a detected breakpoint
    """
    groups = amplitude_ts.groupby("space")
    stmat_out = groups.map(_binseg_single_point, cost_model=cost_model, min_partition_size=min_partition_size)
    return stmat_out


def _binseg_single_point(amplitude_ts: xr.Dataset, cost_model: str, min_partition_size: int) -> xr.Dataset:
    """Compute the BinSeg cost function for breakpoints for a single point.

    Parameters
    ----------
    amplitude_ts: xr.Dataset
      dataset containing the amplitude values of one point.
    cost_model: str
      Which cost model to use
    min_partition_size: int
      minimum size of the partitions

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
            if bkp - breakpoints[-1] >= min_partition_size:
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


def detect_outliers_stm(
    stm: xr.Dataset, db_outlier_detection: bool = True, window_size: int = 15, n_sigma: int = 3
) -> xr.Dataset:
    """Detect outliers based on a hampel filter.

    Parameters
    ----------
    stm: xr.Dataset
      Dataset containing the amplitude values as an amplitude variable, with coordinates `space` and `time`
    db_outlier_detection: bool, optional
      whether or not to do outlier detection in dB. Default and advised True
    window_size: int, optional
      window size of the hampel filter used for detection. Default 15
    n_sigma: int, optional
      number of standard deviations difference required before outlier is detected. Default 3


    Returns
    -------
    xr.Dataset
      The original dataset with a new variable 'outliers': a boolean variable with True indicating an outlier detected
      for that point at that epoch

    """
    for key in ["space", "time", "amplitude"]:
        assert key in stm.keys(), f"Expected {key} in stm but it is missing!"
    outliers = xr.map_blocks(
        _detect_outliers,
        stm["amplitude"],
        kwargs={
            "db_outlier_detection": db_outlier_detection,
            "window_size": window_size,
            "n_sigma": n_sigma,
        },
        template=stm["amplitude"],
    )
    stm = stm.assign({"outliers": (["space", "time"], outliers.data)})

    return stm


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
    if db_outlier_detection:
        amplitude_ts = 10 * np.log10(amplitude_array)
    else: 
        amplitude_ts = amplitude_array

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


def _compute_partition_nad_nmad_amp_stats(amp: xr.DataArray) -> xr.DataArray:
    """Compute the NAD and NMAD and the amplitude statistics on a partition basis.

    Parameters
    ----------
    amp: xr.DataArray
      a DataArray containing the amplitude values of a single partition, loaded into memory

    Returns
    -------
    xr.DataArray
      the same DataArray, with in the coordinates the calculated NAD and NMAD value in the same shape as the partition.

    """
    data = amp.data
    std = np.std(data)
    mean = np.mean(data)
    nad = std / (mean + np.finfo(data.dtype).eps)
    median = np.median(data)
    mad = np.median(np.abs(data - median))
    nmad = mad / (median + np.finfo(data.dtype).eps)

    amp["partition_nad"] = (amp.dims, np.ones_like(data) * nad)
    amp["partition_nmad"] = (amp.dims, np.ones_like(data) * nmad)
    amp["partition_amplitude_sigma"] = (amp.dims, np.ones_like(data) * std)
    amp["partition_amplitude_mean"] = (amp.dims, np.ones_like(data) * mean)
    amp["partition_mad"] = (amp.dims, np.ones_like(data) * mad)
    amp["partition_amplitude_median"] = (amp.dims, np.ones_like(data) * median)
    return amp
