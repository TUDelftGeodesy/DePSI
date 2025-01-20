"""Functions for scatterer selection related operations."""

from datetime import datetime
from typing import Literal

import dask.array as da
import numpy as np
import pyproj
import pytz
import xarray as xr
from scipy.spatial import KDTree

from depsi.point_quality import _detect_outliers, _estimate_breakpoints, _nad_nmad_quality_metrics
from depsi.utils import _npdatetime64_to_datetime, crop_slc_spacetime

# partitioning and outlier detection when requested in ps_selection require a fixed set of keywords to function.
# The following lines are the required keywords against which the input dictionaries will be checked.
REQUIRED_PARTITIONING_KEYS = ["db_partitioning", "search_method", "cost_function", "min_obs_partition"]
REQUIRED_OUTLIER_DETECTION_KEYS = ["db_outlier_detection", "window_size", "n_sigma"]


def ps_selection(
    slcs: xr.Dataset,
    threshold: float,
    method: Literal["nad", "nmad"] = "nad",
    output_chunks: int = 10000,
    mem_persist: bool = False,
    ps_selection_start_date: datetime | str | None = None,
    ps_selection_end_date: datetime | str | int | None = None,
    recalibration_jump_size: int = 10,
    do_rd_coordinate_conversion: bool = False,
    do_partitioning: bool = False,
    partitioning_kwargs: dict | None = None,
    do_outlier_detection: bool = False,
    outlier_detection_kwargs: dict | None = None,
    single_difference_mother: datetime | str = "auto",
) -> xr.Dataset:
    """Select Persistent Scatterers (PS) from an SLC stack, and return a Space-Time Matrix.

    The selection method is defined by `method` and `threshold`.
    The selected pixels will be reshaped to (space, time), where `space` is the number of selected pixels.
    The unselected pixels will be discarded.
    The original `azimuth` and `range` coordinates will be persisted.
    The computed NAD or NMAD will be added to the output dataset as a new variable. It can be persisted in
    memory if `mem_persist` is True.
    The original time axis will be preserved in all cases. If `ps_selection_start_date` and `ps_selection_end_date`
    are provided, the selection will only use the images in the provided time window. However, the full time axis
    will be preserved and returned. In this case, the layers `selection_nmad` / `selection_nad` and `full_ts_nmad` /
    `full_ts_nad` are thus different.

    Parameters
    ----------
    slcs : xr.Dataset
        Input SLC stack. It should have the following dimensions: ("azimuth", "range", "time").
        There should be a `amplitude` variable in the dataset.
    threshold : float
        Threshold value for selection for "nad" / "nmad".
    method : Literal["nad", "nmad"], optional
        Method of selection, by default "nad".
        - "nad": Normalized Amplitude Dispersion
        - "nmad": Normalized median absolute deviation
    output_chunks : int, optional
        Chunk size in the `space` dimension, by default 10000
    mem_persist : bool, optional
        If true persist the NAD or NMAD in memory, by default False.
    ps_selection_start_date : datetime | str | None, optional
      the start date of the time window to be used for the ps_selection, in one of three formats:
      - datetime object
      - str object, formatted as YYYYMMDD
      - None, no cropping in time requested for the ps_selection (default)
    ps_selection_end_date : datetime | str | int | None, optional
      the end date of the time window to be used for the ps_selection, in one of four formats:
      - datetime object
      - str object, formatted as YYYYMMDD
      - int object, which is interpreted as the number of images intended in the crop (including the start date). If
        more images are requested than exist since the start date, all images from start_date until the last image
        are provided.
      - None, no cropping in time requested for the ps_selection (default)
    recalibration_jump_size: int, optional
      the number of images that the recalibration NAD / NMAD variables remains constant. Will start after the
      initialization epoch (if ps_selection_start_date and ps_selection_end_date are not None). Defaults to 10.
    do_rd_coordinate_conversion: bool, optional
      boolean to trigger coordinate conversion from latitude/longitude (WGS84) to RD_X/RD_Y (Rijksdriehoek). This only
      makes sense for AoIs located in the Netherlands. Defaults to False.
    do_partitioning: bool, optional
      boolean to trigger the breakpoint analysis. Defaults to False.
    partitioning_kwargs: dict | None, optional
      the keyword arguments required for the breakpoint analysis. Required if do_breakpoint_analysis is set to True.
      Formatted as a dictionary with required keys:
      - db_partitioning: True or False, whether or not to do partitioning in dB. Advised False
      - search_method: 'pelt' or 'binseg'. Advised 'pelt'
      - cost_function: 'l#' with # replaced by 0-3. Advised 'l2'
      - min_obs_partition: integer. Advised min 0.5 years converted to # images, for Sentinel-1 27 (6 day interval)
    do_outlier_detection: bool, optional
      boolean to trigger the outlier detection. Defaults to False.
    outlier_detection_kwargs: dict | None, optional
      the keyword arguments required for the outlier detection. Required if do_outlier_detection is set to True.
      Formatted as a dictionary with required keys:
      - db_outlier_detection: True or False, whether or not to do outlier detection in dB. Advised True
      - window_size: window size of the hampel filter used for detection. Advised 15
      - n_sigma: number of standard deviations difference required before outlier is detected. Advised 3
    single_difference_mother: datetime | str
      the date to be used as the mother image for the single difference computations, in one of three formats:
      - 'auto' : will detect the mother image in the input SLC dataset, and use that epoch.
      - datetime object
      - str object, formatted as YYYYMMDD

    Returns
    -------
    xr.Dataset
        Selected STM, in form of an xarray.Dataset with dimensions:
        - space ( # PS selected)
        - time ( # epochs of input dataset)
        with coordinates:
        - time: epoch in np.datetime64 format
        - space: index of the PS
        - azimuth: azimuth coordinate of the PS
        - range: range coordinate of the PS
        with attributes:
        - ps_selection_start_date: the epoch of the first image used for the PS selection
        - ps_selection_end_date: the epoch of the last image used for the PS selection
        - ps_sd_mother: the epoch of the mother used for the single differences
        with variables:
        - h2ph (space, time): the height to phase conversion
        - lat (space): latitude of the PS
        - lon (space): longitude of the PS
        - complex (space, time): the complex value of the PS at each epoch
        - amplitude (space, time): the amplitude of the PS at each epoch
        - phase (space, time): the phase of the PS at each epoch
        - selection_nad / selection_nmad (space): the value used for selection of the PS, dependent on method
        - full_ts_nad (space): the Normalized Amplitude Dispersion of the PS
        - full_ts_nmad (space): the Normalized Median Amplitude Dispersion of the PS
        - incremental_nad (space, time): the NAD of all images up to and including that epoch (initialization epochs
            will yield the initialization period NAD if selected). If ps_selection_start_date is not None, images
            before ps_selection_start_date will yield 0
        - incremental_nmad (space, time): the NMAD of all images up to and including that epoch (initialization epochs
            will yield the initialization period NMAD if selected). If ps_selection_start_date is not None, images
            before ps_selection_start_date will yield 0
        - recalibration_nad (space, time): the NAD of all images up to and including the last recalibration epoch
            (dictated by recalibration_jump_size, initialization epochs will yield the initialization period NAD if
            selected). If ps_selection_start_date is not None, images before ps_selection_start_date will yield 0
        - recalibration_nmad (space, time): the NMAD of all images up to and including the last recalibration epoch
            (dictated by recalibration_jump_size, initialization epochs will yield the initialization period NMAD if
            selected). If ps_selection_start_date is not None, images before ps_selection_start_date will yield 0
        - sd_h2ph (space, time): single difference height to phase conversion with respect to single_difference_mother
        - sd_complex (space, time): single difference complex phasor with respect to single_difference_mother
        - sd_amplitude_unnormalized (space, time): single difference complex phasor amplitude to
            single_difference_mother, not normalized
        - sd_phase (space, time): single difference phase with respect to single_difference_mother
        - days_since_first_img (time): number of days since the first image
        - years_since_first_img (time): number of years since the first image (assuming 365.2425 days per year)
        - classification_flag (space): 1 for all selected PS
        If do_rd_coordinate_conversion is set to True:
        - rd_x (space): Rijksdriehoek x coordinate if requested (only recommended in the Netherlands)
        - rd_y (space): Rijksdriehoek y coordinate if requested (only recommended in the Netherlands)
        If do_partitioning is set to True:
        - breakpoints (space, time): boolean array of the breakpoint locations
        - partition_id (space, time): unique identifier for each partition
        - partition_nmad (space, time): NMAD calculated per partition
        - partition_nmad_quality (space, time): NMAD per partition converted to standard deviation using mean + 2 sigma
        - partition_nad (space, time): NAD calculated per partition
        - partition_sd_amplitude_sigma (space, time): standard deviation of the unnormalized single difference
            amplitude per partition
        - partition_sd_amplitude_mean (space, time): mean of the unnormalized single difference
            amplitude per partition
        - partition_sd_amplitude_median (space, time): median of the unnormalized single difference
            amplitude per partition
        - partition_sd_mad (space, time): Median Absolute Deviation of the unnormalized single difference
            amplitude per partition
        If do_outlier_detection is set to True:
        - outliers (space, time): boolean array, where True indicates an outlier detected based on the
            outlier_detection_kwargs and a hampel filter

    Raises
    ------
    NotImplementedError
        Raised when an unsupported method is provided.
    AssertionError
        Raised when:
        - do_partitioning is True but partitioning_kwargs is not provided
        - Keywords are missing in partitioning_kwargs
        - partitioning_kwargs is not a dictionary
        - do_outlier_detection is True but outlier_detection_kwargs is not provided
        - Keywords are missing in outlier_detection_kwargs
        - outlier_detection_kwargs is not a dictionary
    ValueError
        Raised when:
        - single_difference_mother is of an unsupported format
        - the date provided to single_difference_mother is not in the input stack
    """
    if do_partitioning:
        assert partitioning_kwargs is not None, "Breakpoint analysis requested without keyword arguments!"
        assert isinstance(
            partitioning_kwargs, dict
        ), f"breakpoint_kwargs should be dict but is {type(partitioning_kwargs)}"
        assert np.all(
            [key in partitioning_kwargs.keys() for key in REQUIRED_PARTITIONING_KEYS]
        ), f"Keys {REQUIRED_PARTITIONING_KEYS} are required but received {partitioning_kwargs.keys()}!"

    if do_outlier_detection:
        assert outlier_detection_kwargs is not None, "Outlier detection requested without keyword arguments!"
        assert isinstance(
            outlier_detection_kwargs, dict
        ), f"outlier_detection_kwargs should be dict but is {type(outlier_detection_kwargs)}"
        assert np.all(
            [key in outlier_detection_kwargs.keys() for key in REQUIRED_OUTLIER_DETECTION_KEYS]
        ), f"Keys {REQUIRED_OUTLIER_DETECTION_KEYS} are required but received {outlier_detection_kwargs.keys()}!"

    # Make sure there is no temporal chunk
    # since later a block function assumes all temporal data is available in a spatial block
    slcs = slcs.chunk({"time": -1})

    # Calculate selection mask
    match method:
        case "nad":
            if ps_selection_start_date is not None:
                selection_crop = crop_slc_spacetime(
                    slcs, start_date=ps_selection_start_date, end_date=ps_selection_end_date
                )
                nad = xr.map_blocks(
                    _nad_block,
                    selection_crop["amplitude"],
                    template=selection_crop["amplitude"].isel(time=0).drop_vars("time"),
                )
                ps_selection_times = selection_crop["time"].values
            else:
                nad = xr.map_blocks(
                    _nad_block, slcs["amplitude"], template=slcs["amplitude"].isel(time=0).drop_vars("time")
                )
                ps_selection_times = []
            nad = nad.compute() if mem_persist else nad
            slcs = slcs.assign(selection_nad=nad)
            mask = nad < threshold
        case "nmad":
            if ps_selection_start_date is not None:
                selection_crop = crop_slc_spacetime(
                    slcs, start_date=ps_selection_start_date, end_date=ps_selection_end_date
                )
                nmad = xr.map_blocks(
                    _nmad_block,
                    selection_crop["amplitude"],
                    template=selection_crop["amplitude"].isel(time=0).drop_vars("time"),
                )
                ps_selection_times = selection_crop["time"].values
            else:
                nmad = xr.map_blocks(
                    _nmad_block, slcs["amplitude"], template=slcs["amplitude"].isel(time=0).drop_vars("time")
                )
                ps_selection_times = []
            nmad = nmad.compute() if mem_persist else nmad
            slcs = slcs.assign(selection_nmad=nmad)
            mask = nmad < threshold
        case _:
            raise NotImplementedError

    # Get the 1D index on space dimension
    mask_1d = mask.stack(space=("azimuth", "range")).drop_vars(["azimuth", "range", "space"])  # Drop multi-index coords
    index = mask_1d["space"].where(mask_1d.compute(), other=0, drop=True)  # Evaluate the 1D mask to index

    # Reshape from Stack ("azimuth", "range", "time") to Space-Time Matrix  ("space", "time")
    stacked = slcs.stack(space=("azimuth", "range"))

    # Drop multi-index coords for space coordinates
    # This will also azimuth and range coordinates, as they are part of the multi-index coordinates
    stm = stacked.drop_vars(["space", "azimuth", "range"])

    # Assign a continuous index the space dimension
    # Assign azimuth and range back as coordinates
    stm = stm.assign_coords(
        {
            "space": (["space"], range(stm.sizes["space"])),
            "azimuth": (["space"], stacked["azimuth"].values),
            "range": (["space"], stacked["range"].values),
        }
    )  # keep azimuth and range as coordinates

    # Apply selection
    stm_masked = stm.sel(space=index)

    # Re-order the dimensions to community preferred ("space", "time") order
    stm_masked = stm_masked.transpose("space", "time")

    # Rechunk is needed because after apply maksing, the chunksize will be inconsistant
    stm_masked = stm_masked.chunk(
        {
            "space": output_chunks,
            "time": -1,
        }
    )

    # Reset space coordinates
    stm_masked = stm_masked.assign_coords(
        {
            "space": (["space"], range(stm_masked.sizes["space"])),
        }
    )

    # add full timeseries NAD and NMAD
    nad = xr.map_blocks(
        _nad_block, stm_masked["amplitude"], template=stm_masked["amplitude"].isel(time=0).drop_vars("time")
    )
    nmad = xr.map_blocks(
        _nmad_block, stm_masked["amplitude"], template=stm_masked["amplitude"].isel(time=0).drop_vars("time")
    )
    stm_masked = stm_masked.assign({"full_ts_nmad": (["space"], nmad.data)})
    stm_masked = stm_masked.assign({"full_ts_nad": (["space"], nad.data)})
    if ps_selection_start_date is None:
        start_date = _npdatetime64_to_datetime(stm_masked["time"].values[0])
        end_date = _npdatetime64_to_datetime(stm_masked["time"].values[-1])

        stm_masked.attrs["ps_selection_start_date"] = f"{start_date.year}{start_date.month:0>2d}{start_date.day:0>2d}"
        stm_masked.attrs["ps_selection_end_date"] = f"{end_date.year}{end_date.month:0>2d}{end_date.day:0>2d}"
    else:
        start_date = _npdatetime64_to_datetime(ps_selection_times[0])
        end_date = _npdatetime64_to_datetime(ps_selection_times[-1])
        stm_masked.attrs["ps_selection_start_date"] = f"{start_date.year}{start_date.month:0>2d}{start_date.day:0>2d}"
        stm_masked.attrs["ps_selection_end_date"] = f"{end_date.year}{end_date.month:0>2d}{end_date.day:0>2d}"

    # add incremental and recalibration NAD / NMAD
    for loop_method in ["nmad", "nad"]:
        incremental_imgs = []
        recalibration_imgs = []
        recalibration_idx = -1  # start at -1 so that the first addition will trigger a reset of the data layer
        recalibration_data_layer = None
        for date in stm_masked["time"].values:
            if ps_selection_start_date is not None:
                if (
                    date in ps_selection_times
                ):  # only gets triggered in case there is an initialization epoch for the duration
                    # of the initialization epoch
                    start_date = ps_selection_start_date
                    end_date = ps_selection_end_date
                    recalibration_idx = 0
                else:
                    start_date = _npdatetime64_to_datetime(ps_selection_times[0])
                    end_date = _npdatetime64_to_datetime(date)
                    recalibration_idx += 1
                    # add, and do modulo the jump size, so that it will be 0 every time a new image
                    # should be loaded
                    recalibration_idx %= recalibration_jump_size
                    if start_date > end_date:
                        end_date = start_date  # will always be 0 as they are equal
            else:
                start_date = _npdatetime64_to_datetime(stm_masked["time"].values[0])
                end_date = _npdatetime64_to_datetime(date)
                recalibration_idx += 1  # add, and do modulo the jump size, so that it will be 0 every time a new image
                # should be loaded
                recalibration_idx %= recalibration_jump_size

            current_crop = crop_slc_spacetime(stm_masked, start_date=start_date, end_date=end_date)
            match loop_method:
                case "nad":
                    nad = xr.map_blocks(
                        _nad_block,
                        current_crop["amplitude"],
                        template=current_crop["amplitude"].isel(time=0).drop_vars("time"),
                    )
                    incremental_imgs.append(nad)
                    if recalibration_idx == 0:
                        recalibration_data_layer = nad.copy()
                case "nmad":
                    nmad = xr.map_blocks(
                        _nmad_block,
                        current_crop["amplitude"],
                        template=current_crop["amplitude"].isel(time=0).drop_vars("time"),
                    )
                    incremental_imgs.append(nmad)
                    if recalibration_idx == 0:
                        recalibration_data_layer = nmad.copy()

            recalibration_imgs.append(recalibration_data_layer.copy())

        # format the images, and add them to the data array
        incremental_nad_nmad = da.vstack(incremental_imgs).T
        recalibration_nad_nmad = da.vstack(recalibration_imgs).T
        match loop_method:
            case "nad":
                stm_masked = stm_masked.assign({"incremental_nad": (["space", "time"], incremental_nad_nmad)})
                stm_masked = stm_masked.assign({"recalibration_nad": (["space", "time"], recalibration_nad_nmad)})
            case "nmad":
                stm_masked = stm_masked.assign({"incremental_nmad": (["space", "time"], incremental_nad_nmad)})
                stm_masked = stm_masked.assign({"recalibration_nmad": (["space", "time"], recalibration_nad_nmad)})

    # Rechunk is needed because after calculating incremental NAD/NMAD, the chunksize will be inconsistant
    stm_masked_inc = stm_masked.chunk(
        {
            "space": output_chunks,
            "time": -1,
        }
    )

    # Compute the single differences
    # Identify the mother image
    if isinstance(single_difference_mother, datetime):
        format_mother_date = datetime(
            single_difference_mother.year, single_difference_mother.month, single_difference_mother.day, tzinfo=pytz.UTC
        )
        mother_index = [
            idx
            for idx, date in enumerate(stm_masked_inc["time"].values)
            if format_mother_date == _npdatetime64_to_datetime(date)
        ]
    elif isinstance(single_difference_mother, str):
        if single_difference_mother == "auto":
            mother_index = np.where(abs(stm_masked_inc["h2ph"]).sum(axis=0).values == 0)[0]
        elif len(single_difference_mother) == 8:
            format_mother_date = datetime(
                eval(single_difference_mother[:4]),
                eval(single_difference_mother[4:6].lstrip("0")),
                eval(single_difference_mother[6:].lstrip("0")),
                tzinfo=pytz.UTC,
            )
            mother_index = [
                idx
                for idx, date in enumerate(stm_masked_inc["time"].values)
                if format_mother_date == _npdatetime64_to_datetime(date)
            ]
        else:
            raise ValueError(f'Cannot parse {single_difference_mother}, not of type "auto" or "YYYYMMDD"!')
    else:
        raise ValueError(f"Unknown format {type(single_difference_mother)} for single_difference_mother!")
    if len(mother_index) == 0:
        raise ValueError(
            f"Cannot find provided mother date {single_difference_mother}, "
            f"please provide a date that is part of the stack! Possible dates: "
            f"{stm_masked_inc.time.values}"
        )
    sd_mother_index = mother_index[0]  # 0 in case somehow more than 1 image is detected
    # In that case we take the first image that was detected, as this is expected
    sd_mother = _npdatetime64_to_datetime(stm_masked_inc["time"].values[sd_mother_index])

    # Format the single difference mother, and save it to the STM
    sd_mother_formatted = f"{sd_mother.year}{sd_mother.month:0>2d}{sd_mother.day:0>2d}"
    stm_masked_inc.attrs["ps_sd_mother"] = sd_mother_formatted

    # calculate the h2ph single difference (= daughter - mother)
    sd_h2ph = stm_masked_inc["h2ph"] - stm_masked_inc["h2ph"][:, sd_mother_index]
    stm_masked_inc = stm_masked_inc.assign({"sd_h2ph": (["space", "time"], sd_h2ph.data)})

    # calculate the complex single difference, the amplitude, and the phase
    mother_comp = stm_masked_inc["complex"][:, sd_mother_index].conj()
    sd_complex_transposed = stm_masked_inc["complex"].transpose() * mother_comp
    sd_complex = sd_complex_transposed.transpose()
    sd_phase = da.angle(sd_complex)
    sd_amplitude_unnormalized = da.abs(sd_complex)
    stm_masked_inc = stm_masked_inc.assign({"sd_complex": (["space", "time"], sd_complex.data)})
    stm_masked_inc = stm_masked_inc.assign(
        {"sd_amplitude_unnormalized": (["space", "time"], sd_amplitude_unnormalized.data)}
    )
    stm_masked_inc = stm_masked_inc.assign({"sd_phase": (["space", "time"], sd_phase.data)})

    # Add RD coordinates if requested
    if do_rd_coordinate_conversion:
        wgs84 = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:28992", always_xy=True).transform
        # Convert Lat and Lon to RD-coordinates
        rd_x, rd_y = wgs84(stm_masked_inc["lon"], stm_masked_inc["lat"])

        # Add RD coordinates to the dataset
        stm_masked_inc = stm_masked_inc.assign({"rd_x": (["space"], rd_x)})
        stm_masked_inc = stm_masked_inc.assign({"rd_y": (["space"], rd_y)})

    # Add extra time coordinate variables for time intervals since first image
    days = np.array(
        [
            (_npdatetime64_to_datetime(date) - _npdatetime64_to_datetime(stm_masked["time"].values[0])).days
            for date in stm_masked["time"].values
        ]
    )
    stm_masked_inc = stm_masked_inc.assign({"days_since_first_img": (["time"], days)})
    stm_masked_inc = stm_masked_inc.assign({"years_since_first_img": (["time"], days / 365.2425)})

    if do_partitioning:
        breakpoints, partition_identifiers = _estimate_breakpoints(
            stm_masked_inc["amplitude"],
            partitioning_kwargs["db_partitioning"],
            partitioning_kwargs["search_method"],
            partitioning_kwargs["cost_function"],
            partitioning_kwargs["min_obs_partition"],
        )
        stm_masked_inc = stm_masked_inc.assign({"breakpoints": (["space", "time"], breakpoints.data)})
        stm_masked_inc = stm_masked_inc.assign({"partition_id": (["space", "time"], partition_identifiers)})

        # persist the amplitude values to memory to facilitate IO during the groups
        groups = stm_masked_inc["amplitude"]
        groups.data = groups.values
        groups = groups.groupby(stm_masked_inc["partition_id"])

        groups_sd = stm_masked_inc["sd_amplitude_unnormalized"]
        groups_sd.data = groups_sd.values
        groups_sd = groups_sd.groupby(stm_masked_inc["partition_id"])

        # Calculate the partition NAD and NMAD
        partition_stats = groups.map(_compute_partition_nad_nmad_amp_stats)
        partition_nmad = partition_stats["partition_nmad"].data
        partition_nad = partition_stats["partition_nad"].data

        partition_stats_sd = groups_sd.map(_compute_partition_nad_nmad_amp_stats)

        # calculate the quality metrics

        quality_nmad = _nad_nmad_quality_metrics(partition_nmad, "nmad", "2sigma")

        # Save to the STM
        stm_masked_inc = stm_masked_inc.assign({"partition_nmad": (["space", "time"], partition_nmad)})
        stm_masked_inc = stm_masked_inc.assign({"partition_nmad_quality": (["space", "time"], quality_nmad)})
        stm_masked_inc = stm_masked_inc.assign({"partition_nad": (["space", "time"], partition_nad)})
        stm_masked_inc = stm_masked_inc.assign(
            {"partition_sd_amplitude_sigma": (["space", "time"], partition_stats_sd["partition_amp_std"].data)}
        )
        stm_masked_inc = stm_masked_inc.assign(
            {"partition_sd_amplitude_mean": (["space", "time"], partition_stats_sd["partition_amp_mean"].data)}
        )
        stm_masked_inc = stm_masked_inc.assign(
            {"partition_sd_amplitude_median": (["space", "time"], partition_stats_sd["partition_amp_median"].data)}
        )
        stm_masked_inc = stm_masked_inc.assign(
            {"partition_sd_mad": (["space", "time"], partition_stats_sd["partition_amp_mad"].data)}
        )

    if do_outlier_detection:
        outliers = xr.map_blocks(
            _detect_outliers,
            stm_masked_inc["amplitude"],
            kwargs={
                "db_outlier_detection": outlier_detection_kwargs["db_outlier_detection"],
                "window_size": outlier_detection_kwargs["window_size"],
                "n_sigma": outlier_detection_kwargs["n_sigma"],
            },
            template=stm_masked_inc["amplitude"],
        )
        stm_masked_inc = stm_masked_inc.assign({"outliers": (["space", "time"], outliers.data)})

    # Add the classification flag
    stm_masked_inc = stm_masked_inc.assign(
        {"classification_flag": (["space"], np.ones_like(stm_masked_inc.space.values).astype(np.int8))}
    )

    # Compute NAD or NMAD if mem_persist is True
    # This only evaluate a very short task graph, since NAD or NMAD is already in memory
    if mem_persist:
        match method:
            case "nad":
                for key in [
                    "selection_nad",
                    "incremental_nmad",
                    "incremental_nad",
                    "recalibration_nad",
                    "recalibration_nmad",
                ]:
                    stm_masked_inc[key] = stm_masked[key].compute()
            case "nmad":
                for key in [
                    "selection_nmad",
                    "incremental_nmad",
                    "incremental_nad",
                    "recalibration_nad",
                    "recalibration_nmad",
                ]:
                    stm_masked_inc[key] = stm_masked[key].compute()

    return stm_masked_inc


def network_stm_selection(
    stm: xr.Dataset,
    min_dist: int | float,
    include_index: list[int] = None,
    sortby_var: str = "pnt_nmad",
    crs: int | str = "radar",
    x_var: str = "azimuth",
    y_var: str = "range",
    azimuth_spacing: float = None,
    range_spacing: float = None,
):
    """Select a Space-Time Matrix (STM) from a candidate STM for network processing.

    The selection is based on two criteria:
    1. A minimum distance between selected points.
    2. A sorting metric to select better points.

    The candidate STM will be sorted by the sorting metric.
    The selection will be performed iteratively, starting from the best point.
    In each iteration, the best point will be selected, and points within the minimum distance will be removed.
    The process will continue until no points are left in the candidate STM.

    Parameters
    ----------
    stm : xr.Dataset
        candidate Space-Time Matrix (STM).
    min_dist : int | float
        Minimum distance between selected points.
    include_index : list[int], optional
        Index of points in the candidate STM that must be included in the selection, by default None
    sortby_var : str, optional
        Sorting metric for selecting points, by default "pnt_nmad"
    crs : int | str, optional
        EPSG code of Coordinate Reference System of `x_var` and `y_var`, by default "radar".
        If crs is "radar", the distance will be calculated based on radar coordinates, and
        azimuth_spacing and range_spacing must be provided.
    x_var : str, optional
        Data variable name for x coordinate, by default "azimuth"
    y_var : str, optional
        Data variable name for y coordinate, by default "range"
    azimuth_spacing : float, optional
        Azimuth spacing, by default None. Required if crs is "radar".
    range_spacing : float, optional
        Range spacing, by default None. Required if crs is "radar".

    Returns
    -------
    xr.Dataset
        Selected network Space-Time Matrix (STM).

    Raises
    ------
    ValueError
        Raised when `azimuth_spacing` or `range_spacing` is not provided for radar coordinates.
    NotImplementedError
        Raised when an unsupported Coordinate Reference System is provided.
    """
    match crs:
        case "radar":
            if (azimuth_spacing is None) or (range_spacing is None):
                raise ValueError("Azimuth and range spacing must be provided for radar coordinates.")
        case _:
            raise NotImplementedError

    # Get coordinates and sorting metric, load them into memory
    stm_select = None
    stm_remain = stm[[x_var, y_var, sortby_var]].compute()

    # Select the include_index if provided
    if include_index is not None:
        stm_select = stm_remain.isel(space=include_index)

        # Remove points within min_dist of the included points
        coords_include = np.column_stack(
            [stm_select["azimuth"].values * azimuth_spacing, stm_select["range"].values * range_spacing]
        )
        coords_remain = np.column_stack(
            [stm_remain["azimuth"].values * azimuth_spacing, stm_remain["range"].values * range_spacing]
        )
        idx_drop = _idx_within_distance(coords_include, coords_remain, min_dist)
        if idx_drop is not None:
            stm_remain = stm_remain.where(~(stm_remain["space"].isin(idx_drop)), drop=True)

    # Reorder the remaining points by the sorting metric
    stm_remain = stm_remain.sortby(sortby_var)

    # Build a list of the index of selected points
    if stm_select is None:
        space_idx_sel = []
    else:
        space_idx_sel = stm_select["space"].values.tolist()

    while stm_remain.sizes["space"] > 0:
        # Select one point with best sorting metric
        stm_now = stm_remain.isel(space=0)

        # Append the selected point index
        space_idx_sel.append(stm_now["space"].values.tolist())

        # Remove the selected point from the remaining points
        stm_remain = stm_remain.isel(space=slice(1, None)).copy()

        # Remove points in stm_remain within min_dist of stm_now
        coords_remain = np.column_stack(
            [stm_remain["azimuth"].values * azimuth_spacing, stm_remain["range"].values * range_spacing]
        )
        coords_stmnow = np.column_stack(
            [stm_now["azimuth"].values * azimuth_spacing, stm_now["range"].values * range_spacing]
        )
        idx_drop = _idx_within_distance(coords_stmnow, coords_remain, min_dist)
        if idx_drop is not None:
            stm_drop = stm_remain.isel(space=idx_drop)
            stm_remain = stm_remain.where(~(stm_remain["space"].isin(stm_drop["space"])), drop=True)

    # Get the selected points by space index from the original stm
    stm_out = stm.sel(space=space_idx_sel)

    return stm_out


def _nad_block(amp: xr.DataArray) -> xr.DataArray:
    """Compute Normalized Amplitude Dispersion (NAD) for a block of amplitude data.

    Parameters
    ----------
    amp : xr.DataArray
        Amplitude data, with dimensions ("azimuth", "range", "time").
        This can be extracted from an SLC xr.Dataset.

    Returns
    -------
    xr.DataArray
        Normalized Amplitude Dispersion (NAD) data, with dimensions ("azimuth", "range").
    """
    # Compute amplitude dispersion
    # By defalut, the mean and std function from Xarray will skip NaN values
    # However, if there is NaN value in time series, we want to discard the pixel
    # Therefore, we set skipna=False
    # Adding epsilon to avoid zero division
    nad_da = amp.std(dim="time", skipna=False) / (amp.mean(dim="time", skipna=False) + np.finfo(amp.dtype).eps)

    return nad_da


def _nmad_block(amp: xr.DataArray) -> xr.DataArray:
    """Compute Normalized Median Absolute Deviation(NMAD) for a block of amplitude data.

    Parameters
    ----------
    amp : xr.DataArray
        Amplitude data, with dimensions ("azimuth", "range", "time").
        This can be extracted from an SLC xr.Dataset.

    Returns
    -------
    xr.DataArray
        Normalized Median Absolute Dispersion (NMAD) data, with dimensions ("azimuth", "range").
    """
    # Compoute NMAD
    median_amplitude = amp.median(dim="time", skipna=False)
    mad = (np.abs(amp - median_amplitude)).median(dim="time")  # Median Absolute Deviation
    nmad = mad / (median_amplitude + np.finfo(amp.dtype).eps)  # Normalized Median Absolute Deviation

    return nmad


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
    amp["partition_amp_std"] = (amp.dims, np.ones_like(data) * std)
    amp["partition_amp_mean"] = (amp.dims, np.ones_like(data) * mean)
    amp["partition_amp_mad"] = (amp.dims, np.ones_like(data) * mad)
    amp["partition_amp_median"] = (amp.dims, np.ones_like(data) * median)
    return amp


def _idx_within_distance(coords_ref, coords_others, min_dist):
    """Get the index of points in coords_others that are within min_dist of coords_ref.

    Parameters
    ----------
    coords_ref : np.ndarray
        Coordinates of reference points. Shape (n, 2).
    coords_others : np.ndarray
        Coordinates of other points. Shape (m, 2).
    min_dist : int, float
        distance threshold.

    Returns
    -------
    np.ndarray
        Index of points in coords_others that are within `min_dist` of `coords_ref`.
    """
    kd_ref = KDTree(coords_ref)
    kd_others = KDTree(coords_others)
    sdm = kd_ref.sparse_distance_matrix(kd_others, min_dist)
    if len(sdm) > 0:
        idx = np.array(list(sdm.keys()))[:, 1]
        return idx
    else:
        return None
