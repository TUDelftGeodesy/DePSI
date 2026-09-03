"""Functions for adding contextual information to a space-time matrix (STM)."""

from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import stmtools  # noqa: F401  (registers the `.stm` accessor used below)
import xarray as xr


def _load_station_series(data: dict | str | Path, attribute_name: str) -> pd.Series:
    """Load one station's data for `attribute_name` into a date-indexed pandas Series.

    Parameters
    ----------
    data : dict | str | Path
        Either a dict shaped `{"YYYY-MM-DD": {attribute_name: value, ...}, ...}`, or a path to a CSV file
        with a "Date" index column and a column named `attribute_name`.
    attribute_name : str
        The key (dict) or column (CSV) to extract.

    Returns
    -------
    pd.Series
        The station's values for `attribute_name`, indexed by date, sorted by date.
    """
    if isinstance(data, dict):
        try:
            series = pd.Series({pd.Timestamp(date): values[attribute_name] for date, values in data.items()})
        except KeyError as e:
            raise ValueError(
                f"Key '{attribute_name}' was not found for at least one date in the provided data. "
                f"Available keys for that date: {e}"
            ) from e
    else:
        df = pd.read_csv(data, index_col="Date", parse_dates=True)
        if attribute_name not in df.columns:
            raise ValueError(f"Column '{attribute_name}' not found in {data}. Available columns: {list(df.columns)}")
        series = df[attribute_name]

    return series.sort_index()


def add_contextual_information_from_point(
    stm: xr.Dataset,
    selected_attributes: list[str],
    attribute_sources: dict[str, tuple[float, float, dict | str | Path]],
) -> xr.Dataset:
    """Add one or more contextual attributes to an STM from a single station each, via stmtools.

    For each name in `selected_attributes`, looks up where its data comes from in
    `attribute_sources`, reshapes that data into the small `xr.Dataset` `stmtools`' `enrich_from_dataset`
    needs, and calls it directly. Since there's only one location per attribute, its data is broadcast to
    every point in the STM.

    Parameters
    ----------
    stm : xr.Dataset
        Space-time matrix to add the attribute(s) to.
    selected_attributes : list[str]
        Which attributes to add, e.g. `["temperature", "pressure"]`. Each name must be a key in
        `attribute_sources`. Attributes are added in this order.
    attribute_sources : dict[str, tuple[float, float, dict | str | Path]]
        Maps each attribute name to where its data comes from: `lat`/`lon` (WGS84) for the
        station/location it was collected at, and `data`, either a dict shaped
        `{"YYYY-MM-DD": {attribute_name: value, ...}, ...}`, or a path to a CSV file with a "Date" index
        column and a column named `attribute_name`.

    Returns
    -------
    xr.Dataset
        The STM with each of `selected_attributes` added, dims `space`, `time`.

    Examples
    --------
    Add one station time series and broadcast it to all points:

    >>> import numpy as np
    >>> import xarray as xr
    >>> from depsi.contextual_information import add_contextual_information_from_point
    >>>
    >>> time = np.array(["2026-07-01", "2026-07-13", "2026-07-25"], dtype="datetime64[ns]")
    >>> stm = xr.Dataset(
    ...     data_vars={"sd_h2ph": (["space", "time"], np.ones((2, 3)))},
    ...     coords={"lat": ("space", [52.365, 52.366]), "lon": ("space", [4.895, 4.896]), "time": time},
    ... )
    >>> stm_out = add_contextual_information_from_point(
    ...     stm,
    ...     selected_attributes=["temperature"],
    ...     attribute_sources={
    ...         "temperature": (
    ...             52.3676,
    ...             4.9041,
    ...             {
    ...                 "2026-07-01": {"temperature": 18.0},
    ...                 "2026-07-13": {"temperature": 20.0},
    ...                 "2026-07-25": {"temperature": 19.0},
    ...             },
    ...         )
    ...     },
    ... )
    >>> stm_out["temperature"].shape
    (2, 3)
    >>> stm_out["temperature"].values[0].tolist()
    [18.0, 20.0, 19.0]
    """
    stm = stm.copy()
    for attribute_name in selected_attributes:
        if attribute_name not in attribute_sources:
            raise ValueError(
                f"No data source given for '{attribute_name}'. Available: {list(attribute_sources.keys())}"
            )
        source = attribute_sources[attribute_name]
        if not (isinstance(source, tuple) and len(source) == 3):
            raise ValueError(
                f"attribute_sources['{attribute_name}'] must be a (lat, lon, data) tuple, got {source!r}."
            )
        lat, lon, data = source
        series = _load_station_series(data, attribute_name)

        enrichment_ds = xr.Dataset(
            data_vars={attribute_name: (["space", "time"], series.to_numpy()[np.newaxis, :])},
            coords={
                # "space" is deliberately left without its own coordinate: stmtools' enrich_from_dataset
                # only matches on lon/lat (real spatial distance) when "space" isn't itself a coordinate
                # -- if it is, it matches on the raw integer index instead, silently ignoring lon/lat.
                "lat": ("space", [lat]),
                "lon": ("space", [lon]),
                "time": series.index.to_numpy().astype(stm["time"].values.dtype),
            },
        )

        stm = stm.stm.enrich_from_dataset(enrichment_ds, attribute_name)

    return stm


def add_contextual_information_from_grid(
    stm: xr.Dataset,
    selected_attributes: list[str],
    attribute_sources: dict[str, list[tuple[float, float, dict | str | Path]]],
) -> xr.Dataset:
    """Add one or more contextual attributes to an STM from several stations each, via stmtools.

    Like `add_contextual_information_from_point`, but each attribute comes from a list of locations
    instead of one. `stmtools`' `enrich_from_dataset` then does real nearest-neighbor spatial matching.

    All locations for one attribute are reindexed onto the union of their dates before being combined
    (missing dates become NaN), since the small enrichment dataset needs one shared `time` axis across all
    locations.

    Parameters
    ----------
    stm : xr.Dataset
        Space-time matrix to add the attribute(s) to.
    selected_attributes : list[str]
        Which attributes to add. Each name must be a key in `attribute_sources`.
    attribute_sources : dict[str, list[tuple[float, float, dict | str | Path]]]
        Maps each attribute name to a list of locations it's available at -- same `(lat, lon, data)` shape
        per location as `add_contextual_information_from_point`'s `attribute_sources` values.

    Returns
    -------
    xr.Dataset
        The STM with each of `selected_attributes` added, dims `space`, `time`.

    Examples
    --------
    Match each STM point to its nearest station and align differing station dates:

    >>> import numpy as np
    >>> import xarray as xr
    >>> from depsi.contextual_information import add_contextual_information_from_grid
    >>>
    >>> time = np.array(["2026-07-01", "2026-07-13", "2026-07-25"], dtype="datetime64[ns]")
    >>> stm = xr.Dataset(
    ...     data_vars={"sd_h2ph": (["space", "time"], np.ones((2, 3)))},
    ...     coords={"lat": ("space", [52.40, 52.30]), "lon": ("space", [4.90, 4.90]), "time": time},
    ... )
    >>> out = add_contextual_information_from_grid(
    ...     stm,
    ...     selected_attributes=["temperature"],
    ...     attribute_sources={
    ...         "temperature": [
    ...             (
    ...                 52.40,
    ...                 4.90,
    ...                 {
    ...                     "2026-07-01": {"temperature": 14.0},
    ...                     "2026-07-13": {"temperature": 15.0},
    ...                 },
    ...             ),
    ...             (
    ...                 52.30,
    ...                 4.90,
    ...                 {
    ...                     "2026-07-01": {"temperature": 22.0},
    ...                     "2026-07-13": {"temperature": 23.0},
    ...                     "2026-07-25": {"temperature": 24.0},
    ...                 },
    ...             ),
    ...         ]
    ...     },
    ... )
    >>> out["temperature"].values[0].tolist()
    [14.0, 15.0, nan]
    >>> out["temperature"].values[1].tolist()
    [22.0, 23.0, 24.0]
    """
    stm = stm.copy()
    for attribute_name in selected_attributes:
        if attribute_name not in attribute_sources:
            raise ValueError(
                f"No data source given for '{attribute_name}'. Available: {list(attribute_sources.keys())}"
            )
        locations = attribute_sources[attribute_name]
        if not (isinstance(locations, list) and locations):
            raise ValueError(
                f"attribute_sources['{attribute_name}'] must be a non-empty list of (lat, lon, data) "
                f"locations, got {locations!r}."
            )
        for location in locations:
            if not (isinstance(location, tuple) and len(location) == 3):
                raise ValueError(
                    f"Each location in attribute_sources['{attribute_name}'] must be a (lat, lon, data) "
                    f"tuple, got {location!r}."
                )
        series_list = [_load_station_series(data, attribute_name) for _, _, data in locations]

        shared_index = series_list[0].index
        for series in series_list[1:]:
            shared_index = shared_index.union(series.index)
        values = np.stack([series.reindex(shared_index).to_numpy() for series in series_list])

        enrichment_ds = xr.Dataset(
            data_vars={attribute_name: (["space", "time"], values)},
            coords={
                "lat": ("space", [lat for lat, _, _ in locations]),
                "lon": ("space", [lon for _, lon, _ in locations]),
                "time": shared_index.to_numpy().astype(stm["time"].values.dtype),
            },
        )

        stm = stm.stm.enrich_from_dataset(enrichment_ds, attribute_name)

    return stm


def add_contextual_information_from_polygon(
    stm: xr.Dataset,
    selected_attributes: list[str],
    attribute_sources: dict[str, tuple[str | Path | gpd.GeoDataFrame, str]],
) -> xr.Dataset:
    """Add one or more static, polygon-sourced contextual attributes to an STM, via stmtools.

    For each name in `selected_attributes`, looks up its polygon and field name in `attribute_sources`,
    and calls `stm.stm.enrich_from_polygon(polygon, [field])` directly. Unlike
    `add_contextual_information_from_point`/`_from_grid`, the added variable has dims `space`
    only, not `time`: this is for static attributes (e.g. land use, soil type), not time-varying ones.

    If `field` differs from `attribute_name`, the resulting variable is renamed to `attribute_name` after
    enrichment, so the STM ends up with a variable named after the attribute regardless of what it's
    called in the source polygon.

    Unlike `_from_point`/`_from_grid`, `enrich_from_polygon` requires `stm` to be dask-backed (chunked) --
    it reads `stm.chunksizes["space"]` internally, which raises a `KeyError` on a plain numpy-backed STM.
    Chunk first (e.g. `stm.chunk({"space": -1, "time": -1})`) if needed.

    Parameters
    ----------
    stm : xr.Dataset
        Space-time matrix to add the attribute(s) to. Must be dask-backed (chunked).
    selected_attributes : list[str]
        Which attributes to add. Each name must be a key in `attribute_sources`.
    attribute_sources : dict[str, tuple[str | Path | gpd.GeoDataFrame, str]]
        Maps each attribute name to the (multi-)polygon its value comes from (a `geopandas.GeoDataFrame`,
        or a path to a file `geopandas.read_file` can open) and the field name within it to use.

    Returns
    -------
    xr.Dataset
        The STM with each of `selected_attributes` added, dim `space`.
    """
    stm = stm.copy()
    for attribute_name in selected_attributes:
        if attribute_name not in attribute_sources:
            raise ValueError(
                f"No data source given for '{attribute_name}'. Available: {list(attribute_sources.keys())}"
            )
        source = attribute_sources[attribute_name]
        if not (isinstance(source, tuple) and len(source) == 2):
            raise ValueError(
                f"attribute_sources['{attribute_name}'] must be a (polygon, field) tuple, got {source!r}."
            )
        polygon, field = source
        stm = stm.stm.enrich_from_polygon(polygon, [field])
        if field != attribute_name:
            stm = stm.rename({field: attribute_name})

    return stm
