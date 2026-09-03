# Adding contextual information (e.g. temperature)

DePSI can enrich an STM with attributes sourced from outside the SAR processing chain, such as
temperature. `depsi.contextual_information` only handles attaching already-fetched data to an STM --
it has no knowledge of where that data came from, only that it's available as a dict or a file in the
shape it expects (see below). Fetching the raw data itself (e.g. from ERA5 or KNMI) is a separate concern,
handled outside this repo, upstream of whatever calls these functions.

The [demo notebook](../notebooks/demo_get_and_add_contextualinformation.ipynb) still includes an example
of fetching real data via the [`get_weatherdata`](https://github.com/LinaHagenah/Get_WeatherData) package,
purely for illustration -- `get_weatherdata` is not a dependency of `depsi` itself, only of the `demo`
dependency group used to run that notebook.

## Adding attributes directly

There are three functions in `depsi.contextual_information`, one per kind of data source. Each adds
one or more attributes to an STM by name, looking up where each attribute's data comes from and calling
the matching [`stmtools`](https://github.com/TUDelftGeodesy/stmtools) enrichment method directly -- no
enrichment logic of its own.

### One station per attribute

`add_contextual_information_from_point` -- one location's data is broadcast to every point in the STM.
This is the usual case for something like temperature, where one weather station stands in for the whole
area:

```python
from depsi.contextual_information import add_contextual_information_from_point

stm_out = add_contextual_information_from_point(
    stm,
    selected_attributes=["temperature"],
    attribute_sources={
        # attribute name -> (lat, lon, data)
        # data is either a dict shaped {"YYYY-MM-DD": {"temperature": value}, ...},
        # or a path to a CSV with a "Date" index column and a "temperature" column
        "temperature": (52.3676, 4.9041, "path/to/temperature.csv"),
    },
)
```

### Several stations per attribute

`add_contextual_information_from_grid` -- each attribute comes from a list of locations instead of one.
`stmtools` does real nearest-neighbor spatial matching: each STM point gets whichever location is closest
to it. Locations with different date coverage are reindexed onto their shared union of dates first
(missing dates become `NaN`):

```python
from depsi.contextual_information import add_contextual_information_from_grid

stm_out = add_contextual_information_from_grid(
    stm,
    selected_attributes=["temperature"],
    attribute_sources={
        # attribute name -> list of (lat, lon, data), same shape per location as the point case
        "temperature": [
            (52.30, 4.85, "path/to/station_a.csv"),
            (52.45, 5.00, "path/to/station_b.csv"),
        ],
    },
)
```

### A polygon, for static attributes

`add_contextual_information_from_polygon` -- for attributes that don't vary in time (e.g. land use, soil
type), matched via a point-in-polygon join rather than nearest-neighbor. The added variable has dims
`space` only, no `time`. This one also requires `stm` to be dask-backed (chunked) -- `enrich_from_polygon`
reads `stm.chunksizes["space"]` internally, which raises on a plain numpy-backed STM:

```python
from depsi.contextual_information import add_contextual_information_from_polygon

stm_out = add_contextual_information_from_polygon(
    stm.chunk({"space": -1, "time": -1}),
    selected_attributes=["land_use"],
    attribute_sources={
        # attribute name -> (polygon, field)
        # polygon is a geopandas.GeoDataFrame, or a path geopandas.read_file can open
        "land_use": ("path/to/land_use.geojson", "LU_CODE"),
    },
)
```

If the field name in the polygon (`"LU_CODE"` above) differs from the attribute name, the result is
renamed to the attribute name after enrichment.

For all three: only the attributes named in `selected_attributes` are added, even if `attribute_sources`
has more available, so the same sources dict can be reused across calls while controlling what actually
gets added. And the STM must not already have a variable with that attribute's name -- `stmtools` refuses
to overwrite an existing field.

There's no built-in dispatcher for combining multiple kinds in one call -- if you have a mix of
point/grid/polygon-sourced attributes, just call each function for the attributes it applies to, same as
calling any of them alone. The [demo notebook](../notebooks/demo_get_and_add_contextualinformation.ipynb)
has a worked example of combining all three in one script.

## Viewing geometry: incidence angle and crossrange

`incidence_angle` and `crossrange` aren't contextual data sourced externally -- they're viewing-geometry
values computed from the orbit configuration, so they live in a separate module,
`depsi.viewing_geometry`, with two functions instead of one:

```python
from depsi.viewing_geometry import add_local_viewing_geometry, add_cross_range

# adds local_incidence_angle and local_azimuth_angle, from the orbit
stm_out = add_local_viewing_geometry(
    stm,
    orbit_config_file="config/drama/S1_XTI.cfg",
    orbit_res=0.01,
    orbit_mode="IWS",
    orbit="s1_dsc_t037",  # must match the track the STM's SLC stack was processed from
)

# adds sd_cr2ph, computed from sd_h2ph (already on the STM, e.g. from single-differencing) and
# local_incidence_angle (just added above)
stm_out = add_cross_range(stm_out)
```

`add_cross_range` doesn't compute incidence angle itself -- call `add_local_viewing_geometry` first. If
either `sd_h2ph` or `local_incidence_angle` is missing from the STM, it raises a `ValueError` naming which
one and, for `local_incidence_angle`, that you need to call `add_local_viewing_geometry` first.
