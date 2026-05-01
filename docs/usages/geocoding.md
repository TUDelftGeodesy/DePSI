# Geocoding

DePSI provides geocoding functionality to convert between geographic coordinates (latitude, longitude) and radar coordinates (range, azimuth).

This covertion is based on the radar coordinates, as well as radar geometry related metadata, such as orbit state vectors, radar timing parameters, and ellipsoid model. The metadata can be read from the SAR processing software. At present, DePSI geocoding module supports metadata read by `sarxarray` from `DORIS`:

```python
import sarxarray
from depsi.transformations import radar_to_latlonh, validate_geocoding_metadata

# Load Doris metadata with sarxarray
# The metadata is a dictionary with necessary fields for geocoding
metadata = sarxarray.read_metadata('metadata.res', driver="doris5")

# (Optional) validate metadata
metadata_validated = validate_geocoding_metadata(metadata)

# Apply conversion
latlonh = radar_to_latlonh(
    azimuth_coords=stm['azimuth'].values,
    range_coords=stm['range'].values,
    elevation=stm['pnt_height'].values,
    metadata=metadata,
)
```

As shown above, one can use `validate_geocoding_metadata` to check if the metadata contains all necessary fields for geocoding. The function will raise an error if any required field is missing.