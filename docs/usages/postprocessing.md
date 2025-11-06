# Post-processing the Space-Time Matrix

After the time series have been estimated, one may only want a subset of all processed points for e.g. visualisation purposes. 
This module provides the functionality for this.

## Filtering

```python
from depsi.postprocessing import stm_point_filter

dict_filtering = {"key1": {"bounds": (low, high), "return_rejected": True | False}, "key2": ...}
rejected_dict = {}

for key in dict_filtering.keys():
    if dict_filtering[key]["return_rejected"]:
        stm, rejected_dict[key] = stm_point_filter(stm, key, vmin=dict_filtering[key]["bounds"][0], 
                                                   vmax=dict_filtering[key]["bounds"][1], return_removed=True)
    else:
        stm = stm_point_filter(stm, key, vmin=dict_filtering[key]["bounds"][0], 
                                                   vmax=dict_filtering[key]["bounds"][1], return_removed=False)
```

After this operation, the variable `stm` will only contain the points meeting all the filters.