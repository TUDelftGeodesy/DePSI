# xarray_dataset_behavior_template

When to use:
- Dataset or DataArray transforms, including STM and SLC operations.

Required imports:
```python
import numpy as np
import xarray as xr
import pytest
from depsi.<module_name> import <function_name>
```

Minimum scaffold:
```python
def _make_minimal_dataset():
    return xr.Dataset(
        data_vars={
            "<input_var>": (("time", "space"), np.array([[<v11>, <v12>], [<v21>, <v22>]], dtype=float)),
        },
        coords={
            "time": np.array([0, 1]),
            "space": np.array([0, 1]),
        },
        attrs={"<required_attr>": <required_value>},
    )


def test_<function_name>_dataset_behavior():
    ds = _make_minimal_dataset()

    out = <function_name>(ds, <args>)

    # Fast guards
    assert set(["<expected_var_1>", "<expected_var_2>"]).issubset(set(out.data_vars))
    assert out["<expected_var_1>"].shape == (<expected_dim_1>, <expected_dim_2>)

    # Semantic invariant
    assert <invariant_expression>
```

Minimum dataset note:
- Reuse existing fixtures when available, but keep the test self-contained and minimal.
- The 2x2 time/space shape is an intentional default for minimal tests.
- Parameterize the helper (for example with n_time and n_space) only when the target function behavior depends on dimension length, neighborhood/window size, broadcasting, or rank constraints.

Assertion checklist:
- Output key and shape checks as fast guards.
- At least one semantic invariant tied to the function contract.
- Attr and coordinate preservation when required.
