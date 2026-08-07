# Unit Test Template Catalog For DePSI

This catalog provides concrete templates for generating new unit tests in DePSI.
Use this as a template catalog together with SKILL.md.

## Global Rules

- Choose exactly one primary template for each generated test case.
- Keep the input scaffold minimal and explicit.
- Prefer semantic assertions first.
- Key and shape assertions are valid fast guards in DePSI, especially for dataset-transform functions.
- Add one error-path test when the API has meaningful validation behavior.
- Avoid fixture-literal marker checks and opaque large exact-value regressions without clear invariants.

## 1. xarray_dataset_behavior_template

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


## 2. pure_numerical_template

When to use:
- Deterministic scalar or small-array functions with known expected outputs.

Required imports:
```python
import numpy as np
import pytest
from depsi.<module_name> import <function_name>
```

Minimum scaffold:
```python
def test_<function_name>_basic_case():
    # Arrange
    x = np.array([<values>])

    # Act
    out = <function_name>(x, <args>)

    # Assert
    expected = np.array([<expected_values>])
    np.testing.assert_allclose(out, expected, rtol=<rtol>, atol=<atol>)
```

Optional validation path:
```python
def test_<function_name>_invalid_input():
    with pytest.raises(<ExceptionType>):
        <function_name>(<invalid_args>)
```

Assertion checklist:
- Exact equality for integers and categorical outputs.
- allclose for floats with explicit tolerance.
- One edge case near a boundary if meaningful.

## 3. synthetic_recovery_template

When to use:
- Estimation, inversion, or round-trip style behavior with known latent truth.

Required imports:
```python
import numpy as np
import pytest
from depsi.<module_name> import <function_name>
```

Minimum scaffold:
```python
def _make_synthetic_problem():
    true_params = np.array([<p1>, <p2>], dtype=float)
    x = np.array([<x_values>], dtype=float)
    y = <forward_model_expression_using_true_params_and_x>
    return x, y, true_params


def test_<function_name>_recovers_truth():
    x, y, true_params = _make_synthetic_problem()

    est = <function_name>(x, y, <args>)

    np.testing.assert_allclose(est, true_params, rtol=<rtol>, atol=<atol>)
```

Assertion checklist:
- Recovery accuracy with explicit tolerance.
- Stable behavior for minimal but valid synthetic setup.

## 4. validation_error_template

When to use:
- APIs with explicit validation rules for attrs, keys, options, dimensions, or metadata.

Required imports:
```python
import pytest
from depsi.<module_name> import <function_name>
```

Minimum scaffold:
```python
def test_<function_name>_rejects_invalid_<case_name>():
    bad_input = <construct_small_invalid_input>

    with pytest.raises(<ExceptionType>):
        <function_name>(bad_input, <args>)
```

Optional message check:
```python
def test_<function_name>_error_message_contract():
    bad_input = <construct_small_invalid_input>

    with pytest.raises(<ExceptionType>, match="<stable_error_fragment>"):
        <function_name>(bad_input, <args>)
```

Assertion checklist:
- Assert exception type always.
- Assert message only when the text is part of external API behavior.
- Keep one invalid condition per test.


## 5. fixture_io_template

When to use:
- File readers, metadata parsers, and IO adapters requiring realistic fixture inputs.

Required imports:
```python
from pathlib import Path
import pytest
from depsi.<module_name> import <function_name>
```

Minimum scaffold:
```python
def test_<function_name>_fixture_contract():
    fixture_path = Path("tests/data/<fixture_name>")

    out = <function_name>(fixture_path, <args>)

    # Contract-level checks
    assert <contract_field_1_check>
    assert <contract_field_2_check>
```

Fixture policy note:
- Reuse an existing file under `tests/data` when one already matches the parser or reader contract.
- Add a new minimal fixture under `tests/data` only when realistic file content is required and no suitable fixture already exists.
- Do not add a fixture when the behavior can be tested more clearly with synthetic in-memory data.

Optional parser-normalization path:
```python
def test_<function_name>_normalizes_fields():
    fixture_path = Path("tests/data/<fixture_name>")

    out = <function_name>(fixture_path, <args>)

    assert <normalized_rows_or_columns_check>
```

Assertion checklist:
- Assert contract fields, not incidental fixture details.
- Prefer parsed semantic properties over raw marker-string checks.
- Keep fixture scope narrow and stable.

## Quick Selection Guide

- Numerical scalar/array transform: pure_numerical_template
- STM or SLC dataset transform: xarray_dataset_behavior_template
- Inversion or latent-parameter recovery: synthetic_recovery_template
- Input/metadata/options rejection: validation_error_template
- Reader or parser behavior with sample files: fixture_io_template
