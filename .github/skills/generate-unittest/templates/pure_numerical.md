# pure_numerical_template

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
