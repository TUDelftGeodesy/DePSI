# synthetic_recovery_template

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
