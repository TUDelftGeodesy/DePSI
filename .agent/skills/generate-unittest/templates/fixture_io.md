# fixture_io_template

When to use:
- File readers, metadata parsers, and IO adapters requiring realistic fixture inputs.

When to avoid:
- The fixture file is absent for the target function.

Template note:
- When this template is used, always ask the user to provide a fixture file path and name.
- Never generate a fixture file automatically.

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
