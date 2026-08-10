# validation_error_template

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
