---
name: depsi-docs
description: 'Write DePSI documentation per module and function. Use when creating or updating API docs and user-facing usage docs from docstrings, tests, and examples.'
argument-hint: '[module_name] [function_name] [optional_example]' 
user-invocable: true
---

## Purpose

Create or update DePSI documentation for one function at a time, including:

1. API documentation derived from the target function docstring.
2. User documentation demonstrating practical usage.

## Inputs

- `module_name`: DePSI module name (for example: `network`, `io`, `densification`).
- `function_name`: Function name inside `depsi/<module_name>.py`.
- `optiona_example` (optional): Path to a Jupyter notebook (`.ipynb`) or Python script (`.py`) that demonstrates calling the function.

## Repository Documentation Conventions

- Keep documentation aligned with `mkdocs.yml`.
- API docs live under `docs/api/`.
- User docs live under `docs/usages/`.
- Preserve existing Material/MkDocs style and page structure.
- If a new page is created, update `mkdocs.yml` navigation under the correct section (`API Reference` or `Usage`).

## API Documentation Template

When generating API documentation for a function, use the following template as default API documentation for a function:

```md
::: depsi.<module_name>.<function_name>
    options:
      show_root_heading: true
      show_source: true
```

## Required Workflow

1. Resolve and validate target function:
   - Locate `depsi/<module_name>.py`.
   - Confirm `function_name` exists and has a docstring.
   - If missing, stop and return a clear error.

2. Add or update API documentation first:
   - Target API page: `docs/api/<module_name>.md` (or the file mapped for that module in `mkdocs.yml`).
   - Use the API Documentation Template section above as a starting point.
   - Use the function docstring as the canonical source for parameters, returns, and behavior.

3. Gather usage evidence for user documentation:
   - Always inspect the function docstring.
   - Search for unit tests that call the function (prefer `tests/test_<module_name>.py`, then broader `tests/` if needed).
   - If `optiona_example` is provided, inspect and use it.

4. Decide whether user documentation can be added:
   - If **no unit test is found** and **no `optiona_example` is provided**, emit a warning and **stop before adding/updating user documentation**.
   - API documentation from step 2 must be completed before proceeding to user documentation.

5. Add or update user documentation (only when allowed by step 4):
   - Update the relevant file under `docs/usages/` based on the module/topic mapping in `mkdocs.yml`.
   - Include:
     - Purpose and context of the function.
     - A runnable usage snippet.
     - Input/output expectations and key caveats.
   - The usage content must be grounded in:
     1) docstring, 2) unit test behavior (if exists), and 3) provided example (if provided).
   - Do not invent behavior not supported by those sources.

## Output Contract

- Exactly one function documented per invocation.
- API docs updated first.
- User docs updated only when evidence requirements are satisfied.
- Any warning/error must be explicit and actionable (no silent fallback).
