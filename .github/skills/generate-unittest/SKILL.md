---
name: generate-unittest
description: 'Generate new pytest unit tests for DePSI functions using curated template archetypes. Use when adding new functions or improving test coverage.'
argument-hint: '[module_path] [function_name] [optional_preferred_archetype]'
user-invocable: true
---

## Purpose

Generate new pytest unit tests for DePSI functions using curated template archetypes.

## Inputs

- Target module path.
- Target function name.
- Optional preferred archetype if explicitly requested.

## Required References

- unit_test_templates.md

## Repository Conventions

- Source modules are under: `depsi/`
- Unit tests are under: `tests/`
- When asked for tests for module `depsi/<name>.py`, prefer adding/updating `tests/test_<name>.py`.
- Reuse existing pytest style and fixtures from `tests/` before introducing new patterns.

## Primary Archetype Selection

Choose exactly one primary archetype per generated test case:

- `xarray_dataset_behavior_template`: dataset or DataArray transforms (STM/SLC) with contract checks on keys, shapes, attrs, and coordinates.
- `pure_numerical_template`: deterministic scalar or small-array transforms with explicit expected outputs.
- `synthetic_recovery_template`: estimation, inversion, or round-trip behavior that should recover known latent truth.
- `validation_error_template`: explicit rejection behavior for invalid attrs, keys, options, dimensions, or metadata.
- `fixture_io_template`: parsers/readers where realistic files are required to validate contract-level fields. When this archetype is selected, always ask the user to provide a fixture file path.

## Excluded Patterns

Avoid these as default generated-test patterns:

- fixture-literal marker checks (for example, asserting only that a raw marker string exists in a fixture)
- opaque large exact-value regressions copied without clear behavioral invariants

## Workflow

1. Identify the target function contract from code and nearby tests.
2. Select exactly one primary archetype from this skill for the first generated test case. If multiple archetypes are plausible, choose the most contract-relevant one first. If none are plausible, raise an error and stop the process.
3. Instantiate the matching template from unit_test_templates.md.
4. Generate unit test code according to the template, using the following rules:
  1. Keep test inputs minimal and explicit.
  2. Include key and shape assertions as fast guards when testing xarray dataset behavior.
  3. Include one semantic invariant assertion tied to function behavior.
  4. Add one validation/error-path test when the API exposes meaningful validation.
  5. Avoid excluded patterns listed above.
5. If multiple archetypes are plausible in step 2, generate additional test cases by repeating steps 2-4, selecting exactly one archetype per additional test case, until all plausible archetypes are covered.


## Output Contract

Produce test code that:
- uses pytest style already present in the repository,
- can be placed in an existing tests/test_*.py file,
- includes a clear Arrange/Act/Assert structure (comments optional),
- uses explicit tolerances for floating-point comparisons,
- does not over-assert incidental implementation details.

## Quality Checklist

- Primary archetype chosen and stated.
- Inputs are minimal but valid.
- At least one contract-level assertion exists.
- Semantic invariant assertion included.
- Key/shape fast guards included when useful for dataset behavior.
- Validation path included when meaningful.
- No excluded-pattern behavior.
