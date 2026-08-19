# Practical info for developers

## Installation guide

The Python implementation of DePSI is under development. At present you can only install it from the GitHub repository.

We strongly recommend using `uv` to handle the Python environment during development process. Please refer to the [uv documentation](https://docs.astral.sh/uv/) for installation instructions of `uv` 

Before you start, make sure you have access to the correct DePSI repository. You can contribute to the public DePSI repository by forking it to your own GitHub account. If you are a member of the TUDelftGeodesy organization, you can also contribute to the group DePSI repository.

Clone this GitHub repository, then `cd` into the cloned repository.

```bash
cd DePSI
```

Then you can create a development environment simply by synchronizing:

```bash
uv sync --all-groups
```

To activate the environment created by `uv`, run under `DePSI` directory:

```bash
source .venv/bin/activate
```

To add a new package during development, you can use `uv add <package-name>` command. 

```bash
uv add <package-name>
```

In the end, install the pre-commit hooks, which will run the checks before each commit:
```bash
pre-commit install
```

## Linting and formatting

We use `ruff` for linting and formatting. If the pre-commit hooks are installed, the checks will be run automatically before each commit.

To manually run the checks, use the following command in the root directory of the repository:

```bash
ruff check .
```

## Testing

We use `pytest` for testing. All tests are located in the `tests` directory.

To run the tests, use the following command in the root directory of the repository:

```bash
pytest tests
```

The [GitHub Actions](https://github.com/TUDelftGeodesy/DePSI/blob/main/.github/workflows/build.yml) will run the tests automatically for each push and pull-request
on the `main` branch.

### Agent skill for unit-test generation

This repository includes an agent skill for generating DePSI unit tests:

- Skill definition: `.github/skills/generate-unittest/SKILL.md`
- Archetype templates: `.github/skills/generate-unittest/catalog_archetypes.md`

Use this skill in Copilot Chat by explicitly naming the module path and function name, with an optional preferred archetype defined in the Archetype templates file.

Example prompts:

Write unit tests for `form_network` function in `depsi/network.py`.

```text
/generate-unittest depsi/network.py form_network
```

Write unit tests for `form_network` function in `depsi/network.py` and prefer the `xarray_dataset_behavior_template` archetype.

```text
/generate-unittest depsi/network.py form_network xarray_dataset_behavior_template
```

The currently supported archetypes in `catalog_archetypes.md` are:

- `xarray_dataset_behavior_template`: dataset or DataArray transforms (STM/SLC) with contract checks on keys, shapes, attrs, and coordinates. Example: `tests/test_network::TestNetworkFormation::test_form_network_simulated_grid` for `depsi.network.form_network`.

- `pure_numerical_template`: deterministic scalar or small-array functions with known expected outputs. Example: `tests/test_transformations::test_seconds_of_day` for `depsi.transformations.seconds_of_day`.

- `synthetic_recovery_template`: estimation, inversion, or round-trip behavior that should recover known latent truth. Example: `tests/test_arc_estimation::test_periodogram` for `depsi.arc_estimation.periodogram`.

- `validation_error_template`: explicit rejection behavior for invalid attrs, keys, options, dimensions, or metadata. Example: `tests/test_model_estimation::test_estimate_model_params_invalid_model` for `depsi.model_estimation.estimate_model_params`.

- `fixture_io_template`: file readers, metadata parsers, and IO adapters requiring realistic fixture inputs. Example: `tests/test_io::test_read_rcs_csv` for `depsi.io.read_rcs_csv`.

Expected behavior:

- Target function is read from `depsi/<module>.py`.
- Relevant tests are added in `tests/test_<module>.py`.
- Per generated test case, exactly one primary archetype is selected from the archetype templates file.
- Generated tests follow existing pytest style and include meaningful contract-level assertions.

Tip:

- It is preferred to provide the target module path as a context when using the skill. 
- Generate unit tests for one function at a time to keep review and debugging focused.
- When using the "fixture_io_template" archetype, a fixture file path must be provided. The skill will not generate a fixture file automatically.

## Documentation

We use `mkdocs` for documentation. 

To check the documentation at local, use the following command in the root directory of the repository:

```bash
mkdocs serve
```

This will build and render the documentation at a local server. Follow the link provided in the terminal to view the documentation in the browser.

## Parallelization

We use `dask` in many functions for delayed computation and parallelization. Since DePSI operates with Xarray, in most cases, we use Xarray's interface with Dask Arrays, such as `xarray.apply_gufunc` or `xarray.map_blocks` to perform parallel computation. Please refer to the [Xarray Tutorial of Parallelizing Custom Functions](https://tutorial.xarray.dev/advanced/parallel-intro.html) as the best practices for implementing parallelization in DePSI.