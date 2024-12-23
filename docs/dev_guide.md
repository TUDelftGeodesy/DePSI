# Developer Guide

## Installation guide

The Python implementation of DePSI is under development. At present you can only install it from the GitHub repository.

It is assumed that you have `mamba` installed. If not, you can find the installation instructions [here](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html). Other package managers like `conda` or `venv` can be used as well.

Clone this repository and `cd` into it:

```bash
git clone git@github.com:TUDelftGeodesy/DePSI.git
cd DePSI
```

Create a new conda environment (here we give an example name `depsi-dev`) with `mamba`.:

```bash
mamba create -c conda-forge -n depsi-dev python=3.12
```

Here we use Python 3.12 since we aim to support python 3.10 and above.

Activate the environment:

```bash
mamba activate depsi-dev
```

Install this package in development mode, with extra dependencies for development and documentation:

```bash
pip install -e ".[dev,docs]"
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

## Documentation

We use `mkdocs` for documentation. 

To check the documentation at local, use the following command in the root directory of the repository:

```bash
mkdocs serve
```

This will build and render the documentation at a local server. Follow the link provided in the terminal to view the documentation in the browser.