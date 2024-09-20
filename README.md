# PyDePSI (Tentative name)

This is repository is WIP, where we are developing a Python package for inteferometric SAR processing. The software will be inspired by the MATLAB software DePSI, but implemented in Python and include recent developments in the field.

## Installation for development

It is assumed that you have `mamba` installed. If not, you can find the installation instructions [here](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html). Other package managers like `conda` or `venv` can be used as well.

Clone this repository and `cd` into it:

```bash
git clone git@github.com:MotionbyLearning/PyDePSI.git
cd PyDePSI
```

Create a new conda environment (here we give an example name `pydepsi-dev`) with `mamba`.:

```bash
mamba create -c conda-forge -n pydepsi-dev python=3.12
```

Here we use Python 3.12 since we aim to support python 3.10 and above.

Activate the environment:

```bash
mamba activate pydepsi-dev
```

Install this package in development mode:

```bash
pip install -e .[dev,docs]
```

In the end, install the pre-commit hooks:
```bash
pre-commit install
```

## Useful reading material

- [Python packaging user guide](https://packaging.python.org/)
- [Testing in Python](https://docs.kedro.org/en/stable/development/automated_testing.html)
- [Code formatting and linting](https://docs.kedro.org/en/stable/development/linting.html)
