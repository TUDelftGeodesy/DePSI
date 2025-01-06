# DePSI 
This is the WIP repository for the Python version of DePSI *(), a Python package for inteferometric SAR processing. The software is inspired by the MATLAB version DePSI. In this repository, we implement classic DePSI algorithms and new InSAR developments in Python.

## Installation for development

It is assumed that you have `mamba` installed. If not, you can find the installation instructions [here](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html). Other package managers like `conda` or `venv` can be used as well.

Clone this repository 

Then`cd` into it:

```bash
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

Install this package in development mode:

```bash
pip install -e ".[dev,docs]"
```

In the end, install the pre-commit hooks:
```bash
pre-commit install
```

## Useful reading material

- [Python packaging user guide](https://packaging.python.org/)
- [Testing in Python](https://docs.kedro.org/en/stable/development/automated_testing.html)
- [Code formatting and linting](https://docs.kedro.org/en/stable/development/linting.html)

## License

Copyright (c) 2023 - 2025, Netherlands eScience Center & Delft University of Technology

Apache Software License 2.0
