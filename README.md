# DePSI

<table style="width:100%; margin-left:auto; margin-right:auto;">
<tr>
<td width="65%">
  <img src="docs/assets/depsi_snapshot.gif" alt="DePSI Snapshot" style="width:100%; height:auto;">
</td>
<td width="35%" valign="middle">

[![DOI](https://zenodo.org/badge/618377671.svg)](https://doi.org/10.5281/zenodo.19951332)

[![PyPI](https://img.shields.io/pypi/v/depsi.svg?colorB=blue)](https://pypi.python.org/project/depsi/)

[![Build and pytest](https://github.com/TUDelftGeodesy/DePSI/actions/workflows/build.yml/badge.svg)](https://github.com/TUDelftGeodesy/DePSI/actions/workflows/build.yml)

[![License](https://img.shields.io/github/license/TUDelftGeodesy/DePSI)](https://opensource.org/licenses/Apache-2.0)

</td>
</tr>
</table>

**This repository holds a beta version of DePSI. The code is still under development and should still be properly tested and documented. Still, you are invited to have a look at the current code base and send us your suggestions.**

DePSI (van Leijen, 2014) is an open source Python software package for Persistent Scatterer Interferometric SAR (PS-InSAR) processing. It provides a comprehensive suite of functions to identify Persistent Scatterer (PS) points from a time series of interferometric SAR data, and estimate their deformation time series based on customizable assumptions, such as a predefined deformation model.

DePSI was originally implemented in MATLAB in 2014. The Python implementation of DePSI is motivated by acconmondating the classic DePSI algorithm into a modern software framework, enabling easier maintenance and contribution from the InSAR community. The Python version also enables the support of parallel computing and handling large datasets by leveraging `Dask` and `Xarray` libraries.

For the original MATLAB version, which is static without further development, please refer to the [MATLAB branch](https://github.com/TUDelftGeodesy/DePSI/tree/stable) of this repository.

## Documentation

You can find more information about DePSI in the [documentation site](https://tudelftgeodesy.github.io/DePSI_group).

## License

Copyright (c) 2023 - 2025, Netherlands eScience Center & Delft University of Technology

Apache Software License 2.0

## References

[1] [Van Leijen, F.J., 2014. "Persistent scatterer interferometry based on geodetic estimation theory." Delft University of Technology, The Netherlands.](https://repository.tudelft.nl/record/uuid:5dba48d7-ee26-4449-b674-caa8df93e71e)
