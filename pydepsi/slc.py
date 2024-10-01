"""slc.py: Functions for SLC related operations."""

import dask.array as da
import numpy as np


def ifg_to_slc(mother_slc, ifgs):
    """Convert a stack of interferograms to SLCs.

    The conversion will be implemented by conjugated multiplication of the interferograms complex values
    with the complex values of the mother SLC, and then dividing by the squared magnitude of the mother complex.

    Parameters
    ----------
    mother_slc : Xarray.Dataset
        Mother SLC. This Dataset should have three dimensions ('azimuth', 'range', 'time').
        The 'azimuth' and 'range' dimensions should be the same as `ifgs`.
        The 'time' dimension should have size 1.
    ifgs : Xarray.Dataset
        Interferograms. This Dataset should have three dimensions ('azimuth', 'range', 'time').
        The 'azimuth' and 'range' dimensions should be the same as `mother_slc`.

    Returns
    -------
    Xarray.Dataset
        SLCS converted from the interferograms.
    """
    slc_out = ifgs.copy()
    meta_arr = np.array((), dtype=np.complex64)
    slc_complex = da.apply_gufunc(
        _slc_complex_recontruct, "(),()->()", mother_slc["complex"], slc_out["complex"], meta=meta_arr
    )
    slc_out = slc_out.assign({"complex": (("azimuth", "range", "time"), slc_complex)})
    return slc_out


def _slc_complex_recontruct(mother_slc_complex, ifg_complex):
    """Reconstruct the SLC complex values from the mother SLC complex and the interferogram complex.

    The implementation equivalent to: (ifg_complex * mother_slc_complex.conj()) / (np.abs(mother_slc_complex) ** 2).
    """
    return ifg_complex / mother_slc_complex
