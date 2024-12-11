"""test_slc.py"""

import numpy as np
import xarray as xr

from depsi.slc import _slc_complex_recontruct, ifg_to_slc


def test__slc_complex_recontruct():
    mother_slc_complex = xr.DataArray([1 + 1j, 2 + 2j, 3 + 3j])
    ifg_complex = xr.DataArray([1 + 1j, 2 + 2j, 3 + 3j])
    assert np.allclose(_slc_complex_recontruct(mother_slc_complex, ifg_complex), np.array([1 + 0j, 1 + 0j, 1 + 0j]))


def test_ifg_to_slc():
    arr = np.expand_dims(np.array([1 + 1j, 2 + 2j, 3 + 3j]), axis=(1, 2))
    res = np.expand_dims(np.array([1 + 0j, 1 + 0j, 1 + 0j]), axis=(1, 2)).repeat(4, axis=2)
    mother_slc = xr.Dataset({"complex": (("azimuth", "range", "time"), arr)})
    ifgs = xr.Dataset({"complex": (("azimuth", "range", "time"), arr.repeat(4, axis=2))})
    assert np.allclose(ifg_to_slc(mother_slc, ifgs)["complex"].values, res)
