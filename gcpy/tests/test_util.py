#!/usr/bin/env python

from __future__ import print_function

import numpy as np
import pandas as pd
import xarray as xr

from xarray.testing import assert_equal

from gcpy.util import *


def test_convert_lon():
    lons_greater_zero = np.arange(0, 360)
    ones = np.ones_like(lons_greater_zero)
    lons_greater_zero = xr.DataArray(ones, dims=['lon', ],
                                     coords={'lon': lons_greater_zero})

    lons_less_zero = np.arange(-180, 180)
    lons_less_zero = xr.DataArray(ones, dims=['lon', ],
                                  coords={'lon': lons_less_zero})

    # Test atlantic
    assert_equal(convert_lon(lons_greater_zero, 'lon', 'atlantic'),
                 lons_less_zero)

    # Test pacific
    assert_equal(convert_lon(lons_less_zero, 'lon', 'pacific'),
                 lons_greater_zero)


def test_maybe_as_array():
    """ Check that we pass through expected values correctly """

    # Float -> array
    arr1 = 0.
    assert maybe_as_array(arr1) == np.asarray(arr1)

    # List of floats -> array
    arr2 = [0., 1., 2., 3.]
    res2 = maybe_as_array(arr2)
    assert isinstance(res2, np.ndarray)
    np.testing.assert_array_equal(res2, np.asarray(arr2))

    # List of list of floats -> array
    arr3 = [[0., 0.,], [1., 1.], [2., 2.]]
    res3 = maybe_as_array(arr3)
    assert isinstance(res3, np.ndarray)
    np.testing.assert_array_equal(res3, np.asarray(arr3))

    # pandas Series
    arr4 = pd.Series(arr2)
    res4 = maybe_as_array(arr4)
    assert isinstance(res4, pd.Series)

    # xarray DataArray
    arr5 = xr.DataArray(arr2)
    res5 = maybe_as_array(arr5)
    assert isinstance(res5, xr.DataArray)
