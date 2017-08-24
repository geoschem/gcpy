""" Internal utilities for helping to manage xarray and numpy
objects used throughout GCPy """

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
import pandas as pd
import xarray as xr

def maybe_as_array(obj):
    """ Try to cast an object as a numpy array, but respect objects which are
    already ndarray-compatible (e.g. pandas Series and xarray DataArrays) """

    if isinstance(obj, (pd.Series, xr.DataArray)):
        return obj
    else:
        return np.asarray(obj)


