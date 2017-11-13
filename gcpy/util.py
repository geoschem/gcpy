""" Internal utilities for helping to manage xarray and numpy
objects used throughout GCPy """

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
import pandas as pd
import xarray as xr


def convert_lon(data, dim='lon', format='atlantic', neg_dateline=True):
    """ Convert longitudes from -180..180 to 0..360, or vice-versa.

    Parameters
    ----------
    data : DataArray or Dataset
        The container holding the data to be converted; the dimension indicated
        by 'dim' must be associated with this container
    dim : str
        Name of dimension holding the longitude coordinates
    format : {'atlantic', 'pacific'}
        Control whether or not to shift from -180..180 to 0..360 ('pacific') or
        from 0..360 to -180..180 ('atlantic')
    neg_dateline : logical
        If True, then the international dateline is set to -180. instead of 180.

    Returns
    -------
    data, with dimension 'dim' altered according to conversion rule.

    """

    data_copy = data.copy()

    lon = data_copy[dim].values
    new_lon = np.empty_like(lon)

    # Tweak offset for rolling the longitudes later
    offset = 0 if neg_dateline else 1

    if format not in ['atlantic', 'pacific']:
        raise ValueError("Cannot convert longitudes for format '{}'; "
                         "please choose one of 'atlantic' or 'pacific'"
                         .format(format))

    # Create a mask to decide how to mutate the longitude values
    if format == 'atlantic':
        mask = lon >= 180 if neg_dateline else lon > 180

        new_lon[mask] = -(360. - lon[mask])
        new_lon[~mask] = lon[~mask]

        roll_len = len(data[dim])//2 - offset

    elif format == 'pacific':
        mask = lon < 0.

        new_lon[mask] = lon[mask] + 360.
        new_lon[~mask] = lon[~mask]

        roll_len = -len(data[dim])//2 - offset

    # Copy mutated longitude values into copied data container
    data_copy[dim].values = new_lon
    data_copy = data_copy.roll(**{dim: roll_len})

    return data_copy


def create_usa_mask():
    # TODO
    pass


def find_cells_by_country():
    # TODO
    pass


def in_range():
    # TODO
    pass


def maybe_as_array(obj):
    """ Try to cast an object as a numpy array, but respect objects which are
    already ndarray-compatible (e.g. pandas Series and xarray DataArrays) """

    if isinstance(obj, (pd.Series, xr.DataArray)):
        return obj
    else:
        return np.asarray(obj)


def search():
    # TODO
    pass


def strpad():
    # TODO
    pass


def strsci():
    # TODO
    pass