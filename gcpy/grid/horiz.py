""" Horizontal grid definitions and helper functions """

import numpy as np
from numpy import asarray
from .. constants import R_EARTH

def calc_rectilinear_lon_edge(lon_stride, center_at_180):
    """ Compute longitude edge vector for a rectilinear grid.

    Parameters
    ----------
    lon_stride : float
        Stride length in degrees. For example, for a standard GEOS-Chem Classic
        4x5 grid, lon_stride would be 5.

    center_at_180: boolean
        Whether or not the grid should have a cell center at 180 degrees (i.e.
        on the date line). If true, the first grid cell is centered on the date
        line; if false, the first grid edge is on the date line.

    Returns
    -------
    Longitudes of cell edges in degrees East.

    Notes
    -----
    All values are forced to be between [-180,180]. For a grid with N cells in
    each band, N+1 edges will be returned, with the first and last value being
    duplicates.

    Examples
    --------
    >>> from gcpy.grid.horiz import calc_rectilinear_lon_edge
    >>> calc_rectilinear_lon_edge(5.0,true)
    np.array([177.5,-177.5,-172.5,...,177.5])

    See Also
    --------
    [NONE]

    """

    n_lon = np.round(360.0/lon_stride)
    lon_edge = np.linspace(-180.0,180.0,num=n_lon+1)
    if center_at_180:
        lon_edge = lon_edge - (lon_stride/2.0)

    lon_edge[lon_edge<-180.0] = lon_edge[lon_edge<-180] + 360.0
    lon_edge[lon_edge>180.0] = lon_edge[lon_edge>180.0] - 360.0

    return lon_edge

def calc_rectilinear_lat_edge(lat_stride, half_polar_grid):
    """ Compute latitude edge vector for a rectilinear grid.

    """

    if half_polar_grid:
        start_pt = 90.0 + (lat_stride/2.0)
    else:
        start_pt = 90.0

    lat_edge = np.linspace(-1.0*start_pt,start_pt,
        num=1+np.round(2.0*start_pt/lat_stride))

    # Force back onto +/- 90
    lat_edge[lat_edge>90.0] = 90.0
    lat_edge[lat_edge<-90.0] = -90.0

    return lat_edge

def calc_rectilinear_grid_area(lon_edge,lat_edge):
    """ Compute grid cell areas (in m2) for a rectilinear grid.

    """
    
    pass
def calc_delta_lon(lon_edge):
    """ Compute grid cell longitude widths from an edge vector.

    Parameters
    ----------
    lon_edge: float
        Vector of longitude edges, in degrees East.

    Returns
    -------
    Width of each cell, degrees East

    Notes
    -----
    Accounts for looping over the domain.

    Examples
    --------

    #TODO
    """

    n_lon = (lon_edge.size) - 1

    lon_edge = asarray(lon_edge)

    # Set up output array
    lon_delta = np.zeros((n_lon))

    offset = 0.0
    next_lon = lon_edge[0]
    for i_lon in range(0,n_lon):
        last_lon = next_lon
        next_lon = lon_edge[i_lon + 1] + offset
        while next_lon < last_lon:
            offset = offset + 360.0
            next_lon = next_lon + 360.0
        lon_delta[i_lon] = next_lon - last_lon

    return lon_delta
