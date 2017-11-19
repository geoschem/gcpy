""" Horizontal grid definitions and helper functions """

#: definitions

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

    pass

def calc_rectilinear_lat_edge(lat_stride, half_polar_grid):
    """ Compute latitude edge vector for a rectilinear grid.

    """

    pass

def calc_rectilinear_grid_area(lon_edge,lat_edge):
    """ Compute grid cell areas (in m2) for a rectilinear grid.

    """
    
    pass
