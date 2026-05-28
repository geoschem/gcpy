"""
Contains tools for working with cubed-sphere data.

Originally developed by Liam Bindle and Sebastian Eastham.
Included into GCPy by Lizzie Lundgren and Bob Yantosca.

Style updates suggested by the "Pylint" python linter have been adopted.

Examples
--------
Generate a C24 grid and find the indices for a given lat/lon:

.. code-block:: python

   import gcpy
   c24_grid = gcpy.gen_grid(24)            # Generate C24 grid
   lat = 40.0                              # Specify target latitude
   lon = 150.0                             # Specify target longitude
   idx = gcpy.find_index(lat, lon, c24_grid)
   # Returns numpy ndarray [[3],[7],[1]] which can be used to index
   # data (nf=3, Ydim=7, Xdim=1)

   # Check and use to get data
   datafile = '/n/home/GeosChem.SpeciesConc.20190701_0000z.nc4'
   import xarray as xr
   ds = xr.open_dataset(datafile)
   nf    = idx[0,0]
   Ydim  = idx[1,0]
   Xdim  = idx[2,0]
   ds['lats'].isel(nf=nf, Ydim=Ydim, Xdim=Xdim).item()
   # prints 38.711082458496094
   ds['lons'].isel(nf=nf, Ydim=Ydim, Xdim=Xdim).item()
   # prints 151.61871337890625
   ds['SpeciesConcVV_O3'].isel(
       time=0, lev=0, nf=nf, Ydim=Ydim, Xdim=Xdim
   ).item()
   # prints 2.7790051149167994e-08
"""
import numpy as np
import xarray as xr
import gcpy
try:
    import pyproj
    import shapely.ops
    import shapely.geometry
except ImportError as exc:
    raise ImportError(
        "gcpy.cstools needs packages 'pyproj' and 'shapely'!"
    ) from exc

# Constants
RAD_TO_DEG = 180.0 / np.pi
DEG_TO_RAD = np.pi / 180.0


def extract_grid(data):
    """
    Extracts grid information from an xarray Dataset and returns it
    as a cubed-sphere xarray Dataset.

    Parameters
    ----------
    data : xarray.Dataset or xarray.DataArray
        The input dataset.

    Returns
    -------
    data_cs : xarray.Dataset or None
        Cubed-sphere grid extracted from ``data``. Returns ``None``
        if the data is not on a cubed-sphere grid.
    """
    gcpy.util.verify_variable_type(data, (xr.DataArray, xr.Dataset))

    if not is_cubed_sphere(data):
        return None

    cs_res = get_cubed_sphere_res(data)
    return gen_grid(cs_res)


def read_gridspec(gs_obj):
    """
    Reads a GridSpec object and returns an xarray Dataset.

    Parameters
    ----------
    gs_obj : GridSpec
        The GridSpec object to read.

    Returns
    -------
    ds : xarray.Dataset
        The same data as an xarray Dataset object.
    """
    n_cs = gs_obj._tiles[0].area.shape[0]
    lon   = np.zeros((6, n_cs, n_cs))
    lon_b = np.zeros((6, n_cs+1, n_cs+1))
    lat   = np.zeros((6, n_cs, n_cs))
    lat_b = np.zeros((6, n_cs+1, n_cs+1))
    area  = np.zeros((6, n_cs, n_cs))
    for i in range(6):
        tile = gs_obj._tiles[i]
        lon_b[i,...] = tile.supergrid_lons[::2,::2]
        lat_b[i,...] = tile.supergrid_lats[::2,::2]
        lon[i,...]   = tile.supergrid_lons[1::2,1::2]
        lat[i,...]   = tile.supergrid_lats[1::2,1::2]
        area[i,...]  = tile.area[...]

    data = xr.Dataset(
        data_vars={
            "area":  (['nf','Ydim','Xdim'], area),
            "lon":   (['nf','Ydim','Xdim'], lon),
            "lat":   (['nf','Ydim','Xdim'], lat),
            "lon_b": (['nf','Ydim_b','Xdim_b'], lon_b),
            "lat_b": (['nf','Ydim_b','Xdim_b'], lat_b),
        },
        coords={
            "nf":     (['nf'],     list(range(6))),
            "Ydim":   (['Ydim'],   list(range(n_cs))),
            "Xdim":   (['Xdim'],   list(range(n_cs))),
            "Ydim_b": (['Ydim_b'], list(range(n_cs+1))),
            "Xdim_b": (['Xdim_b'], list(range(n_cs+1))),
        },
        attrs={
            "description": f"c{n_cs:d} grid data"
        },
    )
    return data


def face_area(lon_b, lat_b, r_sphere=6.375e6):
    r"""
    Calculates the area of cubed-sphere grid cells on one face.

    Inputs must be in degrees. Edge arrays must be shaped
    ``[N+1 x N+1]``.

    Parameters
    ----------
    lon_b : array_like
        Longitude bounds of grid cell edges in degrees.
    lat_b : array_like
        Latitude bounds of grid cell edges in degrees.
    r_sphere : float, optional
        Radius of the Earth in metres.
        Default: ``6.375e6``

    Returns
    -------
    cs_area : numpy.ndarray
        Surface area in m\ :sup`2` for each grid box of a
        cubed-sphere grid face.
    """
    lon_b_rad = lon_b * DEG_TO_RAD
    lat_b_rad = lat_b * DEG_TO_RAD

    r_sq = r_sphere * r_sphere
    n_cs = lon_b.shape[1] - 1

    cs_area = np.zeros((n_cs, n_cs))
    valid_combo = np.array([[1,2,4],[2,3,1],[3,2,4],[4,1,3]]) - 1

    for i_lon in range(n_cs):
        for i_lat in range(n_cs):
            lon_corner = np.zeros(4)
            lat_corner = np.zeros(4)
            xyz_corner = np.zeros((4,3))
            for i_vert in range(4):
                x_lon = i_lon + (i_vert > 1)
                x_lat = i_lat + (i_vert in (0, 3))
                lon_corner[i_vert] = lon_b_rad[x_lon, x_lat]
                lat_corner[i_vert] = lat_b_rad[x_lon, x_lat]
            for i_vert in range(4):
                xyz_corner[i_vert,:] = ll2xyz(
                    lon_corner[i_vert],
                    lat_corner[i_vert]
                )
            tot_ang = 0.0
            for i_corner in range(4):
                curr_combo = valid_combo[i_corner,:]
                xyz_mini = np.zeros((3,3))
                for i_mini in range(3):
                    xyz_mini[i_mini,:] = xyz_corner[curr_combo[i_mini],:]
                curr_ang = sphere_angle(
                    xyz_mini[0,:],
                    xyz_mini[1,:],
                    xyz_mini[2,:]
                )
                tot_ang += curr_ang
            cs_area[i_lon, i_lat] = r_sq * (tot_ang - (2.0 * np.pi))

    return cs_area


def ll2xyz(lon_pt, lat_pt):
    """
    Converts a lon/lat pair in radians to Cartesian coordinates.

    This function is vectorizable.

    Parameters
    ----------
    lon_pt : float
        Longitude in radians.
    lat_pt : float
        Latitude in radians.

    Returns
    -------
    coords : list of numpy.float64
        Cartesian vector coordinates ``[x_pt, y_pt, z_pt]``
        normalised to the unit sphere.
    """
    x_pt = np.cos(lat_pt) * np.cos(lon_pt)
    y_pt = np.cos(lat_pt) * np.sin(lon_pt)
    z_pt = np.sin(lat_pt)

    return [x_pt, y_pt, z_pt]


def sphere_angle(e_1, e_2, e_3):
    """
    Computes the angle between three points on a sphere.

    Parameters
    ----------
    e_1 : array_like
        Cartesian ``(x, y, z)`` coordinates of the midpoint.
    e_2 : array_like
        Cartesian ``(x, y, z)`` coordinates of a point on one
        side of the midpoint.
    e_3 : array_like
        Cartesian ``(x, y, z)`` coordinates of a point on the
        other side of the midpoint.

    Returns
    -------
    angle : float
        The spherical angle at point ``e_1`` in radians.
    """
    p_vec = np.ones(3)
    q_vec = np.ones(3)
    p_vec[0] = e_1[1]*e_2[2] - e_1[2]*e_2[1]
    p_vec[1] = e_1[2]*e_2[0] - e_1[0]*e_2[2]
    p_vec[2] = e_1[0]*e_2[1] - e_1[1]*e_2[0]

    q_vec[0] = e_1[1]*e_3[2] - e_1[2]*e_3[1]
    q_vec[1] = e_1[2]*e_3[0] - e_1[0]*e_3[2]
    q_vec[2] = e_1[0]*e_3[1] - e_1[1]*e_3[0]

    ddd = np.sum(p_vec * p_vec) * np.sum(q_vec * q_vec)
    if ddd <= 0.0:
        angle = 0.0
    else:
        ddd = np.sum(p_vec * q_vec) / np.sqrt(ddd)
        if np.abs(ddd) > 1.0:
            angle = np.pi / 2.0
        else:
            angle = np.arccos(ddd)

    return angle


def grid_area(
        cs_grid=None,
        cs_res=None,
        stretch_factor=None,
        target_lon=None,
        target_lat=None
):
    r"""
    Returns the surface area in m\ :sup:`2` for each cell of a
    cubed-sphere grid.

    Parameters
    ----------
    cs_grid : dict, optional
        Cubed-sphere grid definition containing keys:
        ``'lat'``, ``'lon'``, ``'lat_b'``, ``'lon_b'``,
        where each value has an extra face dimension of length 6.
        Default: ``None``
    cs_res : int, optional
        Cubed-sphere resolution (number of grid cells along one
        face edge). Default: ``None``
    stretch_factor : float, optional
        Stretching factor for a stretched-grid simulation.
        Default: ``None``
    target_lon : float, optional
        Longitude of the centre of the stretched grid face in
        degrees. Default: ``None``
    target_lat : float, optional
        Latitude of the centre of the stretched grid face in
        degrees. Default: ``None``

    Returns
    -------
    cs_area : numpy.ndarray
        Surface area in m\ :sup:`2` for each cell of the cubed-sphere
        grid. Array shape is ``(6, n, n)`` following the GMAO convention, 
        where ``n`` is the number of cells along a face edge.
    """
    if cs_res is None:
        cs_res = cs_grid['lon_b'].shape[-1] - 1
    elif cs_grid is None:
        cs_grid = gcpy.gen_grid(cs_res, stretch_factor, target_lon, target_lat)
    elif cs_grid is not None and cs_res is not None:
        assert cs_res == cs_grid['lon_b'].shape[-1], \
            'Routine grid_area received inconsistent inputs'

    cs_area = np.zeros((6, cs_res, cs_res))
    for i_face in range(6):
        cs_area[i_face,:,:] = face_area(
            cs_grid['lon_b'][i_face,:,:],
            cs_grid['lat_b'][i_face,:,:]
        )

    return cs_area


def gen_grid(
        n_cs,
        stretch_factor=None,
        target_lon=None,
        target_lat=None
):
    """
    Returns an xarray Dataset specifying a cubed-sphere grid,
    optionally with stretching applied.

    Parameters
    ----------
    n_cs : int
        Number of grid boxes along a single face of the
        cubed-sphere.
    stretch_factor : float, optional
        Stretching factor for a stretched-grid simulation.
        Default: ``None``
    target_lon : float, optional
        Longitude of the centre of the stretched grid face in
        degrees. Default: ``None``
    target_lat : float, optional
        Latitude of the centre of the stretched grid face in
        degrees. Default: ``None``

    Returns
    -------
    grid : xarray.Dataset
        Cubed-sphere grid definition containing variables:
        ``'lat'``, ``'lon'``, ``'lat_b'``, ``'lon_b'``,
        ``'area'``, where each has a face dimension of length 6.
    """
    if stretch_factor is not None:
        cs_temp, _ = gcpy.make_grid_sg(
            n_cs,
            stretch_factor,
            target_lon,
            target_lat
        )
    else:
        cs_temp = gcpy.csgrid_gmao(n_cs)

    return xr.Dataset(
        {
            'nf':     (['nf'],           np.array(range(6))),
            'Ydim':   (['Ydim'],         np.array(range(n_cs))),
            'Xdim':   (['Xdim'],         np.array(range(n_cs))),
            'Ydim_b': (['Ydim_b'],       np.array(range(n_cs+1))),
            'Xdim_b': (['Xdim_b'],       np.array(range(n_cs+1))),
            'lat':    (['nf','Ydim','Xdim'],       cs_temp['lat']),
            'lon':    (['nf','Ydim','Xdim'],       cs_temp['lon']),
            'lat_b':  (['nf','Ydim_b','Xdim_b'],   cs_temp['lat_b']),
            'lon_b':  (['nf','Ydim_b','Xdim_b'],   cs_temp['lon_b']),
            'area':   (['nf','Ydim','Xdim'],       grid_area(cs_temp)),
        }
    )


def corners_to_xy(x_c, y_c):
    """
    Creates XY coordinates for each grid box from corner arrays.

    Developed, tested, and supplied by Liam Bindle.

    Parameters
    ----------
    x_c : numpy.ndarray
        Grid-box corner longitudes. Array shape ``(n+1, n+1)``,
        where ``n`` is the cubed-sphere grid size.
    y_c : numpy.ndarray
        Grid-box corner latitudes. Array shape ``(n+1, n+1)``,
        where ``n`` is the cubed-sphere grid size.

    Returns
    -------
    x_y : numpy.ndarray
        Grid-box Cartesian coordinates. Array shape
        ``(n, n, 5, 2)``, where ``n`` is the cubed-sphere grid
        size and the last two dimensions are the 5 corner points
        (first and last identical) and the XY pair respectively.
    """
    p_0 = slice(0, -1)
    p_1 = slice(1, None)
    boxes_x = np.moveaxis(
        np.array(
            [
                x_c[p_0, p_0],
                x_c[p_1, p_0],
                x_c[p_1, p_1],
                x_c[p_0, p_1],
                x_c[p_0, p_0],
            ]
        ),
        0, -1
    )
    boxes_y = np.moveaxis(
        np.array(
            [
                y_c[p_0, p_0],
                y_c[p_1, p_0],
                y_c[p_1, p_1],
                y_c[p_0, p_1],
                y_c[p_0, p_0],
            ]
        ),
        0, -1
    )
    return np.moveaxis(np.array([boxes_x, boxes_y]), 0, -1)


def central_angle(x_0, y_0, x_1, y_1):
    """
    Returns the central angle (great-circle distance) between two
    points on the sphere.

    This function is vectorizable.
    Developed, tested, and supplied by Liam Bindle.

    Parameters
    ----------
    x_0 : float
        Longitude of the first point in degrees.
    y_0 : float
        Latitude of the first point in degrees.
    x_1 : float
        Longitude of the second point in degrees.
    y_1 : float
        Latitude of the second point in degrees.

    Returns
    -------
    distance : float
        Great-circle distance between ``(x_0, y_0)`` and
        ``(x_1, y_1)`` in degrees.
    """
    x_0 = x_0 * DEG_TO_RAD
    x_1 = x_1 * DEG_TO_RAD
    y_0 = y_0 * DEG_TO_RAD
    y_1 = y_1 * DEG_TO_RAD

    return np.arccos(
        np.sin(y_0) * np.sin(y_1) +
        np.cos(y_0) * np.cos(y_1) *
        np.cos(np.abs(x_0 - x_1))
    ) * RAD_TO_DEG


def find_index_single(
        lat,
        lon,
        x_centers_flat,
        y_centers_flat,
        xy_polygon_defs,
        cs_size,
        latlon_crs,
        jitter_size=0.0
):
    """
    Returns the cubed-sphere grid box indices for a single
    latitude/longitude point.

    Called internally by :func:`find_index`.

    Parameters
    ----------
    lat : float
        Latitude of the query point in degrees.
    lon : float
        Longitude of the query point in degrees.
    x_centers_flat : numpy.ndarray
        Flattened array of cubed-sphere grid cell centre longitudes
        in column-major (Fortran) order.
    y_centers_flat : numpy.ndarray
        Flattened array of cubed-sphere grid cell centre latitudes
        in column-major (Fortran) order.
    xy_polygon_defs : numpy.ndarray
        XY polygon definitions for cubed-sphere grid boxes,
        as returned by :func:`corners_to_xy`.
    cs_size : int
        Number of grid cells along one cubed-sphere face edge.
    latlon_crs : pyproj.Proj
        Projection object from ``pyproj.Proj("+proj=latlon")``.
    jitter_size : float, optional
        If the point cannot be matched to a grid box, the longitude
        is shifted by this distance in metres before retrying.
        Useful for points near the poles.
        Default: ``0.0``

    Returns
    -------
    nf_cs : int
        Index of the cubed-sphere face (0–5).
    ydim_cs : int
        YDim (latitude) index of the matched grid box.
    xdim_cs : int
        XDim (longitude) index of the matched grid box.

    Raises
    ------
    ValueError
        If the point cannot be matched to any grid box and
        ``jitter_size`` is ``0.0``.
    """
    x_find = lon
    y_find = lat
    gnomonic_crs = pyproj.Proj(
        f'+proj=gnom +lat_0={y_find} +lon_0={x_find}'
    )

    # Generate all distances
    distances = central_angle(
        x_find, y_find, x_centers_flat, y_centers_flat
    )
    four_nearest_indexes = np.argpartition(distances, 4)[:4]
    four_nearest_indexes = np.unravel_index(
        four_nearest_indexes, (6, cs_size, cs_size)
    )
    four_nearest_xy = xy_polygon_defs[four_nearest_indexes]
    four_nearest_polygons = [
        shapely.geometry.Polygon(polygon_xy)
        for polygon_xy in four_nearest_xy
    ]

    # Transform to gnomonic projection
    gno_transform = pyproj.Transformer.from_proj(
        latlon_crs, gnomonic_crs, always_xy=True
    ).transform
    four_nearest_polygons_gno = [
        shapely.ops.transform(gno_transform, polygon)
        for polygon in four_nearest_polygons
    ]

    # Figure out which polygon contains the point
    xy_find     = shapely.geometry.Point(x_find, y_find)
    xy_find_gno = shapely.ops.transform(gno_transform, xy_find)
    polygon_contains_point = [
        polygon.contains(xy_find_gno)
        for polygon in four_nearest_polygons_gno
    ]

    # If the point cannot be matched (such as can happen near the poles),
    # move the longitude by the jitter_size (in meters) and try again.
    if np.count_nonzero(polygon_contains_point) == 0:
        if jitter_size > 0.0:
            nf_cs, ydim_cs, xdim_cs = find_index_single(
                y_find,
                x_find + jitter_size,
                x_centers_flat,
                y_centers_flat,
                xy_polygon_defs,
                cs_size,
                latlon_crs,
                jitter_size=0.0
            )
        else:
            msg = f'Point at {x_find:8.2f} E, {y_find:8.2f} N '
            msg += 'could not be matched'
            raise ValueError(msg)

    # The first will be selected, if more than one
    polygon_with_point = np.argmax(polygon_contains_point)
    nf_cs   = four_nearest_indexes[0][polygon_with_point]
    ydim_cs = four_nearest_indexes[1][polygon_with_point]
    xdim_cs = four_nearest_indexes[2][polygon_with_point]

    return nf_cs, ydim_cs, xdim_cs


def find_index(lat, lon, grid, jitter_size=0.0):
    """
    Returns the cubed-sphere grid box indices for one or more
    latitude/longitude coordinates.

    Based on a routine developed, tested, and supplied by
    Liam Bindle.

    Parameters
    ----------
    lat : float or array_like
        Latitude(s) of the query point(s) in degrees.
    lon : float or array_like
        Longitude(s) of the query point(s) in degrees.
    grid : xarray.Dataset
        Cubed-sphere grid definition containing variables:
        ``'lat'``, ``'lon'``, ``'lat_b'``, ``'lon_b'``,
        where each has a face dimension of length 6.
    jitter_size : float, optional
        If a point cannot be matched to a grid box, the longitude
        is shifted by this distance in metres before retrying.
        Useful for points near the poles.
        Default: ``0.0``

    Returns
    -------
    ind : numpy.ndarray
        Array of shape ``(3, N)`` containing ``(nf, YDim, XDim)``
        for each of the ``N`` query points, where:

        - ``nf``   — cubed-sphere face index (0–5)
        - ``YDim`` — latitude index on the matched face
        - ``XDim`` — longitude index on the matched face
    """
    gcpy.util.verify_variable_type(grid, xr.Dataset)

    lon_vec = np.asarray(lon)
    lat_vec = np.asarray(lat)
    n_find  = lon_vec.size

    x_corners      = grid['lon_b'].values
    y_corners      = grid['lat_b'].values
    x_centers      = grid['lon'].values
    y_centers      = grid['lat'].values
    x_centers_flat = x_centers.flatten()
    y_centers_flat = y_centers.flatten()
    cs_size        = x_centers.shape[-1]

    xy_polygon_defs = np.zeros((6, cs_size, cs_size, 5, 2))
    for nf_cs in range(6):
        xy_polygon_defs[nf_cs, ...] = corners_to_xy(
            x_c=x_corners[nf_cs, :, :],
            y_c=y_corners[nf_cs, :, :]
        )
    latlon_crs = pyproj.Proj("+proj=latlon")

    idx = np.full((3, n_find), 0)
    for x_find, y_find, i_find in zip(
            np.nditer(lon_vec), np.nditer(lat_vec), range(n_find)
    ):
        nf_cs, ydim_cs, xdim_cs = find_index_single(
            y_find,
            x_find,
            x_centers_flat,
            y_centers_flat,
            xy_polygon_defs,
            cs_size,
            latlon_crs,
            jitter_size=jitter_size
        )
        idx[:, i_find] = [nf_cs, ydim_cs, xdim_cs]

    return idx


def is_cubed_sphere(data):
    """
    Determines whether data is placed on a cubed-sphere grid.

    Parameters
    ----------
    data : xarray.Dataset or xarray.DataArray
        The input data to be tested.

    Returns
    -------
    bool
        ``True`` if the data is on a cubed-sphere grid,
        ``False`` otherwise.

    Notes
    -----
    A cubed-sphere data file satisfies at least one of:

    1. Has a dimension named ``"nf"`` (GCHP/GEOS diagnostic files).
    2. The lat/lon size ratio is exactly 6 (GCHP/GEOS checkpoints).
    """
    gcpy.util.verify_variable_type(data, (xr.DataArray, xr.Dataset))

    if is_cubed_sphere_diag_grid(data):
        return True
    if is_cubed_sphere_rst_grid(data):
        return True
    return False


def is_cubed_sphere_diag_grid(data):
    """
    Determines whether a dataset has cubed-sphere History
    diagnostic dimensions (i.e. a dimension named ``"nf"``).

    Parameters
    ----------
    data : xarray.DataArray or xarray.Dataset
        The input data.

    Returns
    -------
    bool
        ``True`` if the grid has History diagnostic dimensions
        (contains an ``"nf"`` dimension), ``False`` otherwise.
    """
    return "nf" in data.dims


def is_cubed_sphere_rst_grid(data):
    """
    Determines whether a dataset has cubed-sphere restart file
    dimensions (i.e. lat and lon with ``lat = lon * 6``).

    Parameters
    ----------
    data : xarray.DataArray or xarray.Dataset
        The input data.

    Returns
    -------
    bool
        ``True`` if the grid has restart file dimensions,
        ``False`` otherwise.

    Notes
    -----
    If the ``"lat"`` coordinate is present, the check is performed
    by comparing the lengths of ``lat`` and ``lon``. Otherwise,
    variable names are searched for the prefix ``"SPC_"``.

    .. todo::
        Reconsider this approach if GC-Classic restart variables
        are ever renamed to start with ``SPC_``, or if GCHP internal
        state variables are renamed.
    """
    gcpy.util.verify_variable_type(data, (xr.DataArray, xr.Dataset))

    if "lat" in data.coords:
        return len(data.coords["lat"]) == len(data.coords["lon"]) * 6

    if isinstance(data, xr.Dataset):
        return "SPC_" in data.data_vars.keys()
    return "SPC_" in data.name


def get_cubed_sphere_res(data):
    """
    Returns the number of grid cells along one edge of a
    cubed-sphere grid face.

    For example, returns ``24`` for a C24 grid (which has
    24 × 24 grid cells per face).

    Parameters
    ----------
    data : xarray.DataArray or xarray.Dataset
        The input data.

    Returns
    -------
    cs_res : int
        The cubed-sphere resolution. Returns ``0`` if the data
        is not on a cubed-sphere grid.

    Notes
    -----
    ``dims`` is a dict in Dataset objects but a tuple in DataArray
    objects. Returning the length of the corresponding coords array
    works in both cases.
    """
    gcpy.util.verify_variable_type(data, (xr.DataArray, xr.Dataset))

    if not is_cubed_sphere(data):
        return 0

    if is_cubed_sphere_rst_grid(data):
        return len(data.coords["lon"])
    return len(data.coords["Xdim"])


def is_gchp_lev_positive_down(data):
    """
    Determines whether GCHP data levels are ordered from the top of
    the atmosphere downwards.

    The vertical ordering follows these conventions:

    - Checkpoint files: ``lev:positive = "down"``
    - Emissions collection: ``lev:positive = "down"``
    - All other collections: ``lev:positive = "up"``

    Parameters
    ----------
    data : xarray.DataArray or xarray.Dataset
        The input data.

    Returns
    -------
    bool
        ``True`` if levels are ordered top-of-atmosphere downwards,
        ``False`` if ordered surface upwards.
    """
    gcpy.util.verify_variable_type(data, (xr.DataArray, xr.Dataset))

    if is_cubed_sphere_rst_grid(data):
        return True
    if is_cubed_sphere_diag_grid(data):
        emis_vars = [
            var for var in data.data_vars if var.startswith("Emis")
        ]
        if len(emis_vars) > 0:
            return True
    return False
