"""
Contains tools for working with cubed-sphere data.

Originally developed by Liam Bindle and Sebastian Eastham.
Included into GCPy by Lizzie Lundgren and Bob Yantosca.

Style updates suggested by the "Pylint" python linter have been adopted.

Example:
   import gcpy
   c24_grid=gcpy.gen_grid(24)             # Generate C24 grid
   lat=40.0                               # Specify target latitude
   lon=150.0                              # Specify target longitude
   idx=gcpy.find_index(lat,lon,c24_grid)  # Returns numpy ndarray
                                          #  [[3],[7],[1]] which can be
                                          #  used to index data
                                          #  (nf=3, Ydim=7, Xdim=1)

   # Check and use to get data
   datafile='/n/home/GeosChem.SpeciesConc.20190701_0000z.nc4'
   import xarray ax xr
   ds=xr.open_dataset(datafile)
   nf=idx[0,0]
   Ydim=idx[1,0]
   Xdim=idx[2,0]
   ds['lats'].isel(nf=nf,Ydim=Ydim,Xdim=Xdim].item()
   # prints 38.711082458496094
   ds['lons'].isel(nf=nf,Ydim=Ydim,Xdim=Xdim].item()
   # prints 151.61871337890625
   ds['SpeciesConcVV_O3'].isel(time=0,lev=0,nf=nf,Ydim=Ydim,Xdim=Xdim).item()
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


def extract_grid(
        data
):
    """
    Extracts the grid information from an xarray.Dataset object and
    returns the grid information as a cubed-sphere xarray.Dataset.

    Args:
    -----
    data : xarray.Dataset or xarray.DataArray
        The input dataset

    data_cs: xarray.Dataset or None
        Same data as in argument "ds", but on a cubed-sphere grid
        If the data is not placed on a cubed-sphere grid, then
        this will be returned with the value None.
    """
    gcpy.util.verify_variable_type(data, (xr.DataArray, xr.Dataset))

    if not is_cubed_sphere(data):
        return None

    cs_res = get_cubed_sphere_res(data)
    return gen_grid(cs_res)


def read_gridspec(gs_obj):
    """
    Reads a GridSpec object and returns an xarray.Dataset object.

    Args:
    -----
    gs_obj : GridSpec
        The GridSpec object as input

    Returns:
    --------
    ds : xarray.Dataset
        The same data as an xarray.Dataset object.
    """
    n_cs = gs_obj._tiles[0].area.shape[0]
    lon   = np.zeros((6, n_cs, n_cs))
    lon_b = np.zeros((6, n_cs+1, n_cs+1))
    lat   = np.zeros((6, n_cs, n_cs))
    lat_b = np.zeros((6, n_cs+1, n_cs+1))
    area  = np.zeros((6, n_cs, n_cs))
    for i in range(6):
        tile = gs_obj._tiles[i]
        # lon_b is identical to original definition
        lon_b[i,...] = tile.supergrid_lons[::2,::2]
        # lat_b is dentical to original definition
        lat_b[i,...] = tile.supergrid_lats[::2,::2]
        # lon is NOT identical to original definition
        lon[i,...] = tile.supergrid_lons[1::2,1::2]
        # lat is NOT identical to original definition
        lat[i,...] = tile.supergrid_lats[1::2,1::2]
        area[i,...] = tile.area[...]

    data = xr.Dataset(
        data_vars={
            "area": (['nf','Ydim','Xdim'], area),
            "lon": (['nf','Ydim','Xdim'], lon),
            "lat": (['nf','Ydim','Xdim'], lat),
            "lon_b": (['nf','Ydim_b','Xdim_b'], lon_b),
            "lat_b": (['nf','Ydim_b','Xdim_b'], lat_b),
        },
        coords={
            "nf": (['nf'], list(range(6))),
            "Ydim": (['Ydim'], list(range(n_cs))),
            "Xdim": (['Xdim'], list(range(n_cs))),
            "Ydim_b": (['Ydim_b'], list(range(n_cs+1))),
            "Xdim_b": (['Xdim_b'], list(range(n_cs+1))),
        },
        attrs={
            "description": f"c{n_cs:d} grid data"
        },
    )
    return data


def face_area(
        lon_b,
        lat_b,
        r_sphere=6.375e6
):
    """
    Calculates area of cubed-sphere grid cells on one face.
    Inputs must be in degrees. Edge arrays must be shaped [N+1 x N+1].

    Args:
    -----
    lon_b, lat_b : list of float
        Longitude and latitude bounds (degrees)

    Keyword Args (optional):
    ------------------------
    r_sphere : float
        Radius of Earth (meters).  Default value: 6.375e6

    Returns:
    --------
    cs_area : numpy.ndarray
        Array of surface area (m2) in each grid box of a
        cubed-sphere grid face.
    """

    # Convert inputs to radians
    lon_b_rad = lon_b * DEG_TO_RAD
    lat_b_rad = lat_b * DEG_TO_RAD

    r_sq = r_sphere * r_sphere
    n_cs = lon_b.shape[1] - 1

    # Allocate output array
    cs_area = np.zeros((n_cs,n_cs))

    # Ordering
    valid_combo = np.array([[1,2,4],[2,3,1],[3,2,4],[4,1,3]]) - 1

    for i_lon in range(n_cs):
        for i_lat in range(n_cs):
            lon_corner = np.zeros(4)
            lat_corner = np.zeros(4)
            xyz_corner = np.zeros((4,3))
            for i_vert in range(4):
                x_lon = i_lon + (i_vert > 1)
                x_lat = i_lat + (i_vert == 0 or i_vert == 3)
                lon_corner[i_vert] = lon_b_rad[x_lon,x_lat]
                lat_corner[i_vert] = lat_b_rad[x_lon,x_lat]
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
            cs_area[i_lon,i_lat] = r_sq * (tot_ang - (2.0*np.pi))

    return cs_area


def ll2xyz(
        lon_pt,
        lat_pt
):
    """
    Converts a lon/lat pair (in radians) to Cartesian co-ordinates.
    This is vectorizable.

    Args:
    -----
    lon_pt, lat_pt : float
        Longitude & latitude in radians.

    Returns:
    --------
    [x_pt, y_pt, z_pt] : list of numpy.float64
        Cartesian vector coordinates (X,Y.Z) normalized to
        the unit sphere.
    """
    x_pt = np.cos(lat_pt) * np.cos(lon_pt)
    y_pt = np.cos(lat_pt) * np.sin(lon_pt)
    z_pt = np.sin(lat_pt)

    return [x_pt, y_pt, z_pt]


def sphere_angle(
        e_1,
        e_2,
        e_3
):
    """
    Computes the angle between 3 points on a sphere.

    Args:
    -----
    e_1 : list of float
        (x, y) coordinates at mid point

    e_2, e_3: : list of float
        (x, y) coordinates at points on either side of midpoint

    Returns:
    --------
    angle : float
        The spherical angle at point e_1.
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
        cs_res=None
):
    """
    Return area (m2) for each cell in a cubed-sphere grid

    Args:
    -----
    cs_grid : dict
        Cubed-sphere grid definition as a dict of:
           {'lat'   : lat midpoints,
            'lon'   : lon midpoints,
            'lat_b' : lat edges,
            'lon_b' : lon edges}
        where each value has an extra face dimension of length 6.

    Returns:
    --------
    grid_area : numpy.ndarray
        Surface area (m2) for each cell in the cubed-sphere grid.
        NOTE: Uses GMAO convention, array shape = (6, n, n),
        where n is the number of cells along a face edge.
    """
    # Calculate area on a cubed sphere
    if cs_res is None:
        cs_res = cs_grid['lon_b'].shape[-1] - 1
    elif cs_grid is None:
        cs_grid = gcpy.csgrid_GMAO(cs_res)
    elif cs_grid is not None and cs_res is not None:
        assert cs_res == cs_grid['lon_b'].shape[-1], \
        'Routine grid_area received inconsistent inputs'
    cs_area = np.zeros((6,cs_res,cs_res))
    cs_area[0,:,:] = face_area(
        cs_grid['lon_b'][0,:,:],
        cs_grid['lat_b'][0,:,:]
    )
    for i_face in range(1,6):
        cs_area[i_face,:,:] = cs_area[0,:,:].copy()

    return cs_area


def gen_grid(
        n_cs,
        stretch_factor=None,
        target_lon=None,
        target_lat=None
):
    """
    Returns an xarray.Dataset object specifying a cubed-sphere
    stretched grid.

    Args:
    -----
    n_cs : int
        Number of grid boxes along a single face of the cubed-sphere.

    stretch_factor : int
        Specifies the stretching factor.  Default value: None

    target_lon, target_lat : float
        Specifies the longitude and latitude at the center of the
        cubed-sphere grid face that will be stretched.
        Default values: None, None

    Returns:
    --------
    grid : xarray.Dataset
        Cubed-sphere grid definition containing the variables:
           {'lat'   : lat midpoints,
            'lon'   : lon midpoints,
            'lat_b' : lat edges,
            'lon_b' : lon edges}
        where each value has an extra face dimension of length 6.
    """
    if stretch_factor is not None:
        cs_temp, _ = gcpy.make_grid_SG(
            n_cs,
            stretch_factor,
            target_lon,
            target_lat
        )
    else:
        cs_temp = gcpy.csgrid_GMAO(n_cs)

    return xr.Dataset(
        {'nf':     (['nf'], np.array(range(6))),
         'Ydim':   (['Ydim'], np.array(range(n_cs))),
         'Xdim':   (['Xdim'], np.array(range(n_cs))),
         'Ydim_b': (['Ydim_b'], np.array(range(n_cs+1))),
         'Xdim_b': (['Xdim_b'], np.array(range(n_cs+1))),
         'lat':    (['nf','Ydim','Xdim'], cs_temp['lat']),
         'lon':    (['nf','Ydim','Xdim'], cs_temp['lon']),
	 'lat_b':  (['nf','Ydim_b','Xdim_b'], cs_temp['lat_b']),
	 'lon_b':  (['nf','Ydim_b','Xdim_b'], cs_temp['lon_b']),
         'area':   (['nf','Ydim','Xdim'], grid_area(cs_temp))
        }
    )


def corners_to_xy(
        x_c,
        y_c
):
    """
    Creates xy coordinates for each grid-box
    Developed, tested, and supplied by Liam Bindle.

    Args:
    -----
    x_c : numpy.ndarray
        Grid-box corner longitudes; array shape = (n+1, n+1),
        where n is the cubed-sphere grid size.

    y_c : numpy.ndarray
        Grid-box corner longitudes; array shape = (n+1, n+1),
        where n is the cubed-sphere grid size.

    Returns:
    --------
    x_y : numpy.ndarray
        Grid-box cartesian coordinates; array shape = (n, n, 5),
        where n is the cubed-sphere grid size.
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
                x_c[p_0, p_0]
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
                y_c[p_0, p_0]
            ]
        ),
        0, -1
    )
    return np.moveaxis(
        np.array(
            [boxes_x, boxes_y]
        ), 0, -1
    )


def central_angle(
        x_0,
        y_0,
        x_1,
        y_1):
    """
    Returns the distance (central angle) between cartesian
    coordinates (x_0, y_0) and (x_1, y_1).  This is vectorizable.
    Developed, tested, and supplied by Liam Bindle.

    Args:
    -----
    x_0, y_0 : float
        Longitude and latitude (degrees) of coordinates (x_0, y_0).

    x_1, y_1: float
        Longitude and latitude (degrees) of coordinates (x_1, y_1).

    Returns:
    --------
    distance : float
        Distance (degrees) between (x_0, y_0) and (x_1, y_1).
    """
    x_0 = x_0 * DEG_TO_RAD
    x_1 = x_1 * DEG_TO_RAD
    y_0 = y_0 * DEG_TO_RAD
    y_1 = y_1 * DEG_TO_RAD

    return np.arccos(
        np.sin(y_0) * np.sin(y_1) + \
        np.cos(y_0) * np.cos(y_1) * \
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
    Returns the cubed-sphere grid box corresponding to a given
    latitude and longitude.  Called by routine find_index.

    Args:
    -----
    lat, lon : float or list(float)
        Latitude and longitude (degrees) of the point for which
        cubed-sphere grid indices are desired.

    x_centers_flat, y_centers_flat : float or list(float)
        Flattened (i.e. in Fortran column-major notation) arrays
        of cubed-sphere xDim and yDim values (degrees).

    xy_polygon_defs : float or list(float)
        XY polygon definitions for cubed-sphere grid boxes
        (i.e. the output of function corners_to_xy).

    cs_size : int or list(int)
        Cubed-sphere grid size (i.e. the number of points along
        a face edge).

    latlon_crs : ?
        Ouptut of pyproj.Proj("+proj=latlon")

    jitter_size : float
        If the point cannot be matched to a cubed-sphere grid box,
        then shift longitude by the distance [m] specified in
        jitter_size before doing the lookup once more.  A nonzero
        jitter_size value may be needed when the latitude is close
        to +90 or -90.

    Returns:
    --------
    nf_cs, xdim_cs, ydim_cs : int, list(float), list(float)
        nf_cs is the number of cube-sphere face
        xdim_cs (aka XDim) and ydim_cs (aka YDim) are the longitude
        and latitude arrays for each cell of the cubed-sphere grid.
    """
    # Center on x_find, y_find
    x_find = lon
    y_find = lat
    gnomonic_crs = pyproj.Proj(f'+proj=gnom +lat_0={y_find} +lon_0={x_find}')

    # Generate all distances
    distances = central_angle(
        x_find,
        y_find,
        x_centers_flat,
        y_centers_flat
    )
    four_nearest_indexes = np.argpartition(distances, 4)[:4]

    # Unravel 4 smallest indexes
    four_nearest_indexes = np.unravel_index(
        four_nearest_indexes,
        (6, cs_size, cs_size)
    )
    four_nearest_xy = xy_polygon_defs[four_nearest_indexes]
    four_nearest_polygons = [
        shapely.geometry.Polygon(polygon_xy) for polygon_xy in four_nearest_xy
    ]

    # Transform to gnomonic projection
    gno_transform = pyproj.Transformer.from_proj(
        latlon_crs,
        gnomonic_crs,
        always_xy=True
    ).transform
    four_nearest_polygons_gno = [
        shapely.ops.transform(gno_transform, polygon) \
        for polygon in four_nearest_polygons
    ]

    # Figure out which polygon contains the point
    xy_find = shapely.geometry.Point(x_find, y_find)
    xy_find_gno = shapely.ops.transform(gno_transform, xy_find)
    polygon_contains_point = [
        polygon.contains(xy_find_gno) for polygon in four_nearest_polygons_gno
    ]

    # If the point cannot be matched (such as can happen near the poles),
    # move the longitude by the jitter_size (in meters) and try again.
    if np.count_nonzero(polygon_contains_point) == 0:
        if jitter_size > 0.0:
            nf_cs, ydim_cs, xdim_cs = find_index_single(
                y_find,
                x_find+jitter_size,
                x_centers_flat,
                y_centers_flat,
                xy_polygon_defs,
                cs_size,
                latlon_crs,
                jitter_size=0.0
            )
        else:
            msg = f'Point at {x_find:8.2f} E, {y_find:8.2f} N '
            msg+= 'could not be matched'
            raise ValueError(msg)

    # The first will be selected, if more than one
    polygon_with_point = np.argmax(polygon_contains_point)

    # Get original index
    nf_cs   = four_nearest_indexes[0][polygon_with_point]
    ydim_cs = four_nearest_indexes[1][polygon_with_point]
    xdim_cs = four_nearest_indexes[2][polygon_with_point]

    return nf_cs, ydim_cs, xdim_cs


def find_index(
        lat,
        lon,
        grid,
        jitter_size=0.0
):
    """
    Returns the cubed-sphere grid box indices corresponding to
    given latitude and longitude coordinates.

    Based on a routine developed, tested, and supplied by Liam Bindle.

    Args:
    -----
    lat, lon : float
        Latitude and longitude (degrees) of the point for which
        cubed-sphere indices are desired.

    grid : xarray.Dataset
        Cubed-sphere grid definition with the following variables:
           {'lat'   : lat midpoints,
            'lon'   : lon midpoints,
            'lat_b' : lat edges,
            'lon_b' : lon edges}
        where each value has an extra face dimension of length 6.

    Keyword Args (optional):
    ------------------------
    jitter_size : float
        If the point cannot be matched to a cubed-sphere grid box,
        then shift longitude by the distance [m] specified in
        jitter_size before doing the lookup once more.  A nonzero
        jitter_size value may be needed when the latitude is close
        to +90 or -90.  Default value: 0

    Returns:
    --------
    ind : numpy.ndarray
        Array containing (nf, YDim, XDim), where:
            nf is the number of cubed-sphere faces (= 6)
            YDim is the cubed-sphere longitude index at (lat, lon)
            XDim is the cubed-sphere latitude index at (lat, lon)
    """
    gcpy.util.verify_variable_type(grid, xr.Dataset)

    lon_vec = np.asarray(lon)
    lat_vec = np.asarray(lat)
    n_find = lon_vec.size

    # Get the corners
    x_corners = grid['lon_b'].values
    y_corners = grid['lat_b'].values
    x_centers = grid['lon'].values
    y_centers = grid['lat'].values
    x_centers_flat = x_centers.flatten()
    y_centers_flat = y_centers.flatten()

    cs_size = x_centers.shape[-1]

    # Generate everything that will be reused
    # Get XY polygon definitions for grid boxes
    # 5 (x,y) points defining polygon corners (first and last are same)
    xy_polygon_defs = np.zeros((6, cs_size, cs_size, 5, 2))
    for nf_cs in range(6):
        xy_polygon_defs[nf_cs, ...] = corners_to_xy(
            x_c=x_corners[nf_cs, :, :],
            y_c=y_corners[nf_cs, :, :]
        )
    latlon_crs = pyproj.Proj("+proj=latlon")

    # Find 4 shortest distances to (x_find, y_find)
    idx = np.full((3,n_find), 0)
    for x_find, y_find, i_find in \
        zip(np.nditer(lon_vec), np.nditer(lat_vec), list(range(n_find))):

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
        idx[:,i_find] = [nf_cs, ydim_cs, xdim_cs]

    return idx


def is_cubed_sphere(
        data
):
    """
    Given an xarray Dataset or DataArray object, determines if the
    data is placed on a cubed-sphere grid.

    Args:
    -----
    data : xarray.Dataset or xarray.DataArray
        The input data to be tested

    Returns:
    --------
    is_gchp : bool
        Returns True if data is placed on a cubed-sphere grid,
        and False otherwise.

    Remarks:
    --------
    A cubed-sphere data file has one of the following attributes
    (1) A dimension named "nf" (GCHP/GEOS diagnostic files)
    (2) The lat/lon ratio is exactly 6 (GCHP/GEOS checkpoints)
    """
    gcpy.util.verify_variable_type(data, (xr.DataArray, xr.Dataset))

    if is_cubed_sphere_diag_grid(data):
        return True
    if is_cubed_sphere_rst_grid(data):
        return True
    return False


def is_cubed_sphere_diag_grid(data):
    """
    Determines if a cubed-sphere grid has History file dimensions.
    (i.e. a dimension named "nf", aka number of grid faces).

    Args:
    -----
    data : xarray.DataArray or xarray.Dataset
        The input data.

    Returns:
    --------
    True if the grid has History diagnostic dimensions,
    False otherwise.
    """
    return "nf" in data.dims


def is_cubed_sphere_rst_grid(data):
    """
    Determines if a cubed-sphere grid has restart file dimensions.
    (i.e. lat and lon, with lat = lon*6).

    Args:
    -----
    data : xarray.DataArray or xarray.Dataset
        The input data.

    Returns:
    --------
    True if the grid has restart dimensions, False otherwise.
    """
    gcpy.util.verify_variable_type(data, (xr.DataArray, xr.Dataset))

    # TODO: Rethink this if we ever end up changing the GC-Classic
    # restart variables to start with SPC, or if we ever rename the
    # internal state variables in GCHP. A more robust back-up check
    # could be to see if all the lats and lons are integer, since
    # that will be the case with the GCHP restart file format.

    # NOTE: in DataArray objects, dims is a tuple but not a dict!
    # Comparing the len of the lat & lon coords will work for both.
    if "lat" in data.coords:
         return len(data.coords["lat"]) == len(data.coords["lon"]) * 6

    # Dataset: search data.data_vars for "SPC_"
    # DataArray: search data.name for "SPC_"
    if isinstance(data, xr.Dataset):
        return "SPC_" in data.data_vars.keys()
    return "SPC_" in data.name


def get_cubed_sphere_res(data):
    """
    Given a Dataset or DataArray object, returns the number of
    grid cells along one side of the cubed-sphere grid face
    (e.g. 24 for grid resolution C24, which has 24x25 grid cells
    per face).

    Args:
    -----
    data : xarray.DataArray or xarray.Dataset
        The input data.

    Returns:
    --------
    cs_res : int
        The cubed-sphere resolution.  Will return 0 if the data
        is not placed on a cubed-sphere grid.
    """
    gcpy.util.verify_variable_type(data, (xr.DataArray, xr.Dataset))

    if not is_cubed_sphere(data):
        return 0

    # NOTE: In Dataset objects "dims" is a dict, but in DataArray
    # objects "dims" is a tuple.  Returning the length of the
    # corresponding coords array should work in both cases.
    if is_cubed_sphere_rst_grid(data):
        return len(data.coords["lon"])
    return len(data.coords["Xdim"])


def is_gchp_lev_positive_down(data):
    """
    Determines if GCHP data is arranged vertically from the top of the
    atmosphere downwards or from the surface upwards, according to:

    (1) Checkpoint files:     lev:positive="down
    (2) Emissions collection: lev:positive="down"
    (3) Other collections     lev:positive="up"

    Args:
    -----
    data : xarray.DataArray or xarray.Dataset
       The input data

    Returns:
    --------
    True if the data is arranged from top-of-atm downwards.
    False if the data is arranged from the surface upwards.
    """
    gcpy.util.verify_variable_type(data, (xr.DataArray, xr.Dataset))

    if is_cubed_sphere_rst_grid(data):
        return True
    if is_cubed_sphere_diag_grid(data):
        emis_vars = [var for var in data.data_vars if var.startswith("Emis")]
        if len(emis_vars) > 0:
            return True
    return False
