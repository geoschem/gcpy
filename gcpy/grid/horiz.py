""" Horizontal grid definitions and helper functions """

import numpy as np
from numpy import asarray
from itertools import product

INV_SQRT_3 = 1.0 / np.sqrt(3.0)
ASIN_INV_SQRT_3 = np.arcsin(INV_SQRT_3)

def make_grid_LL(llres, minlon=-180, maxlon=175, minlat=-90, maxlat=90):
    [dlat,dlon] = list(map(float, llres.split('x')))
    #temporary measure until automatic reading of coordiante descriptions is in place
    #maxlon = maxlon + dlon
    #lon_b = np.linspace(minlon - dlon/2, maxlon - dlon/2, int((abs(minlon)+abs(maxlon))/dlon) + 1, endpoint=True)
    #lat_b = np.linspace(minlat - dlat/2, maxlat + dlat/2, 
    #                    int((abs(minlat)+abs(maxlat))/dlat) + 2, endpoint=True).clip(minlat,maxlat)
    lon_b = np.linspace(-180 - dlon/2, 180 - dlon/2, int(360/dlon) + 1, endpoint=True)
    lat_b = np.linspace(-90 - dlat/2, 90 + dlat/2, int(180/dlat) + 2, endpoint=True).clip(-90,90)
    lat = (lat_b[1:] + lat_b[:-1]) / 2
    lon = (lon_b[1:] + lon_b[:-1]) / 2
    llgrid = {'lat': lat, 
              'lon': lon, 
              'lat_b': lat_b, 
              'lon_b': lon_b}
    return llgrid

def make_grid_CS(csres):
    csgrid = csgrid_GMAO(csres)
    csgrid_list = [None]*6
    for i in range(6):
        csgrid_list[i] = {'lat': csgrid['lat'][i], 
                          'lon': csgrid['lon'][i],
                          'lat_b': csgrid['lat_b'][i], 
                          'lon_b': csgrid['lon_b'][i]}
    return [csgrid, csgrid_list]

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

    Parameters
    ----------
    lat_stride : float
        Stride length in degrees. For example, for a standard GEOS-Chem Classic
        4x5 grid, lat_stride would be 4.

    half_polar_grid: boolean
        Whether or not the grid should be "half-polar" (i.e. bands at poles are
        half the size). In either case the grid will start and end at -/+ 90,
        but when half_polar_grid is True, the first and last bands will have a
        width of 1/2 the normal lat_stride.

    Returns
    -------
    Latitudes of cell edges in degrees North.

    Notes
    -----
    All values are forced to be between [-90,90]. For a grid with N cells in
    each band, N+1 edges will be returned, with the first and last value being
    duplicates.

    Examples
    --------
    >>> from gcpy.grid.horiz import calc_rectilinear_lat_edge
    >>> calc_rectilinear_lat_edge(4.0,true)
    np.array([-90,-88,-84,-80,...,84,88,90])

    See Also
    --------
    [NONE]

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

    Parameters
    ----------
    #TODO

    Returns
    -------
    #TODO

    Notes
    -----
    #TODO

    Examples
    --------
    #TODO

    See Also
    --------
    [NONE]

    """
    from .. constants import R_EARTH

    # Convert from km to m
    _radius_earth_m = R_EARTH * 1000.0

    lon_edge = asarray(lon_edge, dtype=float)
    lat_edge = asarray(lat_edge, dtype=float)

    n_lon = (lon_edge.size) - 1
    n_lat = (lat_edge.size) - 1

    grid_area = np.zeros((n_lat,n_lon))

    sfc_area_const = 2.0*np.pi*_radius_earth_m*_radius_earth_m

    # Longitudes loop, so need to be careful
    lon_delta = calc_delta_lon(lon_edge)

    # Convert into weights relative to the total circle
    lon_delta = lon_delta/360.0

    # Precalculate this
    sin_lat_edge = np.sin(np.deg2rad(lat_edge))

    for i_lat in range(0,n_lat):
        sin_diff = sin_lat_edge[i_lat+1] - sin_lat_edge[i_lat]
        grid_area[i_lat,:] = sin_diff * sfc_area_const * lon_delta

    return grid_area

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

def csgrid_GMAO(res, offset=-10):
    """
    Return cubedsphere coordinates with GMAO face orientation

    Parameters
    ----------
    res : cubed-sphere Resolution

    This function was originally written by Jiawei Zhuange and included 
    in package cubedsphere: https://github.com/JiaweiZhuang/cubedsphere

    """

    CS = CSGrid(res, offset=offset)

    lon = CS.lon_center.transpose(2, 0, 1)
    lon_b = CS.lon_edge.transpose(2, 0, 1)
    lat = CS.lat_center.transpose(2, 0, 1)
    lat_b = CS.lat_edge.transpose(2, 0, 1)

    lon[lon < 0] += 360
    lon_b[lon_b < 0] += 360

    for a in [lon, lon_b, lat, lat_b]:

        for tile in [0, 1, 3, 4]:
            a[tile] = a[tile].T
        for tile in [3, 4]:
            a[tile] = np.flip(a[tile], 1)
        for tile in [3, 4, 2, 5]:
            a[tile] = np.flip(a[tile], 0)

        a[2], a[5] = a[5].copy(), a[2].copy()  # swap north&south pole

    return {'lon': lon, 'lat': lat, 'lon_b': lon_b, 'lat_b': lat_b}


class CSGrid(object):
    """Generator for cubed-sphere grid geometries.

    CSGrid computes the latitutde and longitudes of cell centers and edges
    on a cubed-sphere grid, providing a way to retrieve these geometries
    on-the-fly if your model output data does not include them.

    Attributes
    ----------
    {lon,lat}_center : np.ndarray
        lat/lon coordinates for each cell center along the cubed-sphere mesh
    {lon,lat}_edge : np.ndarray
        lat/lon coordinates for the midpoint of the edges separating each
        element on the cubed-sphere mesh.
    xyz_{center,edge} : np.ndarray
        As above, except coordinates are projected into a 3D cartesian space
        with common origin to the original lat/lon coordinate system, assuming
        a unit sphere.

    This class was originally written by Jiawei Zhuange and included 
    in package cubedsphere: https://github.com/JiaweiZhuang/cubedsphere

    """

    def __init__(self, c, offset=None):
        """
        Parameters
        ----------
        c : int
            Number edges along each cubed-sphere edge.
            ======= ====================
               C    Lat/Lon Resolution
            ------- --------------------
               24    4 deg x 5 deg
             48,45   2 deg x 2.5 deg
             96,90   1 deg x 1.25 deg
            192,180  0.5 deg x 0.625 deg
            384,360  0.25 deg x 0.3125 deg
              720    0.12g deg x 0.15625 deg
        offset : float (optional)
            Degrees to offset the first faces' edge in the latitudinal
            direction. If not passed, then the western edge of the first face
            will align with the prime meridian.

       This function was originally written by Jiawei Zhuange and included 
       in package cubedsphere: https://github.com/JiaweiZhuang/cubedsphere

        """
        self.c = c
        self.delta_y = 2. * ASIN_INV_SQRT_3 / c
        self.nx = self.ny = c + 1
        self.offset = offset

        self._initialize()

    def _initialize(self):

        c = self.c
        nx, ny = self.nx, self.ny

        lambda_rad = np.zeros((nx, ny))
        lambda_rad[ 0, :] = 3.*np.pi/4. # West edge
        lambda_rad[-1, :] = 5.*np.pi/4. # East edge

        theta_rad = np.zeros((nx, ny))
        theta_rad[ 0, :] = -ASIN_INV_SQRT_3 + (self.delta_y*np.arange(c+1)) # West edge
        theta_rad[-1, :] = theta_rad[0, :] # East edge

        # Cache the reflection points - our upper-left and lower-right corners
        lonMir1, lonMir2 = lambda_rad[0, 0], lambda_rad[-1, -1]
        latMir1, latMir2 = theta_rad[0, 0], theta_rad[-1, -1]

        xyzMir1 = latlon_to_cartesian(lonMir1, latMir1)
        xyzMir2 = latlon_to_cartesian(lonMir2, latMir2)

        xyzCross = np.cross(xyzMir1, xyzMir2)
        norm = np.sqrt(np.sum(xyzCross**2))
        xyzCross /= norm

        for i in range(1, c):

            lonRef, latRef = lambda_rad[0, i], theta_rad[0, i]
            xyzRef = np.asarray(latlon_to_cartesian(lonRef, latRef, ))

            xyzDot = np.sum(xyzCross*xyzRef)
            xyzImg = xyzRef - (2. * xyzDot * xyzCross)

            xsImg, ysImg, zsImg = xyzImg
            lonImg, latImg = cartesian_to_latlon(xsImg, ysImg, zsImg)

            lambda_rad[i, 0] = lonImg
            lambda_rad[i, -1] = lonImg
            theta_rad[i, 0] = latImg
            theta_rad[i, -1] = -latImg

        pp = np.zeros([3, c+1, c+1])

        # Set the four corners
        # print("CORNERS")
        for i, j in product([0, -1], [0, -1]):
            # print(i, j)
            pp[:, i, j] = latlon_to_cartesian(lambda_rad[i, j], theta_rad[i, j])

        # Map the edges on the sphere back to the cube. Note that all intersections are at x = -rsq3
        # print("EDGES")
        for ij in range(1, c+1):
            # print(ij)
            pp[:, 0, ij] = latlon_to_cartesian(lambda_rad[0, ij], theta_rad[0, ij])
            pp[1, 0, ij] = -pp[1, 0, ij] * INV_SQRT_3 / pp[0, 0, ij]
            pp[2, 0, ij] = -pp[2, 0, ij] * INV_SQRT_3 / pp[0, 0, ij]

            pp[:, ij, 0] = latlon_to_cartesian(lambda_rad[ij, 0], theta_rad[ij, 0])
            pp[1, ij, 0] = -pp[1, ij, 0] * INV_SQRT_3 / pp[0, ij, 0]
            pp[2, ij, 0] = -pp[2, ij, 0] * INV_SQRT_3 / pp[0, ij, 0]

        # # Map interiors
        pp[0, :, :] = -INV_SQRT_3
        # print("INTERIOR")
        for i in range(1, c+1):
            for j in range(1, c+1):
                # Copy y-z face of the cube along j=1
                pp[1, i, j] = pp[1, i, 0]
                # Copy along i=1
                pp[2, i, j] = pp[2, 0, j]

        _pp = pp.copy()
        llr, ttr = vec_cartesian_to_latlon(_pp[0], _pp[1], _pp[2])

        lambda_rad, theta_rad = llr.copy(), ttr.copy()

        # Make grid symmetrical to i = im/2 + 1
        for j in range(1, c+1):
            for i in range(1, c+1):
                # print("({}, {}) -> ({}, {})".format(i, 0, i, j))
                lambda_rad[i, j] = lambda_rad[i, 0]

        for j in range(c+1):
            for i in range(c//2):
                isymm = c - i
                # print(isymm)
                avgPt = 0.5*(lambda_rad[i, j] - lambda_rad[isymm, j])
                # print(lambda_rad[i, j], lambda_rad[isymm, j], avgPt)
                lambda_rad[i, j] = avgPt + np.pi
                lambda_rad[isymm, j] = np.pi - avgPt

                avgPt = 0.5*(theta_rad[i, j] + theta_rad[isymm, j])
                theta_rad[i, j] = avgPt
                theta_rad[isymm, j] = avgPt

        # Make grid symmetrical to j = im/2 + 1
        for j in range(c//2):
            jsymm = c - j
            for i in range(1, c+1):
                avgPt = 0.5*(lambda_rad[i, j] + lambda_rad[i, jsymm])
                lambda_rad[i, j] = avgPt
                lambda_rad[i, jsymm] = avgPt

                avgPt = 0.5*(theta_rad[i, j] - theta_rad[i, jsymm])
                theta_rad[i, j] = avgPt
                theta_rad[i, jsymm] = -avgPt

        # Final correction
        lambda_rad -= np.pi

        llr, ttr = lambda_rad.copy(), theta_rad.copy()

        #######################################################################
        ## MIRROR GRIDS
        #######################################################################

        new_xgrid = np.zeros((c+1, c+1, 6))
        new_ygrid = np.zeros((c+1, c+1, 6))

        xgrid = llr.copy()
        ygrid = ttr.copy()

        new_xgrid[..., 0] = xgrid.copy()
        new_ygrid[..., 0] = ygrid.copy()

        # radius = 6370.0e3
        radius = 1.

        for face in range(1, 6):
            for j in range(c+1):
                for i in range(c+1):
                    x = xgrid[i, j]
                    y = ygrid[i, j]
                    z = radius

                    if face == 1:
                        # Rotate about z only
                        new_xyz = rotate_sphere_3D(x, y, z, -np.pi/2., 'z')

                    elif face == 2:
                        # Rotate about z, then x
                        temp_xyz = rotate_sphere_3D(x, y, z, -np.pi/2., 'z')
                        x, y, z = temp_xyz[:]
                        new_xyz = rotate_sphere_3D(x, y, z, np.pi/2., 'x')

                    elif face == 3:
                        temp_xyz = rotate_sphere_3D(x, y, z, np.pi, 'z')
                        x, y, z = temp_xyz[:]
                        new_xyz = rotate_sphere_3D(x, y, z, np.pi/2., 'x')

                        if ((c % 2) != 0) and (j == c//2 - 1):
                            print(i, j, face)
                            new_xyz[0] = np.pi

                    elif face == 4:
                        temp_xyz = rotate_sphere_3D(x, y, z, np.pi/2., 'z')
                        x, y, z = temp_xyz[:]
                        new_xyz = rotate_sphere_3D(x, y, z,  np.pi/2., 'y')

                    elif face == 5:
                        temp_xyz = rotate_sphere_3D(x, y, z,  np.pi/2., 'y')
                        x, y, z = temp_xyz[:]
                        new_xyz = rotate_sphere_3D(x, y, z, 0., 'z')

                    # print((x, y, z), "\n", new_xyz, "\n" + "--"*40)

                    new_x, new_y, _ = new_xyz
                    new_xgrid[i, j, face] = new_x
                    new_ygrid[i, j, face] = new_y

        lon_edge, lat_edge = new_xgrid.copy(), new_ygrid.copy()

        #######################################################################
        ## CLEANUP GRID
        #######################################################################

        for i, j, f in product(range(c+1), range(c+1), range(6)):
            new_lon = lon_edge[i, j, f]
            if new_lon < 0:
                new_lon+= 2*np.pi
            if np.abs(new_lon) < 1e-10:
                new_lon = 0.
            lon_edge[i, j, f] = new_lon

            if np.abs(lat_edge[i, j, f]) < 1e-10:
                lat_edge[i, j, f] = 0.

        lon_edge_deg = np.rad2deg(lon_edge)
        lat_edge_deg = np.rad2deg(lat_edge)

        #######################################################################
        ## COMPUTE CELL CENTROIDS
        #######################################################################

        lon_ctr = np.zeros((c, c, 6))
        lat_ctr = np.zeros((c, c, 6))
        xyz_ctr = np.zeros((3, c, c, 6))
        xyz_edge = np.zeros((3, c+1, c+1, 6))

        for f in range(6):
            for i in range(c):
                last_x = (i == (c-1))
                for j in range(c):
                    last_y = (j == (c-1))

                    # Get the four corners
                    lat_corner = [lat_edge[  i,   j, f], lat_edge[i+1,   j, f],
                                  lat_edge[i+1, j+1, f], lat_edge[  i, j+1, f]]
                    lon_corner = [lon_edge[  i,   j, f], lon_edge[i+1,   j, f],
                                  lon_edge[i+1, j+1, f], lon_edge[  i, j+1, f]]

                    # Convert from lat-lon back to cartesian
                    xyz_corner = np.asarray(vec_latlon_to_cartesian(lon_corner, lat_corner))

                    # Store the edge information
                    xyz_edge[:, i, j, f] = xyz_corner[:, 0]
                    if last_x:
                        xyz_edge[:, i+1, j, f] = xyz_corner[:, 1]
                    if last_x or last_y:
                        xyz_edge[:, i+1, j+1, f] = xyz_corner[:, 2]
                    if last_y:
                        xyz_edge[:, i, j+1, f] = xyz_corner[:, 3]

                    e_mid = np.sum(xyz_corner, axis=1)
                    e_abs = np.sqrt(np.sum(e_mid * e_mid))
                    if e_abs > 0:
                        e_mid = e_mid / e_abs

                    xyz_ctr[:, i, j, f] = e_mid
                    _lon, _lat = cartesian_to_latlon(*e_mid)
                    lon_ctr[i, j, f] = _lon
                    lat_ctr[i, j, f] = _lat

        lon_ctr_deg = np.rad2deg(lon_ctr)
        lat_ctr_deg = np.rad2deg(lat_ctr)

        if self.offset is not None:
            lon_edge_deg += self.offset
            lon_ctr_deg += self.offset

        #######################################################################
        ## CACHE
        #######################################################################

        self.lon_center = lon_ctr_deg
        self.lat_center = lat_ctr_deg

        self.lon_edge = lon_edge_deg
        self.lat_edge = lat_edge_deg

        self.xyz_center = xyz_ctr
        self.xyz_edge = xyz_edge

def latlon_to_cartesian(lon, lat):
    """ Convert latitude/longitude coordinates along the unit sphere to cartesian
    coordinates defined by a vector pointing from the sphere's center to its
    surface.

    This function was originally written by Jiawei Zhuange and included 
    in package cubedsphere: https://github.com/JiaweiZhuang/cubedsphere

    """

    x = np.cos(lat) * np.cos(lon)
    y = np.cos(lat) * np.sin(lon)
    z = np.sin(lat)

    return x, y, z
vec_latlon_to_cartesian = np.vectorize(latlon_to_cartesian)


def cartesian_to_latlon(x, y, z, ret_xyz=False):
    """ Convert a cartesian coordinate to latitude/longitude coordinates.

    Optionally return the original cartesian coordinate as a tuple.

    This function was originally written by Jiawei Zhuange and included 
    in package cubedsphere: https://github.com/JiaweiZhuang/cubedsphere

    """

    xyz = np.array([x, y, z])
    vector_length = np.sqrt(np.sum(xyz*xyz, axis=0))
    xyz /= vector_length
    x, y, z = xyz

    if (np.abs(x) + np.abs(y)) < 1e-20:
        lon = 0.
    else:
        lon = np.arctan2(y, x)
    if lon < 0.:
        lon += 2*np.pi

    lat = np.arcsin(z)
    # If not normalizing vector, take lat = np.arcsin(z/vector_length)

    if ret_xyz:
        return lon, lat, xyz
    else:
        return lon, lat
vec_cartesian_to_latlon = np.vectorize(cartesian_to_latlon)


def spherical_to_cartesian(theta, phi, r=1):
    """ Convert spherical coordinates in the form (theta, phi[, r]) to
    cartesian, with the origin at the center of the original spherical
    coordinate system.

    This function was originally written by Jiawei Zhuange and included 
    in package cubedsphere: https://github.com/JiaweiZhuang/cubedsphere

    """
    x = r * np.cos(phi) * np.cos(theta)
    y = r * np.cos(phi) * np.sin(theta)
    z = r * np.sin(phi)
    return x, y, z
vec_spherical_to_cartesian = np.vectorize(spherical_to_cartesian)


def cartesian_to_spherical(x, y, z):
    """ Convert cartesian coordinates to spherical in the form
    (theta, phi[, r]) with the origin remaining at the center of the
    original spherical coordinate system.

    This function was originally written by Jiawei Zhuange and included 
    in package cubedsphere: https://github.com/JiaweiZhuang/cubedsphere

    """
    r = np.sqrt(x**2 + y**2 + z**2)
    #theta = np.arccos(z / r)
    theta = np.arctan2(y, x)
    phi = np.arctan2(z, np.sqrt(x**2 + y**2))

    # if np.abs(x) < 1e-16:
    #     phi = np.pi
    # else:
    #     phi = np.arctan(y / x)
    return theta, phi, r
vec_cartesian_to_spherical = np.vectorize(cartesian_to_spherical)


def rotate_sphere_3D(theta, phi, r, rot_ang, rot_axis='x'):
    """ Rotate a spherical coordinate in the form (theta, phi[, r])
    about the indicating axis, 'rot_axis'.

    This method accomplishes the rotation by projecting to a
    cartesian coordinate system and performing a solid body rotation
    around the requested axis.

    This function was originally written by Jiawei Zhuange and included 
    in package cubedsphere: https://github.com/JiaweiZhuang/cubedsphere

    """
    cos_ang = np.cos(rot_ang)
    sin_ang = np.sin(rot_ang)

    x, y, z = spherical_to_cartesian(theta, phi, r)
    if rot_axis == 'x':
        x_new = x
        y_new = cos_ang*y + sin_ang*z
        z_new = -sin_ang*y + cos_ang*z
    elif rot_axis == 'y':
        x_new = cos_ang*x - sin_ang*z
        y_new = y
        z_new = sin_ang*x + cos_ang*z
    elif rot_axis == 'z':
        x_new = cos_ang*x + sin_ang*y
        y_new = -sin_ang*x + cos_ang*y
        z_new = z

    theta_new, phi_new, r_new = cartesian_to_spherical(x_new, y_new, z_new)

    return theta_new, phi_new, r_new
