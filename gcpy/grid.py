import numpy as np
import xarray as xr
from numpy import asarray
import scipy.sparse
from itertools import product
from .util import get_shape_of_data
from .grid_stretching_transforms import scs_transform
from .constants import R_EARTH_m


def get_troposphere_mask(ds):
    """
    Returns a mask array for picking out the tropospheric grid boxes.

    Args:
        ds: xarray Dataset
            Dataset containing certain met field variables (i.e.
            Met_TropLev, Met_BXHEIGHT).

    Returns:
        tropmask: numpy ndarray
            Tropospheric mask.  False denotes grid boxes that are
            in the troposphere and True in the stratosphere
            (as per Python masking logic).
    """

    # ==================================================================
    # Initialization
    # ==================================================================

    # Make sure ds is an xarray Dataset object
    if not isinstance(ds, xr.Dataset):
        raise TypeError("The ds argument must be an xarray Dataset!")

    # Make sure certain variables are found
    if "Met_BXHEIGHT" not in ds.data_vars.keys():
        raise ValueError("Met_BXHEIGHT could not be found!")
    if "Met_TropLev" not in ds.data_vars.keys():
        raise ValueError("Met_TropLev could not be found!")

    # Mask of tropospheric grid boxes in the Ref dataset
    shape = get_shape_of_data(np.squeeze(ds["Met_BXHEIGHT"]))

    # Determine if this is GCHP data
    is_gchp = "nf" in ds["Met_BXHEIGHT"].dims

    # ==================================================================
    # Create the mask arrays for the troposphere
    #
    # Convert the Met_TropLev DataArray objects to numpy ndarrays of
    # integer.  Also subtract 1 to convert from Fortran to Python
    # array index notation.
    # ==================================================================

    multi_time_slices = (is_gchp and len(shape) == 5) or \
                        (not is_gchp and len(shape) == 4)

    if multi_time_slices:
        # --------------------------------------------------------------
        # GCC: There are multiple time slices
        # --------------------------------------------------------------

        # Create the tropmask array with dims
        # (time, lev, nf*lat*lon) for GCHP, or
        # (time, lev, lat*lon   ) for GCC
        tropmask = np.ones((shape[0], shape[1],
                            np.prod(np.array(shape[2:]))), bool)

        # Loop over each time
        for t in range(tropmask.shape[0]):

            # Pick the tropopause level and make a 1-D array
            values = ds["Met_TropLev"].isel(time=t).values
            lev = np.int_(np.squeeze(values) - 1)
            lev_1d = lev.flatten()

            # Create the tropospheric mask array
            for x in range(tropmask.shape[2]):
                tropmask[t, 0: lev_1d[x], x] = False

    else:
        # --------------------------------------------------------------
        # There is only one time slice
        # --------------------------------------------------------------

        # Create the tropmask array with dims (lev, lat*lon)
        tropmask = np.ones((shape[0], np.prod(np.array(shape[1:]))), bool)

        # Pick the tropopause level and make a 1-D array
        values = ds["Met_TropLev"].values
        lev = np.int_(np.squeeze(values) - 1)
        lev_1d = lev.flatten()

        # Create the tropospheric mask array
        for x in range(tropmask.shape[1]):
            tropmask[0: lev_1d[x], x] = False

    # Reshape into the same shape as Met_BxHeight
    return tropmask.reshape(shape)


def get_input_res(data):
    """
    Returns resolution of dataset passed to compare_single_level or compare_zonal_means

    Args:
        data: xarray Dataset
            Input GEOS-Chem dataset

    Returns:
        res: str or int
            Lat/lon res of the form 'latresxlonres' or cubed-sphere resolution
        gridtype: str
            'll' for lat/lon or 'cs' for cubed-sphere

    """
    vdims = data.dims
    if "lat" in vdims and "lon" in vdims:
        lat = data["lat"].values
        lon = data["lon"].values
        if lat.size / 6 == lon.size:
            return lon.size, "cs"
        else:
            lat.sort()
            lon.sort()
            # use increment of second and third coordinates
            # to avoid polar mischief
            lat_res = np.abs(lat[2] - lat[1])
            lon_res = np.abs(lon[2] - lon[1])
            return str(lat_res) + "x" + str(lon_res), "ll"

    else:
        #print("grid is cs: ", vdims)
        # GCHP data using MAPL v1.0.0+ has dims time, lev, nf, Ydim, and Xdim
        if isinstance(data.dims, tuple):
            return len(data["Xdim"].values), "cs"
        else:
            return data.dims["Xdim"], "cs"


def call_make_grid(res, gridtype, in_extent=[-180, 180, -90, 90],
                   out_extent=[-180, 180, -90, 90], sg_params=[1, 170, -90]):
    """
    Create a mask with NaN values removed from an input array

    Args:
        res: str or int
            Resolution of grid (format 'latxlon' or csres)
        gridtype: str
            'll' for lat/lon or 'cs' for cubed-sphere

    Keyword Args (optional):
        in_extent: list[float, float, float, float]
            Describes minimum and maximum latitude and longitude of input data
            in the format [minlon, maxlon, minlat, maxlat]
            Default value: [-180, 180, -90, 90]
        out_extent: list[float, float, float, float]
            Desired minimum and maximum latitude and longitude of output grid
            in the format [minlon, maxlon, minlat, maxlat]
            Default value: [-180, 180, -90, 90]
        sg_params: list[float, float, float] (stretch_factor, target_longitude, target_latitude)
            Desired stretched-grid parameters in the format
            [stretch_factor, target_longitude, target_latitude].
            Will trigger stretched-grid creation if not default values.
            Default value: [1, 170, -90] (no stretching)

    Returns:
        [grid, grid_list]: list(dict, list(dict))
            Returns the created grid. 
            grid_list is a list of grids if gridtype is 'cs', else it is None
    """

    # call appropriate make_grid function and return new grid
    if gridtype == "ll":
        return [make_grid_LL(res, in_extent, out_extent), None]
    elif sg_params == [1, 170, -90]:
        # standard CS
        return make_grid_CS(res)
    else:
        return make_grid_SG(res, *sg_params)


def get_grid_extents(data, edges=True):
    """
    Get min and max lat and lon from an input GEOS-Chem xarray dataset or grid dict

    Args:
        data: xarray Dataset or dict
            A GEOS-Chem dataset or a grid dict
        edges (optional): bool
            Whether grid extents should use cell edges instead of centers
            Default value: True

    Returns:
        minlon: float
            Minimum longitude of data grid
        maxlon: float
            Maximum longitude of data grid
        minlat: float
            Minimum latitude of data grid
        maxlat: float
            Maximum latitude of data grid
    """

    if isinstance(data, dict):
        if "lon_b" in data and edges:
            return np.min(
                data["lon_b"]), np.max(
                data["lon_b"]), np.min(
                data["lat_b"]), np.max(
                data["lat_b"])
        elif not edges:
            return np.min(
                data["lon"]), np.max(
                data["lon"]), np.min(
                data["lat"]), np.max(
                data["lat"])
        else:
            return -180, 180, -90, 90
    elif "lat" in data.dims and "lon" in data.dims:
        lat = data["lat"].values
        lon = data["lon"].values
        if lat.size / 6 == lon.size:
            # No extents for CS plots right now
            return -180, 180, -90, 90
        else:
            lat = np.sort(lat)
            minlat = np.min(lat)
            if abs(abs(lat[1]) - abs(lat[0])
                   ) != abs(abs(lat[2]) - abs(lat[1])):
                #pole is cutoff
                minlat = minlat - 1
            maxlat = np.max(lat)
            if abs(abs(lat[-1]) - abs(lat[-2])
                   ) != abs(abs(lat[-2]) - abs(lat[-3])):
                maxlat = maxlat + 1
            # add longitude res to max longitude
            lon = np.sort(lon)
            minlon = np.min(lon)
            maxlon = np.max(lon) + abs(abs(lon[-1]) - abs(lon[-2]))
            return minlon, maxlon, minlat, maxlat
    else:
        # GCHP data using MAPL v1.0.0+ has dims time, lev, nf, Ydim, and Xdim
        return -180, 180, -90, 90


def get_vert_grid(dataset, AP=[], BP=[]):
    """
    Determine vertical grid of input dataset

    Args:
        dataset: xarray Dataset
            A GEOS-Chem output dataset

    Keyword Args (optional):
        AP: list-like type
            Hybrid grid parameter A in hPa
            Default value: []
        BP: list-like type
            Hybrid grid parameter B (unitless)
            Default value: []

    Returns:
        p_edge: numpy array
            Edge pressure values for vertical grid
        p_mid: numpy array
            Midpoint pressure values for vertical grid
        nlev: int
            Number of levels in vertical grid
    """

    if dataset.sizes["lev"] in (72, 73):
        return GEOS_72L_grid.p_edge(), GEOS_72L_grid.p_mid(), 72
    elif dataset.sizes["lev"] in (47, 48):
        return GEOS_47L_grid.p_edge(), GEOS_47L_grid.p_mid(), 47
    elif AP == [] or BP == []:
        if dataset.sizes["lev"] == 1:
            AP = [1, 1]
            BP = [1]
            new_grid = vert_grid(AP, BP)
            return new_grid.p_edge(), new_grid.p_mid(), np.size(AP)
        else:
            raise ValueError(
                "Only 72/73 or 47/48 level vertical grids are automatically determined" +
                "from input dataset by get_vert_grid(), please pass grid parameters AP and BP" +
                "as keyword arguments")
    else:
        new_grid = vert_grid(AP, BP)
        return new_grid.p_edge(), new_grid.p_mid(), np.size(AP)


def get_pressure_indices(pedge, pres_range):
    """
    Get indices where edge pressure values are within a given pressure range

    Args:
        pedge: numpy array
            A GEOS-Chem output dataset
        pres_range: list(float, float)
            Contains minimum and maximum pressure

    Returns:
        numpy array
            Indices where edge pressure values are within a given pressure range
    """

    return np.where(
        (pedge <= np.max(pres_range)) & (
            pedge >= np.min(pres_range)))[0]


def pad_pressure_edges(pedge_ind, max_ind, pmid_len):
    """
    Add outer indices to edge pressure index list

    Args:
        pedge_ind: list
            List of edge pressure indices
        max_ind: int
            Maximum index
        pmid_len: int
            Length of pmid which should not be exceeded by indices

    Returns:
        pedge_ind: list
            List of edge pressure indices, possibly with new minimum and maximum indices
    """

    if max_ind > pmid_len:
        # don't overstep array bounds for full array
        max_ind = max_ind - 1
    if min(pedge_ind) != 0:
        pedge_ind = np.append(min(pedge_ind) - 1, pedge_ind)
    if max(pedge_ind) != max_ind:
        pedge_ind = np.append(pedge_ind, max(pedge_ind) + 1)
    return pedge_ind


def get_ind_of_pres(dataset, pres):
    """
    Get index of pressure level that contains the requested pressure value.

    Args:
        dataset: xarray Dataset
            GEOS-Chem dataset
        pres: int or float
            Desired pressure value

    Returns:
        index: int
            Index of level in dataset that corresponds to requested pressure

    """
    pedge, pmid, _ = get_vert_grid(dataset)
    converted_dataset = convert_lev_to_pres(dataset, pmid, pedge)
    return np.argmin(np.abs(converted_dataset['lev'] - pres).values)


def convert_lev_to_pres(dataset, pmid, pedge, lev_type='pmid'):
    """
    Convert lev dimension to pressure in a GEOS-Chem dataset

    Args:
        dataset: xarray Dataset
            GEOS-Chem dataset
        pmid: np.array
            Midpoint pressure values
        pedge: np.array
            Edge pressure values
        lev_type (optional): str
            Denote whether lev is 'pedge' or 'pmid' if grid is not 72/73 or 47/48 levels
            Default value: 'pmid'

    Returns:
        dataset: xarray Dataset
            Input dataset with "lev" dimension values replaced with pressure values
    """

    if dataset.sizes["lev"] in (72, 47):
        dataset["lev"] = pmid
    elif dataset.sizes["lev"] in (73, 48):
        dataset["lev"] = pedge
    elif lev_type == 'pmid':
        print('Warning: Assuming levels correspond with midpoint pressures')
        dataset["lev"] = pmid
    else:
        dataset["lev"] = pedge
    dataset["lev"].attrs["unit"] = "hPa"
    dataset["lev"].attrs["long_name"] = "level pressure"
    return dataset


class vert_grid:
    def __init__(self, AP=None, BP=None, p_sfc=1013.25):
        if (len(AP) != len(BP)) or (AP is None):
            # Throw error?
            print('Inconsistent vertical grid specification')
        self.AP = np.array(AP)
        self.BP = np.array(BP)
        self.p_sfc = p_sfc

    def p_edge(self):
        # Calculate pressure edges using eta coordinate
        return self.AP + self.BP * self.p_sfc

    def p_mid(self):
        p_edge = self.p_edge()
        return (p_edge[1:] + p_edge[:-1]) / 2.0


# Standard vertical grids
_GEOS_72L_AP = np.array([0.000000e+00,
                         4.804826e-02,
                         6.593752e+00,
                         1.313480e+01,
                         1.961311e+01,
                         2.609201e+01,
                         3.257081e+01,
                         3.898201e+01,
                         4.533901e+01,
                         5.169611e+01,
                         5.805321e+01,
                         6.436264e+01,
                         7.062198e+01,
                         7.883422e+01,
                         8.909992e+01,
                         9.936521e+01,
                         1.091817e+02,
                         1.189586e+02,
                         1.286959e+02,
                         1.429100e+02,
                         1.562600e+02,
                         1.696090e+02,
                         1.816190e+02,
                         1.930970e+02,
                         2.032590e+02,
                         2.121500e+02,
                         2.187760e+02,
                         2.238980e+02,
                         2.243630e+02,
                         2.168650e+02,
                         2.011920e+02,
                         1.769300e+02,
                         1.503930e+02,
                         1.278370e+02,
                         1.086630e+02,
                         9.236572e+01,
                         7.851231e+01,
                         6.660341e+01,
                         5.638791e+01,
                         4.764391e+01,
                         4.017541e+01,
                         3.381001e+01,
                         2.836781e+01,
                         2.373041e+01,
                         1.979160e+01,
                         1.645710e+01,
                         1.364340e+01,
                         1.127690e+01,
                         9.292942e+00,
                         7.619842e+00,
                         6.216801e+00,
                         5.046801e+00,
                         4.076571e+00,
                         3.276431e+00,
                         2.620211e+00,
                         2.084970e+00,
                         1.650790e+00,
                         1.300510e+00,
                         1.019440e+00,
                         7.951341e-01,
                         6.167791e-01,
                         4.758061e-01,
                         3.650411e-01,
                         2.785261e-01,
                         2.113490e-01,
                         1.594950e-01,
                         1.197030e-01,
                         8.934502e-02,
                         6.600001e-02,
                         4.758501e-02,
                         3.270000e-02,
                         2.000000e-02,
                         1.000000e-02])

_GEOS_72L_BP = np.array([1.000000e+00,
                         9.849520e-01,
                         9.634060e-01,
                         9.418650e-01,
                         9.203870e-01,
                         8.989080e-01,
                         8.774290e-01,
                         8.560180e-01,
                         8.346609e-01,
                         8.133039e-01,
                         7.919469e-01,
                         7.706375e-01,
                         7.493782e-01,
                         7.211660e-01,
                         6.858999e-01,
                         6.506349e-01,
                         6.158184e-01,
                         5.810415e-01,
                         5.463042e-01,
                         4.945902e-01,
                         4.437402e-01,
                         3.928911e-01,
                         3.433811e-01,
                         2.944031e-01,
                         2.467411e-01,
                         2.003501e-01,
                         1.562241e-01,
                         1.136021e-01,
                         6.372006e-02,
                         2.801004e-02,
                         6.960025e-03,
                         8.175413e-09,
                         0.000000e+00,
                         0.000000e+00,
                         0.000000e+00,
                         0.000000e+00,
                         0.000000e+00,
                         0.000000e+00,
                         0.000000e+00,
                         0.000000e+00,
                         0.000000e+00,
                         0.000000e+00,
                         0.000000e+00,
                         0.000000e+00,
                         0.000000e+00,
                         0.000000e+00,
                         0.000000e+00,
                         0.000000e+00,
                         0.000000e+00,
                         0.000000e+00,
                         0.000000e+00,
                         0.000000e+00,
                         0.000000e+00,
                         0.000000e+00,
                         0.000000e+00,
                         0.000000e+00,
                         0.000000e+00,
                         0.000000e+00,
                         0.000000e+00,
                         0.000000e+00,
                         0.000000e+00,
                         0.000000e+00,
                         0.000000e+00,
                         0.000000e+00,
                         0.000000e+00,
                         0.000000e+00,
                         0.000000e+00,
                         0.000000e+00,
                         0.000000e+00,
                         0.000000e+00,
                         0.000000e+00,
                         0.000000e+00,
                         0.000000e+00])

GEOS_72L_grid = vert_grid(_GEOS_72L_AP, _GEOS_72L_BP)

# Reduced grid
_GEOS_47L_AP = np.zeros(48)
_GEOS_47L_BP = np.zeros(48)

# Fill in the values for the surface
_GEOS_47L_AP[0] = _GEOS_72L_AP[0]
_GEOS_47L_BP[0] = _GEOS_72L_BP[0]

# Build the GEOS 72-layer to 47-layer mapping matrix at the same time
_xmat_i = np.zeros((72))
_xmat_j = np.zeros((72))
_xmat_s = np.zeros((72))

# Index here is the 1-indexed layer number
for _i_lev in range(1, 37):
    # Map from 1-indexing to 0-indexing
    _x_lev = _i_lev - 1
    # Sparse matrix for regridding
    # Below layer 37, it's 1:1
    _xct = _x_lev
    _xmat_i[_xct] = _x_lev
    _xmat_j[_xct] = _x_lev
    _xmat_s[_xct] = 1.0
    # Copy over the pressure edge for the top of the grid cell
    _GEOS_47L_AP[_i_lev] = _GEOS_72L_AP[_i_lev]
    _GEOS_47L_BP[_i_lev] = _GEOS_72L_BP[_i_lev]

# Now deal with the lumped layers
_skip_size_vec = [2, 4]
_number_lumped = [4, 7]

# Initialize
_i_lev = 36
_i_lev_72 = 36
for _lump_seg in range(2):
    _skip_size = _skip_size_vec[_lump_seg]
    # 1-indexed starting point in the 47-layer grid
    _first_lev_47 = _i_lev + 1
    _first_lev_72 = _i_lev_72 + 1

    # Loop over the coarse vertical levels (47-layer grid)
    for _i_lev_offset in range(_number_lumped[_lump_seg]):
        # i_lev is the index for the current level on the 47-level grid
        _i_lev = _first_lev_47 + _i_lev_offset
        # Map from 1-indexing to 0-indexing
        _x_lev = _i_lev - 1

        # Get the 1-indexed location of the last layer in the 72-layer grid
        # which is below the start of the current lumping region
        _i_lev_72_base = _first_lev_72 + (_i_lev_offset * _skip_size) - 1

        # Get the 1-indexed location of the uppermost level in the 72-layer
        # grid which is within the target layer on the 47-layer grid
        _i_lev_72 = _i_lev_72_base + _skip_size

        # Do the pressure edges first
        # These are the 0-indexed locations of the upper edge for the
        # target layers in 47- and 72-layer grids
        _GEOS_47L_AP[_i_lev] = _GEOS_72L_AP[_i_lev_72]
        _GEOS_47L_BP[_i_lev] = _GEOS_72L_BP[_i_lev_72]

        # Get the total pressure delta across the layer on the lumped grid
        # We are within the fixed pressure levels so don't need to account
        # for variations in surface pressure
        _dp_total = _GEOS_47L_AP[_i_lev - 1] - _GEOS_47L_AP[_i_lev]

        # Now figure out the mapping
        for _i_lev_offset_72 in range(_skip_size):
            # Source layer in the 72 layer grid (0-indexed)
            _x_lev_72 = _i_lev_72_base + _i_lev_offset_72
            _xct = _x_lev_72
            _xmat_i[_xct] = _x_lev_72
            # Target in the 47 layer grid
            _xmat_j[_xct] = _x_lev

            # Proportion of 72-layer grid cell, by pressure, within expanded
            # layer
            _xmat_s[_xct] = (_GEOS_72L_AP[_x_lev_72] -
                             _GEOS_72L_AP[_x_lev_72 + 1]) / _dp_total
    _start_pt = _i_lev

# Do last entry separately (no layer to go with it)
_xmat_72to47 = scipy.sparse.coo_matrix(
    (_xmat_s, (_xmat_i, _xmat_j)), shape=(72, 47))

GEOS_47L_grid = vert_grid(_GEOS_47L_AP, _GEOS_47L_BP)

# CAM 26-layer grid
_CAM_26L_AP = np.flip(np.array([219.4067, 489.5209, 988.2418, 1805.201,
                                2983.724, 4462.334, 6160.587, 7851.243,
                                7731.271, 7590.131, 7424.086, 7228.744,
                                6998.933, 6728.574, 6410.509, 6036.322,
                                5596.111, 5078.225, 4468.96, 3752.191,
                                2908.949, 2084.739, 1334.443, 708.499,
                                252.136, 0., 0.]), axis=0) * 0.01
_CAM_26L_BP = np.flip(np.array([0., 0., 0., 0.,
                                0., 0., 0., 0.,
                                0.01505309, 0.03276228, 0.05359622, 0.07810627,
                                0.1069411, 0.14086370, 0.180772, 0.227722,
                                0.2829562, 0.3479364, 0.4243822, 0.5143168,
                                0.6201202, 0.7235355, 0.8176768, 0.8962153,
                                0.9534761, 0.9851122, 1.]), axis=0)

CAM_26L_grid = vert_grid(_CAM_26L_AP, _CAM_26L_BP)


def make_grid_LL(llres, in_extent=[-180, 180, -90, 90], out_extent=[]):
    """
    Creates a lat/lon grid description.

    Args:
        llres: str
            lat/lon resolution in 'latxlon' format (e.g. '4x5')

    Keyword Args (optional):
        in_extent: list[float, float, float, float]
            Describes minimum and maximum latitude and longitude of initial grid
            in the format [minlon, maxlon, minlat, maxlat]
            Default value: [-180, 180, -90, 90]
        out_extent: list[float, float, float, float]
            Describes minimum and maximum latitude and longitude of target grid
            in the format [minlon, maxlon, minlat, maxlat]. Needed when intending
            to use grid to trim extent of input data
            Default value: [] (assumes value of in_extent)

    Returns:
        llgrid: dict
            dict grid description of format {'lat'   : lat midpoints,
                                             'lon'   : lon midpoints,
                                             'lat_b' : lat edges,
                                             'lon_b' : lon edges}
    """

    # get initial bounds of grid
    [minlon, maxlon, minlat, maxlat] = in_extent
    [dlat, dlon] = list(map(float, llres.split('x')))
    lon_b = np.linspace(minlon - dlon / 2, maxlon - dlon /
                        2, int((maxlon - minlon) / dlon) + 1)
    lat_b = np.linspace(minlat - dlat / 2, maxlat + dlat / 2,
                        int((maxlat - minlat) / dlat) + 2)
    if minlat <= -90:
        lat_b = lat_b.clip(-90, None)
    if maxlat >= 90:
        lat_b = lat_b.clip(None, 90)
    lat = (lat_b[1:] + lat_b[:-1]) / 2
    lon = (lon_b[1:] + lon_b[:-1]) / 2

    # trim grid bounds when your desired extent is not the same as your
    # initial grid extent
    if out_extent == []:
        out_extent = in_extent
    if out_extent != in_extent:
        [minlon, maxlon, minlat, maxlat] = out_extent
        minlon_ind = np.nonzero(lon >= minlon)
        maxlon_ind = np.nonzero(lon <= maxlon)
        lon_inds = np.intersect1d(minlon_ind, maxlon_ind)
        lon = lon[lon_inds]
        # make sure to get edges of grid correctly
        lon_inds = np.append(lon_inds, np.max(lon_inds) + 1)
        lon_b = lon_b[lon_inds]

        minlat_ind = np.nonzero(lat >= minlat)
        maxlat_ind = np.nonzero(lat <= maxlat)
        lat_inds = np.intersect1d(minlat_ind, maxlat_ind)
        lat = lat[lat_inds]
        # make sure to get edges of grid correctly
        lat_inds = np.append(lat_inds, np.max(lat_inds) + 1)
        lat_b = lat_b[lat_inds]

    llgrid = {'lat': lat,
              'lon': lon,
              'lat_b': lat_b,
              'lon_b': lon_b}
    return llgrid


def make_grid_CS(csres):
    """
    Creates a cubed-sphere grid description.

    Args:
        csres: int
            cubed-sphere resolution of target grid

    Returns:
        [csgrid, csgrid_list]: list[dict, list[dict]]
            csgrid is a dict of format {'lat'   : lat midpoints,
                                        'lon'   : lon midpoints,
                                        'lat_b' : lat edges,
                                        'lon_b' : lon edges}
            where each value has an extra face dimension of length 6.
            csgrid_list is a list of dicts separated by face index
    """

    csgrid = csgrid_GMAO(csres)
    csgrid_list = [None] * 6
    for i in range(6):
        csgrid_list[i] = {'lat': csgrid['lat'][i],
                          'lon': csgrid['lon'][i],
                          'lat_b': csgrid['lat_b'][i],
                          'lon_b': csgrid['lon_b'][i]}

    return [csgrid, csgrid_list]


def make_grid_SG(csres, stretch_factor, target_lon, target_lat):
    """
    Creates a stretched-grid grid description.

    Args:
        csres: int
            cubed-sphere resolution of target grid
        stretch_factor: float
            stretch factor of target grid
        target_lon: float
            target stretching longitude of target grid
        target_lon: float
            target stretching latitude of target grid

    Returns:
        [csgrid, csgrid_list]: list[dict, list[dict]]
            csgrid is a dict of format {'lat'   : lat midpoints,
                                        'lon'   : lon midpoints,
                                        'lat_b' : lat edges,
                                        'lon_b' : lon edges}
            where each value has an extra face dimension of length 6.
            csgrid_list is a list of dicts separated by face index
    """

    csgrid = csgrid_GMAO(csres, offset=0)
    csgrid_list = [None] * 6
    for i in range(6):
        lat = csgrid['lat'][i].flatten()
        lon = csgrid['lon'][i].flatten()
        lon, lat = scs_transform(
            lon, lat, stretch_factor, target_lon, target_lat)
        lat = lat.reshape((csres, csres))
        lon = lon.reshape((csres, csres))
        lat_b = csgrid['lat_b'][i].flatten()
        lon_b = csgrid['lon_b'][i].flatten()
        lon_b, lat_b = scs_transform(
            lon_b, lat_b, stretch_factor, target_lon, target_lat)
        lat_b = lat_b.reshape((csres + 1, csres + 1))
        lon_b = lon_b.reshape((csres + 1, csres + 1))
        csgrid_list[i] = {'lat': lat,
                          'lon': lon,
                          'lat_b': lat_b,
                          'lon_b': lon_b}
    for i in range(6):
        csgrid['lat'][i] = csgrid_list[i]['lat']
        csgrid['lon'][i] = csgrid_list[i]['lon']
        csgrid['lat_b'][i] = csgrid_list[i]['lat_b']
        csgrid['lon_b'][i] = csgrid_list[i]['lon_b']
    return [csgrid, csgrid_list]


def calc_rectilinear_lon_edge(lon_stride, center_at_180):
    """ Compute longitude edge vector for a rectilinear grid.
    Parameters
    ----------
    lon_stride: float
        Stride length in degrees. For example, for a standard GEOS-Chem Classic
        4x5 grid, lon_stride would be 5.
    center_at_180: bool
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

    n_lon = np.round(360.0 / lon_stride)
    lon_edge = np.linspace(-180.0, 180.0, num=n_lon + 1)
    if center_at_180:
        lon_edge = lon_edge - (lon_stride / 2.0)

    lon_edge[lon_edge < -180.0] = lon_edge[lon_edge < -180] + 360.0
    lon_edge[lon_edge > 180.0] = lon_edge[lon_edge > 180.0] - 360.0

    return lon_edge


def calc_rectilinear_lat_edge(lat_stride, half_polar_grid):
    """ Compute latitude edge vector for a rectilinear grid.
    Parameters
    ----------
    lat_stride: float
        Stride length in degrees. For example, for a standard GEOS-Chem Classic
        4x5 grid, lat_stride would be 4.
    half_polar_grid: bool
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
        start_pt = 90.0 + (lat_stride / 2.0)
    else:
        start_pt = 90.0

    lat_edge = np.linspace(-1.0 * start_pt, start_pt,
                           num=1 + np.round(2.0 * start_pt / lat_stride))

    # Force back onto +/- 90
    lat_edge[lat_edge > 90.0] = 90.0
    lat_edge[lat_edge < -90.0] = -90.0

    return lat_edge


def calc_rectilinear_grid_area(lon_edge, lat_edge):
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

    # Convert from km to m
    _radius_earth_m = R_EARTH_m

    lon_edge = asarray(lon_edge, dtype=float)
    lat_edge = asarray(lat_edge, dtype=float)

    n_lon = (lon_edge.size) - 1
    n_lat = (lat_edge.size) - 1

    grid_area = np.zeros((n_lat, n_lon))

    sfc_area_const = 2.0 * np.pi * _radius_earth_m * _radius_earth_m

    # Longitudes loop, so need to be careful
    lon_delta = calc_delta_lon(lon_edge)

    # Convert into weights relative to the total circle
    lon_delta = lon_delta / 360.0

    # Precalculate this
    sin_lat_edge = np.sin(np.deg2rad(lat_edge))

    for i_lat in range(0, n_lat):
        sin_diff = sin_lat_edge[i_lat + 1] - sin_lat_edge[i_lat]
        grid_area[i_lat, :] = sin_diff * sfc_area_const * lon_delta

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
    for i_lon in range(0, n_lon):
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
    res: cubed-sphere Resolution
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


_INV_SQRT_3 = 1.0 / np.sqrt(3.0)
_ASIN_INV_SQRT_3 = np.arcsin(_INV_SQRT_3)


class CSGrid(object):
    """Generator for cubed-sphere grid geometries.
    CSGrid computes the latitutde and longitudes of cell centers and edges
    on a cubed-sphere grid, providing a way to retrieve these geometries
    on-the-fly if your model output data does not include them.
    Attributes
    ----------
    {lon,lat}_center: np.ndarray
        lat/lon coordinates for each cell center along the cubed-sphere mesh
    {lon,lat}_edge: np.ndarray
        lat/lon coordinates for the midpoint of the edges separating each
        element on the cubed-sphere mesh.
    xyz_{center,edge}: np.ndarray
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
        c: int
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
        offset: float (optional)
            Degrees to offset the first faces' edge in the latitudinal
            direction. If not passed, then the western edge of the first face
            will align with the prime meridian.
       This function was originally written by Jiawei Zhuange and included
       in package cubedsphere: https://github.com/JiaweiZhuang/cubedsphere
        """
        self.c = c
        self.delta_y = 2. * _ASIN_INV_SQRT_3 / c
        self.nx = self.ny = c + 1
        self.offset = offset

        self._initialize()

    def _initialize(self):

        c = self.c
        nx, ny = self.nx, self.ny

        lambda_rad = np.zeros((nx, ny))
        lambda_rad[0, :] = 3. * np.pi / 4.  # West edge
        lambda_rad[-1, :] = 5. * np.pi / 4.  # East edge

        theta_rad = np.zeros((nx, ny))
        theta_rad[0, :] = -_ASIN_INV_SQRT_3 + \
            (self.delta_y * np.arange(c + 1))  # West edge
        theta_rad[-1, :] = theta_rad[0, :]  # East edge

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

            xyzDot = np.sum(xyzCross * xyzRef)
            xyzImg = xyzRef - (2. * xyzDot * xyzCross)

            xsImg, ysImg, zsImg = xyzImg
            lonImg, latImg = cartesian_to_latlon(xsImg, ysImg, zsImg)

            lambda_rad[i, 0] = lonImg
            lambda_rad[i, -1] = lonImg
            theta_rad[i, 0] = latImg
            theta_rad[i, -1] = -latImg

        pp = np.zeros([3, c + 1, c + 1])

        # Set the four corners
        # print("CORNERS")
        for i, j in product([0, -1], [0, -1]):
            # print(i, j)
            pp[:, i, j] = latlon_to_cartesian(
                lambda_rad[i, j], theta_rad[i, j])

        # Map the edges on the sphere back to the cube. 
        #Note that all intersections are at x = -rsq3
        # print("EDGES")
        for ij in range(1, c + 1):
            # print(ij)
            pp[:, 0, ij] = latlon_to_cartesian(
                lambda_rad[0, ij], theta_rad[0, ij])
            pp[1, 0, ij] = -pp[1, 0, ij] * _INV_SQRT_3 / pp[0, 0, ij]
            pp[2, 0, ij] = -pp[2, 0, ij] * _INV_SQRT_3 / pp[0, 0, ij]

            pp[:, ij, 0] = latlon_to_cartesian(
                lambda_rad[ij, 0], theta_rad[ij, 0])
            pp[1, ij, 0] = -pp[1, ij, 0] * _INV_SQRT_3 / pp[0, ij, 0]
            pp[2, ij, 0] = -pp[2, ij, 0] * _INV_SQRT_3 / pp[0, ij, 0]

        # # Map interiors
        pp[0, :, :] = -_INV_SQRT_3
        # print("INTERIOR")
        for i in range(1, c + 1):
            for j in range(1, c + 1):
                # Copy y-z face of the cube along j=1
                pp[1, i, j] = pp[1, i, 0]
                # Copy along i=1
                pp[2, i, j] = pp[2, 0, j]

        _pp = pp.copy()
        llr, ttr = vec_cartesian_to_latlon(_pp[0], _pp[1], _pp[2])

        lambda_rad, theta_rad = llr.copy(), ttr.copy()

        # Make grid symmetrical to i = im/2 + 1
        for j in range(1, c + 1):
            for i in range(1, c + 1):
                # print("({}, {}) -> ({}, {})".format(i, 0, i, j))
                lambda_rad[i, j] = lambda_rad[i, 0]

        for j in range(c + 1):
            for i in range(c // 2):
                isymm = c - i
                # print(isymm)
                avgPt = 0.5 * (lambda_rad[i, j] - lambda_rad[isymm, j])
                # print(lambda_rad[i, j], lambda_rad[isymm, j], avgPt)
                lambda_rad[i, j] = avgPt + np.pi
                lambda_rad[isymm, j] = np.pi - avgPt

                avgPt = 0.5 * (theta_rad[i, j] + theta_rad[isymm, j])
                theta_rad[i, j] = avgPt
                theta_rad[isymm, j] = avgPt

        # Make grid symmetrical to j = im/2 + 1
        for j in range(c // 2):
            jsymm = c - j
            for i in range(1, c + 1):
                avgPt = 0.5 * (lambda_rad[i, j] + lambda_rad[i, jsymm])
                lambda_rad[i, j] = avgPt
                lambda_rad[i, jsymm] = avgPt

                avgPt = 0.5 * (theta_rad[i, j] - theta_rad[i, jsymm])
                theta_rad[i, j] = avgPt
                theta_rad[i, jsymm] = -avgPt

        # Final correction
        lambda_rad -= np.pi

        llr, ttr = lambda_rad.copy(), theta_rad.copy()

        #######################################################################
        # MIRROR GRIDS
        #######################################################################

        new_xgrid = np.zeros((c + 1, c + 1, 6))
        new_ygrid = np.zeros((c + 1, c + 1, 6))

        xgrid = llr.copy()
        ygrid = ttr.copy()

        new_xgrid[..., 0] = xgrid.copy()
        new_ygrid[..., 0] = ygrid.copy()

        # radius = 6370.0e3
        radius = 1.

        for face in range(1, 6):
            for j in range(c + 1):
                for i in range(c + 1):
                    x = xgrid[i, j]
                    y = ygrid[i, j]
                    z = radius

                    if face == 1:
                        # Rotate about z only
                        new_xyz = rotate_sphere_3D(x, y, z, -np.pi / 2., 'z')

                    elif face == 2:
                        # Rotate about z, then x
                        temp_xyz = rotate_sphere_3D(x, y, z, -np.pi / 2., 'z')
                        x, y, z = temp_xyz[:]
                        new_xyz = rotate_sphere_3D(x, y, z, np.pi / 2., 'x')

                    elif face == 3:
                        temp_xyz = rotate_sphere_3D(x, y, z, np.pi, 'z')
                        x, y, z = temp_xyz[:]
                        new_xyz = rotate_sphere_3D(x, y, z, np.pi / 2., 'x')

                        if ((c % 2) != 0) and (j == c // 2 - 1):
                            print(i, j, face)
                            new_xyz = (np.pi, *new_xyz)

                    elif face == 4:
                        temp_xyz = rotate_sphere_3D(x, y, z, np.pi / 2., 'z')
                        x, y, z = temp_xyz[:]
                        new_xyz = rotate_sphere_3D(x, y, z, np.pi / 2., 'y')

                    elif face == 5:
                        temp_xyz = rotate_sphere_3D(x, y, z, np.pi / 2., 'y')
                        x, y, z = temp_xyz[:]
                        new_xyz = rotate_sphere_3D(x, y, z, 0., 'z')

                    # print((x, y, z), "\n", new_xyz, "\n" + "--"*40)

                    new_x, new_y, _ = new_xyz
                    new_xgrid[i, j, face] = new_x
                    new_ygrid[i, j, face] = new_y

        lon_edge, lat_edge = new_xgrid.copy(), new_ygrid.copy()

        #######################################################################
        # CLEANUP GRID
        #######################################################################

        for i, j, f in product(range(c + 1), range(c + 1), range(6)):
            new_lon = lon_edge[i, j, f]
            if new_lon < 0:
                new_lon += 2 * np.pi
            if np.abs(new_lon) < 1e-10:
                new_lon = 0.
            lon_edge[i, j, f] = new_lon

            if np.abs(lat_edge[i, j, f]) < 1e-10:
                lat_edge[i, j, f] = 0.

        lon_edge_deg = np.rad2deg(lon_edge)
        lat_edge_deg = np.rad2deg(lat_edge)

        #######################################################################
        # COMPUTE CELL CENTROIDS
        #######################################################################

        lon_ctr = np.zeros((c, c, 6))
        lat_ctr = np.zeros((c, c, 6))
        xyz_ctr = np.zeros((3, c, c, 6))
        xyz_edge = np.zeros((3, c + 1, c + 1, 6))

        for f in range(6):
            for i in range(c):
                last_x = (i == (c - 1))
                for j in range(c):
                    last_y = (j == (c - 1))

                    # Get the four corners
                    lat_corner = [
                        lat_edge[i, j, f],
                        lat_edge[i + 1, j, f],
                        lat_edge[i + 1, j + 1, f],
                        lat_edge[i, j + 1, f]]
                    lon_corner = [
                        lon_edge[i, j, f],
                        lon_edge[i + 1, j, f],
                        lon_edge[i + 1, j + 1, f],
                        lon_edge[i, j + 1, f]]

                    # Convert from lat-lon back to cartesian
                    xyz_corner = np.asarray(
                        vec_latlon_to_cartesian(
                            lon_corner, lat_corner))

                    # Store the edge information
                    xyz_edge[:, i, j, f] = xyz_corner[:, 0]
                    if last_x:
                        xyz_edge[:, i + 1, j, f] = xyz_corner[:, 1]
                    if last_x or last_y:
                        xyz_edge[:, i + 1, j + 1, f] = xyz_corner[:, 2]
                    if last_y:
                        xyz_edge[:, i, j + 1, f] = xyz_corner[:, 3]

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
        # CACHE
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
    vector_length = np.sqrt(np.sum(xyz * xyz, axis=0))
    xyz /= vector_length
    x, y, z = xyz

    if (np.abs(x) + np.abs(y)) < 1e-20:
        lon = 0.
    else:
        lon = np.arctan2(y, x)
    if lon < 0.:
        lon += 2 * np.pi

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
        y_new = cos_ang * y + sin_ang * z
        z_new = -sin_ang * y + cos_ang * z
    elif rot_axis == 'y':
        x_new = cos_ang * x - sin_ang * z
        y_new = y
        z_new = sin_ang * x + cos_ang * z
    elif rot_axis == 'z':
        x_new = cos_ang * x + sin_ang * y
        y_new = -sin_ang * x + cos_ang * y
        z_new = z

    theta_new, phi_new, r_new = cartesian_to_spherical(x_new, y_new, z_new)

    return theta_new, phi_new, r_new
