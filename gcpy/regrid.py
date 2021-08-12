''' Functions for creating xesmf regridder objects '''

import os
import xesmf as xe
from .grid import make_grid_LL, make_grid_CS, make_grid_SG, get_input_res, call_make_grid, \
    get_grid_extents, get_vert_grid
import hashlib
import numpy as np
import xarray as xr
import pandas as pd
import scipy.sparse
import warnings

def make_regridder_L2L(
        llres_in, llres_out, weightsdir='.', reuse_weights=False,
        in_extent=[-180, 180, -90, 90],
        out_extent=[-180, 180, -90, 90]):
    """
    Create an xESMF regridder between two lat/lon grids

    Args:
        llres_in: str
            Resolution of input grid in format 'latxlon', e.g. '4x5'
        llres_out: str
            Resolution of output grid in format 'latxlon', e.g. '4x5'

    Keyword Args (optional):
        weightsdir: str
            Directory in which to create xESMF regridder NetCDF files
            Default value: '.'
        reuse_weights: bool
            Set this flag to True to reuse existing xESMF regridder NetCDF files
            Default value: False
        in_extent: list[float, float, float, float]
            Describes minimum and maximum latitude and longitude of input grid
            in the format [minlon, maxlon, minlat, maxlat]
            Default value: [-180, 180, -90, 90]
        out_extent: list[float, float, float, float]
            Desired minimum and maximum latitude and longitude of output grid
            in the format [minlon, maxlon, minlat, maxlat]
            Default value: [-180, 180, -90, 90]

    Returns:
        regridder: xESMF regridder
            regridder object between the two specified grids
    """

    llgrid_in = make_grid_LL(llres_in, in_extent, out_extent)
    llgrid_out = make_grid_LL(llres_out, out_extent)
    if in_extent == [-180, 180, -90,
                     90] and out_extent == [-180, 180, -90, 90]:
        weightsfile = os.path.join(
            weightsdir, 'conservative_{}_{}.nc'.format(
                llres_in, llres_out))
    else:
        in_extent_str = str(in_extent).replace(
            '[', '').replace(
            ']', '').replace(
            ', ', 'x')
        out_extent_str = str(out_extent).replace(
            '[', '').replace(
            ']', '').replace(
            ', ', 'x')
        weightsfile = os.path.join(
            weightsdir, 'conservative_{}_{}_{}_{}.nc'.format(
                llres_in, llres_out, in_extent_str, out_extent_str))
        
    if not os.path.isfile(weightsfile) and reuse_weights:
        #prevent error with more recent versions of xesmf
        reuse_weights=False
        
    try:
        regridder = xe.Regridder(
            llgrid_in,
            llgrid_out,
            method='conservative',
            filename=weightsfile,
            reuse_weights=reuse_weights)
    except BaseException:
        regridder = xe.Regridder(
            llgrid_in,
            llgrid_out,
            method='conservative',
            filename=weightsfile,
            reuse_weights=reuse_weights)
    return regridder


def make_regridder_C2L(csres_in, llres_out, weightsdir='.',
                       reuse_weights=True, sg_params=[1, 170, -90]):
    """
    Create an xESMF regridder from a cubed-sphere to lat/lon grid

    Args:
        csres_in: int
            Cubed-sphere resolution of input grid
        llres_out: str
            Resolution of output grid in format 'latxlon', e.g. '4x5'

    Keyword Args (optional):
        weightsdir: str
            Directory in which to create xESMF regridder NetCDF files
            Default value: '.'
        reuse_weights: bool
            Set this flag to True to reuse existing xESMF regridder NetCDF files
            Default value: False
        sg_params: list[float, float, float] (stretch_factor, target_longitude, target_latitude)
            Input grid stretched-grid parameters in the format
            [stretch_factor, target_longitude, target_latitude].
            Will trigger stretched-grid creation if not default values.
            Default value: [1, 170, -90] (no stretching)

    Returns:
        regridder_list: list[6 xESMF regridders]
            list of regridder objects (one per cubed-sphere face) between the two specified grids
    """
    [sf_in, tlon_in, tlat_in] = sg_params
    if sg_params == [1, 170, -90]:
        _, csgrid_list = make_grid_CS(csres_in)
    else:
        _, csgrid_list = make_grid_SG(
            csres_in, stretch_factor=sg_params[0],
            target_lon=sg_params[1],
            target_lat=sg_params[2])
    llgrid = make_grid_LL(llres_out)
    regridder_list = []
    for i in range(6):
        if sg_params == [1, 170, -90]:
            weightsfile = os.path.join(
                weightsdir, 'conservative_c{}_{}_{}.nc'.format(
                    str(csres_in), llres_out, str(i)))
        else:
            weights_fname = f'conservative_sg{sg_hash(csres_in, sf_in, tlat_in, tlon_in)}_ll{llres_out}_F{i}.nc'
            weightsfile = os.path.join(weightsdir, weights_fname)
            
        if not os.path.isfile(weightsfile) and reuse_weights:
            #prevent error with more recent versions of xesmf
            reuse_weights=False

        try:
            regridder = xe.Regridder(
                csgrid_list[i],
                llgrid,
                method='conservative',
                filename=weightsfile,
                reuse_weights=reuse_weights)
        except BaseException:
            regridder = xe.Regridder(
                csgrid_list[i],
                llgrid,
                method='conservative',
                filename=weightsfile,
                reuse_weights=reuse_weights)
        regridder_list.append(regridder)
    return regridder_list


def make_regridder_S2S(
        csres_in,
        csres_out,
        sf_in=1,
        tlon_in=170,
        tlat_in=-90,
        sf_out=1,
        tlon_out=170,
        tlat_out=-90,
        weightsdir='.',
        verbose=True):
    """
    Create an xESMF regridder from a cubed-sphere / stretched-grid grid
    to another cubed-sphere / stretched-grid grid.
    Stretched-grid params of 1, 170, -90 indicate no stretching.

    Args:
        csres_in: int
            Cubed-sphere resolution of input grid
        csres_out: int
            Cubed-sphere resolution of output grid

    Keyword Args (optional):
        sf_in: float
            Stretched-grid factor of input grid
            Default value: 1
        tlon_in: float
            Target longitude for stretching in input grid
            Default value: 170
        tlat_in: float
            Target longitude for stretching in input grid
            Default value: -90
        sf_out: float
            Stretched-grid factor of output grid
            Default value: 1
        tlon_out: float
            Target longitude for stretchingg in output grid
            Default value: 170
        tlat_out: float
            Target longitude for stretching in output grid
            Default value: -90
        weightsdir: str
            Directory in which to create xESMF regridder NetCDF files
            Default value: '.'
        verbose: bool
            Set this flag to True to enable printing when output faces do not
            intersect input faces when regridding
            Default value: True

    Returns:
        regridder_list: list[6 xESMF regridders]
            list of regridder objects (one per cubed-sphere face) between the two specified grids
    """

    _, igrid_list = make_grid_SG(
        csres_in, stretch_factor=sf_in, target_lat=tlat_in, target_lon=tlon_in)
    _, ogrid_list = make_grid_SG(
        csres_out, stretch_factor=sf_out, target_lat=tlat_out,
        target_lon=tlon_out)
    regridder_list = []
    for o_face in range(6):
        regridder_list.append({})
        for i_face in range(6):
            weights_fname = f'conservative_sg{sg_hash(csres_in, sf_in, tlat_in, tlon_in)}_F{i_face}_sg{sg_hash(csres_out, sf_out, tlat_out, tlon_out)}_F{o_face}.nc'
            weightsfile = os.path.join(weightsdir, weights_fname)
            reuse_weights = os.path.exists(weightsfile)
            if not os.path.isfile(weightsfile) and reuse_weights:
                #prevent error with more recent versions of xesmf
                reuse_weights=False

            try:
                regridder = xe.Regridder(igrid_list[i_face],
                                         ogrid_list[o_face],
                                         method='conservative',
                                         filename=weightsfile,
                                         reuse_weights=reuse_weights)
                regridder_list[-1][i_face] = regridder
            except ValueError:
                if verbose:
                    print(f"iface {i_face} doesn't intersect oface {o_face}")

    return regridder_list


def make_regridder_L2S(llres_in, csres_out, weightsdir='.',
                       reuse_weights=True, sg_params=[1, 170, -90]):
    """
    Create an xESMF regridder from a lat/lon to a cubed-sphere grid

    Args:
        llres_in: str
            Resolution of input grid in format 'latxlon', e.g. '4x5'
        csres_out: int
            Cubed-sphere resolution of output grid

    Keyword Args (optional):
        weightsdir: str
            Directory in which to create xESMF regridder NetCDF files
            Default value: '.'
        reuse_weights: bool
            Set this flag to True to reuse existing xESMF regridder NetCDF files
            Default value: False
        sg_params: list[float, float, float] (stretch_factor, target_longitude, target_latitude)
            Output grid stretched-grid parameters in the format
            [stretch_factor, target_longitude, target_latitude].
            Will trigger stretched-grid creation if not default values.
            Default value: [1, 170, -90] (no stretching)

    Returns:
        regridder_list: list[6 xESMF regridders]
            list of regridder objects (one per cubed-sphere face) between the two specified grids
    """

    llgrid = make_grid_LL(llres_in)
    if sg_params == [1, 170, -90]:
        _, csgrid_list = make_grid_CS(csres_out)
    else:
        _, csgrid_list = make_grid_SG(
            csres_out, stretch_factor=sg_params[0],
            target_lon=sg_params[1],
            target_lat=sg_params[2])

    regridder_list = []
    for i in range(6):
        if sg_params == [1, 170, -90]:
            weightsfile = os.path.join(
                weightsdir, 'conservative_{}_c{}_{}.nc'.format(
                    llres_in, str(csres_out), str(i)))
        else:
            weights_fname = f'conservative_ll{llres_in}_sg{sg_hash(csres_out, *sg_params)}_F{i}.nc'
            weightsfile = os.path.join(weightsdir, weights_fname)
        if not os.path.isfile(weightsfile) and reuse_weights:
            #prevent error with more recent versions of xesmf
            reuse_weights=False

        try:
            regridder = xe.Regridder(
                llgrid,
                csgrid_list[i],
                method='conservative',
                filename=weightsfile,
                reuse_weights=reuse_weights)
        except BaseException:
            regridder = xe.Regridder(
                llgrid,
                csgrid_list[i],
                method='conservative',
                filename=weightsfile,
                reuse_weights=reuse_weights)
        regridder_list.append(regridder)
    return regridder_list


def create_regridders(
        refds, devds, weightsdir='.', reuse_weights=True, cmpres=None,
        zm=False, sg_ref_params=[1, 170, -90],
        sg_dev_params=[1, 170, -90]):
    """
    Internal function used for creating regridders between two datasets.
    Follows decision logic needed for plotting functions.
    Originally code from compare_single_level and compare_zonal_mean.

    Args:
        refds: xarray Dataset
            Input dataset
        devds: xarray Dataset
            Output dataset

    Keyword Args (optional):
        weightsdir: str
            Directory in which to create xESMF regridder NetCDF files
            Default value: '.'
        reuse_weights: bool
            Set this flag to True to reuse existing xESMF regridder NetCDF files
            Default value: False
        cmpres: int or str
            Specific target resolution for comparison grid used in difference and ratio plots
            Default value: None (will follow logic chain below)
        zm: bool
            Set this flag to True if regridders will be used in zonal mean plotting
            Default value: False
        sg_ref_params: list[float, float, float] (stretch_factor, target_longitude, target_latitude)
            Ref grid stretched-grid parameters in the format
            [stretch_factor, target_longitude, target_latitude].
            Default value: [1, 170, -90] (no stretching)
        sg_dev_params: list[float, float, float] (stretch_factor, target_longitude, target_latitude)
            Dev grid stretched-grid parameters in the format
            [stretch_factor, target_longitude, target_latitude].
            Default value: [1, 170, -90] (no stretching)

    Returns:
        list of many different quantities needed for regridding in plotting functions
            refres, devres, cmpres: bool
                 Resolution of a dataset grid
            refgridtype, devgridtype, cmpgridtype: str
                 Gridtype of a dataset ('ll' or 'cs')
            regridref, regriddev, regridany: bool
                 Whether to regrid a dataset
            refgrid, devgrid, cmpgrid: dict
                 Grid definition of a dataset
            refregridder, devregridder: xESMF regridder
                 Regridder object between refgrid or devgrid and cmpgrid
                 (will be None if input grid is not lat/lon)
            refregridder_list, devregridder_list: list[6 xESMF regridders]
                 List of regridder objects for each face between refgrid or devgrid and cmpgrid
                 (will be None if input grid is not cubed-sphere)
    """

    # Take two lat/lon or cubed-sphere xarray datasets and regrid them if
    # needed
    refres, refgridtype = get_input_res(refds)
    devres, devgridtype = get_input_res(devds)

    refminlon, refmaxlon, refminlat, refmaxlat = get_grid_extents(refds)
    devminlon, devmaxlon, devminlat, devmaxlat = get_grid_extents(devds)
    # choose smallest overlapping area for comparison
    cmpminlon = max(x for x in [refminlon, devminlon] if x is not None)
    cmpmaxlon = min(x for x in [refmaxlon, devmaxlon] if x is not None)
    cmpminlat = max(x for x in [refminlat, devminlat] if x is not None)
    cmpmaxlat = min(x for x in [refmaxlat, devmaxlat] if x is not None)

    ref_extent = [refminlon, refmaxlon, refminlat, refmaxlat]
    cmp_extent = [cmpminlon, cmpmaxlon, cmpminlat, cmpmaxlat]
    dev_extent = [devminlon, devmaxlon, devminlat, devmaxlat]
    # ==================================================================
    # Determine comparison grid resolution and type (if not passed)
    # ==================================================================

    # If no cmpres is passed then choose highest resolution between ref and dev.
    # If one dataset is lat-lon and the other is cubed-sphere, and no comparison
    # grid resolution is passed, then default to 1x1.25. If both cubed-sphere and
    # plotting zonal mean, over-ride to be 1x1.25 lat-lon with a warning
    sg_cmp_params = [1, 170, -90]
    if cmpres is None:
        if refres == devres and refgridtype == "ll":
            cmpres = refres
            cmpgridtype = refgridtype
        elif refgridtype == "ll" and devgridtype == "ll":
            cmpres = min([refres, devres])
            cmpgridtype = refgridtype
        elif refgridtype == "cs" and devgridtype == "cs":
            if zm:
                print(
                    "Warning: zonal mean comparison must be lat-lon. Defaulting to 1x1.25")
                cmpres = '1x1.25'
                cmpgridtype = "ll"
            elif sg_ref_params != [] or sg_dev_params != []:
                # pick ref grid when a stretched-grid and non-stretched-grid
                # are passed
                cmpres = refres
                cmpgridtype = "cs"
                sg_cmp_params = sg_ref_params
            else:
                # pick higher resolution CS grid out of two standard
                # cubed-sphere grids
                cmpres = max([refres, devres])
                cmpgridtype = "cs"
        elif refgridtype == "ll" and float(refres.split('x')[0]) < 1 and float(refres.split('x')[1]) < 1.25:
            cmpres = refres
            cmpgridtype = "ll"
        elif devgridtype == "ll" and float(devres.split('x')[0]) < 1 and float(devres.split('x')[1]) < 1.25:
            cmpres = devres
            cmpgridtype = "ll"
        else:
            # default to 1x1.25 lat/lon grid for mixed CS and LL grids
            cmpres = "1x1.25"
            cmpgridtype = "ll"
    elif "x" in cmpres:
        cmpgridtype = "ll"
    elif zm:
        print("Warning: zonal mean comparison must be lat-lon. Defaulting to 1x1.25")
        cmpres = '1x1.25'
        cmpgridtype = "ll"
    else:
        # cubed-sphere cmpres
        if isinstance(cmpres, list):
            # stretched-grid resolution
            # first element is cubed-sphere resolution, rest are sg params
            sg_cmp_params = cmpres[1:]
            cmpres = int(cmpres[0])
        else:
            cmpres = int(cmpres)  # must cast to integer for cubed-sphere
        cmpgridtype = "cs"

    # Determine what, if any, need regridding.
    regridref = refres != cmpres or sg_ref_params != sg_cmp_params
    regriddev = devres != cmpres or sg_dev_params != sg_cmp_params
    regridany = regridref or regriddev
    # ==================================================================
    # Make grids (ref, dev, and comparison)
    # ==================================================================
    [refgrid, _] = call_make_grid(
        refres, refgridtype, ref_extent, cmp_extent, sg_ref_params)

    [devgrid, _] = call_make_grid(
        devres, devgridtype, dev_extent, cmp_extent, sg_dev_params)

    [cmpgrid, _] = call_make_grid(
        cmpres, cmpgridtype, cmp_extent, cmp_extent, sg_cmp_params)

    # =================================================================
    # Make regridders, if applicable
    # =================================================================
    refregridder = None
    refregridder_list = None
    devregridder = None
    devregridder_list = None
    if regridref:
        if refgridtype == "ll":
            if cmpgridtype == "cs":
                refregridder_list = make_regridder_L2S(
                    refres,
                    cmpres,
                    weightsdir=weightsdir,
                    reuse_weights=reuse_weights,
                    sg_params=sg_cmp_params)
            else:
                refregridder = make_regridder_L2L(
                    refres,
                    cmpres,
                    weightsdir=weightsdir,
                    reuse_weights=reuse_weights,
                    in_extent=ref_extent)
        else:
            if cmpgridtype == "cs":
                refregridder_list = make_regridder_S2S(
                    refres,
                    cmpres,
                    *sg_ref_params,
                    *sg_cmp_params,
                    weightsdir=weightsdir,
                    verbose=False)
            else:
                refregridder_list = make_regridder_C2L(
                    refres,
                    cmpres,
                    weightsdir=weightsdir,
                    reuse_weights=reuse_weights,
                    sg_params=sg_ref_params)
    if regriddev:
        if devgridtype == "ll":
            if cmpgridtype == "cs":
                devregridder_list = make_regridder_L2S(
                    devres,
                    cmpres,
                    weightsdir=weightsdir,
                    reuse_weights=reuse_weights,
                    sg_params=sg_cmp_params)
            else:
                devregridder = make_regridder_L2L(
                    devres,
                    cmpres,
                    weightsdir=weightsdir,
                    reuse_weights=reuse_weights,
                    in_extent=dev_extent)
        else:
            if cmpgridtype == "cs":
                devregridder_list = make_regridder_S2S(
                    devres,
                    cmpres,
                    *sg_dev_params,
                    *sg_cmp_params,
                    weightsdir=weightsdir,
                    verbose=False)
            else:
                devregridder_list = make_regridder_C2L(
                    devres,
                    cmpres,
                    weightsdir=weightsdir,
                    reuse_weights=reuse_weights,
                    sg_params=sg_dev_params)

    return [
        refres,
        refgridtype,
        devres,
        devgridtype,
        cmpres,
        cmpgridtype,
        regridref,
        regriddev,
        regridany,
        refgrid,
        devgrid,
        cmpgrid,
        refregridder,
        devregridder,
        refregridder_list,
        devregridder_list]


def regrid_comparison_data(
        data,
        res,
        regrid,
        regridder,
        regridder_list,
        global_cmp_grid,
        gridtype,
        cmpgridtype,
        cmpminlat_ind=0,
        cmpmaxlat_ind=-2,
        cmpminlon_ind=0,
        cmpmaxlon_ind=-2,
        nlev=1):
    """
    Regrid comparison datasets to cubed-sphere (including stretched-grid) or lat/lon format.

    Args:
        data: xarray DataArray
            DataArray containing a GEOS-Chem data variable
        res: int
            Cubed-sphere resolution for comparison grid
        regrid: bool
            Set to true to regrid dataset
        regridder: xESMF regridder
            Regridder between the original data grid and the comparison grid
        regridder_list: list(xESMF regridder)
            List of regridders for cubed-sphere data
        global_cmp_grid: xarray DataArray
            Comparison grid
        gridtype: str
            Type of input data grid (either 'll' or 'cs')
        cmpgridtype: str
            Type of input data grid (either 'll' or 'cs')

    Keyword Args (optional):
        cmpminlat_ind: int
            Index of minimum latitude extent for comparison grid
            Default value: 0
        cmpmaxlat_ind: int
            Index (minus 1) of maximum latitude extent for comparison grid
            Default value: -2
        cmpminlon_ind: int
            Index of minimum longitude extent for comparison grid
            Default value: 0
        cmpmaxlon_ind: int
            Index (minus 1) of maximum longitude extent for comparison grid
            Default value: -2
        nlev: int
            Number of levels of input grid and comparison grid
            Default value: 1

    Returns:
        data: xarray DataArray
            Original DataArray regridded to comparison grid (including resolution and extent changes)
    """

    if regrid:
        if gridtype == "ll":
            if cmpgridtype == "ll":
                # regrid ll to ll
                new_data = regridder(data)
            elif cmpgridtype == "cs":
                # ll to CS
                new_data = np.zeros([nlev, 6, res, res]).squeeze()
                for j in range(6):
                    new_data[j, ...] = regridder_list[j](data)
            if nlev == 1:
                # limit to extent of cmpgrid
                new_data=new_data[cmpminlat_ind:cmpmaxlat_ind +
                                  1, cmpminlon_ind:cmpmaxlon_ind + 1].squeeze()
            return new_data
        elif cmpgridtype == "ll":
            # CS to ll
            if nlev == 1:
                new_data = np.zeros([global_cmp_grid['lat'].size,
                                     global_cmp_grid['lon'].size])
                data_reshaped = data.data.reshape(6, res, res)
            else:
                new_data = np.zeros([nlev, global_cmp_grid['lat'].size,
                                     global_cmp_grid['lon'].size])
                data_reshaped = data.data.reshape(
                    nlev, 6, res, res).swapaxes(0, 1)
            for j in range(6):
                regridder = regridder_list[j]
                new_data = new_data + regridder(data_reshaped[j])
            if nlev == 1:
                # limit to extent of cmpgrid
                new_data=new_data[cmpminlat_ind:cmpmaxlat_ind +
                                  1, cmpminlon_ind:cmpmaxlon_ind + 1].squeeze()
            return new_data
        elif cmpgridtype == "cs":
            # CS to CS
            # Reformat dimensions to T, Z, F, Y, X
            if 'Xdim' in data.dims:
                data_format = 'diagnostic'
            else:
                data_format = 'checkpoint'
            new_data = reformat_dims(
                data, format=data_format, towards_common=True)
            # Transpose to T, Z, F, Y, X
            if len(new_data.dims) == 5:
                new_data = new_data.transpose('T', 'Z', 'F', 'Y', 'X')
            elif len(new_data.dims) == 4:
                # no time
                new_data = new_data.transpose('Z', 'F', 'Y', 'X')
            elif len(new_data.dims) == 3:
                # no time or vertical
                new_data = new_data.transpose('F', 'Y', 'X')
            # For each output face, sum regridded input faces
            oface_datasets = []
            for oface in range(6):
                oface_regridded = []
                for iface, regridder in regridder_list[oface].items():
                    ds_iface = new_data.isel(F=iface)
                    if 'F' in ds_iface.coords:
                        ds_iface = ds_iface.drop('F')
                    oface_regridded.append(
                        regridder(ds_iface, keep_attrs=True))
                oface_regridded = xr.concat(
                    oface_regridded, dim='intersecting_ifaces').sum(
                    'intersecting_ifaces', keep_attrs=True)
                oface_datasets.append(oface_regridded)
            new_data = xr.concat(oface_datasets, dim='F')

            new_data = new_data.rename({
                'y': 'Y',
                'x': 'X',
            })
            # lat, lon are from xESMF which we don't want
            new_data = new_data.drop(['lat', 'lon'])

            # reformat dimensions to previous format
            new_data = reformat_dims(
                new_data, format=data_format, towards_common=False)
            return new_data
    else:
        return data


def reformat_dims(ds, format, towards_common):
    """
    Reformat dimensions of a cubed-sphere / stretched-grid grid between different GCHP formats

    Args:
        ds: xarray Dataset
             Dataset to be reformatted
        format: str
             Format from or to which to reformat ('checkpoint' or 'diagnostic')
        towards_common: bool
             Set this flag to True to move towards a common dimension format

    Returns:
        ds: xarray Dataset
             Original dataset with reformatted dimensions
    """
    def unravel_checkpoint_lat(ds_in):
        if isinstance(ds_in, xr.Dataset):
            cs_res = ds_in.dims['lon']
            assert cs_res == ds_in.dims['lat'] // 6
        else:
            cs_res = ds_in['lon'].size
            assert cs_res == ds_in['lat'].size // 6
        mi = pd.MultiIndex.from_product([
            np.linspace(1, 6, 6),
            np.linspace(1, cs_res, cs_res)
        ])
        ds_in = ds_in.assign_coords({'lat': mi})
        ds_in = ds_in.unstack('lat')
        return ds_in

    def ravel_checkpoint_lat(ds_out):
        if isinstance(ds, xr.Dataset):
            cs_res = ds_out.dims['lon']
        else:
            cs_res = ds_out['lon'].size
        ds_out = ds_out.stack(lat=['lat_level_0', 'lat_level_1'])
        ds_out = ds_out.assign_coords({
            'lat': np.linspace(1, 6 * cs_res, 6 * cs_res)
        })
        return ds_out

    dim_formats = {
        'checkpoint': {
            'unravel': [unravel_checkpoint_lat],
            'ravel': [ravel_checkpoint_lat],
            'rename': {
                'lon': 'X',
                'lat_level_0': 'F',
                'lat_level_1': 'Y',
                'time': 'T',
                'lev': 'Z',
            },
            'transpose': ('time', 'lev', 'lat', 'lon')
        },
        'diagnostic': {
            'rename': {
                'nf': 'F',
                'lev': 'Z',
                'Xdim': 'X',
                'Ydim': 'Y',
                'time': 'T',
            },
            'transpose': ('time', 'lev', 'nf', 'Ydim', 'Xdim')
        }
    }
    if towards_common:
        # Unravel dimensions
        for unravel_callback in dim_formats[format].get('unravel', []):
            ds = unravel_callback(ds)

        # Rename dimensions
        ds = ds.rename(dim_formats[format].get('rename', {}))

        return ds
    else:
        # Reverse rename
        ds = ds.rename(
            {v: k for k, v in dim_formats[format].get('rename', {}).items()})

        # Ravel dimensions
        for ravel_callback in dim_formats[format].get('ravel', []):
            ds = ravel_callback(ds)

        # Transpose
        if len(ds.dims) == 5 or (len(ds.dims) == 4 and 'lev' in list(
                ds.dims) and 'time' in list(ds.dims)):
            # full dim dataset
            ds = ds.transpose(*dim_formats[format].get('transpose', []))
        elif len(ds.dims) == 4:
            # single time
            ds = ds.transpose(*dim_formats[format].get('transpose', [])[1:])
        elif len(ds.dims) == 3:
            # single level / time
            ds = ds.transpose(*dim_formats[format].get('transpose', [])[2:])
        return ds


def sg_hash(
        cs_res,
        stretch_factor: float,
        target_lat: float,
        target_lon: float):
    return hashlib.sha1(
        'cs={cs_res},sf={stretch_factor:.5f},tx={target_lon:.5f},ty={target_lat:.5f}'.format(
            stretch_factor=stretch_factor,
            target_lat=target_lat,
            target_lon=target_lon,
            cs_res=cs_res).encode()).hexdigest()[
        :7]

def regrid_vertical_datasets(ref, dev, target_grid_choice='ref', ref_vert_params=[[],[]],
                             dev_vert_params=[[],[]], target_vert_params=[[],[]]):
    """
    Perform complete vertical regridding of GEOS-Chem datasets to 
    the vertical grid of one of the datasets or an entirely different 
    vertical grid.
    
    Args:
        ref: xarray.Dataset
            First dataset
        dev: xarray.Dataset
            Second dataset
        target_grid_choice (optional): str
            Will regrid to the chosen dataset among the two datasets
            unless target_vert_params is provided
            Default value: 'ref'
        ref_vert_params (optional): list(list, list) of list-like types
            Hybrid grid parameter A in hPa and B (unitless) in [AP, BP] format. 
            Needed if ref grid is not 47 or 72 levels
            Default value: [[], []]
        dev_vert_params (optional): list(list, list) of list-like types
            Hybrid grid parameter A in hPa and B (unitless) in [AP, BP] format. 
            Needed if dev grid is not 47 or 72 levels
            Default value: [[], []]
        target_vert_params (optional): list(list, list) of list-like types
            Hybrid grid parameter A in hPa and B (unitless) in [AP, BP] format. 
            Will override target_grid_choice as target grid
            Default value: [[], []]
    Returns:
        new_ref: xarray.Dataset
            First dataset, possibly regridded to a new vertical grid
        new_dev: xarray.Dataset
            Second dataset, possibly regridded to a new vertical grid
    """

    # Get mid-point pressure and edge pressures for this grid
    ref_pedge, ref_pmid, _ = get_vert_grid(ref, *ref_vert_params)
    dev_pedge, dev_pmid, _ = get_vert_grid(dev, *dev_vert_params)
    
    new_ref, new_dev = ref, dev
    
    if len(ref_pedge) != len(dev_pedge) or target_vert_params != [[],[]]:
        if target_vert_params != [[],[]]:
            #use a specific target grid for regridding if passed
            target_grid = vert_grid(*target_vert_params)
            target_pedge, target_pmid = target_grid.p_edge(), target_grid.p_mid()        
        elif target_grid_choice == 'ref':
            target_pedge, target_pmid = ref_pedge, ref_pmid
        else:
            target_pedge, target_pmid = dev_pedge, dev_pmid
        
        def regrid_one_vertical_dataset(ds, ds_pedge, target_pedge, target_pmid):
            new_ds = ds
            if len(ds_pedge) != len(target_pedge):
                #regrid all 3D (plus possible time dimension) variables
                xmat_ds = gen_xmat(ds_pedge, target_pedge)
                regrid_variables = [v for v in ds.data_vars if (("lat" in ds[v].dims or "Xdim" in ds[v].dims)
                                                                 and ("lon" in ds[v].dims or "Ydim" in ds[v].dims)
                                                                 and ("lev" in ds[v].dims))]
                new_ds = xr.Dataset()
                #currently drop data vars that have lev but don't also have x and y coordinates
                for v in (set(ds.data_vars)-set(regrid_variables)): 
                    if 'lev' not in ds[v].dims:
                        new_ds[v] = ds[v]
                new_ds.attrs = ds.attrs
                for v in regrid_variables:
                    if "time" in ds[v].dims:
                        new_ds_temp = []
                        for time in range(len(ds[v].time)):
                            new_ds_v = regrid_vertical(ds[v].isel(time=time), xmat_ds, target_pmid)
                            new_ds_temp.append(new_ds_v.expand_dims("time"))
                        new_ds[v] = xr.concat(new_ds_temp, "time")
                    else:
                        new_ds[v] = regrid_vertical(ds[v], xmat, target_pmid)
            return new_ds
        
        new_ref = regrid_one_vertical_dataset(ref, ref_pedge, target_pedge, target_pmid)
        new_dev = regrid_one_vertical_dataset(dev, dev_pedge, target_pedge, target_pmid)

    return new_ref, new_dev

def regrid_vertical(src_data_3D, xmat_regrid, target_levs=[]):
    """
    Performs vertical regridding using a sparse regridding matrix
    This function was originally written by Sebastian Eastham and included
    in package gcgridobj: https://github.com/sdeastham/gcgridobj

    Args:
        src_data_3D: xarray DataArray or numpy array
            Data to be regridded
        xmat_regrid: sparse scipy coordinate matrix
            Regridding matrix from input data grid to target grid
        target_levs (optional): list
            Values for Z coordinate of returned data (if returned data is of type xr.DataArray)
            Default value: []

    Returns:
        out_data: xarray DataArray or numpy array
            Data regridded to target grid
    """

    # Assumes that the FIRST dimension of the input data is vertical
    nlev_in = src_data_3D.shape[0]
    if xmat_regrid.shape[1] == nlev_in:
        # Current regridding matrix is for the reverse regrid
        # Rescale matrix to get the contributions right
        # Warning: this assumes that the same vertical range is covered
        warnings.warn(
            'Using inverted regridding matrix. This may cause incorrect extrapolation')
        xmat_renorm = xmat_regrid.transpose().toarray()
        for ilev in range(xmat_renorm.shape[1]):
            norm_fac = np.sum(xmat_renorm[:, ilev])
            if np.abs(norm_fac) < 1.0e-20:
                norm_fac = 1.0
            xmat_renorm[:, ilev] /= norm_fac

        xmat_renorm = scipy.sparse.coo_matrix(xmat_renorm)
    elif xmat_regrid.shape[0] == nlev_in:
        # Matrix correctly dimensioned
        xmat_renorm = xmat_regrid.copy()
    else:
        print(src_data_3D, xmat_regrid.shape)
        raise ValueError('Regridding matrix not correctly sized')

    nlev_out = xmat_renorm.shape[1]
    out_shape = [nlev_out] + list(src_data_3D.shape[1:])
    n_other = np.product(src_data_3D.shape[1:])
    temp_data = np.zeros((nlev_out, n_other))
    in_data = np.reshape(np.array(src_data_3D), (nlev_in, n_other))
    for ix in range(n_other):
        in_data_vec = np.matrix(in_data[:, ix])
        temp_data[:, ix] = in_data_vec * xmat_renorm
    out_data = np.reshape(temp_data, out_shape)

    # Transfer over old / create new coordinates for xarray DataArrays
    if isinstance(src_data_3D, xr.DataArray):
        new_coords = {
            coord: src_data_3D.coords[coord].data
            for coord in src_data_3D.coords if coord != 'lev'}
        if target_levs == []:
            new_coords['lev'] = np.arange(
                1, out_data.shape[0], out_data.shape[0])
        else:
            new_coords['lev'] = target_levs
        # GCHP-specific
        if 'lats' in src_data_3D.coords:
            new_coords['lats'] = (
                ('lat', 'lon'), src_data_3D.coords['lats'].data)
        if 'lons' in src_data_3D.coords:
            new_coords['lons'] = (
                ('lat', 'lon'), src_data_3D.coords['lons'].data)
        out_data = xr.DataArray(out_data,
                                dims=tuple([dim for dim in src_data_3D.dims]),
                                coords=new_coords,
                                attrs=src_data_3D.attrs)

    return out_data


def gen_xmat(p_edge_from, p_edge_to):
    """
    Generates regridding matrix from one vertical grid to another.
    This function was originally written by Sebastian Eastham and included
    in package gcgridobj: https://github.com/sdeastham/gcgridobj

    Args:
        p_edge_from: numpy array
            Edge pressures of the input grid
        p_edge_to: numpy array
            Edge pressures of the target grid

    Returns:
        xmat: sparse scipy coordinate matrix
            Regridding matrix from input grid to target grid
    """
    n_from = len(p_edge_from) - 1
    n_to = len(p_edge_to) - 1

    # Guess - max number of entries?
    n_max = max(n_to, n_from) * 5

    # Index being mapped from
    xmat_i = np.zeros(n_max)
    # Index being mapped to
    xmat_j = np.zeros(n_max)
    # Weights
    xmat_s = np.zeros(n_max)

    # Find the first output box which has any commonality with the input box
    first_from = 0
    i_to = 0
    if p_edge_from[0] > p_edge_to[0]:
        # "From" grid starts at lower altitude (higher pressure)
        while p_edge_to[0] < p_edge_from[first_from + 1]:
            first_from += 1
    else:
        # "To" grid starts at lower altitude (higher pressure)
        while p_edge_to[i_to + 1] > p_edge_from[0]:
            i_to += 1

    frac_to_total = 0.0

    i_weight = 0
    for i_from in range(first_from, n_from):
        p_base_from = p_edge_from[i_from]
        p_top_from = p_edge_from[i_from + 1]

        # Climb the "to" pressures until you intersect with this box
        while i_to < n_to and p_base_from <= p_edge_to[i_to + 1]:
            i_to += 1
            frac_to_total = 0.0

        # Now, loop over output layers as long as there is any overlap,
        # i.e. as long as the base of the "to" layer is below the
        # top of the "from" layer
        last_box = False

        while p_edge_to[i_to] >= p_top_from and not last_box and not i_to >= n_to:
            p_base_common = min(p_base_from, p_edge_to[i_to])
            p_top_common = max(p_top_from, p_edge_to[i_to + 1])
            # Fraction of target box
            frac_to = (p_base_common - p_top_common) / \
                (p_edge_to[i_to] - p_edge_to[i_to + 1])

            xmat_i[i_weight] = i_from
            xmat_j[i_weight] = i_to
            xmat_s[i_weight] = frac_to

            i_weight += 1
            last_box = p_edge_to[i_to + 1] <= p_top_from
            if not last_box:
                i_to += 1

    return scipy.sparse.coo_matrix(
        (xmat_s[: i_weight],
         (xmat_i[: i_weight],
          xmat_j[: i_weight])),
        shape=(n_from, n_to))
