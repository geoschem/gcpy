''' Functions for creating xesmf regridder objects '''

import os
import xesmf as xe
from .grid import make_grid_LL, make_grid_CS, get_input_res, call_make_grid, get_grid_extents
import numpy as np

def make_regridder_L2L( llres_in, llres_out, weightsdir='.', reuse_weights=False,
                        in_extent=[-180,180,-90,90], out_extent=[-180,180,-90,90]):
    [in_minlon, in_maxlon, in_minlat, in_maxlat] = in_extent
    [out_minlon, out_maxlon, out_minlat, out_maxlat] = out_extent
    llgrid_in = make_grid_LL(llres_in, in_extent, out_extent)
    llgrid_out = make_grid_LL(llres_out, out_extent)
    if in_extent == [-180,180,-90,90] and out_extent == [-180,180,-90,90]:
        weightsfile = os.path.join(weightsdir,'conservative_{}_{}.nc'.format(llres_in, llres_out))
    else:
        in_extent_str = str(in_extent).replace('[', '').replace(']','').replace(', ', 'x')
        out_extent_str = str(out_extent).replace('[', '').replace(']','').replace(', ', 'x')
        weightsfile = os.path.join(weightsdir,'conservative_{}_{}_{}_{}.nc'.format(llres_in, llres_out, 
                                                                                   in_extent_str, out_extent_str))
    regridder = xe.Regridder(llgrid_in, llgrid_out, method='conservative', filename=weightsfile, reuse_weights=reuse_weights)
    return regridder

def make_regridder_C2L( csres_in, llres_out, weightsdir='.', reuse_weights=True ):
    csgrid, csgrid_list = make_grid_CS(csres_in)
    llgrid = make_grid_LL(llres_out)
    regridder_list = []
    for i in range(6):
        weightsfile = os.path.join(weightsdir, 'conservative_c{}_{}_{}.nc'.format(str(csres_in), llres_out, str(i)))
        regridder = xe.Regridder(csgrid_list[i], llgrid, method='conservative', filename=weightsfile, reuse_weights=reuse_weights)
        regridder_list.append(regridder)
    return regridder_list

def create_regridders(refds, devds, weightsdir='.', reuse_weights=True, cmpres=None, zm=False):
    #Take two lat/lon or cubed-sphere xarray datasets and regrid them if needed
    refres, refgridtype = get_input_res(refds)
    devres, devgridtype = get_input_res(devds)
    
    refminlon, refmaxlon, refminlat, refmaxlat = get_grid_extents(refds)
    devminlon, devmaxlon, devminlat, devmaxlat = get_grid_extents(devds)
    #choose smallest overlapping area for comparison 
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
    
    if cmpres == None:
        if refres == devres and refgridtype == "ll":
            cmpres = refres
            cmpgridtype = refgridtype
        elif refgridtype == "ll" and devgridtype == "ll":
            cmpres = min([refres, devres])
            cmpgridtype = refgridtype
        elif refgridtype == "cs" and devgridtype == "cs":
            # CS to CS regridding is not enabled yet, so default to 1x1.25
            # cmpres = max([refres, devres])
            # cmpgridtype = 'cs'
            if refres==devres and not zm:
                cmpres=int(refres)
                cmpgridtype="cs"
            else:
                cmpres = "1x1.25"
                cmpgridtype = "ll"
        elif refgridtype == "ll" and float(refres.split('x')[0])<1 and float(refres.split('x')[1])<1.25:
            cmpres = refres
            cmpgridtype = "ll"
        elif devgridtype == "ll" and float(devres.split('x')[0])<1 and float(devres.split('x')[1])<1.25:
            cmpres = devres
            cmpgridtype = "ll"
        else:
            cmpres = "1x1.25"
            cmpgridtype = "ll"

    elif "x" in cmpres:
        cmpgridtype = "ll"
    else:
        if zm:
            print("Warning: zonal mean comparison must be lat-lon. Defaulting to 1x1.25")
            cmpres='1x1.25'
        else:
            cmpgridtype = "cs"
            cmpres = int(cmpres)  # must cast to integer for cubed-sphere
        

    # Determine what, if any, need regridding.
    regridref = refres != cmpres
    regriddev = devres != cmpres
    regridany = regridref or regriddev
    # ==================================================================
    # Make grids (ref, dev, and comparison)
    # ==================================================================
    [refgrid, regrid_list] = call_make_grid(refres, refgridtype, zm, False, ref_extent, cmp_extent)
    
    [devgrid, devgrid_list] = call_make_grid(devres, devgridtype, zm, False, dev_extent, cmp_extent)

    [cmpgrid, cmpgrid_list] = call_make_grid(cmpres, cmpgridtype, zm, True, cmp_extent, cmp_extent)
    
    # =================================================================
    # Make regridders, if applicable
    # TODO: Make CS to CS regridders
    # =================================================================


    msg = "CS to CS regridding is not yet implemented in gcpy. " \
        + "Ref and dev cubed sphere grids must be the same resolution, " \
        + "or pass cmpres to compare_single_level as a lat-lon grid resolution."
    refregridder = None
    refregridder_list = None
    devregridder = None
    devregridder_list = None
    if regridref:
        if refgridtype == "ll":
            refregridder = make_regridder_L2L(
                refres, cmpres, weightsdir=weightsdir, reuse_weights=reuse_weights,
                in_extent = ref_extent, out_extent = cmp_extent
            )
        else:
            if cmpgridtype == "cs":
                raise ValueError(msg)
            else:
                refregridder_list = make_regridder_C2L(
                    refres, cmpres, weightsdir=weightsdir, reuse_weights=reuse_weights
                )
    if regriddev:
        if devgridtype == "ll":
            devregridder = make_regridder_L2L(
                devres, cmpres, weightsdir=weightsdir, reuse_weights=reuse_weights,
                in_extent = dev_extent, out_extent = cmp_extent
            )
        else:
            if cmpgridtype == "cs":
                raise ValueError(msg)
            else:
                devregridder_list = make_regridder_C2L(
                    devres, cmpres, weightsdir=weightsdir, reuse_weights=reuse_weights
                )


    return [refres, refgridtype, devres, devgridtype, cmpres, cmpgridtype,
    regridref, regriddev, regridany, refgrid, devgrid, cmpgrid, refregridder, 
    devregridder, refregridder_list, devregridder_list]

def dict_72_to_47():
    """
    Get 72L-to-47L conversion dict which stores weights from 72 levels to 47 levels
    Returns:
    --------
        conv_dict : dict {72L (int) : (47L (int), weight (int))}
              Mapping of 72L to 47L
    """

    conv_dict = {L72 : (xmat_72to47.col[L72], xmat_72to47.data[L72]) for L72 in xmat_72to47.row}
    return conv_dict

def reduce_72_to_47(DataArray, conv_dict, pmid_ind_72, pmid_ind_47):
    """
    Reduce 72 level DataArray to 47 level-equivalent for full or restricted
    pressure ranges.
    Args:
    -----
        DataArray : xarray DataArray
            72 level DataArray
    
        conv_dict : dict
            Mapping of 72L to 47L
        pmid_ind_72 : list(int)
            List of midpoint indices for 72L grid
        pmid_ind_47 : list(int)
            List of midpoint indices for 47L grid
        
    Returns:
    --------
         xarray DataArray
            DataArray now reduced to 47L grid
    """

    #reduce 72 level DataArray to 47 level-equivalent
    #This function works for both full and restricted pressure ranges
    new_shape = list(DataArray.data.shape)
    #assumes first dim is level
    new_shape[0] = len(pmid_ind_47)
    reduced_offset = min(pmid_ind_47)
    reduced_data = np.zeros(new_shape)
    for i in range(0, len(pmid_ind_72)):
        lev = pmid_ind_72[i]
        reduced_data[conv_dict[lev][0]-reduced_offset] = reduced_data[conv_dict[lev][0]-reduced_offset] + DataArray.data[i]*conv_dict[lev][1]
    new_coords = {coord : DataArray.coords[coord].data for coord in DataArray.coords if coord != 'lev'}
    new_coords['lev'] = np.array(pmid_ind_47)
    #GCHP-specific
    if 'lats' in DataArray.coords:
        new_coords['lats'] = (('lon', 'lat'), DataArray.coords['lats'].data)
    if 'lons' in DataArray.coords:
        new_coords['lons'] = (('lon', 'lat'), DataArray.coords['lons'].data)
    return xr.DataArray(reduced_data, dims=tuple([dim for dim in DataArray.dims]),
                        coords = new_coords, attrs = DataArray.attrs)

def regrid_comparison_data(data, res, regrid, regridder, regridder_list, global_cmp_grid, gridtype,
                           cmpminlat_ind=0, cmpmaxlat_ind=-2, cmpminlon_ind=0, cmpmaxlon_ind=-2, nlev=1):
    """
    Regrid comparison datasets to lat/lon format.
    Args:
    -----
        data : xarray DataArray
            DataArray containing a GEOS-Chem data variable
        res : int
            Cubed-sphere resolution for comparison grid
        
        regrid : boolean
            Set to true to regrid dataset
        regridder : xESMF regridder
            Regridder between the original data grid and the comparison grid
     
        regridder_list : list(xESMFW regridder)
            List of regridders for cubed-sphere data
        global_cmp_grid : xarray DataArray
            Comparison grid
    
        gridtype : str
            Type of input data grid (either 'll' or 'cs')
        
        cmpminlat_ind : int
            Index of minimum latitude extent for comparison grid
        cmpmaxlat_ind : int
            Index (minus 1) of maximum latitude extent for comparison grid
        cmpminlon_ind : int
            Index of minimum longitude extent for comparison grid
        cmpmaxlon_ind : int
            Index (minus 1) of maximum longitude extent for comparison grid
        nlev : int
            Number of levels of input grid and comparison grid
    Returns:
    --------
        data : xarray DataArray
            Original DataArray regridded to comparison grid (including resolution and extent changes)
    """

    if regrid:
        if gridtype == "ll":
            # regrid ll to ll
            return regridder(data)
        else:
            if nlev is 1:
                new_data = np.zeros([global_cmp_grid['lat'].size,
                                     global_cmp_grid['lon'].size])
                data_reshaped = data.data.reshape(6, res, res)
            else:
                new_data = np.zeros([nlev, global_cmp_grid['lat'].size,
                                     global_cmp_grid['lon'].size])
                data_reshaped = data.data.reshape(nlev, 6, res, res).swapaxes(0, 1)
            for j in range(6):
                regridder = regridder_list[j]
                new_data += regridder(data_reshaped[j])
            if nlev is 1:
                #limit to extent of cmpgrid
                return new_data[cmpminlat_ind:cmpmaxlat_ind+1,cmpminlon_ind:cmpmaxlon_ind+1].squeeze()
            else:
                return new_data
    else:
        return data
