''' Functions for creating xesmf regridder objects '''

import os
import xesmf as xe
from .horiz import make_grid_LL, make_grid_CS
from ..core import get_input_res, call_make_grid, get_grid_extents
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
    [refgrid, regrid_list] = call_make_grid(refres, refgridtype, True, False, ref_extent, cmp_extent)
    
    [devgrid, devgrid_list] = call_make_grid(devres, devgridtype, True, False, dev_extent, cmp_extent)

    [cmpgrid, cmpgrid_list] = call_make_grid(cmpres, cmpgridtype, True, True, cmp_extent, cmp_extent)
    
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

