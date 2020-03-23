''' Functions for creating xesmf regridder objects '''

import os
import xesmf as xe
from .horiz import make_grid_LL, make_grid_CS
from ..core import get_input_res, call_make_grid, get_grid_extents


def make_regridder_L2L( llres_in, llres_out, weightsdir='.', reuse_weights=False,
                        minlon=-180, maxlon=180, minlat=-90, maxlat=90 ):
    llgrid_in = make_grid_LL(llres_in, minlon, maxlon, minlat, maxlat)
    llgrid_out = make_grid_LL(llres_out, minlon, maxlon, minlat, maxlat)
    weightsfile = os.path.join(weightsdir,'conservative_{}_{}.nc'.format(llres_in, llres_out))
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

    [refgrid, regrid_list] = call_make_grid(refres, refgridtype, True, False, 
                                            minlon=refminlon, maxlon=refmaxlon,
                                            minlat=refminlat, maxlat=refmaxlat)
    [devgrid, devgrid_list] = call_make_grid(devres, devgridtype, True, False,
                                            minlon=devminlon, maxlon=devmaxlon,
                                            minlat=devminlat, maxlat=devmaxlat)
    [cmpgrid, cmpgrid_list] = call_make_grid(cmpres, cmpgridtype, True, True)

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
                refres, cmpres, weightsdir=weightsdir, reuse_weights=reuse_weights
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
                devres, cmpres, weightsdir=weightsdir, reuse_weights=reuse_weights
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

