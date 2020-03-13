''' Functions for creating xesmf regridder objects '''

import os
import xesmf as xe
from .horiz import make_grid_LL, make_grid_CS

def make_regridder_L2L( llres_in, llres_out, weightsdir='.', reuse_weights=False,
minlon=-180, maxlon=180, minlat=-90, maxlat=90 ):
    llgrid_in = make_grid_LL(llres_in, minlon, maxlon, minlat, maxlat)
    llgrid_out = make_grid_LL(llres_out, minlon, maxlon, minlat, maxlat)
    weightsfile = os.path.join(weightsdir,'conservative_{}_{}.nc'.format(llres_in, llres_out))
    regridder = xe.Regridder(llgrid_in, llgrid_out, method='conservative', filename=weightsfile, reuse_weights=reuse_weights)
    return regridder

def make_regridder_C2L( csres_in, llres_out, weightsdir='.', reuse_weights=False ):
    csgrid, csgrid_list = make_grid_CS(csres_in)
    llgrid = make_grid_LL(llres_out)
    regridder_list = []
    for i in range(6):
        weightsfile = os.path.join(weightsdir, 'conservative_c{}_{}_{}.nc'.format(str(csres_in), llres_out, str(i)))
        regridder = xe.Regridder(csgrid_list[i], llgrid, method='conservative', filename=weightsfile, reuse_weights=reuse_weights)
        regridder_list.append(regridder)
    return regridder_list
