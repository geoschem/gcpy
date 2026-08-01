"""
Unit tests for gcpy.regrid, focused on the cubed-sphere to lat/lon
conservative-regridding path implicated in GitHub issue #330 (zonal
mean plots of constant fields showed color striping caused by tiny
numerical noise introduced during regridding).
"""
import numpy as np
import xarray as xr

from gcpy.grid import make_grid_ll
from gcpy.regrid import make_regridder_cs2ll, regrid_comparison_data

# Small cubed-sphere resolution so that regridder weights can be
# computed quickly within a unit test.  Use the same comparison grid
# (1x1.25) as the sample data in GitHub issue #330, since the
# magnitude of the pre-fix regridding noise is grid-dependent.
CSRES = 24
LLRES = "1x1.25"
CONST_VALUE = 100.00064087


def test_constant_cs_field_regrids_to_constant(tmp_path):
    """
    Regridding a perfectly constant cubed-sphere field to a lat/lon
    grid should return a field that is constant to within a very
    tight relative tolerance.  Before the fix in
    regrid_comparison_data() (renormalizing the summed per-face
    conservative-regridding weights), this could differ by as much
    as ~1e-6 relative, which produced visible "striping" once
    plotted (GitHub issue #330).
    """
    llgrid = make_grid_ll(LLRES)
    regridder_list = make_regridder_cs2ll(
        CSRES,
        LLRES,
        weightsdir=str(tmp_path),
        reuse_weights=False,
    )
    
    data = xr.DataArray(np.full((6 * CSRES, CSRES), CONST_VALUE))
    result = regrid_comparison_data(
        data,
        CSRES,
        True,
        None,
        regridder_list,
        llgrid,
        "cs",
        "ll",
        nlev=1,
    )
    
    finite = result[np.isfinite(result)]
    assert np.allclose(finite, CONST_VALUE, rtol=1e-14, atol=1e-13)
