"""
Unit tests for gcpy.grid, focused on the "pole cutoff" latitude grid
handling in get_grid_extents() and make_grid_ll() implicated in
GitHub issue #425 (zonal mean plots of sub-1-degree-resolution
lat/lon grids, e.g. HEMCO diagnostics, failed with a pcolormesh
shape-mismatch error).
"""
import numpy as np
import xarray as xr

from gcpy.grid import get_grid_extents, make_grid_ll
from gcpy.plot.single_panel import single_panel

# Reproduce the pole-cutoff 0.5x0.625 grid described in the issue:
# lat centers [-89.875, -89.5, ..., 89.5, 89.875] (a narrower boundary
# band than the 0.5 degree interior spacing), 361 centers total.
LLRES = "0.5x0.625"
_INTERIOR_LAT = np.round(np.arange(-89.5, 89.5 + 1e-6, 0.5), 8)
POLE_CUTOFF_LAT = np.concatenate(([-89.875], _INTERIOR_LAT, [89.875]))
POLE_CUTOFF_LON = np.round(np.arange(-180, 180, 0.625), 8)


def test_get_grid_extents_snaps_pole_cutoff_lat_to_90():
    """
    get_grid_extents() should snap a detected pole-cutoff latitude
    band to the true pole (+/-90) rather than an arbitrary
    resolution-independent offset (GitHub issue #425).
    """
    ds = xr.Dataset(
        coords={"lat": POLE_CUTOFF_LAT, "lon": POLE_CUTOFF_LON})

    minlon, maxlon, minlat, maxlat = get_grid_extents(ds)

    assert (minlon, maxlon, minlat, maxlat) == (-180, 180, -90, 90)


def test_make_grid_ll_pole_cutoff_extent_has_no_duplicate_edges():
    """
    make_grid_ll() should return the correct number of latitude edges
    (362, i.e. 361 cell centers) for a 0.5x0.625 grid, with no
    duplicate/degenerate edges at the poles. Before the fix, an
    overshot pole-cutoff extent of (-90.875, 90.875) produced 365
    edges with triplicate values at each pole (GitHub issue #425).
    """
    extent = list(get_grid_extents(
        xr.Dataset(coords={"lat": POLE_CUTOFF_LAT, "lon": POLE_CUTOFF_LON})
    ))

    llgrid = make_grid_ll(LLRES, extent, extent)

    assert llgrid["lat_b"].shape == (362,)
    assert llgrid["lat"].shape == (361,)
    assert np.unique(llgrid["lat_b"]).size == llgrid["lat_b"].size


def test_single_panel_zonal_mean_pole_cutoff_grid():
    """
    single_panel(..., plot_type="zonal_mean") should not raise a
    pcolormesh shape-mismatch error when plotting data on a
    pole-cutoff lat/lon grid (GitHub issue #425). Uses a 47-level
    grid so get_vert_grid() resolves pressure edges without requiring
    explicit AP/BP parameters.
    """
    nlev = 47
    data = xr.DataArray(
        np.random.rand(nlev, POLE_CUTOFF_LAT.size, POLE_CUTOFF_LON.size),
        dims=["lev", "lat", "lon"],
        coords={
            "lev": np.arange(1, nlev + 1),
            "lat": POLE_CUTOFF_LAT,
            "lon": POLE_CUTOFF_LON,
        },
    )

    plot = single_panel(data, plot_type="zonal_mean")

    assert plot is not None
