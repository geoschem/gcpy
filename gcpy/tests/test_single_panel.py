"""
Unit tests for gcpy/plot/single_panel.py.

Covers the near-constant / noise tolerances a standalone panel applies
when it builds its own norm (norm=None).  These were not reaching
normalize_colors, so a standalone panel did not render the way the
corresponding panel of a six-panel plot does, which the comment at
that call site says is the intent.
"""
import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from gcpy.plot.core import CONSTANT_REL_TOL, constant_rel_tol
from gcpy.plot.single_panel import single_panel

NLEV, NLAT = 6, 8
LAT_B = np.linspace(-90.0, 90.0, NLAT + 1)
PEDGE = np.linspace(1000.0, 100.0, NLEV + 1)
PEDGE_IND = np.arange(NLEV + 1)

# Magnitude of a PassiveTracer field [ppb], and a Dev - Ref spread
# that is negligible against it (3.5e-4 / 100.06 = 3.5e-6 relative)
DATA_SCALE = 100.06
NEGLIGIBLE_DIFF = 3.5e-4


def _zonal_mean_norm(values, **kwargs):
    """Returns the (vmin, vmax) of the norm single_panel builds."""
    darr = xr.DataArray(values, dims=("lev", "lat"))
    plt.figure()
    try:
        plot = single_panel(
            darr,
            plot_type="zonal_mean",
            grid={"lat_b": LAT_B},
            gridtype="ll",
            extent=(-180.0, 180.0, -90.0, 90.0),
            pedge=PEDGE,
            pedge_ind=PEDGE_IND,
            yaxis_units="level",
            add_cb=False,
            **kwargs,
        )
        return float(plot.norm.vmin), float(plot.norm.vmax)
    finally:
        plt.close("all")


def _difference_field():
    arr = np.full((NLEV, NLAT), NEGLIGIBLE_DIFF, dtype=np.float32)
    arr[0, 0] = -NEGLIGIBLE_DIFF
    return arr


def test_difference_panel_collapses_when_given_a_data_scale():
    # With the magnitude of the data it was differenced from, the
    # panel can tell that its spread is only noise.  (-1, 1) is the
    # collapsed sentinel normalize_colors returns for a difference.
    norm = _zonal_mean_norm(
        _difference_field(),
        use_cmap_RdBu=True,
        data_scale=DATA_SCALE,
    )
    assert norm == (-1.0, 1.0)


def test_difference_panel_keeps_its_range_without_a_data_scale():
    # Without it there is no basis for calling the spread negligible,
    # so the panel must keep its real range rather than guess.
    vmin, vmax = _zonal_mean_norm(
        _difference_field(),
        use_cmap_RdBu=True,
    )
    assert (vmin, vmax) != (-1.0, 1.0)
    assert np.isclose(vmax, NEGLIGIBLE_DIFF, rtol=1e-6)


def test_ref_dev_panel_uses_the_precision_of_its_own_data():
    # A real*8 field carrying 1e-7 of real relative structure sits
    # below the real*4 tolerance but above the real*8 one, so it must
    # survive.  Before rel_tol was threaded through, this panel got
    # the real*4 fallback and collapsed to a blank 0-1 scale.
    rel = 1.0e-7
    assert constant_rel_tol(np.zeros(1, np.float64)) < rel < CONSTANT_REL_TOL

    base = 1.0
    arr = np.linspace(base * (1.0 - rel), base, NLEV * NLAT,
                      dtype=np.float64).reshape(NLEV, NLAT)
    vmin, vmax = _zonal_mean_norm(arr)
    assert (vmin, vmax) != (0.0, 1.0), "collapsed a real*8 field"
    assert np.isclose(vmax, base)
