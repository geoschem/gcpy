"""
Unit tests covering the plotting-side "near-constant" tolerance guard
added for GitHub issue #330 (zonal mean plots of constant fields
showed color striping because a colorbar range computed from
noisy-but-nominally-constant data was not collapsed to a flat scale).

Also covers three follow-on fixes for the same over-eager tolerance
collapsing real, non-constant, small-magnitude fields to a blank/
white panel: PR #439 (single-level plots, e.g. HEMCO emissions),
zonal-mean Ref/Dev panels (e.g. TransportTracers species like Rn222),
and zonal-mean difference panels (e.g. EmisCO_Aircraft), where the
absolute tolerance is now scaled to the magnitude of the Ref and Dev
data instead of being a fixed number in unknown units.
"""
import numpy as np
import pytest
from matplotlib import colors

from gcpy.util import get_nan_mask, is_nearly_constant
import matplotlib.pyplot as plt

from gcpy.plot.core import CONSTANT_REL_TOL, NOISE_REL_TOL, \
    constant_rel_tol, diff_is_negligible, mask_meaningless_ratio, \
    noise_atol, normalize_colors
from gcpy.plot.six_plot import (
    colorbar_for_all_zero_or_nan,
    vmin_vmax_for_absdiff_plots,
    colorbar_ticks_and_format,
    colorbar_for_constant_field,
    compute_norm_for_plot,
    ref_dev_data_scale,
    ref_equals_dev,
)

# Magnitude of the regridding noise reported in GitHub issue #330
# (~5e-9 relative on a field with a value of ~100).
NOISY_LO = 100.00064076
NOISY_HI = 100.00064126

# Order-of-magnitude of a real (non-noise) HEMCO emissions flux, e.g.
# kg/m2/s, that is small in absolute terms but spatially varying.
EMIS_VMIN = 0.0
EMIS_VMAX = 3.0e-10

# Order-of-magnitude of a real (non-noise) TransportTracers zonal-mean
# field (e.g. SpeciesConcVV_Rn222 in ppb, or WetLossLS_Strat), which is
# small enough in absolute terms to trip a fixed absolute tolerance
# tuned for issue #330's ~100-magnitude test case.
TRACER_VMIN = 0.0
TRACER_VMAX = 5.0e-11

# Order-of-magnitude of a real (non-noise) zonal-mean emissions field
# and of its Dev - Ref difference: EmisCO_Aircraft (kg/m2/s), whose
# differences are ~10% of the field but ~1e-14 in absolute terms.
AIRCRAFT_SCALE = 2.0e-13
AIRCRAFT_DIFF = 1.0e-14

# Spread of the absdiff panel produced by the issue #330 regridding
# noise, i.e. two independently-noisy copies of a field of ~100.
NOISY_DIFF = NOISY_HI - NOISY_LO          # 5.0e-7

# Granularity of a float32 field of magnitude ~100.  GEOS-Chem writes
# float32 diagnostics, so Dev - Ref for a genuinely identical field
# can be this large in absolute terms.  Note that the old fixed
# atol=1e-7 was far too small to recognize this as noise.
FLOAT32_DIFF = 1.5e-5


def test_is_nearly_constant_true_for_regridding_noise():
    assert is_nearly_constant([NOISY_LO, NOISY_HI])


def test_is_nearly_constant_false_for_real_difference():
    assert not is_nearly_constant([1.0, 1.2])


def test_is_nearly_constant_true_for_all_nan():
    assert is_nearly_constant([np.nan, np.nan])


def test_is_nearly_constant_with_target():
    assert is_nearly_constant([1.0 + 1e-7, 1.0 - 1e-7, np.nan], target=1.0)
    assert not is_nearly_constant([1.5, 1.6], target=1.0)


def test_normalize_colors_collapses_near_constant_range():
    norm = normalize_colors(NOISY_LO, NOISY_HI)
    assert isinstance(norm, colors.Normalize)
    assert (norm.vmin, norm.vmax) == (0.0, 1.0)


def test_normalize_colors_collapses_near_constant_difference_range():
    # A difference panel collapses only when its spread is negligible
    # *relative to the Ref/Dev data it came from*, so the caller has
    # to supply that magnitude.  Here it is issue #330's ~100.
    norm = normalize_colors(
        -1e-11, 1e-11, is_difference=True, data_scale=NOISY_HI
    )
    assert (norm.vmin, norm.vmax) == (-1.0, 1.0)


def test_normalize_colors_does_not_collapse_real_range():
    norm = normalize_colors(1.0, 1.2)
    assert (norm.vmin, norm.vmax) == (1.0, 1.2)


def test_ref_equals_dev_tolerant_to_noise():
    noisy_ratio = np.array([1.0 + 1e-7, 1.0 - 1e-7, np.nan])
    assert ref_equals_dev(noisy_ratio)


def test_ref_equals_dev_false_for_real_difference():
    real_ratio = np.array([1.0, 1.5, np.nan])
    assert not ref_equals_dev(real_ratio)


def test_normalize_colors_collapses_near_constant_ratio_range_around_one():
    # Regression test: a near-constant ratio range (Ref ~= Dev) used
    # to fall through to the 0-centered difference collapse, which
    # placed the real data value (~1.0) at the extreme edge of the
    # colorbar instead of the middle, and displaced the "Ref and Dev
    # equal" ticklabel (always positioned at data value 1.0) to the
    # far right instead of the center.
    norm = normalize_colors(
        0.9999999, 1.0000001,
        is_difference=True, log_color_scale=True, ratio_log=True,
        is_ratio=True,
    )
    assert np.isclose(norm(1.0), 0.5)


def test_normalize_colors_default_preserves_small_emissions_like_range():
    # A small but real (non-noise) HEMCO flux range must not collapse
    # to a flat/blank color scale, even with the tolerance guard
    # enabled (the default), since Ref/Dev panels no longer use a
    # fixed absolute tolerance (see comments in normalize_colors()).
    norm = normalize_colors(EMIS_VMIN, EMIS_VMAX)
    assert (norm.vmin, norm.vmax) == (EMIS_VMIN, EMIS_VMAX)


def test_normalize_colors_no_tolerance_preserves_small_real_range():
    # This is the fix contract: single-level plots pass
    # use_tolerance=False, so a real (non-noise) small-magnitude range
    # like an emissions flux is preserved instead of being collapsed.
    norm = normalize_colors(EMIS_VMIN, EMIS_VMAX, use_tolerance=False)
    assert (norm.vmin, norm.vmax) == (EMIS_VMIN, EMIS_VMAX)


def test_normalize_colors_no_tolerance_still_collapses_exact_zero():
    # Even with use_tolerance=False, a genuinely all-zero field must
    # still collapse to a flat scale (guards against e.g. LogNorm
    # errors on all-zero/all-NaN data for single-level plots).
    norm = normalize_colors(0.0, 0.0, use_tolerance=False)
    assert (norm.vmin, norm.vmax) == (0.0, 1.0)

    norm_nan = normalize_colors(np.nan, np.nan, use_tolerance=False)
    assert (norm_nan.vmin, norm_nan.vmax) == (0.0, 1.0)


def test_normalize_colors_no_tolerance_does_not_collapse_regridding_noise():
    # Without the tolerance guard, single-level plots no longer treat
    # tiny CS->LL regridding noise as "exactly constant" -- this is the
    # deliberate, narrower scope tradeoff: that protection remains for
    # zonal-mean plots only (see the compute_norm_for_plot tests below).
    norm = normalize_colors(NOISY_LO, NOISY_HI, use_tolerance=False)
    assert (norm.vmin, norm.vmax) == (NOISY_LO, NOISY_HI)


def test_compute_norm_for_plot_single_level_ref_preserves_small_emissions_range():
    # End-to-end regression test for the reported bug: a single-level
    # "ref" subplot (e.g. Total_Emissions.pdf) with a small but real
    # vmin/vmax must not collapse to a blank/flat panel.
    plot_val, norm = compute_norm_for_plot(
        np.array([EMIS_VMIN, EMIS_VMAX]),
        EMIS_VMIN,
        EMIS_VMAX,
        "ref",
        plot_type="single_level",
    )
    assert (norm.vmin, norm.vmax) == (EMIS_VMIN, EMIS_VMAX)


def test_compute_norm_for_plot_zonal_mean_ref_still_collapses_noise():
    # Zonal-mean plots must keep today's tolerance-based behavior, so
    # the original #330 striping fix is unaffected by this change.
    plot_val, norm = compute_norm_for_plot(
        np.array([NOISY_LO, NOISY_HI]),
        NOISY_LO,
        NOISY_HI,
        "ref",
        plot_type="zonal_mean",
    )
    assert (norm.vmin, norm.vmax) == (0.0, 1.0)


def test_normalize_colors_zonal_mean_ref_preserves_small_real_range():
    # Regression test for the blank/white Ref & Dev zonal-mean panels
    # reported for TransportTracers fields (e.g. SpeciesConcVV_Rn222,
    # WetLossLS_Strat): a small but real (non-noise) range must not be
    # collapsed just because use_tolerance=True (the zonal-mean default).
    norm = normalize_colors(TRACER_VMIN, TRACER_VMAX)
    assert (norm.vmin, norm.vmax) == (TRACER_VMIN, TRACER_VMAX)


def test_compute_norm_for_plot_zonal_mean_ref_preserves_small_tracer_range():
    # Same fix, but end-to-end through compute_norm_for_plot(), the way
    # it's actually called for a zonal-mean "ref" subplot (e.g.
    # SpeciesConcVV_Rn222_ZonalMean.pdf, WetLossLS_Strat_ZonalMean.pdf).
    plot_val, norm = compute_norm_for_plot(
        np.array([TRACER_VMIN, TRACER_VMAX]),
        TRACER_VMIN,
        TRACER_VMAX,
        "ref",
        plot_type="zonal_mean",
    )
    assert (norm.vmin, norm.vmax) == (TRACER_VMIN, TRACER_VMAX)


# ----------------------------------------------------------------------
# Difference panels: the absolute tolerance is scaled to the magnitude
# of the Ref & Dev data being differenced, rather than being a fixed
# number in the (unknown) units of the field.
# ----------------------------------------------------------------------

def test_noise_atol_scales_with_data_scale():
    assert noise_atol(100.0) == NOISE_REL_TOL * 100.0
    assert noise_atol(2.0e-13) == NOISE_REL_TOL * 2.0e-13
    # Sign of the scale is irrelevant; only its magnitude matters
    assert noise_atol(-100.0) == NOISE_REL_TOL * 100.0


@pytest.mark.parametrize("data_scale", [None, np.nan, np.inf, -np.inf])
def test_noise_atol_returns_zero_for_unusable_scale(data_scale):
    # No usable scale means "collapse only on exact equality", which
    # errs towards showing the data rather than hiding it.
    assert noise_atol(data_scale) == 0.0


def test_noise_atol_returns_zero_for_zero_scale():
    assert noise_atol(0.0) == 0.0


def test_normalize_colors_difference_preserves_small_magnitude_real_diff():
    # Headline regression test for the reported bug: EmisCO_Aircraft
    # zonal-mean differences are ~1e-14 in absolute terms but ~5% of
    # the field itself, and were being collapsed to a blank panel by
    # the old fixed atol=1e-7.
    norm = normalize_colors(
        -AIRCRAFT_DIFF, AIRCRAFT_DIFF,
        is_difference=True, data_scale=AIRCRAFT_SCALE,
    )
    assert (norm.vmin, norm.vmax) == (-AIRCRAFT_DIFF, AIRCRAFT_DIFF)


def test_normalize_colors_difference_collapses_regridding_noise():
    # Defense in depth for issue #330: the same spread, on a field
    # 12 orders of magnitude larger, is noise and must collapse.
    norm = normalize_colors(
        -NOISY_DIFF, NOISY_DIFF,
        is_difference=True, data_scale=NOISY_HI,
    )
    assert (norm.vmin, norm.vmax) == (-1.0, 1.0)


def test_normalize_colors_difference_collapses_float32_roundoff():
    # New coverage that the old fixed atol=1e-7 did not provide: on a
    # large-magnitude field, plain float32 round-off exceeds 1e-7 and
    # so used to be plotted as if it were signal.
    norm = normalize_colors(
        -FLOAT32_DIFF, FLOAT32_DIFF,
        is_difference=True, data_scale=100.0,
    )
    assert (norm.vmin, norm.vmax) == (-1.0, 1.0)


def test_normalize_colors_difference_without_data_scale_uses_rtol_only():
    # With no scale supplied there is no way to know that a bare
    # +/-1e-11 is noise, so the safe default is to plot it.
    norm = normalize_colors(-1e-11, 1e-11, is_difference=True)
    assert (norm.vmin, norm.vmax) == (-1e-11, 1e-11)


def test_normalize_colors_ref_dev_ignores_data_scale():
    # Ref/Dev panels use the relative tolerance only, no matter what
    # scale is passed (guards the fix from commit 96b0200).
    norm = normalize_colors(
        TRACER_VMIN, TRACER_VMAX, data_scale=1.0
    )
    assert (norm.vmin, norm.vmax) == (TRACER_VMIN, TRACER_VMAX)


def test_normalize_colors_ratio_ignores_data_scale():
    # Ratio panels are dimensionless and anchored at 1; a physical
    # data_scale must not be able to blank them out.
    norm = normalize_colors(
        0.5, 2.0,
        is_difference=True, log_color_scale=True, ratio_log=True,
        is_ratio=True, data_scale=1.0e6,
    )
    assert np.isclose(norm(0.5), 0.0)
    assert np.isclose(norm(2.0), 1.0)


# ----------------------------------------------------------------------
# End-to-end through compute_norm_for_plot, the way six_plot calls it
# ----------------------------------------------------------------------

@pytest.mark.parametrize("subplot", ["dyn_absdiff", "res_absdiff"])
def test_compute_norm_for_plot_zm_absdiff_preserves_small_real_diff(subplot):
    _, norm = compute_norm_for_plot(
        np.array([-AIRCRAFT_DIFF, AIRCRAFT_DIFF]),
        -AIRCRAFT_DIFF,
        AIRCRAFT_DIFF,
        subplot,
        plot_type="zonal_mean",
        data_scale=AIRCRAFT_SCALE,
    )
    assert (norm.vmin, norm.vmax) == (-AIRCRAFT_DIFF, AIRCRAFT_DIFF)


@pytest.mark.parametrize("subplot", ["dyn_absdiff", "res_absdiff"])
def test_compute_norm_for_plot_zm_absdiff_collapses_noise(subplot):
    _, norm = compute_norm_for_plot(
        np.array([-NOISY_DIFF, NOISY_DIFF]),
        -NOISY_DIFF,
        NOISY_DIFF,
        subplot,
        plot_type="zonal_mean",
        data_scale=NOISY_HI,
    )
    assert (norm.vmin, norm.vmax) == (-1.0, 1.0)


def test_compute_norm_for_plot_ref_does_not_receive_data_scale():
    # compute_norm_for_plot must not forward data_scale to the Ref/Dev
    # branch; a huge scale must not blank a small-magnitude panel.
    _, norm = compute_norm_for_plot(
        np.array([TRACER_VMIN, TRACER_VMAX]),
        TRACER_VMIN,
        TRACER_VMAX,
        "ref",
        plot_type="zonal_mean",
        data_scale=1.0e6,
    )
    assert (norm.vmin, norm.vmax) == (TRACER_VMIN, TRACER_VMAX)


def test_compute_norm_for_plot_single_level_absdiff_ignores_data_scale():
    # Single-level plots use exact equality (use_tolerance=False), so
    # data_scale must have no effect at all there.
    _, norm = compute_norm_for_plot(
        np.array([-AIRCRAFT_DIFF, AIRCRAFT_DIFF]),
        -AIRCRAFT_DIFF,
        AIRCRAFT_DIFF,
        "dyn_absdiff",
        plot_type="single_level",
        data_scale=1.0e6,
    )
    assert (norm.vmin, norm.vmax) == (-AIRCRAFT_DIFF, AIRCRAFT_DIFF)


def test_compute_norm_for_plot_diff_of_diffs_row3_preserves_real_signal():
    # In a diff-of-diffs plot, row 3 is named "dyn_absdiff"/"res_absdiff"
    # but holds a difference of *fractional* differences, which is
    # dimensionless.  six_plot passes data_scale=1.0 for that row; a
    # real O(0.1) signal there must survive.
    _, norm = compute_norm_for_plot(
        np.array([-0.1, 0.1]),
        -0.1,
        0.1,
        "dyn_absdiff",
        plot_type="zonal_mean",
        data_scale=1.0,
    )
    assert (norm.vmin, norm.vmax) == (-0.1, 0.1)


# ----------------------------------------------------------------------
# ref_dev_data_scale
# ----------------------------------------------------------------------

@pytest.mark.parametrize(
    "vmins, vmaxs, expected",
    [
        ([1.0, 2.0, 1.0], [3.0, 4.0, 4.0], 4.0),
        ([-5.0, 0.0, -5.0], [0.0, 1.0, 1.0], 5.0),
        ([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 0.0),
        # compare_zonal_mean builds vmin_both/vmax_both with np.min &
        # np.max, not their NaN-safe forms, so an all-NaN Ref poisons
        # index 2.  Scanning all six entries recovers Dev's scale.
        ([np.nan, 1.0, np.nan], [np.nan, 2.0, np.nan], 2.0),
    ],
)
def test_ref_dev_data_scale(vmins, vmaxs, expected):
    assert ref_dev_data_scale(vmins, vmaxs) == expected


def test_ref_dev_data_scale_all_nan_returns_none():
    assert ref_dev_data_scale([np.nan] * 3, [np.nan] * 3) is None


# ----------------------------------------------------------------------
# Coherence between the difference row and the ratio row
# ----------------------------------------------------------------------

@pytest.mark.parametrize("rel_diff", [1e-9, 1e-8, 1e-7, 1e-6, 5e-6, 1e-5,
                                      1.2e-5, 2e-5, 5e-5, 1e-4, 1e-3,
                                      1e-2])
def test_negligible_absdiff_implies_ref_equals_dev(rel_diff):
    """
    Whenever the difference panel is blanked as "negligible", the ratio
    panel must agree and report "Ref and Dev equal throughout domain".

    The reverse need not hold (the ratio panel is the more forgiving of
    the two), but this direction must: it would be contradictory to
    hide the differences while the ratio panel still plots them.
    ref_equals_dev now takes its rtol from NOISE_REL_TOL, so the two
    cannot drift apart; the invariant holds for any value <= 2e-5.
    """
    scale = 100.0
    vmax = rel_diff * scale

    absdiff_is_negligible = is_nearly_constant(
        [-vmax, vmax], atol=noise_atol(scale)
    )
    ratio_says_equal = ref_equals_dev(np.array([1.0 - rel_diff,
                                                1.0 + rel_diff]))

    if absdiff_is_negligible:
        assert ratio_says_equal, (
            f"difference panel blanked at rel_diff={rel_diff:g} while "
            "the ratio panel still plots a difference"
        )


# ----------------------------------------------------------------------
# Sparse fields (zero over most of the domain, e.g. aircraft emissions)
# ----------------------------------------------------------------------

def _absdiff_cbar_label(vmin, vmax, data_scale, subplot="res_absdiff"):
    """Render one absdiff colorbar and return its tick label(s)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from gcpy.plot.core import normalize_colors

    norm = normalize_colors(vmin, vmax, is_difference=True,
                            data_scale=data_scale)
    fig, axes = plt.subplots()
    mesh = axes.pcolormesh(np.zeros((4, 4)), cmap="RdBu_r", norm=norm)
    cbar = fig.colorbar(mesh, ax=axes, orientation="horizontal")
    cbar.mappable.set_norm(norm)
    cbar = colorbar_ticks_and_format(
        np.zeros((4, 4)), cbar, vmin, vmax, subplot,
        data_scale=data_scale, use_tolerance=True,
    )
    fig.canvas.draw()
    labels = [t.get_text() for t in cbar.ax.get_xticklabels()]
    plt.close(fig)
    return labels


def test_restricted_panel_of_sparse_field_is_not_called_negligible():
    # A sparse field (aircraft emissions are zero over most of the
    # domain) makes the 5th and 95th percentiles of Dev - Ref both
    # zero, collapsing the restricted-range panel.  That is NOT the
    # same as the differences being negligible -- the dynamic-range
    # panel beside it still shows real structure -- so it must not
    # claim they are.
    labels = _absdiff_cbar_label(-0.0, 0.0, AIRCRAFT_SCALE)
    assert labels == ["Zero within the 5th-95th percentile range"]


def test_genuinely_negligible_diff_is_still_labeled_negligible():
    labels = _absdiff_cbar_label(-NOISY_DIFF, NOISY_DIFF, NOISY_HI,
                                 subplot="dyn_absdiff")
    assert labels == ["Differences negligible throughout domain"]


def test_real_diff_panel_keeps_its_numeric_colorbar():
    labels = _absdiff_cbar_label(-AIRCRAFT_DIFF, AIRCRAFT_DIFF,
                                 AIRCRAFT_SCALE, subplot="dyn_absdiff")
    assert "Differences negligible throughout domain" not in labels
    assert "Zero within the 5th-95th percentile range" not in labels


# ----------------------------------------------------------------------
# Ratio panels: noise divided by noise, and the collapsed-colorbar anchor
# ----------------------------------------------------------------------

# Magnitude of a field whose real values are ~1e-8, and the residue
# that CS->LL regridding leaves where the field is really zero.
FIELD_SCALE = 1.0e-8
RESIDUE = 1.0e-25


def test_mask_meaningless_ratio_masks_noise_over_noise():
    # Both Ref and Dev are regridding residue: their ratio is
    # arbitrary and used to paint the panel solid red.
    ref = np.array([RESIDUE, 2 * RESIDUE])
    dev = np.array([5 * RESIDUE, RESIDUE])
    out = mask_meaningless_ratio(np.abs(dev) / np.abs(ref),
                                 ref, dev, FIELD_SCALE)
    assert np.isnan(out).all()


def test_mask_meaningless_ratio_keeps_real_data():
    ref = np.array([2.0e-9, 4.0e-9])
    dev = np.array([2.2e-9, 3.6e-9])
    out = mask_meaningless_ratio(np.abs(dev) / np.abs(ref),
                                 ref, dev, FIELD_SCALE)
    assert not np.isnan(out).any()
    assert np.allclose(out, [1.1, 0.9])


def test_mask_meaningless_ratio_keeps_nothing_to_something():
    # Ref negligible but Dev real is a genuine change, not noise, so
    # its large ratio must survive.
    ref = np.array([RESIDUE])
    dev = np.array([2.0e-9])
    out = mask_meaningless_ratio(np.abs(dev) / np.abs(ref),
                                 ref, dev, FIELD_SCALE)
    assert not np.isnan(out).any()


def test_mask_meaningless_ratio_without_scale_is_a_no_op():
    ref = np.array([RESIDUE, 2.0e-9])
    dev = np.array([5 * RESIDUE, 2.2e-9])
    frac = np.abs(dev) / np.abs(ref)
    out = mask_meaningless_ratio(frac, ref, dev, None)
    assert np.array_equal(out, frac)


@pytest.mark.parametrize("subplot, expected", [
    ("dyn_ratio", 1.0),    # MidpointLogNorm spans [0.5, 2.0] about 1.0
    ("res_ratio", 1.0),
    ("dyn_absdiff", 0.0),  # Normalize spans [-1, 1] about 0.0
    ("res_absdiff", 0.0),
])
def test_collapsed_colorbar_tick_matches_the_panel_anchor(subplot, expected):
    # Regression test: a tick at 0.0 on a ratio panel falls outside
    # its [0.5, 2.0] norm, which stretched the colorbar axes and left
    # it rendering completely blank -- no gradient, ticks or label.
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from gcpy.plot.core import normalize_colors

    is_ratio = "ratio" in subplot
    norm = normalize_colors(
        np.nan, np.nan, is_difference=True,
        log_color_scale=is_ratio, ratio_log=is_ratio, is_ratio=is_ratio,
    )
    fig, axes = plt.subplots()
    mesh = axes.pcolormesh(np.full((4, 4), np.nan), cmap="RdBu_r", norm=norm)
    cbar = fig.colorbar(mesh, ax=axes, orientation="horizontal")
    cbar.mappable.set_norm(norm)
    cbar = colorbar_for_all_zero_or_nan(cbar, subplot, all_nan=True)
    fig.canvas.draw()
    ticks = list(cbar.get_ticks())
    xlim = cbar.ax.get_xlim()
    plt.close(fig)

    assert ticks == [expected]
    # The tick must lie inside the norm, so the axes are not stretched
    assert xlim == (norm.vmin, norm.vmax)


def test_res_absdiff_range_ignores_nans():
    # Regression test: np.percentile returns NaN if the array holds a
    # single NaN, which made vmin/vmax NaN.  is_nearly_constant treats
    # an all-non-finite range as constant, so the panel collapsed to a
    # flat -1..1 scale and got labeled "negligible" -- while still
    # drawing the real, now fully saturated, differences underneath.
    data = np.array([[-50.0, 50.0, np.nan], [10.0, -10.0, 20.0]])

    vmin, vmax = vmin_vmax_for_absdiff_plots(data, "res_absdiff", (1, 1))

    assert np.isfinite(vmin) and np.isfinite(vmax)
    assert vmax > 0.0 and vmin == -vmax
    # and it must therefore NOT be mistaken for a constant panel
    norm = normalize_colors(vmin, vmax, is_difference=True, data_scale=50.0)
    assert (norm.vmin, norm.vmax) != (-1.0, 1.0)


# ----------------------------------------------------------------------
# Single-level parity: the same protections, wired to the same places
# ----------------------------------------------------------------------

def _latlon_dataset(values, name="EmisCO_Total"):
    """Wraps a (lat, lon) array as a one-level GEOS-Chem-like Dataset."""
    import xarray as xr
    from gcpy.grid import make_grid_ll

    grid = make_grid_ll("4x5")
    dset = xr.Dataset(
        {name: (("time", "lev", "lat", "lon"), values[None, None, ...])},
        coords={"time": [0], "lev": [1.0],
                "lat": np.asarray(grid["lat"]),
                "lon": np.asarray(grid["lon"])},
    )
    dset[name].attrs.update(units="kg/m2/s", long_name=name)
    return dset


def test_single_level_ratio_masks_noise_over_noise():
    # compare_single_level computes Dev/Ref the same way
    # compare_zonal_mean does, so it needs the same protection: half
    # this map is real emissions and half is regridding residue.
    import matplotlib
    matplotlib.use("Agg")
    import sys
    from gcpy.plot.compare_single_level import compare_single_level

    rng = np.random.default_rng(2)
    nlat, nlon, half = 46, 72, 36
    ref = np.zeros((nlat, nlon))
    dev = np.zeros((nlat, nlon))
    ref[:, :half] = 2.0e-9 * rng.random((nlat, half))
    dev[:, :half] = ref[:, :half] * 1.05
    ref[:, half:] = 1.0e-25 * rng.random((nlat, half))
    dev[:, half:] = 1.0e-25 * rng.random((nlat, half))

    seen = {}
    six = sys.modules["gcpy.plot.six_plot"]
    original = six.compute_norm_for_plot

    def spy(plot_val, vmin, vmax, subplot, **kwargs):
        out, norm = original(plot_val, vmin, vmax, subplot, **kwargs)
        if "ratio" in subplot:
            bad = (np.ma.getmaskarray(out)
                   | np.isnan(np.ma.filled(out, 0.0)))
            seen[subplot] = (bad[:, :half], bad[:, half:])
        return out, norm

    six.compute_norm_for_plot = spy
    try:
        compare_single_level(
            _latlon_dataset(ref), "Ref", _latlon_dataset(dev), "Dev",
            varlist=["EmisCO_Total"],
        )
    finally:
        six.compute_norm_for_plot = original

    assert seen, "no ratio panel was drawn"
    for subplot, (real_half, noise_half) in seen.items():
        assert not real_half.any(), f"{subplot} masked real data"
        assert noise_half.all(), f"{subplot} left noise unmasked"


def test_single_panel_fallback_norm_respects_plot_type():
    # single_panel builds its own norm when none is supplied.  That
    # fallback used to omit use_tolerance, which defaults to True, so
    # a standalone single-level panel took the near-constant tolerance
    # path that six_plot deliberately keeps it out of (GitHub #439).
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import xarray as xr
    from gcpy.grid import make_grid_ll
    from gcpy.plot.single_panel import single_panel

    grid = make_grid_ll("4x5")
    lat, lon = np.asarray(grid["lat"]), np.asarray(grid["lon"])
    # Constant to within the relative tolerance, but far from zero
    const = xr.DataArray(
        np.full((lat.size, lon.size), NOISY_LO),
        dims=("lat", "lon"), coords={"lat": lat, "lon": lon},
    )

    plt.figure()
    try:
        plot = single_panel(const, plot_type="single_level", grid=grid)
        vmin, vmax = plot.norm.vmin, plot.norm.vmax
    finally:
        plt.close("all")

    # Must not collapse onto the 0..1 sentinel: the field is ~100
    assert (vmin, vmax) != (0.0, 1.0)
    assert vmin > 1.0 and vmax > 1.0


def test_ref_equals_dev_survives_the_nan_mask():
    # six_plot passes the ratio panel's data through get_nan_mask
    # before handing it to colorbar_ticks_and_format, so ref_equals_dev
    # sees a masked array.  Reading through that mask exposed the fill
    # value and the "Ref and Dev equal throughout domain" label stopped
    # appearing on any sparse field whose Ref and Dev were identical.
    ratio = np.array([1.0, 1.0, np.nan, 1.0, np.nan])

    assert ref_equals_dev(ratio)
    assert ref_equals_dev(get_nan_mask(ratio))


# ----------------------------------------------------------------------
# Ratio panels with nothing to show: say which side vanished
# ----------------------------------------------------------------------

@pytest.mark.parametrize(
    "subplot, all_zero, all_nan, ref_is_zero, dev_is_zero, expected",
    [
        # Ref zero, Dev real: Dev/Ref is a division by zero everywhere
        ("dyn_ratio", False, True, True, False,
         "Ref is zero throughout domain"),
        ("res_ratio", False, True, True, False,
         "Ref is zero throughout domain"),
        # Dev zero, Ref real: the ratio is identically zero.  Note this
        # case arrives via all_zero, not all_nan.
        ("dyn_ratio", True, False, False, True,
         "Dev is zero throughout domain"),
        ("res_ratio", True, False, False, True,
         "Dev is zero throughout domain"),
        # Both zero is 0/0, which really is undefined -- left alone
        ("dyn_ratio", False, True, True, True,
         "Undefined throughout domain"),
        # Neither zero, but the data is all NaN for some other reason
        ("dyn_ratio", False, True, False, False,
         "Undefined throughout domain"),
        # Non-ratio panels must be untouched by the new wording
        ("dyn_absdiff", True, False, True, False, "Zero throughout domain"),
        ("ref", True, False, True, False, "Zero throughout domain"),
        ("dev", False, True, True, False, "Undefined throughout domain"),
    ],
)
def test_ratio_colorbar_names_the_side_that_is_zero(
        subplot, all_zero, all_nan, ref_is_zero, dev_is_zero, expected
):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots()
    mesh = axes.pcolormesh(np.full((4, 4), np.nan), cmap="RdBu_r")
    cbar = fig.colorbar(mesh, ax=axes, orientation="horizontal")
    cbar = colorbar_for_all_zero_or_nan(
        cbar, subplot, all_nan=all_nan,
        ref_is_zero=ref_is_zero, dev_is_zero=dev_is_zero,
    )
    fig.canvas.draw()
    labels = [t.get_text() for t in cbar.ax.get_xticklabels()]
    plt.close(fig)

    assert labels == [expected]


# ======================================================================
# Ref/Dev "is this field constant" tolerance, scaled to the data's
# precision.  Regression tests for the PassiveTracer restart (100 ppb)
# whose zonal-mean panels collapsed onto a blank 0-1 ppb colorbar.
# ======================================================================

# A PassiveTracer restart initialized to 1e-7 v/v, i.e. 100 ppb, with
# the real structure it carries: a relative spread of 7.5e-6, which a
# 1e-5 tolerance swallowed.
TRACER_RST_VMIN = 99.99925
TRACER_RST_VMAX = 100.00000

# The issue #330 CS->LL regridding noise, in relative terms
REGRID_NOISE_REL = (NOISY_HI - NOISY_LO) / NOISY_HI


def test_constant_rel_tol_defaults_to_real4_when_precision_unknown():
    assert constant_rel_tol() == CONSTANT_REL_TOL
    assert constant_rel_tol(None) == CONSTANT_REL_TOL
    assert constant_rel_tol([1.0, 2.0]) == CONSTANT_REL_TOL


def test_constant_rel_tol_is_looser_for_real4_than_real8():
    tol32 = constant_rel_tol(np.zeros(4, dtype=np.float32))
    tol64 = constant_rel_tol(np.zeros(4, dtype=np.float64))
    assert tol32 > tol64, "real*4 round-off needs the looser tolerance"


def test_constant_rel_tol_real4_straddles_roundoff_and_real_signal():
    # Must sit above real*4 round-off but below the real structure in
    # a PassiveTracer restart, or one of the two gets misclassified.
    tol = constant_rel_tol(np.zeros(4, dtype=np.float32))
    eps32 = float(np.finfo(np.float32).eps)
    signal = (TRACER_RST_VMAX - TRACER_RST_VMIN) / TRACER_RST_VMAX
    assert eps32 < tol < signal


def test_constant_rel_tol_floored_at_regridding_noise_for_real8():
    # Interpolation error is a property of the regridding scheme, not
    # of the dtype, so real*8 must not drop below the floor and let
    # issue #330's striping back in.
    tol = constant_rel_tol(np.zeros(4, dtype=np.float64))
    assert tol >= REGRID_NOISE_REL


def test_normalize_colors_preserves_real_tracer_restart_range():
    # The reported bug: these panels must keep their true 100 ppb
    # range, not collapse to a dimensionless 0-1 scale.
    norm = normalize_colors(
        TRACER_RST_VMIN,
        TRACER_RST_VMAX,
        rel_tol=constant_rel_tol(np.zeros(4, dtype=np.float32)),
    )
    assert (norm.vmin, norm.vmax) == (TRACER_RST_VMIN, TRACER_RST_VMAX)


def test_normalize_colors_still_collapses_regridding_noise():
    # ...while the issue #330 noise it was built for still collapses.
    norm = normalize_colors(
        NOISY_LO,
        NOISY_HI,
        rel_tol=constant_rel_tol(np.zeros(4, dtype=np.float32)),
    )
    assert (norm.vmin, norm.vmax) == (0.0, 1.0)


def test_normalize_colors_rel_tol_is_scale_invariant():
    # atol is 0 on this path, so the same relative spread must be
    # judged identically whether the field is ppb (~1e2) or an
    # aircraft emission flux (~1e-13).
    rel = 7.5e-6
    tol = constant_rel_tol(np.zeros(4, dtype=np.float32))
    for mag in (1.0e2, 1.0e0, 1.0e-6, 1.0e-13):
        norm = normalize_colors(mag * (1.0 - rel), mag, rel_tol=tol)
        assert (norm.vmin, norm.vmax) != (0.0, 1.0), f"collapsed at {mag:g}"


def test_compute_norm_for_plot_zonal_mean_keeps_tracer_restart_range():
    # End-to-end through the call the zonal-mean Ref panel actually
    # makes, with the field supplied so its precision is used.
    plot_val = np.array(
        [TRACER_RST_VMIN, TRACER_RST_VMAX],
        dtype=np.float32,
    )
    _, norm = compute_norm_for_plot(
        plot_val,
        TRACER_RST_VMIN,
        TRACER_RST_VMAX,
        "ref",
        plot_type="zonal_mean",
    )
    assert (norm.vmin, norm.vmax) == (TRACER_RST_VMIN, TRACER_RST_VMAX)


def test_colorbar_for_constant_field_names_the_value():
    # A collapsed panel must say what the constant is.  Leaving the
    # numeric ticks of the dimensionless norm on a colorbar labeled in
    # the field's units is what made a 100 ppb field read as 1 ppb.
    _, axes = plt.subplots()
    mappable = axes.imshow(np.full((2, 2), 100.0))
    cbar = plt.colorbar(mappable, ax=axes)
    cbar = colorbar_for_constant_field(cbar, 100.0, 100.0)
    labels = [tick.get_text() for tick in cbar.ax.get_yticklabels()]
    plt.close("all")
    assert labels == ["Constant at 100 throughout domain"]


def test_colorbar_ticks_and_format_labels_constant_ref_panel():
    # The missing branch: a Ref panel that normalize_colors collapsed
    # had no label path and fell through to generic numeric ticks.
    _, axes = plt.subplots()
    mappable = axes.imshow(np.full((2, 2), 0.5))
    cbar = plt.colorbar(mappable, ax=axes)
    cbar = colorbar_ticks_and_format(
        np.full((2, 2), NOISY_LO, dtype=np.float32),
        cbar,
        NOISY_LO,
        NOISY_HI,
        "ref",
        use_tolerance=True,
    )
    labels = [tick.get_text() for tick in cbar.ax.get_yticklabels()]
    plt.close("all")
    assert len(labels) == 1
    assert labels[0].startswith("Constant at ")


# ======================================================================
# Agreement between the difference row and the ratio row.  Both are
# governed by NOISE_REL_TOL, but the difference panels are symmetric
# about zero, so testing their full span measured 2 * max|Dev - Ref|
# and applied an effective tolerance of NOISE_REL_TOL / 2.  A real
# 7.5e-6 relative difference (a TransportTracers PassiveTracer restart,
# 14.7.0 vs 14.8.0) fell in the resulting factor-of-two gap: the
# difference row plotted the offset while the ratio row beside it
# announced "Ref and Dev equal throughout domain".
# ======================================================================

# Magnitude of a PassiveTracer restart, in ppb
DIFF_SCALE = 100.0

# Relative Dev - Ref difference between the 14.7.0 and 14.8.0 restarts
TRACER_REL_DIFF = 7.5e-6


def _rows_agree(rel_diff, data_scale=DIFF_SCALE):
    """Returns (difference row says negligible, ratio row says equal)."""
    half = rel_diff * data_scale
    norm = normalize_colors(
        -half, half, is_difference=True, data_scale=data_scale
    )
    diff_negligible = (norm.vmin, norm.vmax) == (-1.0, 1.0)
    ratio_equal = bool(ref_equals_dev(np.full(8, 1.0 + rel_diff)))
    return diff_negligible, ratio_equal


def test_diff_is_negligible_tests_half_the_span():
    # The criterion is on max|Dev - Ref|, which is half of a symmetric
    # panel's span -- not on the span itself.
    scale = DIFF_SCALE
    atol = noise_atol(scale)
    assert diff_is_negligible(-atol, atol, scale)
    assert not diff_is_negligible(-2.0 * atol, 2.0 * atol, scale)


def test_diff_is_negligible_requires_flat_panel_without_data_scale():
    # No usable data_scale means no tolerance to apply (cf. noise_atol)
    assert diff_is_negligible(0.0, 0.0, None)
    assert not diff_is_negligible(-1.0, 1.0, None)


def test_diff_and_ratio_rows_agree_on_the_tracer_restart():
    # The reported contradiction: these two must not disagree.
    diff_negligible, ratio_equal = _rows_agree(TRACER_REL_DIFF)
    assert diff_negligible == ratio_equal, (
        "difference and ratio rows disagree about whether Ref and Dev "
        "differ"
    )
    # ...and both should report the difference, since it is real
    assert not diff_negligible
    assert not ratio_equal


def test_diff_and_ratio_rows_agree_across_magnitudes():
    # Sweep the band the two rows used to straddle, plus decades to
    # either side, at magnitudes spanning ppb to aircraft-emission flux.
    for scale in (1.0e2, 1.0e0, 1.0e-13):
        for rel_diff in np.logspace(-8, -3, 60):
            diff_negligible, ratio_equal = _rows_agree(rel_diff, scale)
            assert diff_negligible == ratio_equal, (
                f"rows disagree at scale={scale:g}, "
                f"rel_diff={rel_diff:.3e}"
            )


def test_noise_rel_tol_is_the_threshold_both_rows_use():
    # Straddle NOISE_REL_TOL itself and confirm both rows switch there.
    below = _rows_agree(NOISE_REL_TOL * 0.5)
    above = _rows_agree(NOISE_REL_TOL * 2.0)
    assert below == (True, True), "both rows should call this noise"
    assert above == (False, False), "both rows should call this signal"


def test_ratio_panel_collapse_uses_the_shared_constant():
    # The ratio panel's own color-scale collapse must move with the
    # label from ref_equals_dev, rather than tracking
    # is_nearly_constant's separate default.
    lo, hi = 1.0, 1.0 + NOISE_REL_TOL * 0.5
    norm = normalize_colors(lo, hi, is_difference=True, is_ratio=True)
    assert (norm.vmin, norm.vmax) == (0.5, 2.0)

    lo, hi = 1.0, 1.0 + NOISE_REL_TOL * 4.0
    norm = normalize_colors(lo, hi, is_difference=True, is_ratio=True)
    assert (norm.vmin, norm.vmax) != (0.5, 2.0)
