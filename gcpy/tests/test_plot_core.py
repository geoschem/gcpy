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

from gcpy.util import is_nearly_constant
from gcpy.plot.core import NOISE_REL_TOL, noise_atol, normalize_colors
from gcpy.plot.six_plot import (
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
