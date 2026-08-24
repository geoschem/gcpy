"""
Unit tests covering the plotting-side "near-constant" tolerance guard
added for GitHub issue #330 (zonal mean plots of constant fields
showed color striping because a colorbar range computed from
noisy-but-nominally-constant data was not collapsed to a flat scale).

Also covers two follow-on fixes for the same over-eager tolerance
collapsing real, non-constant, small-magnitude fields to a blank/
white panel: PR #439 (single-level plots, e.g. HEMCO emissions) and
zonal-mean Ref/Dev panels (e.g. TransportTracers species like Rn222).
"""
import numpy as np
from matplotlib import colors

from gcpy.util import is_nearly_constant
from gcpy.plot.core import normalize_colors
from gcpy.plot.six_plot import ref_equals_dev, compute_norm_for_plot

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
    norm = normalize_colors(-1e-11, 1e-11, is_difference=True)
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
