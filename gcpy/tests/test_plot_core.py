"""
Unit tests covering the plotting-side "near-constant" tolerance guard
added for GitHub issue #330 (zonal mean plots of constant fields
showed color striping because a colorbar range computed from
noisy-but-nominally-constant data was not collapsed to a flat scale).
"""
import numpy as np
from matplotlib import colors

from gcpy.util import is_nearly_constant
from gcpy.plot.core import normalize_colors
from gcpy.plot.six_plot import ref_equals_dev

# Magnitude of the regridding noise reported in GitHub issue #330
# (~5e-9 relative on a field with a value of ~100).
NOISY_LO = 100.00064076
NOISY_HI = 100.00064126


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
