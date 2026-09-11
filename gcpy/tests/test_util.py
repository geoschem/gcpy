"""
Unit tests for helper routines in gcpy/util.py.
"""
import warnings

import numpy as np
import pytest

from gcpy.util import get_nan_mask, is_nearly_constant, \
    warn_if_flip_levels_mismatch


def test_get_nan_mask_masks_the_nans():
    # Regression test: the mask used to be built from the original
    # array rather than from the filled one.  The fill value is chosen
    # to lie above every element of the input, so the test never
    # matched and the routine masked nothing at all.
    data = np.array([1.0, 2.0, np.nan, 4.0])
    out = get_nan_mask(data)

    assert np.array_equal(np.ma.getmaskarray(out),
                          [False, False, True, False])
    assert not np.isnan(np.ma.filled(out, 0.0)).any()
    assert np.array_equal(out.compressed(), [1.0, 2.0, 4.0])


def test_get_nan_mask_leaves_finite_data_alone():
    data = np.array([1.0, 2.0, 3.0])
    out = get_nan_mask(data)

    assert not np.ma.getmaskarray(out).any()
    assert np.array_equal(np.ma.filled(out, 0.0), data)


def test_get_nan_mask_handles_all_nan():
    out = get_nan_mask(np.array([np.nan, np.nan]))

    assert np.ma.getmaskarray(out).all()


def test_is_nearly_constant_respects_a_mask():
    # A masked array's underlying data still holds whatever sits under
    # the mask.  get_nan_mask deliberately fills with a value far above
    # the data, so reading through the mask would make a constant array
    # look wildly variable.
    masked = get_nan_mask(np.array([1.0, 1.0, np.nan, 1.0]))

    assert is_nearly_constant(masked)
    assert is_nearly_constant(masked, target=1.0)


def test_is_nearly_constant_still_detects_real_variation_when_masked():
    masked = get_nan_mask(np.array([0.5, 1.0, np.nan, 2.0]))

    assert not is_nearly_constant(masked)
    assert not is_nearly_constant(masked, target=1.0)


# ----------------------------------------------------------------------
# flip_levels guard
# ----------------------------------------------------------------------

@pytest.mark.parametrize("flip_ref, flip_dev", [(False, False), (True, True)])
def test_no_warning_when_flip_levels_agree(flip_ref, flip_dev):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        warned = warn_if_flip_levels_mismatch(flip_ref, flip_dev)

    assert warned is False
    assert caught == []


@pytest.mark.parametrize("flip_ref, flip_dev, flipped, other", [
    (True, False, "Ref", "Dev"),
    (False, True, "Dev", "Ref"),
])
def test_warns_when_only_one_side_is_flipped(flip_ref, flip_dev, flipped, other):
    # Comparing opposite ends of the vertical grid renders one panel
    # entirely zero for a surface-only field, which reads as a real
    # change rather than as a configuration mistake.
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        warned = warn_if_flip_levels_mismatch(flip_ref, flip_dev)

    assert warned is True
    assert len(caught) == 1
    assert issubclass(caught[0].category, UserWarning)
    message = str(caught[0].message)
    assert f"levels of {flipped} will be flipped" in message
    assert f"of {other} will not" in message


def test_flip_levels_guard_accepts_non_bool_truthiness():
    # compare_diags reads these straight out of YAML, so they can
    # arrive as anything truthy rather than as real bools.
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        warned = warn_if_flip_levels_mismatch(1, 0)

    assert warned is True
    assert len(caught) == 1
