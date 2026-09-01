"""
Unit tests for helper routines in gcpy/util.py.
"""
import numpy as np

from gcpy.util import get_nan_mask


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
