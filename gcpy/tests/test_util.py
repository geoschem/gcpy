"""
Unit tests for helper routines in gcpy/util.py.
"""
import warnings

import numpy as np
import pytest
import xarray as xr

from gcpy.util import get_nan_mask, is_nearly_constant, \
    warn_if_flip_levels_mismatch, normalize_mapl3_vertical_dims


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


# ----------------------------------------------------------------------
# MAPL3 vertical dimensions ("lev" for level centers, "edge" for edges)
# ----------------------------------------------------------------------

def _make_mapl3_dataset():
    """MAPL3-style cubed-sphere dataset mixing lev (72) and edge (73)."""
    dims_lev = ("time", "lev", "nf", "Ydim", "Xdim")
    dims_edge = ("time", "edge", "nf", "Ydim", "Xdim")
    shape_lev = (1, 72, 6, 2, 2)
    shape_edge = (1, 73, 6, 2, 2)
    return xr.Dataset(
        {
            "MFX": (dims_lev,
                    np.arange(np.prod(shape_lev), dtype=float)
                    .reshape(shape_lev)),
            "PLE0": (dims_edge,
                     np.arange(np.prod(shape_edge), dtype=float)
                     .reshape(shape_edge)),
            "AREA": (("nf", "Ydim", "Xdim"), np.ones((6, 2, 2))),
        },
        coords={
            "lev": np.arange(1, 73, dtype=float),
            "edge": np.arange(1, 74, dtype=float),
        },
    )


def test_normalize_mapl3_leaves_datasets_without_edge_alone():
    dset = _make_mapl3_dataset().drop_dims("edge")

    assert normalize_mapl3_vertical_dims(dset, ["MFX"]) is dset


def test_normalize_mapl3_edge_variables_become_lev():
    dset = _make_mapl3_dataset()
    out = normalize_mapl3_vertical_dims(dset, ["PLE0"])

    assert "edge" not in out.dims
    assert out.sizes["lev"] == 73
    assert "MFX" not in out.data_vars
    assert "AREA" in out.data_vars
    assert out["PLE0"].dims[1] == "lev"
    # Level order must be exactly as in the file (no flipping)
    assert np.array_equal(out["PLE0"].values, dset["PLE0"].values)
    assert np.array_equal(out["lev"].values, dset["edge"].values)


def test_normalize_mapl3_lev_variables_drop_edge():
    dset = _make_mapl3_dataset()
    out = normalize_mapl3_vertical_dims(dset, ["MFX"])

    assert "edge" not in out.dims
    assert out.sizes["lev"] == 72
    assert "PLE0" not in out.data_vars
    assert "AREA" in out.data_vars
    assert np.array_equal(out["MFX"].values, dset["MFX"].values)


def test_normalize_mapl3_rejects_mixed_lev_and_edge():
    with pytest.raises(ValueError, match="separate calls"):
        normalize_mapl3_vertical_dims(_make_mapl3_dataset(), ["MFX", "PLE0"])


def test_normalize_mapl3_without_varlist_skips_edge_with_warning():
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        out = normalize_mapl3_vertical_dims(_make_mapl3_dataset())

    assert "edge" not in out.dims
    assert "PLE0" not in out.data_vars
    assert "MFX" in out.data_vars
    assert len(caught) == 1
    assert "PLE0" in str(caught[0].message)
