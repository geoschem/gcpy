"""
Common variables and functions used by modules in gcpy.plot.
"""
from os import path
import warnings
from matplotlib import colors
import numpy as np
from gcpy.util import is_nearly_constant

# ======================================================================
# Constants used for plotting
# ======================================================================

# Save warnings format to undo overwriting built into pypdf
_warning_format = warnings.showwarning

# Current directory
_plot_dir = path.dirname(__file__)

# Colormap definitions
_rgb_WhGrYlRd = np.genfromtxt(
    path.join(_plot_dir, 'colormaps', 'WhGrYlRd.txt'),
    delimiter=' '
)
WhGrYlRd = colors.ListedColormap(_rgb_WhGrYlRd / 255.0)

# Use a style sheet to control plot attributes
gcpy_style = path.join(_plot_dir, "gcpy_plot_style")

# Relative tolerance for deciding when a Ref vs. Dev difference is
# numerical noise rather than real signal.  For difference panels it
# is scaled by the magnitude of the Ref & Dev data (see function
# gcpy.plot.core.noise_atol), so it does not depend on the units of
# the field.  The same value serves as the relative tolerance in
# function gcpy.plot.six_plot.ref_equals_dev, so that the difference
# and ratio rows will agree about when Ref and Dev differ.
#
# The criterion is on the largest |Dev - Ref| in the panel, i.e. the
# rows agree when max|Dev - Ref| <= NOISE_REL_TOL * data_scale.
NOISE_REL_TOL = 5.0e-6

# Absolute tolerance for comparing dimensionless ratio data against
# 1.0 (see function gcpy.plot.six_plot.ref_equals_dev).
RATIO_ABS_TOL = 1.0e-10

# Number of ULPs (units in the last place) of the data's own
# precision that a Ref or Dev field may span and still count as
# constant across the domain.  See gcpy.plot.six_plot.constant_rel_tol.
CONSTANT_TOL_ULPS = 8

# Relative tolerance for numerical noise caused by regridding
REGRID_NOISE_REL_TOL = 1.0e-8

# Tolerance that function gcpy.plot.core.constant_rel_tol
# returns for float32 (aka REAL*4) data, which is the precision
# that most of GEOS-Chem's output fields have.
CONSTANT_REL_TOL = max(
    CONSTANT_TOL_ULPS * float(np.finfo(np.float32).eps),
    REGRID_NOISE_REL_TOL,
)


def constant_rel_tol(data=None):
    """
    Returns the relative tolerance for deciding that a Ref or Dev
    field is constant across the domain, and so should get a flat
    color scale instead of one stretched across numerical noise
    (see GitHub issue #330).  The tolerance is scaled to the precision
    the data is carried at and floored at the regridding noise level;
    see CONSTANT_TOL_ULPS and REGRID_NOISE_REL_TOL.

    Parameters
    ----------
    data : xarray.DataArray or numpy.ndarray, optional
        The field being plotted.  If it is missing, or its precision
        cannot be determined, CONSTANT_REL_TOL value is returned,
        since most GEOS-Chem output fields have float32 (aka REAL*4)
        precision.
        Default value: None

    Returns
    -------
    rtol : float
        Relative tolerance value to pass to the function
        :func:`gcpy.util.is_nearly_constant`.

    Notes
    -----
    The on-disk dtype (precision) is preferred over the in-memory one.
    An array read from a float32 (aka REAL*4) file and promoted to
    float64 by regridding still carries only float32 information.
    Judging it by float64 round-off would mistake that file's round-off
    for real signal.
    """
    dtype = None

    # Prefer the on-disk dtype recorded by xarray at read time
    encoding = getattr(data, "encoding", None)
    if isinstance(encoding, dict):
        dtype = encoding.get("dtype", None)

    # Otherwise fall back to the in-memory dtype
    if dtype is None:
        dtype = getattr(data, "dtype", None)

    # Tested explicitly: np.finfo(None) quietly reports float64
    # rather than raising, which would silently apply real*8
    # precision to data whose dtype we never actually learned.
    if dtype is None:
        return CONSTANT_REL_TOL

    try:
        eps = float(np.finfo(dtype).eps)
    except (TypeError, ValueError):
        # Not a float dtype (or no dtype at all): assume real*4
        return CONSTANT_REL_TOL

    return max(CONSTANT_TOL_ULPS * eps, REGRID_NOISE_REL_TOL)


def six_panel_subplot_names(diff_of_diffs):
    """
    Returns the names of the subplots for the 6-panel plots.

    Parameters
    ----------
    diff_of_diffs : bool
        Indicates if this is a diff-of-diffs benchmark (True)
        or not (False).  Ratio plots are only included if
        diff_of_diffs is False.

    Returns
    -------
    subplots : list of str
        List of names of each of the subplots in the 6-panel plot.
    """
    if diff_of_diffs:
        return ["ref", "dev",
                "dyn_absdiff", "res_absdiff",
                "dyn_absdiff", "res_absdiff"]

    return ["ref", "dev",
            "dyn_absdiff", "res_absdiff",
            "dyn_ratio", "res_ratio",
    ]


def noise_atol(data_scale):
    """
    Returns the absolute tolerance for deciding whether a difference
    panel holds real signal or only numerical noise.  Scaling the
    tolerance by the magnitude of the data lets one criterion serve
    both a field of magnitude ~100 and one of magnitude ~1e-13.

    Parameters
    ----------
    data_scale : float or None
        Magnitude of the Ref and Dev data (i.e. the largest absolute
        value in either field).

    Returns
    -------
    atol : float
        Absolute tolerance to pass to :func:`gcpy.util.is_nearly_constant`.
        A missing or non-finite data_scale yields 0.0, so that the panel
        collapses to a flat color scale only on exact equality.
    """
    if data_scale is None:
        return 0.0

    scale = abs(float(data_scale))
    if not np.isfinite(scale):
        return 0.0

    return NOISE_REL_TOL * scale


def diff_is_negligible(vmin, vmax, data_scale):
    """
    Returns whether a difference panel holds only numerical noise
    rather than real signal, so that its color scale should be
    collapsed to a flat color scale and its colorbar labeled as such.

    Parameters
    ----------
    vmin, vmax : float
        Min and max of the difference panel's data range.  These are
        symmetric about zero (see six_plot.vmin_vmax_for_absdiff_plots).
    data_scale : float or None
        Magnitude of the Ref and Dev data (see
        :func:`gcpy.plot.six_plot.ref_dev_data_scale`).

    Returns
    -------
    is_negligible : bool
        Whether the panel's differences are negligible compared to the
        magnitude of the data that was differenced.

    Notes
    -----
    The test is on half the span, not the whole of it.  Because the
    range is symmetric about zero, half the span is the largest
    |Dev - Ref| in the panel.
    """
    atol = noise_atol(data_scale)

    # Without a usable data_scale there is no tolerance to apply, so
    # collapse only a panel that is exactly flat (cf. noise_atol).
    if atol <= 0.0:
        return float(vmin) == float(vmax)

    return abs(0.5 * (float(vmax) - float(vmin))) <= atol


def mask_meaningless_ratio(fracdiff, ref, dev, data_scale):
    """
    Masks the cells of a Dev/Ref ratio array in which both Ref and Dev
    are negligible compared to the magnitude of the field.  This avoids
    presenting numerical noise as if it were signal.

    Parameters
    ----------
    fracdiff : numpy array
        The Dev/Ref ratio array.
    ref, dev : numpy array
        The Ref and Dev data that were divided, in physical units.
    data_scale : float or None
        Magnitude of the Ref and Dev data (see
        six_plot.ref_dev_data_scale).  If None, nothing is masked.

    Returns
    -------
    fracdiff : numpy array
        The ratio array, with the noise-only cells set to NaN.

    Notes
    -----
    Only cells where *both* Ref and Dev are negligible are masked.
    A cell where Ref is negligible but Dev is not represents a real
    change from nothing to something, and its large ratio is meaningful,
    so it is left alone.
    """
    atol = noise_atol(data_scale)
    if atol <= 0.0:
        return fracdiff

    is_noise = (np.abs(ref) <= atol) & (np.abs(dev) <= atol)

    return np.where(is_noise, np.nan, fracdiff)


def normalize_colors(
        vmin,
        vmax,
        is_difference=False,
        log_color_scale=False,
        ratio_log=False,
        is_ratio=False,
        use_tolerance=True,
        data_scale=None,
        rel_tol=None,
):
    """
    Normalizes a data range to the colormap range used by matplotlib
    functions. For log-color scales, special handling is done to prevent
    taking the log of data that is all zeroes.

    Parameters
    ----------
    vmin : float
        Minimum value of the data range.
    vmax : float
        Maximum value of the data range.
    is_difference : bool, optional
        Set this switch to denote that we are using a difference
        color scale (i.e. with zero in the middle of the range).
        Default value: False
    log_color_scale : bool, optional
        Logical flag to denote that we are using a logarithmic
        color scale instead of a linear color scale.
        Default value: False
    ratio_log : bool, optional
        Indicates whether we are using log scaling for ratio plots
        (True) or not (False).
        Default value: False
    is_ratio : bool, optional
        Set this switch to denote that we are using a ratio color
        scale (i.e. with 1, not 0, in the middle of the range).
        Default value: False
    use_tolerance : bool, optional
        Logical flag that determines when a difference or ratio panel
        will be "collapsed" to a flat color scale.  Set to True to
        collapse when vmin/vmax is "nearly" constant.  Set to False to
        collapse only upon an exact match (vmin==vmax, vmin==vmax==0,
        or vmin==vmax==NaN).
        Default value: True
    data_scale : float, optional
        Magnitude of the Ref and Dev data that were differenced to
        produce this panel.  Used only for difference panels, to scale
        the tolerance that decides whether the panel holds real signal
        or only numerical noise.  If None, such a panel collapses only
        when vmin and vmax are exactly equal.
        Default value: None
    rel_tol : float, optional
        Relative tolerance for deciding that a Ref or Dev panel is
        constant across the domain.  Callers should pass
        constant_rel_tol(data) so that the tolerance matches the
        precision the field is carried at (float32 vs float64).  If
        None, the float32 value CONSTANT_REL_TOL is used.
        Default value: None

    Returns
    -------
    norm : matplotlib.colors.Normalize
        The normalized matplotlib color range, stored in
        a matplotlib Normalize object.

    Notes
    -----
    For log color scales, we will use a range of 3 orders of
    magnitude (i.e. from vmax/1e3 to vmax).
    """

    # Define class for logarithmic non-symmetric color scheme
    class MidpointLogNorm(colors.LogNorm):
        """
        Class for logarithmic non-symmetric color scheme
        """
        def __init__(
                self,
                vmin=None,
                vmax=None,
                midpoint=None,
                clip=False
        ):
            super().__init__(vmin, vmax, clip)
            self.midpoint = midpoint

        def __call__(self, value, clip=None):
            result, _ = self.process_value(value)
            x_val = [
                np.log(self.vmin),
                np.log(self.midpoint),
                np.log(self.vmax)
            ]
            y_val = [0, 0.5, 1]
            return np.ma.array(
                np.interp(np.log(value), x_val, y_val),
                mask=result.mask,
                copy=False
            )

    # If use_tolerance is True, then allow a data field to be called
    # "constant" if it is either exactly equal everywhere or if it
    # falls within a certain tolerance of being equal everywhere.
    # This will let us know when to collapse to a flat color scale,
    # avoiding spurious striping artifacts as in GitHub issue #330.
    if use_tolerance:

        # Ratio plots: Anchored at 1.0
        if is_ratio:
            is_constant = is_nearly_constant(
                [vmin, vmax],
                rtol=NOISE_REL_TOL,
                atol=RATIO_ABS_TOL,
                target=1.0,
            )

        # Difference plots: Symmetric around zero
        elif is_difference:
            is_constant = diff_is_negligible(vmin, vmax, data_scale)

        # Ref and Dev have no natural anchor, so only use the
        # relative tolerance in making this determination.
        else:
            if rel_tol is None:
                rel_tol = CONSTANT_REL_TOL
            is_constant = is_nearly_constant(
                [vmin, vmax],
                rtol=rel_tol,
                atol=0.0,
            )

    # If use_tolerance is False, then the data field can only be
    # called "constant" if it is either zero or undefined (NaN)
    # everywhere.
    else:
        is_constant = (
            (vmin == 0 and vmax == 0)
            or (np.isnan(vmin) and np.isnan(vmax))
        )

    # For the case where the data is zero everywhere (vmin==vmax==0)
    # or undefined everywhere (vmin==vmax==NaN), normalize the data range
    # so that the color corresponding to zero (white) will be placed in
    # the middle of the colorbar, where we will add a single tick.
    if is_constant:

        # Ratio plots: Place "white" at 1.0
        if is_ratio:
            return MidpointLogNorm(vmin=0.5, vmax=2.0, midpoint=1.0)

        # Difference plots: Place "white" at 0.0
        if is_difference:
            return colors.Normalize(vmin=-1.0, vmax=1.0)

        # Ref and Dev plots: Place "white" at 0.5
        return colors.Normalize(vmin=0.0, vmax=1.0)

    # For log color scales, assume a range 3 orders of magnitude
    # below the maximum value.  Otherwise use a linear scale.
    if log_color_scale and not ratio_log:
        return colors.LogNorm(vmin=vmax/1e3, vmax=vmax)
    if log_color_scale:
        return MidpointLogNorm(vmin=vmin, vmax=vmax, midpoint=1)

    # For linear color scales: normalize to the data range
    return colors.Normalize(vmin=vmin, vmax=vmax)


def text_to_data_units(ax, text):
    """
    Computes the width of a label in data units.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Matplotlib Axes object.
    text : matplotlib.text.Text
        Text object that is being plotted.

    Returns
    -------
    width_in_data_units : float
        Length of the text label in data units.
    """

    # Get the extent of the text as a Bbox object
    bbox = text.get_window_extent()

    # Convert Bbox width from pixels to data units
    inv = ax.transData.inverted()
    pixel_to_data = inv.transform([[bbox.x0, bbox.y0], [bbox.x1, bbox.y0]])
    width_in_data_units = abs(pixel_to_data[1][0] - pixel_to_data[0][0])

    return width_in_data_units
