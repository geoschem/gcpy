"""
Common variables and functions used by modules in gcpy.plot.
"""
import logging
from os import path
import warnings
from matplotlib import colors
import numpy as np

from gcpy.util import is_nearly_constant

# Save warnings format to undo overwriting built into pypdf
_warning_format = warnings.showwarning

# Silence benign "findfont: Failed to find font weight ..." messages.
# These come from Matplotlib's font_manager logger (not the warnings
# module), and only indicate that the installed font family lacks an
# exact-weight face for a style setting below (e.g. "medium");
# Matplotlib already falls back to the closest available weight.
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)

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

# Relative tolerance for deciding that a Ref vs. Dev difference is
# numerical noise rather than real signal.  For difference panels it
# is scaled by the magnitude of the Ref & Dev data (see noise_atol),
# so it does not depend on the units of the field.  The same value
# serves as the rtol of six_plot.ref_equals_dev, so that the
# difference and ratio rows agree about when Ref and Dev differ.
NOISE_REL_TOL = 1.0e-5

# Absolute tolerance for comparing dimensionless ratio data against
# 1.0 (see six_plot.ref_equals_dev).  Ratios are anchored at 1, so
# NOISE_REL_TOL dominates; this is only a floor for values near zero.
RATIO_ABS_TOL = 1.0e-10


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
        Absolute tolerance to pass to is_nearly_constant.  A missing
        or non-finite data_scale yields 0.0, so that the panel
        collapses only on exact equality.
    """
    if data_scale is None:
        return 0.0

    scale = abs(float(data_scale))
    if not np.isfinite(scale):
        return 0.0

    return NOISE_REL_TOL * scale


def mask_meaningless_ratio(fracdiff, ref, dev, data_scale):
    """
    Masks the cells of a Dev/Ref ratio array in which both Ref and Dev
    are negligible compared to the magnitude of the field, so that the
    ratio panel does not present numerical noise as though it were
    signal.

    Where a field is really zero, regridding (and float32 round-off)
    leaves behind a tiny residue instead of an exact zero.  Dividing
    one such residue by another gives a ratio that is arbitrary and
    routinely saturates the color scale, painting whole regions solid
    red next to regions of exact zeros left gray.

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
    Only cells where *both* Ref and Dev are negligible are masked.  A
    cell where Ref is negligible but Dev is not represents a real
    change from nothing to something, and its large ratio is
    meaningful, so it is left alone.
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
        Whether to treat vmin/vmax as "nearly constant" using
        is_nearly_constant()'s relative/absolute tolerance (True),
        or to require an exact match (vmin == vmax, or both zero, or
        both NaN) before collapsing to a flat color scale (False).
        The tolerance-based check exists to avoid color striping on
        zonal-mean plots of constant fields with tiny CS->LL
        regridding noise (see GitHub issue #330).  Single-level plots
        pass False (see GitHub issue #439).
        Default value: True
    data_scale : float, optional
        Magnitude of the Ref and Dev data that were differenced to
        produce this panel.  Used only for difference panels, to scale
        the tolerance that decides whether the panel holds real signal
        or only numerical noise.  If None, such a panel collapses only
        when vmin and vmax are exactly equal.
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

    if use_tolerance:
        # Tolerating tiny numerical noise (rather than requiring exact
        # equality) keeps it from stretching the colorbar into visible
        # color "striping" (see GitHub issue #330).  Each panel type is
        # anchored to a different value, so each needs its own atol.
        if is_ratio:
            # Anchored at 1 and dimensionless: rtol alone suffices
            is_constant = is_nearly_constant([vmin, vmax])
        elif is_difference:
            # Symmetric about zero (see vmin_vmax_for_absdiff_plots),
            # so rtol can never fire.  Scale atol to the data: a fixed
            # one, being in unknown units, blanked out small-magnitude
            # fields such as EmisCO_Aircraft (~1e-13 kg/m2/s).
            is_constant = is_nearly_constant(
                [vmin, vmax],
                atol=noise_atol(data_scale)
            )
        else:
            # Ref/Dev have no natural anchor: rtol only.  Still catches
            # the issue #330 noise, which is a relative error, without
            # mistaking small-but-real fields for constant.
            is_constant = is_nearly_constant([vmin, vmax], atol=0.0)
    else:
        is_constant = (
            (vmin == 0 and vmax == 0)
            or (np.isnan(vmin) and np.isnan(vmax))
        )

    if is_constant:
        # If the data is zero (or nearly constant) everywhere, or
        # undefined everywhere (vmin=vmax=NaN), then normalize the
        # data range so that the color corresponding to zero (white)
        # will be placed in the middle of the colorbar, where we will
        # add a single tick.
        #
        # Ratio data is centered on 1, not 0, so it needs its own
        # branch: using the 0-centered Normalize below would place
        # the near-1 "Ref and Dev are equal" data at the extreme edge
        # of the colorbar instead of in the middle.
        if is_ratio:
            return MidpointLogNorm(vmin=0.5, vmax=2.0, midpoint=1.0)
        if is_difference:
            return colors.Normalize(vmin=-1.0, vmax=1.0)
        return colors.Normalize(vmin=0.0, vmax=1.0)

    # For log color scales, assume a range 3 orders of magnitude
    # below the maximum value.  Otherwise use a linear scale.
    if log_color_scale and not ratio_log:
        return colors.LogNorm(vmin=vmax / 1e3, vmax=vmax)
    if log_color_scale:
        return MidpointLogNorm(vmin=vmin, vmax=vmax, midpoint=1)

    # For linear color scales: Normalize between min & max
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
