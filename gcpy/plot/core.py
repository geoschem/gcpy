"""
Common variables and functions used by modules in gcpy.plot.
"""
from os import path
import warnings
from matplotlib import colors
import numpy as np

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


def six_panel_subplot_names(diff_of_diffs):
    """
    Returns the names of the subplots for the 6-panel plots.

    Args:
    -----
    diff_of_diffs : bool
        Indicates if this is a diff-of-diffs benchmark (True)
        or not (False),  Ratio plots are only included if
        diff_of_diffs is False.

    Returns:
    --------
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


def normalize_colors(
        vmin,
        vmax,
        is_difference=False,
        log_color_scale=False,
        ratio_log=False
):
    """
    Normalizes a data range to the colormap range used by matplotlib
    functions. For log-color scales, special handling is done to prevent
    taking the log of data that is all zeroes.

    Args:
        vmin: float
            Minimum value of the data range.
        vmax: float
            Maximum value of the data range.

    Keyword Args (optional):
        is_difference: bool
            Set this switch to denote that we are using a difference
            color scale (i.e. with zero in the middle of the range).
            Default value: False
        log_color_scale: bool
            Logical flag to denote that we are using a logarithmic
            color scale instead of a linear color scale.
            Default value: False
        ratio_log : bool
            Indicates whether we are using log scaling for ratio plots
            (True) or not (False).
            Default value: False

    Returns:
        norm: matplotlib Norm
            The normalized matplotlib color range, stored in
            a matplotlib Norm object.

    Remarks:
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

    # Absolute value of v
    abs_vmin = abs(vmin)
    abs_vmax = abs(vmax)

    if (abs_vmin == 0 and abs_vmax == 0) or \
       (np.isnan(vmin) and np.isnan(vmax)):
        # If the data is zero everywhere (vmin=vmax=0) or undefined
        # everywhere (vmin=vmax=NaN), then normalize the data range
        # so that the color corresponding to zero (white) will be
        # placed in the middle of the colorbar, where we will
        # add a single tick.
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
