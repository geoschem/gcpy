"""
Creates a six-panel comparison plot.

Row 1: Model output (Ref version, Dev version)
Row 2: Abs difference (dynamic range and restricted range)
Row 3: Ratio (dynamic range and restricted range)

NOTE: For diff-of-diffs comparisons, Row 3 (Ratio) is replaced
by Fractional Difference (dynamic range and restricted range).

Also contains several helper routines that were split off
from the gcpy/plot.py.
"""
from matplotlib import ticker
import matplotlib.pyplot as plt
import numpy as np
from dask.array import Array as DaskArray
import xarray as xr
import cartopy.crs as ccrs
from gcpy.util import get_nan_mask, verify_variable_type
from gcpy.plot.core import gcpy_style, normalize_colors
from gcpy.plot.single_panel import single_panel

# Suppress numpy divide by zero warnings to prevent output spam
np.seterr(divide="ignore", invalid="ignore")

# Use a style sheet to control plot attributes
plt.style.use(gcpy_style)


def six_plot(
        subplot,
        all_zero,
        all_nan,
        plot_val,
        grid,
        axes,
        rowcol,
        title,
        comap,
        unit,
        extent,
        masked_data,
        other_all_nan,
        gridtype,
        vmins,
        vmaxs,
        use_cmap_RdBu,
        match_cbar,
        verbose,
        log_color_scale,
        pedge=np.full((1, 1), -1),
        pedge_ind=np.full((1, 1), -1),
        log_yaxis=False,
        xtick_positions=None,
        xticklabels=None,
        plot_type="single_level",
        ratio_log=False,
        proj=ccrs.PlateCarree(),
        ll_plot_func='imshow',
        **extra_plot_args
):
    """
    Plotting function to be called from compare_single_level or
    compare_zonal_mean. Primarily exists to eliminate code redundancy
    in the prior listed functions and has not been tested separately.

    Args:
    -----
    subplot: str
        Type of plot to create (ref, dev, absolute difference or
        fractional difference).
    all_zero: bool
        Set this flag to True if the data to be plotted consist
        only of zeros.
    all_nan: bool
        Set this flag to True if the data to be plotted consist
        only of NaNs.
    plot_val: xarray.DataArray, numpy.ndarray, or dask.array.Array
        Single data variable to plot.
    grid: dict
        Dictionary mapping plot_val to plottable coordinates
    axes: matplotlib.axes
        Axes object to plot information. Will create a new axes
        if none is passed.
    rowcol: tuple
        Subplot position in overall Figure.
    title: str
        Title to print on axes
    comap: matplotlib Colormap
        Colormap for plotting data values.
    unit: str
        Units of plotted data.
    extent: tuple (minlon, maxlon, minlat, maxlat)
        Describes minimum and maximum latitude and longitude of
        input data.
    masked_data: numpy array
        Masked area for cubed-sphere plotting.
    other_all_nan: bool
        Set this flag to True if plotting ref/dev and the other
        of ref/dev is all nan.
    gridtype: str
        "ll" for lat/lon or "cs" for cubed-sphere.
    vmins: list of float
        list of length 3 of minimum ref value, dev value,
        and absdiff value.
    vmaxs: list of float
        list of length 3 of maximum ref value, dev value,
        and absdiff value.
    use_cmap_RdBu: bool
        Set this flag to True to use a blue-white-red colormap
    match_cbar: bool
        Set this flag to True if you are plotting with the
        same colorbar for ref and dev.
    verbose: bool
        Set this flag to True to enable informative printout.
    log_color_scale: bool
        Set this flag to True to enable log-scale colormapping.

    Keyword Args (optional):
    ------------------------
    pedge: numpy array
        Edge pressures of grid cells in data to be plotted.
        Default value: np.full((1,1), -1)
    pedge_ind: numpy array
        Indices where edge pressure values are within a given
        pressure range.
        Default value: np.full((1,1), -1)
    log_yaxis: bool
        Set this flag to True to enable log scaling of pressure
        in zonal mean plots.
        Default value: False
     xtick_positions: list of float
         Locations of lat/lon or lon ticks on plot.
         Default value: None
     xticklabels: list of str
         Labels for lat/lon ticks.
         Default value: None
     plot_type: str
         Type of plot, either "single_level" or "zonal"mean".
         Default value: "single_level"
     ratio_log: bool
         Set this flag to True to enable log scaling for ratio plots
         Default value: False
     proj: cartopy projection
         Projection for plotting data.
         Default value: ccrs.PlateCarree()
     ll_plot_func: str
         Function to use for lat/lon single level plotting with
         possible values 'imshow' and 'pcolormesh'. imshow is much
         faster but is slightly displaced when plotting from dateline
         to dateline and/or pole to pole.
         Default value: 'imshow'
     extra_plot_args: various
         Any extra keyword arguments are passed through the
         plotting functions to be used in calls to pcolormesh() (CS)
         or imshow() (Lat/Lon).
    """
    verify_variable_type(plot_val, (np.ndarray, xr.DataArray, DaskArray))

    # Compute the min & max values
    vmin, vmax = compute_vmin_vmax_for_plot(
        plot_val,
        vmins,
        vmaxs,
        subplot,
        rowcol,
        all_zero=all_zero,
        all_nan=all_nan,
        other_all_nan=other_all_nan,
        match_cbar=match_cbar,
        use_cmap_RdBu=use_cmap_RdBu,
        verbose=verbose,
    )

    # Compute the norm object (i.e. put the colorscale on a
    # range of 0..1, which are matplotlib color coordinates)
    # (also remove NaNs in data for ratio plots)
    plot_val, norm = compute_norm_for_plot(
        plot_val,
        vmin,
        vmax,
        subplot,
        use_cmap_RdBu=use_cmap_RdBu,
        log_color_scale=log_color_scale,
        ratio_log=ratio_log
    )

    # Create one of the 6 subplots
    plot = single_panel(
        plot_val,
        axes,
        plot_type,
        grid,
        gridtype,
        title,
        comap,
        norm,
        unit,
        extent,
        masked_data,
        use_cmap_RdBu,
        log_color_scale,
        add_cb=False,
        pedge=pedge,
        pedge_ind=pedge_ind,
        log_yaxis=log_yaxis,
        xtick_positions=xtick_positions,
        xticklabels=xticklabels,
        proj=proj,
        ll_plot_func=ll_plot_func,
        **extra_plot_args)

    # Control how close to the plot the colorbar will go
    pad = 0.15
    if "single_level" in plot_type:
        pad = 0.025

    # Define the colorbar for the plot
    cbar = plt.colorbar(
        plot,
        ax=axes,
        orientation="horizontal",
        norm=norm,
        pad=pad
    )
    cbar.mappable.set_norm(norm)
    cbar = colorbar_ticks_and_format(
        plot_val,
        cbar,
        vmin,
        vmax,
        subplot,
        all_zero=all_zero,
        all_nan=all_nan,
        use_cmap_RdBu=use_cmap_RdBu,
        log_color_scale=log_color_scale,
    )
    cbar.set_label(unit)


def verbose_print(verbose, rowcol, vmin, vmax):
    """
    Routine to print the vmin & vmax values for each subplot.

    Args:
    -----
    verbose : bool
       Toggles informative prrintout on (True) or off (False).
    rowcol : int
       Subplot index.
    vmin, vmax : float
       Minimum and maximum of data range.
    """
    if verbose:
        print(f"Subplot ({rowcol}) vmin, vmax: {vmin}, {vmax}")


def compute_vmin_vmax_for_plot(
        plot_val,
        vmins,
        vmaxs,
        subplot,
        rowcol,
        all_zero=False,
        all_nan=False,
        other_all_nan=False,
        match_cbar=False,
        use_cmap_RdBu=False,
        verbose=False
):
    """
    Computes the min & max values for a subplot of a six-panel plot.

    Args:
    -----
    plot_val: xarray.DataArray, numpy.ndarray, or dask.array.Array
        Single data variable to plot.
    subplot: str
        Subplot name (see routine six_panel_subplot_names)
    vmins: list of float
        [minimum ref value, minimum dev value, absdiff value]
    vmaxs: list of float
        [maximum ref value, maximum dev value, absdiff value]

    Keyword Arguments (optional):
    -----------------------------
    all_zero: bool
        Indicates if the data consists of all zeros (True)
        or not (False)
    all_nan: bool
        Indicates if the data consists of all NaN values (True)
        or not (False)
    other_all_nan: bool
        Indicates if plotting ref/dev and the other of ref/dev contains
        all NaN values (True) or not (False).
    match_cbar: bool
        Toggles using the same colorbar for ref and dev on (True)
        or off (False).
    use_cmap_RdBu: bool
        Toggles a blue-white-red colormap on (True) or off (False).
    verbose: bool
        Toggles informative printout on (True) or off (False).

    Returns:
    --------
    vmin, vmax : float
        Min and max values for this subplot of a 6-panel plot
    """
    # ==================================================================
    # Get min and max values for Ref or Dev subplots
    # ==================================================================
    if subplot in ("ref", "dev"):
        return vmin_vmax_for_ref_dev_plots(
            subplot,
            rowcol,
            vmins,
            vmaxs,
            all_zero=all_zero,
            all_nan=all_nan,
            other_all_nan=other_all_nan,
            match_cbar=match_cbar,
            use_cmap_RdBu=use_cmap_RdBu,
            verbose=verbose
        )

    # ==================================================================
    # Get min and max values for Absdiff and Ratio subplots
    # ==================================================================

    # First check if all data is zero or NaN
    if all_zero:
        verbose_print(verbose, rowcol, 0, 0)
        return 0, 0
    if all_nan:
        verbose_print(verbose, rowcol, np.nan, np.nan)
        return np.nan, np.nan

    # Absdiff
    if subplot in ("dyn_absdiff", "res_absdiff"):
        return vmin_vmax_for_absdiff_plots(
            plot_val,
            subplot,
            rowcol,
            verbose=verbose
        )

    # Ratio
    if subplot in ("dyn_ratio", "res_ratio"):
        return vmin_vmax_for_ratio_plots(
            plot_val,
            subplot,
            rowcol,
            verbose=verbose
        )

    # Make sure the function returns a value.  This will avoid
    # an "inconsistent-return-statements" warning from Pylint.
    return None


def vmin_vmax_for_ref_dev_plots(
        subplot,
        rowcol,
        vmins,
        vmaxs,
        all_zero=False,
        all_nan=False,
        other_all_nan=False,
        match_cbar=False,
        use_cmap_RdBu=False,
        verbose=False,
):
    """
    Returns the vmin and vmax values for the "Ref" or "Dev"
    subplots of a six-panel plot.

    Args:
    -----
    subplot: str
        Subplot name (see routine six_panel_subplot_names).
    rowcol : int
        Subplot index.
    vmins: list of float
        [minimum ref value, minimum dev value, absdiff value]
    vmaxs: list of float
        [maximum ref value, maximum dev value, absdiff value]

    Keyword Arguments (optional):
    -----------------------------
    all_zero: bool
        Indicates if the data consists of all zeros (True)
        or not (False).
    all_nan: bool
        Indicates if the data consists of all NaN values (True)
        or not (False).
    other_all_nan: bool
        Indicates if plotting ref/dev and the other of ref/dev contains
        all NaN values (True) or not (False).
    match_cbar: bool
        Toggles using the same colorbar for ref and dev on (True)
        or off (False).
    use_cmap_RdBu: bool
        Toggles a blue-white-red colormap on (True) or off (False).
    verbose: bool
        Toggles informative printout on (True) or off (False).

    Returns:
    --------
    vmin, vmax : float
        Min and max values to plot.
    """
    #---------------------------------------------------------------
    # Data is all zero or Nan
    #---------------------------------------------------------------
    if all_zero or all_nan:
        [vmin, vmax] = [vmins[1], vmaxs[1]]
        if subplot == "ref":
            [vmin, vmax] = [vmins[0], vmaxs[0]]
        verbose_print(verbose, rowcol, vmin, vmax)
        return vmin, vmax

    #---------------------------------------------------------------
    # We are using a difference colormap (diff of diffs)
    #---------------------------------------------------------------
    if use_cmap_RdBu:

        # Ref supblot, diff-of-diffs
        if subplot in "ref":
            vmax = max([np.abs(vmins[0]), np.abs(vmaxs[0])])
            if match_cbar and not other_all_nan:
                vmax = max([np.abs(vmins[2]), np.abs(vmaxs[2])])
            verbose_print(verbose, rowcol, -vmax, vmax)
            return -vmax, vmax

        # Dev subplot, diff-of-diffs
        vmax = max([np.abs(vmins[1]), np.abs(vmaxs[1])])
        if match_cbar and not other_all_nan:
            vmax = max([np.abs(vmins[2]), np.abs(vmaxs[2])])
        verbose_print(verbose, rowcol, -vmax, vmax)
        return -vmax, vmax

    #---------------------------------------------------------------
    # We are using a gradient colormap
    #---------------------------------------------------------------

    # Ref subplot
    if subplot in "ref":
        [vmin, vmax] = [vmins[0], vmaxs[0]]
        if match_cbar and not other_all_nan:
            [vmin, vmax] = [vmins[2], vmaxs[2]]
        verbose_print(verbose, rowcol, vmin, vmax)
        return vmin, vmax

    # Dev subplot
    [vmin, vmax] = [vmins[1], vmaxs[1]]
    if match_cbar and not other_all_nan:
        [vmin, vmax] = [vmins[2], vmaxs[2]]
    verbose_print(verbose, rowcol, vmin, vmax)
    return vmin, vmax


def vmin_vmax_for_absdiff_plots(
        plot_val,
        subplot,
        rowcol,
        verbose=False,
):
    """
    Returns the vmin and vmax values for the "Absolute Difference
    (dynamic range)" or "Absolute Difference (restricted range)"
    subplots of a of a six-panel plot.

    Args:
    -----
    plot_val: xarray.DataArray, numpy.ndarray, or dask.array.Array
        Single data variable of GEOS-Chem output to plot.
    subplot: str
        Subplot name (see routine six_panel_subplot_names).
    rowcol : int
        Subplot index.

    Keyword Arguments (optional):
    -----------------------------
    verbose: bool
        Toggles informative printout on (True) or off (False).

    Returns:
    --------
    vmin, vmax : float
        Min and max values to plot.
    """
    # Absdiff (dynamic range) subplot: min & max (excluding NaNs)
    if subplot in "dyn_absdiff":
        vmax = max(
            [np.abs(np.nanmin(plot_val)), np.abs(np.nanmax(plot_val))]
        )
        verbose_print(verbose, rowcol, -vmax, vmax)
        return -vmax, vmax

    # Absdiff (restricted range) subplot
    if subplot in "res_absdiff":
        [pct5, pct95] = [
            np.percentile(plot_val, 5),
            np.percentile(plot_val, 95),
        ]
        vmax = np.max([np.abs(pct5), np.abs(pct95)])
        verbose_print(verbose, rowcol, -vmax, vmax)
        return -vmax, vmax

    # Make sure the function returns a value.  This will avoid
    # an "inconsistent-return-statements" warning from Pylint.
    return None


def vmin_vmax_for_ratio_plots(
        plot_val,
        subplot,
        rowcol,
        verbose=False,
):
    """
    Returns the vmin and vmax values for the "Ratio (dynamic range)"
    or "Ratio (restricted range) subplot of a six-panel plot.

    Args:
    -----
    plot_val: xarray.DataArray, numpy.ndarray, or dask.array.Array
        Single data variable to plot.
    subplot: str
        Subplot name (see routine six_panel_subplot_names).
    rowcol : int
        Subplot index.

    Keyword Arguments (optional):
    -----------------------------
    verbose: bool
        Toggles informative printout on (True) or off (False).

    Returns:
    --------
    vmin, vmax : float
        Min and max values to plot.
    """
    # Ratio (dynamic range) subplot)
    if subplot in "dyn_ratio":
        vmin = np.min(
            [np.abs(np.nanmin(plot_val)), np.abs(np.nanmax(plot_val))]
        )
        if np.abs(vmin) > 0.0:                     # If vmin > 0, compute
            vmax = 1.0 / vmin                      # vmax as its reciprocal
        else:
            vmax = np.abs(np.nanmax(plot_val))     # Otherwise compute vmin
            vmin = 1.0 / vmax                      # as reciprocal of vmax
        if vmin > vmax:
            vmin, vmax = vmax, vmin                # Swap values if needed
        verbose_print(verbose, rowcol, vmin, vmax)
        return vmin, vmax

    # Ratio (restricted range) subplot
    verbose_print(verbose, rowcol, 0.5, 2.0)
    return 0.5, 2.0


def compute_norm_for_plot(
        plot_val,
        vmin,
        vmax,
        subplot,
        use_cmap_RdBu=False,
        log_color_scale=False,
        ratio_log=False,
):
    """
    Normalize colors (put into range [0..1] for matplotlib methods).

    Args:
    -----
    plot_val: xarray.DataArray, numpy.ndarray, or dask.array.Array
        Single data variable GEOS-Chem output to plot
    vmin, vmax : float
        Min and max value for this subplot of a 6-panel plot.
    subplot: str
        Subplot name (see routine six_panel_subplot_names)

    Keyword Arguments (optional):
    -----------------------------
    use_cmap_RdBu: bool
        Toggles a blue-white-red colormap on (True) or off (False).
    log_color_scale : bool
        Toggles a logarithmic color scale on (True) or off (False).
    ratio_log : bool
        Toggles log scaling for ratio plots on (True) or not (False).
    verbose: bool
        Toggles informative printout on (True) or off (False).

    Returns:
    --------
    vmin, vmax : float
        Min and max values for this subplot of a 6-panel plot
    """
    # ==================================================================
    # Ref and Dev subplots
    # ==================================================================
    if subplot in ("ref", "dev"):
        return plot_val, normalize_colors(
            vmin,
            vmax,
            is_difference=use_cmap_RdBu,
            log_color_scale=log_color_scale,
            ratio_log=ratio_log
        )

    # ==================================================================
    # Absdiff (dynamic & restricted range) subplots
    # ==================================================================
    if subplot in ("dyn_absdiff", "res_absdiff"):
        return plot_val, normalize_colors(
            vmin,
            vmax,
            is_difference=True
        )

    # ==================================================================
    # Ratio (dynamic & restricted range) subplots
    # Remove NaNs for compatibility with color normalization
    # ==================================================================
    plot_val = get_nan_mask(plot_val)
    return plot_val, normalize_colors(
        vmin,
        vmax,
        is_difference=True,
        log_color_scale=True,
        ratio_log=ratio_log
    )


def colorbar_ticks_and_format(
        plot_val,
        cbar,
        vmin,
        vmax,
        subplot,
        all_zero=False,
        all_nan=False,
        use_cmap_RdBu=False,
        log_color_scale=False,
):
    """
    Adjusts colorbar tick placement and label formatting style
    for a subplot of a 6-panel plot.  Called from routine six_plot.

    Args:
    -----
    plot_val: xarray.DataArray, numpy.ndarray, or dask.array.Array
        Single data variable to plot.
    cbar : matplotlib.colorbar.Colorbar
        The input colorbar.
    vmin, vmax : float
        Min and max of the data range to plot.
    subplot: str
        Subplot name (see routine six_panel_subplot_names).

    Keyword Arguments (optional):
    -----------------------------
    all_zero: bool
        Indicates if the data consists of all zeros (True)
        or not (False).
    all_nan: bool
        Indicates if the data consists of all NaN values (True)
        or not (False).
    use_cmap_RdBu: bool
        Toggles a blue-white-red colormap on (True) or off (False).
    log_color_scale : bool
        Toggles a logarithmic color scale on (True) or off (False).

    Returns:
    --------
    cbar : matplotlib.colorbar.Colorbar
        The modified colorbar.
    """
    # ==================================================================
    # Data is all zero or NaN:
    # Place a single tick with an appropriate label in the middle.
    # For RdBu colortables this goes at 0.0; otherwise at 0.5.
    # ==================================================================
    if all_zero or all_nan:
        return colorbar_for_all_zero_or_nan(
            cbar,
            subplot,
            all_nan=all_nan,
            use_cmap_RdBu=use_cmap_RdBu,
        )

    # ==================================================================
    # Data is plottable: Pick the locations and format of tick
    # labels depending the subplot and the colormap that is used.
    # ==================================================================

    #-------------------------------------------------------------------
    # Ref and Dev subplots, log scale
    #-------------------------------------------------------------------
    if subplot in ("ref", "dev") and log_color_scale:
        cbar.formatter = ticker.LogFormatter(base=10)
        cbar.minorticks_off()
        return cbar

    #-------------------------------------------------------------------
    # Ratio (dynamic and restricted range) subplots):
    #-------------------------------------------------------------------
    if subplot in ("dyn_ratio", "res_ratio"):

        def ref_equals_dev(array):
            """
            Internal routine to check that returns true if all elements
            of Ref/Dev are equal to 1 or NaN (aka missing value).
            This is needed to be able to add a ticklabel stating
            that Ref & Dev are equal throughout the domain.
            """
            uniq = np.unique(array)
            if len(uniq) == 2:
                return np.any(np.isin(uniq, [1.0])) and np.any(np.isnan(uniq))
            return np.all(np.isin(uniq, [1.0]))

        # When Ref == Dev
        if ref_equals_dev(plot_val):
            return colorbar_for_ref_equals_dev(cbar)

        # Dynamic range ratio subplot
        if subplot in "dyn_ratio":
            return colorbar_for_dyn_ratio_plots(cbar, vmin, vmax)

        # Restricted range ratio subplot
        return colorbar_for_res_ratio_plots(cbar)

    #-------------------------------------------------------------------
    # For the following subplots:
    # (1) Ref & Dev, with non-log color scales
    # (2) Absdiff (dynamic range)
    # (3) Absdiff (restricted range)
    #-------------------------------------------------------------------

    # For data ranges between 0.1 and 100:
    if 0.1 < (vmax - vmin) < 100.0:
        return colorbar_for_small_data_range(
            cbar,
            vmin,
            vmax,
            diff_cmap=(use_cmap_RdBu or "absdiff" in subplot)
        )

    # For larger data ranges, automatically find good tick locations
    # (but not too many that the labels smush together)
    cbar.locator = ticker.MaxNLocator(nbins=4)
    cbar.minorticks_off()
    return cbar


def colorbar_for_all_zero_or_nan(
        cbar,
        subplot,
        all_nan=False,
        use_cmap_RdBu=False,
):
    """
    Formats a colorbar object for the case when Ref or Dev
    contains either all zeroes or all NaNs.

    Args:
    -----
    cbar : matplotlib.colorbar.Colorbar
        The input colorbar.
    subplot : str
        Name of this subplot of a 6-panel plot.

    Keyword Args (optional):
    ------------------------
    all_nan : bool
        Indicates that the data array contains all NaN values (True)
        or not (False).
    use_cmap_RdBu : bool
        Indicates that we are using a difference colortable (True)
        or not (False).

    Returns:
    --------
    cbar : matplotlib.colorbar.Colorbar
        The modified colorbar
    """
    pos = [0.0]
    if subplot in ("ref", "dev"):
        if not use_cmap_RdBu:
            pos = [0.5]
    labels = ["Zero throughout domain"]
    if all_nan:
        labels = ["Undefined throughout domain"]
    cbar.set_ticks(pos, labels=labels)
    cbar.minorticks_off()
    return cbar


def colorbar_for_ref_equals_dev(cbar):
    """
    Formats a colorbar object for the case when Ref and Dev
    are equal throughout the domain.

    Args:
    -----
    cbar : matplotlib.colorbar.Colorbar
        The input colorbar.

    Returns:
    --------
    cbar : matplotlib.colorbar.Colorbar
        The modified colorbar.
    """
    pos = [1.0]
    cbar.set_ticks(
        pos,
        labels=["Ref and Dev equal throughout domain"]
    )
    cbar.minorticks_off()
    return cbar


def colorbar_for_dyn_ratio_plots(
        cbar,
        vmin,
        vmax
):
    """
    Formats a colorbar object for the "dynamic range ratio"
    subplot of a six-panel plot.

    Args:
    -----
    cbar : matplotlib.colorbar.Colorbar
        The input colorbar.
    vmin, vmax : float
        Min and max of the data range.

    Returns:
    --------
    cbar : matplotlib.colorbar.Colorbar
        The modified colorbar.
    """
    # If the ratio is in the range 0.999 and 1.001, then
    # place tickmarks at [vmin, 1, vmax].  This should help
    # to avoid the tick labels from running together.
    if vmin > 0.999 and vmax < 1.001:
        pos = [vmin, 1.0, vmax]
        cbar.set_ticks(pos)
        cbar.formatter = ticker.ScalarFormatter()
        cbar.formatter.set_useOffset(False)
        cbar.minorticks_off()
        return cbar

    # If the ratio is in the range 0.1 .. 10.0, then place
    # tickmarks [vmin, avg(vmin,1), 1, avg(vmax,1), vmax].
    # This should be good enough for most cases.  Perhaps
    # think about implementing a better method later on.
    if vmin > 0.1 and vmax < 10.0:
        pos = [vmin, (vmin+1.0)/2.0, 1.0, (vmax+1.0)/2.0, vmax]
        cbar.set_ticks(pos)
        cbar.formatter = ticker.ScalarFormatter()
        cbar.formatter.set_useOffset(False)
        cbar.minorticks_off()
        return cbar

    # Use LogLocator and LogFormatter for larger data ranges
    cbar.locator = ticker.LogLocator(base=10, subs='all')
    cbar.formatter = ticker.LogFormatter(base=10)
    cbar.minorticks_off()
    return cbar


def colorbar_for_res_ratio_plots(cbar):
    """
    Formats a colorbar object for the "restricted range ratio"
    subplot of a six-panel plot.

    Args:
    -----
    cbar : matplotlib.colorbar.Colorbar
        The input colorbar.

    Returns:
    --------
    cbar : matplotlib.colorbar.Colorbar
        The modified colorbar.
    """
    # Use fixed ticks and ScalarFormatter
    pos = [0.5, 0.75, 1.0, 1.5, 2.0]
    cbar.set_ticks(pos)
    cbar.formatter = ticker.ScalarFormatter()
    cbar.minorticks_off()
    return cbar


def colorbar_for_small_data_range(
        cbar,
        vmin,
        vmax,
        diff_cmap=False,
):
    """
    Formats a colorbar object for data that falls within the range
    of 0.1 .. 100.

    Args:
    -----
    cbar : matplotlib.colorbar.Colorbar
        The input colorbar.
    vmin, vmax : float
        Min and max of the data range.
    diff_cmap : bool
        Indicates that we are using a diverging colortable (True)
        or not (False).

    Returns:
    --------
    cbar : matplotlib.colorbar.Colorbar
        The modified colorbar.
    """
    # If using a difference colormap (e.g. for absdiff),
    # then place ticks symmetrically around zero.
    if diff_cmap:
        pos = [vmin, vmin/2.0, 0.0, vmax/2.0, vmax]
        cbar.set_ticks(pos)
        cbar.formatter = ticker.ScalarFormatter()
        cbar.formatter.set_useOffset(False)
        cbar.minorticks_off()
        return cbar

    # Otherwise place ticks symmetrically along the data range
    vrange = vmax - vmin
    pos = [vmin, vmin+vrange*0.25, vmin+vrange*0.5, vmin+vrange*0.75, vmax]
    cbar.set_ticks(pos)
    cbar.formatter = ticker.ScalarFormatter()
    cbar.formatter.set_useOffset(False)
    cbar.minorticks_off()
    return cbar
