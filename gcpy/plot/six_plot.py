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
from gcpy.util import get_nan_mask, is_nearly_constant, verify_variable_type
from gcpy.plot.core import NOISE_REL_TOL, RATIO_ABS_TOL, \
    constant_rel_tol, gcpy_style, noise_atol, normalize_colors
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
        yaxis_units="pressure",
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

    Parameters
    ----------
    subplot : str
        Type of plot to create (ref, dev, absolute difference or
        fractional difference).
    all_zero : bool
        Set this flag to True if the data to be plotted consist
        only of zeros.
    all_nan : bool
        Set this flag to True if the data to be plotted consist
        only of NaNs.
    plot_val : xarray.DataArray or numpy.ndarray or dask.array.Array
        Single data variable to plot.
    grid : dict
        Dictionary mapping plot_val to plottable coordinates.
    axes : matplotlib.axes.Axes
        Axes object to plot information. Will create a new axes
        if none is passed.
    rowcol : tuple
        Subplot position in overall Figure.
    title : str
        Title to print on axes.
    comap : matplotlib.colors.Colormap
        Colormap for plotting data values.
    unit : str
        Units of plotted data.
    extent : tuple of float
        Describes minimum and maximum latitude and longitude of
        input data in the form (minlon, maxlon, minlat, maxlat).
    masked_data : numpy.ndarray
        Masked area for cubed-sphere plotting.
    other_all_nan : bool
        Set this flag to True if plotting ref/dev and the other
        of ref/dev is all nan.
    gridtype : str
        "ll" for lat/lon or "cs" for cubed-sphere.
    vmins : list of float
        List of length 3 of minimum ref value, dev value,
        and minimum of both (the last for use with match_cbar=True).
        Also scales the difference panels' noise tolerance (see
        routine ref_dev_data_scale).
    vmaxs : list of float
        List of length 3 of maximum ref value, dev value,
        and maximum of both (the last for use with match_cbar=True).
        See the note on vmins above.
    use_cmap_RdBu : bool
        Set this flag to True to use a blue-white-red colormap.
    match_cbar : bool
        Set this flag to True if you are plotting with the
        same colorbar for ref and dev.
    verbose : bool
        Set this flag to True to enable informative printout.
    log_color_scale : bool
        Set this flag to True to enable log-scale colormapping.
    pedge : numpy.ndarray, optional
        Edge pressures of grid cells in data to be plotted.
        Default value: np.full((1,1), -1)
    pedge_ind : numpy.ndarray, optional
        Indices where edge pressure values are within a given
        pressure range.
        Default value: np.full((1,1), -1)
    log_yaxis : bool, optional
        Set this flag to True to enable log scaling of pressure
        in zonal mean plots.
        Default value: False
    yaxis_units : str, optional
        Units to use for the Y-axis of zonal mean plots. Either
        "pressure" (hPa) or "level" (model vertical level index).
        Default value: "pressure"
    xtick_positions : list of float, optional
        Locations of lat/lon or lon ticks on plot.
        Default value: None
    xticklabels : list of str, optional
        Labels for lat/lon ticks.
        Default value: None
    plot_type : str, optional
        Type of plot, either "single_level" or "zonal_mean".
        Default value: "single_level"
    ratio_log : bool, optional
        Set this flag to True to enable log scaling for ratio plots.
        Default value: False
    proj : cartopy.crs.Projection, optional
        Projection for plotting data.
        Default value: ccrs.PlateCarree()
    ll_plot_func : str, optional
        Function to use for lat/lon single level plotting with
        possible values 'imshow' and 'pcolormesh'. imshow is much
        faster but is slightly displaced when plotting from dateline
        to dateline and/or pole to pole.
        Default value: 'imshow'
    **extra_plot_args
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

    # Magnitude of this panel's data, used to scale its noise
    # tolerance.  Row 2 is in the physical units of Ref & Dev; row 3
    # is a ratio (or, for diff-of-diffs, a difference of fractional
    # differences), which is dimensionless.  Key it off the row index,
    # not the subplot name: six_panel_subplot_names() names rows 2 and
    # 3 identically for diff-of-diffs plots.
    data_scale = None
    if rowcol[0] == 1:
        data_scale = ref_dev_data_scale(vmins, vmaxs)
    elif rowcol[0] == 2:
        data_scale = 1.0

    # Whether Ref or Dev is zero everywhere.  vmins & vmaxs hold the
    # Ref, Dev and combined ranges, so a field whose min and max are
    # both zero is zero throughout the domain.  That is what makes a
    # Dev/Ref ratio undefined (Ref zero) or identically zero (Dev zero).
    ref_is_zero = vmins[0] == 0 and vmaxs[0] == 0
    dev_is_zero = vmins[1] == 0 and vmaxs[1] == 0

    # Compute the norm object (i.e. put the colorscale on a
    # range of 0..1, which are matplotlib color coordinates)
    # (also remove NaNs in data for ratio plots)
    plot_val, norm = compute_norm_for_plot(
        plot_val,
        vmin,
        vmax,
        subplot,
        plot_type=plot_type,
        use_cmap_RdBu=use_cmap_RdBu,
        log_color_scale=log_color_scale,
        ratio_log=ratio_log,
        data_scale=data_scale,
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
        yaxis_units=yaxis_units,
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
        data_scale=data_scale,
        use_tolerance="zonal_mean" in plot_type,
        ref_is_zero=ref_is_zero,
        dev_is_zero=dev_is_zero,
    )
    cbar.set_label(unit)


def verbose_print(verbose, rowcol, vmin, vmax):
    """
    Routine to print the vmin & vmax values for each subplot.

    Parameters
    ----------
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

    Parameters
    ----------
    plot_val : xarray.DataArray or numpy.ndarray or dask.array.Array
        Single data variable to plot.
    vmins : list of float
        [minimum ref value, minimum dev value, absdiff value].
    vmaxs : list of float
        [maximum ref value, maximum dev value, absdiff value].
    subplot : str
        Subplot name (see routine six_panel_subplot_names).
    rowcol : int
        Subplot index.
    all_zero : bool, optional
        Indicates if the data consists of all zeros (True)
        or not (False).
        Default value: False
    all_nan : bool, optional
        Indicates if the data consists of all NaN values (True)
        or not (False).
        Default value: False
    other_all_nan : bool, optional
        Indicates if plotting ref/dev and the other of ref/dev contains
        all NaN values (True) or not (False).
        Default value: False
    match_cbar : bool, optional
        Toggles using the same colorbar for ref and dev on (True)
        or off (False).
        Default value: False
    use_cmap_RdBu : bool, optional
        Toggles a blue-white-red colormap on (True) or off (False).
        Default value: False
    verbose : bool, optional
        Toggles informative printout on (True) or off (False).
        Default value: False

    Returns
    -------
    vmin, vmax : float
        Min and max values for this subplot of a 6-panel plot.
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

    Parameters
    ----------
    subplot : str
        Subplot name (see routine six_panel_subplot_names).
    rowcol : int
        Subplot index.
    vmins : list of float
        [minimum ref value, minimum dev value, absdiff value].
    vmaxs : list of float
        [maximum ref value, maximum dev value, absdiff value].
    all_zero : bool, optional
        Indicates if the data consists of all zeros (True)
        or not (False).
        Default value: False
    all_nan : bool, optional
        Indicates if the data consists of all NaN values (True)
        or not (False).
        Default value: False
    other_all_nan : bool, optional
        Indicates if plotting ref/dev and the other of ref/dev contains
        all NaN values (True) or not (False).
        Default value: False
    match_cbar : bool, optional
        Toggles using the same colorbar for ref and dev on (True)
        or off (False).
        Default value: False
    use_cmap_RdBu : bool, optional
        Toggles a blue-white-red colormap on (True) or off (False).
        Default value: False
    verbose : bool, optional
        Toggles informative printout on (True) or off (False).
        Default value: False

    Returns
    -------
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

        # Ref subplot, diff-of-diffs
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

    Parameters
    ----------
    plot_val : xarray.DataArray or numpy.ndarray or dask.array.Array
        Single data variable of GEOS-Chem output to plot.
    subplot : str
        Subplot name (see routine six_panel_subplot_names).
    rowcol : int
        Subplot index.
    verbose : bool, optional
        Toggles informative printout on (True) or off (False).
        Default value: False

    Returns
    -------
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
    # NOTE: Use the NaN-safe percentile, to match the dynamic-range
    # branch above.  np.percentile returns NaN if the array holds even
    # one NaN, which collapses the panel to a flat color scale while
    # it still draws the real (saturated) data underneath.
    if subplot in "res_absdiff":
        [pct5, pct95] = [
            np.nanpercentile(plot_val, 5),
            np.nanpercentile(plot_val, 95),
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

    Parameters
    ----------
    plot_val : xarray.DataArray or numpy.ndarray or dask.array.Array
        Single data variable to plot.
    subplot : str
        Subplot name (see routine six_panel_subplot_names).
    rowcol : int
        Subplot index.
    verbose : bool, optional
        Toggles informative printout on (True) or off (False).
        Default value: False

    Returns
    -------
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


def ref_dev_data_scale(vmins, vmaxs):
    """
    Returns the magnitude of the Ref and Dev data of a six-panel plot,
    which scales the difference panels' noise tolerance (see
    gcpy.plot.core.noise_atol).

    Parameters
    ----------
    vmins : list of float
        [minimum ref value, minimum dev value, minimum of both].
    vmaxs : list of float
        [maximum ref value, maximum dev value, maximum of both].

    Returns
    -------
    data_scale : float or None
        Largest absolute value among the finite entries, or None if
        there are none (e.g. Ref and Dev are both all NaN).

    Notes
    -----
    All six entries are scanned, not just the "both" pair, because
    compare_zonal_mean builds those with np.min & np.max rather than
    their NaN-safe forms, so an all-NaN Ref would poison them.
    """
    vals = np.abs(np.asarray(list(vmins) + list(vmaxs), dtype=float))
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return None

    return float(vals.max())


def compute_norm_for_plot(
        plot_val,
        vmin,
        vmax,
        subplot,
        plot_type="single_level",
        use_cmap_RdBu=False,
        log_color_scale=False,
        ratio_log=False,
        data_scale=None,
):
    """
    Normalize colors (put into range [0..1] for matplotlib methods).

    Parameters
    ----------
    plot_val : xarray.DataArray or numpy.ndarray or dask.array.Array
        Single data variable GEOS-Chem output to plot.
    vmin, vmax : float
        Min and max value for this subplot of a 6-panel plot.
    subplot : str
        Subplot name (see routine six_panel_subplot_names).
    plot_type : str, optional
        Either "single_level" or "zonal_mean".  Only zonal-mean plots
        use a tolerance (rather than exact equality) to decide if
        vmin/vmax are "nearly constant", since that is the only plot
        type where CS->LL regridding noise on constant fields has been
        observed to cause color striping (see GitHub issue #330).
        Default value: "single_level"
    use_cmap_RdBu : bool, optional
        Toggles a blue-white-red colormap on (True) or off (False).
        Default value: False
    log_color_scale : bool, optional
        Toggles a logarithmic color scale on (True) or off (False).
        Default value: False
    ratio_log : bool, optional
        Toggles log scaling for ratio plots on (True) or not (False).
        Default value: False
    data_scale : float, optional
        Magnitude of the Ref and Dev data (see ref_dev_data_scale).
        Used only by the difference panels, to scale the tolerance
        below which their spread counts as noise instead of signal.
        Default value: None

    Returns
    -------
    vmin, vmax : float
        Min and max values for this subplot of a 6-panel plot.
    """
    use_tolerance = "zonal_mean" in plot_type

    # ==================================================================
    # Ref and Dev subplots
    # ==================================================================
    if subplot in ("ref", "dev"):
        return plot_val, normalize_colors(
            vmin,
            vmax,
            is_difference=use_cmap_RdBu,
            log_color_scale=log_color_scale,
            ratio_log=ratio_log,
            use_tolerance=use_tolerance,
            rel_tol=constant_rel_tol(plot_val),
        )

    # ==================================================================
    # Absdiff (dynamic & restricted range) subplots
    # ==================================================================
    if subplot in ("dyn_absdiff", "res_absdiff"):
        return plot_val, normalize_colors(
            vmin,
            vmax,
            is_difference=True,
            use_tolerance=use_tolerance,
            data_scale=data_scale,
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
        ratio_log=ratio_log,
        is_ratio=True,
        use_tolerance=use_tolerance,
    )


def ref_equals_dev(array, rtol=NOISE_REL_TOL, atol=RATIO_ABS_TOL):
    """
    Returns True if all finite elements of a Ref/Dev ratio array are
    within tolerance of 1.0 (i.e. Ref is essentially equal to Dev
    throughout the domain).  This is needed to be able to add a
    ticklabel stating that Ref & Dev are equal throughout the domain.

    A tolerance (rather than exact equality) is used so that tiny
    numerical noise (e.g. from regridding, see GitHub issue #330)
    does not prevent this special case from being detected.

    Parameters
    ----------
    array : numpy array
        Ref/Dev ratio data (may contain NaNs).
    rtol : float, optional
        Relative tolerance.
        Default value: gcpy.plot.core.NOISE_REL_TOL
    atol : float, optional
        Absolute tolerance.
        Default value: gcpy.plot.core.RATIO_ABS_TOL

    Returns
    -------
    is_equal : bool
        Whether Ref and Dev are essentially equal throughout the
        domain.
    """
    return is_nearly_constant(array, rtol=rtol, atol=atol, target=1.0)


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
        data_scale=None,
        use_tolerance=False,
        ref_is_zero=False,
        dev_is_zero=False,
):
    """
    Adjusts colorbar tick placement and label formatting style
    for a subplot of a 6-panel plot.  Called from routine six_plot.

    Parameters
    ----------
    plot_val : xarray.DataArray or numpy.ndarray or dask.array.Array
        Single data variable to plot.
    cbar : matplotlib.colorbar.Colorbar
        The input colorbar.
    vmin, vmax : float
        Min and max of the data range to plot.
    subplot : str
        Subplot name (see routine six_panel_subplot_names).
    all_zero : bool, optional
        Indicates if the data consists of all zeros (True)
        or not (False).
        Default value: False
    all_nan : bool, optional
        Indicates if the data consists of all NaN values (True)
        or not (False).
        Default value: False
    use_cmap_RdBu : bool, optional
        Toggles a blue-white-red colormap on (True) or off (False).
        Default value: False
    log_color_scale : bool, optional
        Toggles a logarithmic color scale on (True) or off (False).
        Default value: False
    data_scale : float, optional
        Magnitude of the Ref and Dev data (see ref_dev_data_scale).
        Used to recognize a difference panel that was collapsed for
        holding only noise, so that we can label it as such.
        Default value: None
    use_tolerance : bool, optional
        Whether normalize_colors was called with its near-constant
        tolerance enabled (i.e. this is a zonal-mean plot).
        Default value: False
    ref_is_zero, dev_is_zero : bool, optional
        Whether the Ref (resp. Dev) data are zero throughout the
        domain.  Used to say why a ratio panel has nothing to show.
        Default value: False

    Returns
    -------
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
            ref_is_zero=ref_is_zero,
            dev_is_zero=dev_is_zero,
        )

    # ==================================================================
    # Data is plottable: Pick the locations and format of tick
    # labels depending the subplot and the colormap that is used.
    # ==================================================================

    #-------------------------------------------------------------------
    # Ref & Dev subplots that normalize_colors collapsed for being
    # constant across the domain.  Without this they keep numeric
    # ticks taken from a dimensionless 0-1 (or -1..1) norm while the
    # colorbar is labeled in the field's units, so the panel reads as
    # though the field topped out at 1 -- e.g. a TransportTracers
    # PassiveTracer restart at 100 ppb appearing to max out at 1 ppb.
    #
    # Checked ahead of the log-scale branch below because
    # normalize_colors likewise settles "is constant" before it
    # considers log scaling, and returns a linear collapsed norm
    # either way.
    #-------------------------------------------------------------------
    if (
            subplot in ("ref", "dev")
            and use_tolerance
            and is_nearly_constant(
                [vmin, vmax],
                rtol=constant_rel_tol(plot_val),
                atol=0.0,
            )
    ):
        return colorbar_for_constant_field(
            cbar,
            vmin,
            vmax,
            use_cmap_RdBu=use_cmap_RdBu,
        )

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

        # When Ref == Dev
        if ref_equals_dev(plot_val):
            return colorbar_for_ref_equals_dev(cbar)

        # Dynamic range ratio subplot
        if subplot in "dyn_ratio":
            return colorbar_for_dyn_ratio_plots(cbar, vmin, vmax)

        # Restricted range ratio subplot
        return colorbar_for_res_ratio_plots(cbar)

    #-------------------------------------------------------------------
    # Absdiff subplots that normalize_colors collapsed for holding
    # only noise.  Without this they show a blank plot over a -1..1
    # colorbar that has nothing to do with the data.
    #-------------------------------------------------------------------
    if (
            "absdiff" in subplot
            and use_tolerance
            and is_nearly_constant(
                [vmin, vmax],
                atol=noise_atol(data_scale)
            )
    ):
        # A restricted-range panel of a sparse field (e.g. aircraft
        # emissions, which are zero over most of the domain) collapses
        # because its 5th and 95th percentiles are both zero, not
        # because the differences are negligible -- the dynamic-range
        # panel beside it may well show real structure.  Say so.
        if vmin == 0 and vmax == 0:
            return colorbar_for_flat_restricted_range(cbar)
        return colorbar_for_negligible_diff(cbar)

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
        ref_is_zero=False,
        dev_is_zero=False,
):
    """
    Formats a colorbar object for the case when Ref or Dev
    contains either all zeroes or all NaNs.

    Parameters
    ----------
    cbar : matplotlib.colorbar.Colorbar
        The input colorbar.
    subplot : str
        Name of this subplot of a 6-panel plot.
    all_nan : bool, optional
        Indicates that the data array contains all NaN values (True)
        or not (False).
        Default value: False
    use_cmap_RdBu : bool, optional
        Indicates that we are using a difference colortable (True)
        or not (False).
        Default value: False
    ref_is_zero, dev_is_zero : bool, optional
        Whether the Ref (resp. Dev) data are zero throughout the
        domain.  Used to say why a ratio panel has nothing to show.
        Default value: False

    Returns
    -------
    cbar : matplotlib.colorbar.Colorbar
        The modified colorbar.
    """
    # Place the tick at the value the panel's color scale is anchored
    # to.  Ratio panels are anchored at 1, and their norm spans
    # [0.5, 2.0]; a tick at 0.0 falls outside that range, which
    # stretches the colorbar axes and leaves it blank.
    pos = [0.0]
    if subplot in ("ref", "dev"):
        if not use_cmap_RdBu:
            pos = [0.5]
    elif subplot in ("dyn_ratio", "res_ratio"):
        pos = [1.0]
    labels = ["Zero throughout domain"]
    if all_nan:
        labels = ["Undefined throughout domain"]

    # A ratio panel with nothing to show is far more informative if it
    # says which side vanished.  Dev/Ref is undefined when Ref is zero
    # and identically zero when Dev is zero; when both are zero it is
    # 0/0, which really is just undefined, so that case is left alone.
    if subplot in ("dyn_ratio", "res_ratio") and ref_is_zero != dev_is_zero:
        if ref_is_zero:
            labels = ["Ref is zero throughout domain"]
        else:
            labels = ["Dev is zero throughout domain"]

    cbar.set_ticks(pos, labels=labels)
    cbar.minorticks_off()
    return cbar


def colorbar_for_ref_equals_dev(cbar):
    """
    Formats a colorbar object for the case when Ref and Dev
    are equal throughout the domain.

    Parameters
    ----------
    cbar : matplotlib.colorbar.Colorbar
        The input colorbar.

    Returns
    -------
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


def colorbar_for_constant_field(cbar, vmin, vmax, use_cmap_RdBu=False):
    """
    Formats a colorbar object for a Ref or Dev subplot whose data is
    constant across the domain (to within constant_rel_tol), and whose
    color scale normalize_colors therefore collapsed to a flat range.

    That collapsed range is dimensionless -- [0, 1], or [-1, 1] for a
    difference colormap -- so leaving numeric ticks on it labels the
    colorbar in the field's units with values the field never takes.
    A single tick naming the constant is both honest and more useful.

    Parameters
    ----------
    cbar : matplotlib.colorbar.Colorbar
        The input colorbar.
    vmin, vmax : float
        Min and max of the data range, in the field's own units
        (i.e. before normalize_colors collapsed them).
    use_cmap_RdBu : bool, optional
        Whether this panel uses a blue-white-red difference colormap
        (True) or not (False).  Sets which anchor the tick goes on.
        Default value: False

    Returns
    -------
    cbar : matplotlib.colorbar.Colorbar
        The modified colorbar.
    """
    # Match the tick to the anchor of the collapsed norm that
    # normalize_colors returned: [-1, 1] for a difference colormap,
    # otherwise [0, 1].
    pos = [0.0] if use_cmap_RdBu else [0.5]

    # vmin and vmax agree to within constant_rel_tol, so either one
    # names the constant; the midpoint avoids favoring an endpoint.
    value = 0.5 * (float(vmin) + float(vmax))
    if np.isfinite(value):
        label = f"Constant at {value:.6g} throughout domain"
    else:
        label = "Constant throughout domain"

    cbar.set_ticks(pos, labels=[label])
    cbar.minorticks_off()
    return cbar


def colorbar_for_negligible_diff(cbar):
    """
    Formats a colorbar object for an "absolute difference" subplot in
    which the Dev - Ref spread is negligible compared to the Ref and
    Dev data themselves, and so was collapsed to a flat color range.

    Parameters
    ----------
    cbar : matplotlib.colorbar.Colorbar
        The input colorbar.

    Returns
    -------
    cbar : matplotlib.colorbar.Colorbar
        The modified colorbar.
    """
    pos = [0.0]
    cbar.set_ticks(
        pos,
        labels=["Differences negligible throughout domain"]
    )
    cbar.minorticks_off()
    return cbar


def colorbar_for_flat_restricted_range(cbar):
    """
    Formats a colorbar object for a "restricted range" subplot whose
    5th and 95th percentiles are both zero, which happens when a field
    is zero over most of the domain (e.g. aircraft emissions).

    Parameters
    ----------
    cbar : matplotlib.colorbar.Colorbar
        The input colorbar.

    Returns
    -------
    cbar : matplotlib.colorbar.Colorbar
        The modified colorbar.
    """
    pos = [0.0]
    cbar.set_ticks(
        pos,
        labels=["Zero within the 5th-95th percentile range"]
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

    Parameters
    ----------
    cbar : matplotlib.colorbar.Colorbar
        The input colorbar.
    vmin, vmax : float
        Min and max of the data range.

    Returns
    -------
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

    Parameters
    ----------
    cbar : matplotlib.colorbar.Colorbar
        The input colorbar.

    Returns
    -------
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

    Parameters
    ----------
    cbar : matplotlib.colorbar.Colorbar
        The input colorbar.
    vmin, vmax : float
        Min and max of the data range.
    diff_cmap : bool, optional
        Indicates that we are using a diverging colortable (True)
        or not (False).
        Default value: False

    Returns
    -------
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
