"""
Specific utilities for creating plots from GEOS-Chem benchmark simulations.
"""

import os
import shutil
import yaml
from yaml import load as yaml_load_file
import copy
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.geoaxes import GeoAxes  # for assertion
import matplotlib as mpl
import matplotlib.ticker as mticker
from matplotlib.backends.backend_pdf import PdfPages
from PyPDF2 import PdfFileWriter, PdfFileReader, PdfFileMerger
from .plot import WhGrYlRd
from .grid.horiz import make_grid_LL, make_grid_CS
from .grid.regrid import make_regridder_C2L, make_regridder_L2L
from .grid.gc_vertical import GEOS_72L_grid
from . import core
from .core import gcplot, call_make_grid, get_input_res, all_zero_or_nan, get_grid_extents
from .units import convert_units
from .constants import skip_these_vars
from joblib import Parallel, delayed, cpu_count, parallel_backend
from multiprocessing import current_process
import warnings

# Save warnings format to undo overwriting built into PyPDF2
warning_format = warnings.showwarning

# Suppress numpy divide by zero warnings to prevent output spam
np.seterr(divide="ignore", invalid="ignore")

# YAML files
aod_spc = "aod_species.yml"
spc_categories = "benchmark_categories.yml"
emission_spc = "emission_species.yml"
emission_inv = "emission_inventories.yml"

#class benchmark_plotter

def sixplot(
    subplot,
    all_zero,
    all_nan,
    plot_val,
    grid,
    ax,
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
    pedge=None,
    pedge_ind=0,
    log_yaxis=False,
    xtick_positions=[],
    xticklabels=[],
    plot_type="single_level"
):

    """
    Plotting function to be called from compare_single_level or compare_zonal_mean.
    Can also be called on its own
    WBD MOVE TO core.py and RENAME?
    Args:
    -----
    
    ^subplot : str
       Type of plot to create (ref, dev, absolute difference or fractional difference)
    
    ^all_zero : boolean
       Set this flag to True if the data to be plotted consist only of zeros
    
    ^all_nan : boolean
       Set this flag to True if the data to be plotted consist only of NaNs

    *plot_val : xarray DataArray
       Single variable GEOS-Chem output values to plot
    
    grid : dict
       Dictionary mapping plot_val to plottable coordinates

    WBD VVVVV
    ax : matplotlib axes 
       Axes object to plot information. Will create a new axes if none is passed.

    ^rowcol : tuple
       Subplot position in overall Figure WBD DELETE?

    title : str
       Title to print on axes

    comap : matplotlib Colormap
       Colormap for plotting data values    

    unit : str
       Units of plotted data
    
    extent : tuple (minlon, maxlon, minlat, maxlat) WBD SHOULD BE KEYWORD?
       Describes minimum and maximum latitude and longitude of input data 

    masked_data : numpy array
       Masked area for cubed-sphere plotting
    
    #Need to modify this name
    ^other_all_nan : boolean
        Set this flag to True if plotting ref/dev and the other of ref/dev is all nan

    gridtype : str
       "ll" for lat/lon or "cs" for cubed-sphere
    
    ^vmins: list of float
       list of length 3 of minimum ref value, dev value, and absdiff value

    ^vmaxs: list of float
       list of length 3 of maximum ref value, dev value, and absdiff value

    use_cmap_RdBu : boolean
       Set this flag to True to use a blue-white-red colormap 

    ^match_cbar : boolean
       Set this flag to True if you are plotting with the same colorbar for ref and dev

    verbose : boolean
       Set this flag to True to enable informative printout.    

    log_color_scale : boolean
       Set this flag to True to enable log-scale colormapping

    pedge : 
       Edge pressures of grid cells in data to be plotted

    pedge_ind : int
       Index of edge pressure values within pressure range  in data to be plotted 

    log_yaxis : boolean
       Set this flag to True to enable log scaling of pressure in zonal mean plots

    xtick_positions : list of float
       Locations of lat/lon or lon ticks on plot

    xtick_labels: list of str
       Labels for lat/lon ticks
    """

    # Set min and max of the data range
    if subplot in ("ref", "dev"):
        if all_zero or all_nan:
            if subplot is "ref":
                [vmin, vmax] = [vmins[0], vmaxs[0]]
            else:
                [vmin, vmax] = [vmins[1], vmaxs[1]]
        elif use_cmap_RdBu:
            if subplot is "ref":
                if match_cbar and (not other_all_nan):
                    absmax = max([np.abs(vmins[2]), np.abs(vmaxs[2])])
                else:
                    absmax = max([np.abs(vmins[0]), np.abs(vmaxs[0])])
            else:
                if match_cbar and (not other_all_nan):
                    absmax = max([np.abs(vmins[2]), np.abs(vmaxs[2])])
                else:
                    absmax = max([np.abs(vmins[1]), np.abs(vmaxs[1])])
        else:
            if subplot is "ref":
                if match_cbar and (not other_all_nan):
                    [vmin, vmax] = [vmins[2], vmaxs[2]]
                else:
                    [vmin, vmax] = [vmins[0], vmaxs[0]]
            else:
                if match_cbar and (not other_all_nan):
                    [vmin, vmax] = [vmins[2], vmaxs[2]]
                else:
                    [vmin, vmax] = [vmins[1], vmaxs[1]]
    else:
        if all_zero:
            [vmin, vmax] = [0, 0]
        elif all_nan:
            [vmin, vmax] = [np.nan, np.nan]
        else:
            if subplot is "dyn_abs_diff":
                # Min and max of abs. diff, excluding NaNs
                diffabsmax = max(
                    [np.abs(np.nanmin(plot_val)), np.abs(np.nanmax(plot_val))]
                )
                [vmin, vmax] = [-diffabsmax, diffabsmax]
            elif subplot is "res_abs_diff":
                [pct5, pct95] = [
                    np.percentile(plot_val, 5),
                    np.percentile(plot_val, 95),
                ]
                abspctmax = np.max([np.abs(pct5), np.abs(pct95)])
                [vmin, vmax] = [-abspctmax, abspctmax]
            elif subplot is "dyn_frac_diff":
                fracdiffabsmax = np.max(
                    [np.abs(np.nanmin(plot_val)), np.abs(np.nanmax(plot_val))]
                )
                [vmin, vmax] = [-fracdiffabsmax, fracdiffabsmax]
            else:
                [vmin, vmax] = [-2, 2]
    if verbose:
        print("Subplot ({}) vmin, vmax: {}, {}".format(rowcol, vmin, vmax))

    # Normalize colors (put into range [0..1] for matplotlib methods)
    if subplot in ("ref", "dev"):
        norm = core.normalize_colors(
            vmin, vmax, is_difference=use_cmap_RdBu,
            log_color_scale=log_color_scale
        )
    else:
        norm = core.normalize_colors(vmin, vmax, is_difference=True)
        
    #Create plot
    plot = gcplot(plot_val, ax, plot_type, grid, gridtype, title, comap,
                  norm, unit, extent, masked_data, use_cmap_RdBu, log_color_scale,
                  add_cb=False, pedge=pedge, pedge_ind=pedge_ind, log_yaxis=log_yaxis,
                  xtick_positions=xtick_positions, xticklabels=xticklabels)

    # Define the colorbar for the plot
    cb = plt.colorbar(plot, ax=ax, orientation="horizontal", pad=0.10)
    cb.mappable.set_norm(norm)
    if all_zero or all_nan:
        if subplot in ("ref", "dev"):
            if use_cmap_RdBu:
                cb.set_ticks([0.0])
            else:
                cb.set_ticks([0.5])
        else:
            cb.set_ticks([0.0])
        if all_nan:
            cb.set_ticklabels(["Undefined throughout domain"])
        else:
            cb.set_ticklabels(["Zero throughout domain"])
    else:
        if subplot in ("ref", "dev") and log_color_scale:
            cb.formatter = mticker.LogFormatter(base=10)
        else:
            if (vmax - vmin) < 0.1 or (vmax - vmin) > 100:
                cb.locator = mticker.MaxNLocator(nbins=4)

    cb.update_ticks()
    cb.set_label(unit)


def compare_single_level(
    refdata,
    refstr,
    devdata,
    devstr,
    varlist=None,
    ilev=0,
    itime=0,
    weightsdir=None,
    pdfname="",
    cmpres=None,
    match_cbar=True,
    normalize_by_area=False,
    enforce_units=True,
    flip_ref=False,
    flip_dev=False,
    use_cmap_RdBu=False,
    verbose=False,
    log_color_scale=False,
    extra_title_txt=None,
    plot_extent = [-1000, -1000, -1000, -1000],
    n_job=-1,
    sigdiff_list=[],
):
    """
    Create single-level 3x2 comparison map plots for variables common 
    in two xarray Datasets. Optionally save to PDF. 

    Args:
    -----
        refdata : xarray dataset
            Dataset used as reference in comparison

        refstr  : str
            String description for reference data to be used in plots
     
        devdata : xarray dataset
            Dataset used as development in comparison

        devstr  : str
            String description for development data to be used in plots
 
    Keyword Args (optional):
    ------------------------
        varlist : list of strings
            List of xarray dataset variable names to make plots for
            Default value: None (will compare all common variables)

        ilev : integer
            Dataset level dimension index using 0-based system
            Default value: 0   

        itime : integer
            Dataset time dimension index using 0-based system
            Default value: 0

        weightsdir : str
            Directory path for storing regridding weights
            Default value: None (will create/store weights in 
            current directory)

        pdfname : str
            File path to save plots as PDF
            Default value: Empty string (will not create PDF)

        cmpres : str
            String description of grid resolution at which 
            to compare datasets
            Default value: None (will compare at highest resolution 
            of ref and dev)

        match_cbar : boolean
            Set this flag to True if you wish to use the same colorbar
            bounds for the Ref and Dev plots.
            Default value: True

        normalize_by_area : boolean
            Set this flag to True if you wish to normalize the Ref and Dev
            raw data by grid area.
            Default value: False

        enforce_units : boolean
            Set this flag to True to force an error if Ref and Dev
            variables have different units.
            Default value: True

        flip_ref : boolean
            Set this flag to True to flip the vertical dimension of 
            3D variables in the Ref dataset.
            Default value: False

        flip_dev : boolean
            Set this flag to True to flip the vertical dimension of
            3D variables in the Dev dataset.
            Default value: False

        use_cmap_RdBu : boolean
            Set this flag to True to use a blue-white-red colormap 
            for plotting the raw data in both the Ref and Dev datasets.
            Default value: False

        verbose : boolean
            Set this flag to True to enable informative printout.
            Default value: False

        log_color_scale: boolean         
            Set this flag to True to plot data (not diffs)
            on a log color scale.
            Default value: False

        extra_title_txt : str
            Specifies extra text (e.g. a date string such as "Jan2016")
            for the top-of-plot title.
            Default value: None

        sigdiff_list: list of str
            Returns a list of all quantities having significant 
            differences (where |max(fractional difference)| > 0.1).
            Default value: []

    Returns:
    --------
        Nothing

    Example:
    --------
        >>> import matplotlib.pyplot as plt
        >>> import xarray as xr
        >>> from gcpy import benchmark
        >>> refds = xr.open_dataset('path/to/ref.nc4')
        >>> devds = xr.open_dataset('path/to/dev.nc4')
        >>> varlist = ['SpeciesConc_O3', 'SpeciesConc_NO']
        >>> benchmark.compare_single_level( refds, '12.3.2', devds, 'bug fix', varlist=varlist )
        >>> plt.show()
    """
    warnings.showwarning = warning_format
    # TODO: refactor this function and zonal mean plot function.
    # There is a lot of overlap and repeated code that could be abstracted.
    # Error check arguments
    if not isinstance(refdata, xr.Dataset):
        raise TypeError("The refdata argument must be an xarray Dataset!")

    if not isinstance(devdata, xr.Dataset):
        raise TypeError("The devdata argument must be an xarray Dataset!")

    # If no varlist is passed, plot all (surface only for 3D)
    if varlist == None:
        quiet = not verbose
        vardict = core.compare_varnames(refdata, devdata, quiet=quiet)
        varlist = vardict["commonvars3D"] + vardict["commonvars2D"]
        print("Plotting all common variables")
    n_var = len(varlist)

    # If no weightsdir is passed, set to current directory in case it is needed
    if weightsdir == None:
        weightsdir = "."

    # If no pdf name passed, then do not save to PDF
    savepdf = True
    if pdfname == "":
        savepdf = False
    # Cleanup previous temp PDFs
    for i in range(n_var):
        try:
            os.remove(pdfname + "BENCHMARKFIGCREATION.pdf" + str(i))
        except:
            continue
    # =================================================================
    # Determine input grid resolutions and types
    # =================================================================

    # GCC output and GCHP output using pre-v1.0.0 MAPL have lat and lon dims

    refres, refgridtype = get_input_res(refdata)
    devres, devgridtype = get_input_res(devdata)

    # =================================================================
    # Determine comparison grid resolution and type (if not passed)
    # =================================================================

    # If no cmpres is passed then choose highest resolution between ref and dev.
    # If one dataset is lat-lon and the other cubed sphere, and no
    # comparison grid resolution is passed, then default to 1x1.25
    if cmpres == None:
        if refres == devres and refgridtype == "ll":
            cmpres = refres
            cmpgridtype = refgridtype
        elif refgridtype == "ll" and devgridtype == "ll":
            cmpres = min([refres, devres])
            cmpgridtype = "ll"
        elif refgridtype == "cs" and devgridtype == "cs":
            # CS to CS regridding is not enabled yet, so default to 1x1.25
            # cmpres = max([refres, devres])
            # cmpgridtype = 'cs'
            cmpres = "1x1.25"
            cmpgridtype = "ll"
        else:
            cmpres = "1x1.25"
            cmpgridtype = "ll"
    elif "x" in cmpres:
        cmpgridtype = "ll"
    else:
        cmpgridtype = "cs"
        cmpres = int(cmpres)  # must cast to integer for cubed-sphere

    # Determine what, if any, need regridding.
    regridref = refres != cmpres
    regriddev = devres != cmpres
    regridany = regridref or regriddev

    # =================================================================
    # Get grid extents, only regridding to extents of each grid
    # =================================================================

    refminlon, refmaxlon, refminlat, refmaxlat = get_grid_extents(refdata)
    devminlon, devmaxlon, devminlat, devmaxlat = get_grid_extents(devdata)
    cmpminlon = min(x for x in [refminlon, devminlon] if x is not None)
    cmpmaxlon = max(x for x in [refmaxlon, devmaxlon] if x is not None)
    cmpminlat = min(x for x in [refminlat, devminlat] if x is not None)
    cmpmaxlat = max(x for x in [refmaxlat, devmaxlat] if x is not None)
    
    # =================================================================
    # Make grids (ref, dev, and comparison)
    # =================================================================

    [refgrid, regrid_list] = call_make_grid(refres, refgridtype, False, False,
                                            refminlon, refmaxlon, refminlat, refmaxlat)
    [devgrid, devgrid_list] = call_make_grid(devres, devgridtype, False, False,
                                             devminlon, devmaxlon, devminlat, devmaxlat)
    [cmpgrid, cmpgrid_list] = call_make_grid(cmpres, cmpgridtype, False, True,
                                             cmpminlon, cmpmaxlon, cmpminlat, cmpmaxlat)

    # =================================================================
    # Make regridders, if applicable
    # TODO: Make CS to CS regridders
    # =================================================================
    #print(refgrid)
    #print(len(refgrid['lat']), len(refgrid['lon']))
    if regridref:
        if refgridtype == "ll":
            refregridder = make_regridder_L2L(
                refres, cmpres, weightsdir=weightsdir, reuse_weights=True,
                #minlon=refminlon, maxlon=refmaxlon, minlat=refminlat, maxlat=refmaxlat
            )
        else:
            if cmpgridtype == "cs":
                print(
                    "ERROR: CS to CS regridding is not yet implemented in gcpy. Ref and dev cubed sphere grids must be the same resolution, or pass cmpres to compare_single_level as a lat-lon grid resolution."
                )
                return
            else:
                refregridder_list = make_regridder_C2L(
                    refres, cmpres, weightsdir=weightsdir, reuse_weights=True
                )
    if regriddev:
        if devgridtype == "ll":
            devregridder = make_regridder_L2L(
                devres, cmpres, weightsdir=weightsdir, reuse_weights=True,
                #minlon=devminlon, maxlon=devmaxlon, minlat=devminlat, maxlat=devmaxlat
            )
        else:
            if cmpgridtype == "cs":
                print(
                    "ERROR: CS to CS regridding is not yet implemented in gcpy. Ref and dev cubed sphere grids must be the same resolution, or pass cmpres to compare_single_level as a lat-lon grid resolution."
                )
                return
            else:
                print(devres, cmpres)
                devregridder_list = make_regridder_C2L(
                    devres, cmpres, weightsdir=weightsdir, reuse_weights=True
                )

    # =================================================================
    # Get lat/lon extents, if applicable
    # =================================================================
    refminlon, refmaxlon, refminlat, refmaxlat = get_grid_extents(refgrid)
    devminlon, devmaxlon, devminlat, devmaxlat = get_grid_extents(devgrid)
    cmpminlon, cmpmaxlon, cmpminlat, cmpmaxlat = get_grid_extents(cmpgrid)    

    ds_refs = [None] * n_var
    ds_devs = [None] * n_var
    for i in range(n_var):
        # ==============================================================
        # Slice the data, allowing for the
        # possibility of no time dimension (bpch)
        # ==============================================================
        varname = varlist[i]
        units_ref, units_dev = check_units(refdata, devdata, varname)
        # Ref
        vdims = refdata[varname].dims
        if "time" in vdims and "lev" in vdims:
            if flip_ref:
                ds_refs[i] = refdata[varname].isel(time=itime, lev=71 - ilev)
            else:
                ds_refs[i] = refdata[varname].isel(time=itime, lev=ilev)
        elif "time" not in vdims and "lev" in vdims:
            if flip_ref:
                ds_refs[i] = refdata[varname].isel(lev=71 - ilev)
            else:
                ds_refs[i] = refdata[varname].isel(lev=ilev)
        elif "time" in vdims and "lev" not in vdims:
            ds_refs[i] = refdata[varname].isel(time=itime)
        else:
            ds_refs[i] = refdata[varname]

        # Dev
        vdims = devdata[varname].dims
        if "time" in vdims and "lev" in vdims:
            if flip_dev:
                ds_devs[i] = devdata[varname].isel(time=itime, lev=71 - ilev)
            else:
                ds_devs[i] = devdata[varname].isel(time=itime, lev=ilev)
        elif "time" not in vdims and "lev" in vdims:
            if flip_dev:
                ds_devs[i] = devdata[varname].isel(lev=71 - ilev)
            else:
                ds_devs[i] = devdata[varname].isel(lev=ilev)
        elif "time" in vdims and "lev" not in vdims:
            ds_devs[i] = devdata[varname].isel(time=itime)
        else:
            ds_devs[i] = devdata[varname]

        # ==============================================================
        # Reshape cubed sphere data if using MAPL v1.0.0+
        # TODO: update function to expect data in this format
        # ==============================================================
        ds_refs[i] = reshape_MAPL_CS(ds_refs[i], refdata[varname].dims)
        ds_devs[i] = reshape_MAPL_CS(ds_devs[i], devdata[varname].dims)

    ds_ref_cmps = [None] * n_var
    ds_dev_cmps = [None] * n_var
    for i in range(n_var):
        ds_ref = ds_refs[i]
        ds_dev = ds_devs[i]
        # Ref
        if regridref:
            if refgridtype == "ll":
                # regrid ll to ll
                ds_ref_cmps[i] = refregridder(ds_ref)
            else:
                # regrid cs to ll
                ds_ref_cmps[i] = np.zeros([cmpgrid["lat"].size, cmpgrid["lon"].size])
                ds_ref_reshaped = ds_ref.data.reshape(6, refres, refres)
                for j in range(6):
                    regridder = refregridder_list[j]
                    ds_ref_cmps[i] += regridder(ds_ref_reshaped[j])
        else:
            ds_ref_cmps[i] = ds_ref

        # Dev
        if regriddev:
            if devgridtype == "ll":
                # regrid ll to ll
                ds_dev_cmps[i] = devregridder(ds_dev)
            else:
                # regrid cs to ll
                ds_dev_cmps[i] = np.zeros([cmpgrid["lat"].size, cmpgrid["lon"].size])
                ds_dev_reshaped = ds_dev.data.reshape(6, devres, devres)
                for j in range(6):
                    regridder = devregridder_list[j]
                    ds_dev_cmps[i] += regridder(ds_dev_reshaped[j])
        else:
            ds_dev_cmps[i] = ds_dev

    # =================================================================
    # Create pdf if saving to file
    # =================================================================

    if savepdf:
        print("Creating {} for {} variables".format(pdfname, n_var))
        pdf = PdfPages(pdfname)
        pdf.close()

    # =================================================================
    # Loop over variables
    # =================================================================

    print_units_warning = True

    # This loop is written as a function so it can be called in parallel
    def createfig(ivar):

        # Suppress harmless run-time warnings (mostly about underflow)
        warnings.filterwarnings('ignore', category=RuntimeWarning)
        warnings.filterwarnings('ignore', category=UserWarning)

        if savepdf:
            print("{} ".format(ivar), end="")
        varname = varlist[ivar]
        varndim_ref = refdata[varname].ndim
        varndim_dev = devdata[varname].ndim

        # Convert mol/mol units to ppb and ensure ref and dev units match
        units_ref, units_dev = check_units(refdata, devdata, varname)

        ds_ref = ds_refs[ivar]
        ds_dev = ds_devs[ivar]

        # ==============================================================
        # Area normalization, if any
        # ==============================================================

        # if normalizing by area, adjust units to be per m2,
        # and adjust title string
        units = units_ref
        subtitle_extra = ""
        varndim = varndim_ref

        # if regridding then normalization by area may be necessary
        # depending on units. Either pass normalize_by_area=True to
        # normalize all, or include units that should always be normalized
        # by area below. GEOS-Chem Classic output files include area and so
        # do not need to be passed.
        exclude_list = ["WetLossConvFrac", "Prod_", "Loss_"]
        if regridany and (units in ("kg", "kgC") or normalize_by_area):
            if not any(s in varname for s in exclude_list):
                if (
                    "AREAM2" in refdata.data_vars.keys()
                    and "AREAM2" in devdata.data_vars.keys()
                ):
                    ds_ref.values = ds_ref.values / refdata["AREAM2"].values
                    ds_dev.values = ds_dev.values / devdata["AREAM2"].values
                else:
                    print(
                        "ERROR: Variables AREAM2 needed for area normalization missing from dataset"
                    )
                    return
                units = "{}/m2".format(units)
                units_ref = units
                units_dev = units
                subtitle_extra = ", Normalized by Area"

        # ==============================================================
        # Get comparison data sets, regridding input slices if needed
        # ==============================================================

        # Reshape ref/dev cubed sphere data, if any
        ds_ref_reshaped = None
        if refgridtype == "cs":
            ds_ref_reshaped = ds_ref.data.reshape(6, refres, refres)
        ds_dev_reshaped = None
        if devgridtype == "cs":
            ds_dev_reshaped = ds_dev.data.reshape(6, devres, devres)

        ds_ref_cmp = ds_ref_cmps[ivar]
        ds_dev_cmp = ds_dev_cmps[ivar]

        # Reshape comparison cubed sphere data, if any
        if cmpgridtype == "cs":
            ds_ref_cmp_reshaped = ds_ref_cmp.data.reshape(6, cmpres, cmpres)
            ds_dev_cmp_reshaped = ds_dev_cmp.data.reshape(6, cmpres, cmpres)

        # ==============================================================
        # Get min and max values for use in the colorbars
        # ==============================================================

        # WBD MOVE TO FUNCTION

        # Ref
        vmin_ref = float(ds_ref.data.min())
        vmax_ref = float(ds_ref.data.max())

        # Dev
        vmin_dev = float(ds_dev.data.min())
        vmax_dev = float(ds_dev.data.max())

        # Comparison
        if cmpgridtype == "cs":
            vmin_ref_cmp = float(ds_ref_cmp.data.min())
            vmax_ref_cmp = float(ds_ref_cmp.data.max())
            vmin_dev_cmp = float(ds_dev_cmp.data.min())
            vmax_dev_cmp = float(ds_dev_cmp.data.max())
            vmin_cmp = np.min([vmin_ref_cmp, vmin_dev_cmp])
            vmax_cmp = np.max([vmax_ref_cmp, vmax_dev_cmp])
        else:
            vmin_cmp = np.min([ds_ref_cmp.min(), ds_dev_cmp.min()])
            vmax_cmp = np.max([ds_ref_cmp.max(), ds_dev_cmp.max()])

        # Get overall min & max
        vmin_abs = np.min([vmin_ref, vmin_dev, vmin_cmp])
        vmax_abs = np.max([vmax_ref, vmax_dev, vmax_cmp])
        if match_cbar:
            [vmin, vmax] = [vmin_abs, vmax_abs]

        if verbose:
            print("vmin_ref: {}".format(vmin_ref))
            print("vmax_ref: {}".format(vmax_ref))
            print("vmin_dev: {}".format(vmin_dev))
            print("vmax_dev: {}".format(vmax_dev))
            if cmpgridtype == "cs":
                print("vmin_ref_cmp: {}".format(vmin_ref_cmp))
                print("vmax_ref_cmp: {}".format(vmax_ref_cmp))
                print("vmin_dev_cmp: {}".format(vmin_dev_cmp))
                print("vmax_dev_cmp: {}".format(vmax_dev_cmp))
            print("vmin_cmp: {}".format(vmin_cmp))
            print("vmax_cmp: {}".format(vmax_cmp))
            print("vmin_abs: {}".format(vmin_abs))
            print("vmax_abs: {}".format(vmax_abs))

        # ==============================================================
        # Test if Ref and/or Dev contain all zeroes or all NaNs.
        # This will have implications as to how we set min and max
        # values for the color ranges below.
        # ==============================================================

        ref_is_all_zero, ref_is_all_nan = all_zero_or_nan(ds_ref.values)
        dev_is_all_zero, dev_is_all_nan = all_zero_or_nan(ds_dev.values)

        # ==============================================================
        # Calculate absolute difference
        # ==============================================================
        if cmpgridtype == "ll":
            absdiff = np.array(ds_dev_cmp) - np.array(ds_ref_cmp)
        else:
            absdiff = ds_dev_cmp_reshaped - ds_ref_cmp_reshaped

        # Test if the abs. diff. is zero everywhere or NaN everywhere
        absdiff_is_all_zero, absdiff_is_all_nan = all_zero_or_nan(absdiff)

        # For cubed-sphere, take special care to avoid a spurious
        # boundary line, as described here: https://stackoverflow.com/questions/46527456/preventing-spurious-horizontal-lines-for-ungridded-pcolormesh-data
        if cmpgridtype == "cs":
            absdiff = np.ma.masked_where(np.abs(cmpgrid["lon"] - 180) < 2, absdiff)

        # ==============================================================
        # Calculate fractional difference, set divides by zero to NaN
        # ==============================================================
        if cmpgridtype == "ll":
            fracdiff = (np.array(ds_dev_cmp) - np.array(ds_ref_cmp)) / np.array(
                ds_ref_cmp
            )
        else:
            fracdiff = (ds_dev_cmp_reshaped - ds_ref_cmp_reshaped) / ds_ref_cmp_reshaped

        # Replace Infinity values with NaN
        fracdiff = np.where(fracdiff == np.inf, np.nan, fracdiff)

        # Test if the frac. diff. is zero everywhere or NaN everywhere
        fracdiff_is_all_zero = absdiff_is_all_zero or (
            not np.any(fracdiff) and not ref_is_all_zero
        )
        fracdiff_is_all_nan = np.isnan(fracdiff).all() or ref_is_all_zero

        # Absolute max value of fracdiff, excluding NaNs
        fracdiffabsmax = max([np.abs(np.nanmin(fracdiff)), np.abs(np.nanmax(fracdiff))])

        # For cubed-sphere, take special care to avoid a spurious
        # boundary line, as described here: https://stackoverflow.com/questions/46527456/preventing-spurious-horizontal-lines-for-ungridded-pcolormesh-data
        if cmpgridtype == "cs":
            fracdiff = np.ma.masked_where(np.abs(cmpgrid["lon"] - 180) < 2, fracdiff)

        # ==============================================================
        # Create 3x2 figure
        # ==============================================================

        # Create figures and axes objects
        # Also define the map projection that will be shown
        figs, ((ax0, ax1), (ax2, ax3), (ax4, ax5)) = plt.subplots(
            3, 2, figsize=[12, 14], subplot_kw={"projection": ccrs.PlateCarree()}
        )
        # Ensure subplots don't overlap when invoking plt.show()
        if not savepdf:
            plt.subplots_adjust(hspace=0.4)
        # Give the figure a title
        offset = 0.96
        fontsize = 25
        if "lev" in ds_ref.dims and "lev" in ds_dev.dims:
            if ilev == 0:
                levstr = "Surface"
            elif ilev == 22:
                levstr = "500 hPa"
            else:
                levstr = "Level " + str(ilev - 1)
            if extra_title_txt is not None:
                figs.suptitle(
                    "{}, {} ({})".format(varname, levstr, extra_title_txt),
                    fontsize=fontsize,
                    y=offset,
                )
            else:
                figs.suptitle(
                    "{}, {}".format(varname, levstr), fontsize=fontsize, y=offset
                )
        elif (
            "lat" in ds_ref.dims
            and "lat" in ds_dev.dims
            and "lon" in ds_ref.dims
            and "lon" in ds_dev.dims
        ):
            if extra_title_txt is not None:
                figs.suptitle(
                    "{} ({})".format(varname, extra_title_txt),
                    fontsize=fontsize,
                    y=offset,
                )
            else:
                figs.suptitle("{}".format(varname), fontsize=fontsize, y=offset)
        else:
            print("Incorrect dimensions for {}!".format(varname))

        # ==============================================================
        # Set colormaps for data plots
        #
        # Use shallow copy (copy.copy() to create color map objects,
        # in order to avoid set_bad() from being applied to the base
        # color table. See: https://docs.python.org/3/library/copy.html
        # ==============================================================

        # Colormaps for 1st row (Ref and Dev)
        if use_cmap_RdBu:
            cmap_toprow_nongray = copy.copy(mpl.cm.RdBu_r)
            cmap_toprow_gray = copy.copy(mpl.cm.RdBu_r)
        else:
            cmap_toprow_nongray = copy.copy(WhGrYlRd)
            cmap_toprow_gray = copy.copy(WhGrYlRd)
        cmap_toprow_gray.set_bad(color="gray")

        if refgridtype == "ll":
            if ref_is_all_nan:
                ref_cmap = cmap_toprow_gray
            else:
                ref_cmap = cmap_toprow_nongray

            if dev_is_all_nan:
                dev_cmap = cmap_toprow_gray
            else:
                dev_cmap = cmap_toprow_nongray

        # Colormaps for 2nd row (Abs. Diff.) and 3rd row (Frac. Diff,)
        cmap_nongray = copy.copy(mpl.cm.RdBu_r)
        cmap_gray = copy.copy(mpl.cm.RdBu_r)
        cmap_gray.set_bad(color="gray")

        # ==============================================================
        # Set plot bounds for non cubed-sphere plotting
        # ==============================================================
        if cmpgridtype == "ll":
            ref_extent = (refminlon, refmaxlon, refminlat, refmaxlat)
            dev_extent = (devminlon, devmaxlon, devminlat, devmaxlat)
            cmp_extent = (cmpminlon, cmpmaxlon, cmpminlat, cmpmaxlat)

        # ==============================================================
        # Set titles for plots
        # ==============================================================

        if refgridtype == "ll":
            ref_title = "{} (Ref){}\n{}".format(refstr, subtitle_extra, refres)
            dev_title = "{} (Dev){}\n{}".format(devstr, subtitle_extra, devres)
        else:
            ref_title = "{} (Ref){}\nc{}".format(refstr, subtitle_extra, refres)
            dev_title = "{} (Dev){}\nc{}".format(devstr, subtitle_extra, devres)

        if regridany:
            absdiff_dynam_title = "Difference ({})\nDev - Ref, Dynamic Range".format(
                cmpres
            )
            absdiff_fixed_title = "Difference ({})\nDev - Ref, Restricted Range [5%,95%]".format(
                cmpres
            )
            fracdiff_dynam_title = "Fractional Difference ({})\n(Dev-Ref)/Ref, Dynamic Range".format(
                cmpres
            )
            fracdiff_fixed_title = "Fractional Difference ({})\n(Dev-Ref)/Ref, Fixed Range".format(
                cmpres
            )
        else:
            absdiff_dynam_title = "Difference\nDev - Ref, Dynamic Range"
            absdiff_fixed_title = "Difference\nDev - Ref, Restricted Range [5%,95%]"
            fracdiff_dynam_title = "Fractional Difference\n(Dev-Ref)/Ref, Dynamic Range"
            fracdiff_fixed_title = "Fractional Difference\n(Dev-Ref)/Ref, Fixed Range"

        # ==============================================================
        # Bundle variables for 6 parallel plotting calls
        # 0 = Ref                 1 = Dev
        # 2 = Dynamic abs diff    3 = Restricted abs diff
        # 4 = Dynamic frac diff   5 = Restricted frac diff
        # ==============================================================

        subplots = [
            "ref",
            "dev",
            "dyn_abs_diff",
            "res_abs_diff",
            "dyn_frac_diff",
            "frac_abs_diff",
        ]

        all_zeros = [
            ref_is_all_zero,
            dev_is_all_zero,
            absdiff_is_all_zero,
            absdiff_is_all_zero,
            fracdiff_is_all_zero,
            fracdiff_is_all_zero,
        ]

        all_nans = [
            ref_is_all_nan,
            dev_is_all_nan,
            absdiff_is_all_nan,
            absdiff_is_all_nan,
            fracdiff_is_all_nan,
            fracdiff_is_all_nan,
        ]
        if not -1000 in plot_extent:
            extents = [plot_extent, plot_extent,
                       plot_extent, plot_extent,
                       plot_extent, plot_extent]
        else:
            extents = [ref_extent, dev_extent, 
                       cmp_extent, cmp_extent, 
                       cmp_extent, cmp_extent]

        plot_vals = [ds_ref, ds_dev, absdiff, absdiff, fracdiff, fracdiff]
        grids = [refgrid, devgrid, cmpgrid, cmpgrid, cmpgrid, cmpgrid]
        axs = [ax0, ax1, ax2, ax3, ax4, ax5]
        rowcols = [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1)]
        titles = [
            ref_title,
            dev_title,
            absdiff_dynam_title,
            absdiff_fixed_title,
            fracdiff_dynam_title,
            fracdiff_fixed_title,
        ]

        if refgridtype == "ll":
            cmaps = [ref_cmap, dev_cmap, cmap_gray, cmap_gray, cmap_gray, cmap_gray]
        else:
            cmaps = [
                cmap_toprow_nongray,
                cmap_toprow_nongray,
                cmap_nongray,
                cmap_nongray,
                cmap_nongray,
                cmap_nongray,
            ]

        ref_masked = None
        dev_masked = None
        if refgridtype == "cs":
            ref_masked = np.ma.masked_where(
                np.abs(refgrid["lon"] - 180) < 2, ds_ref_reshaped
            )
        if devgridtype == "cs":
            dev_masked = np.ma.masked_where(
                np.abs(devgrid["lon"] - 180) < 2, ds_dev_reshaped
            )
        masked = [ref_masked, dev_masked, absdiff, absdiff, fracdiff, fracdiff]

        gridtypes = [
            refgridtype,
            devgridtype,
            cmpgridtype,
            cmpgridtype,
            cmpgridtype,
            cmpgridtype,
        ]

        unit_list = [units_ref, units_dev, units, units, "unitless", "unitless"]

        other_all_nans = [dev_is_all_nan, ref_is_all_nan, False, False, False, False]

        mins = [vmin_ref, vmin_dev, vmin_abs]
        maxs = [vmax_ref, vmax_dev, vmax_abs]

        # Plot
        for i in range(6):
            sixplot(
                subplots[i],
                all_zeros[i],
                all_nans[i],
                plot_vals[i],
                grids[i],
                axs[i],
                rowcols[i],
                titles[i],
                cmaps[i],
                unit_list[i],
                extents[i],
                masked[i],
                other_all_nans[i],
                gridtypes[i],
                mins,
                maxs,
                use_cmap_RdBu,
                match_cbar,
                verbose,
                log_color_scale,
                plot_type="single_level"
            )

        # ==============================================================
        # Update the list of variables with significant differences.
        # Criterion: abs(max(fracdiff)) > 0.1
        # Do not include NaNs in the criterion, because these indicate
        # places where fracdiff could not be computed (div-by-zero).
        # ==============================================================
        if np.abs(np.nanmax(fracdiff)) > 0.1:
            sigdiff_list.append(varname)

        # ==============================================================
        # Add this page of 6-panel plots to a PDF file
        # ==============================================================
        if savepdf:
            pdf = PdfPages(pdfname + "BENCHMARKFIGCREATION.pdf" + str(ivar))
            pdf.savefig(figs)
            pdf.close()
            plt.close(figs)

    #do not attempt nested thread parallelization due to issues with matplotlib
    if current_process().name != "MainProcess":
        n_job = 1
    if savepdf:
        Parallel(n_jobs = n_job) (delayed(createfig)(i) for i in range(n_var))
    else:
        #disable parallel plotting to allow interactive figure plotting
        for i in range(n_var):
            createfig(i)

    # ==================================================================
    # Finish
    # ==================================================================
    if savepdf:
        print("Closed PDF")
        merge = PdfFileMerger()
        for i in range(n_var):
            merge.append(pdfname + "BENCHMARKFIGCREATION.pdf" + str(i))
            os.remove(pdfname + "BENCHMARKFIGCREATION.pdf" + str(i))
        merge.write(pdfname)
        merge.close()
        warnings.showwarning = warning_format


def compare_zonal_mean(
    refdata,
    refstr,
    devdata,
    devstr,
    varlist=None,
    itime=0,
    weightsdir=None,
    pdfname="",
    cmpres=None,
    match_cbar=True,
    pres_range=[0, 2000],
    normalize_by_area=False,
    enforce_units=True,
    flip_ref=False,
    flip_dev=False,
    use_cmap_RdBu=False,
    verbose=False,
    log_color_scale=False,
    log_yaxis=False,
    extra_title_txt=None,
    n_job=-1,
    sigdiff_list=[],
):

    """
    Create single-level 3x2 comparison zonal-mean plots for variables
    common in two xarray Daatasets. Optionally save to PDF. 

    Args:
    -----
        refdata : xarray dataset
            Dataset used as reference in comparison

        refstr  : str
            String description for reference data to be used in plots
     
        devdata : xarray dataset
            Dataset used as development in comparison

        devstr  : str
            String description for development data to be used in plots
 
    Keyword Args (optional):
    ------------------------
        varlist : list of strings
            List of xarray dataset variable names to make plots for
            Default value: None (will compare all common 3D variables)

        itime : integer
            Dataset time dimension index using 0-based system
            Default value: 0

        weightsdir : str
            Directory path for storing regridding weights
            Default value: None (will create/store weights in 
            current directory)

        pdfname : str
            File path to save plots as PDF
            Default value: Empty string (will not create PDF)

        cmpres : str
            String description of grid resolution at which
            to compare datasets
            Default value: None (will compare at highest resolution 
            of Ref and Dev)

        match_cbar : boolean
            Set this flag to True to use same the colorbar bounds
            for both Ref and Dev plots.
            Default value: True

        pres_range : list of two integers
            Pressure range of levels to plot [hPa]. The vertical axis will
            span the outer pressure edges of levels that contain pres_range 
            endpoints.
            Default value: [0,2000]

        normalize_by_area : boolean
            Set this flag to True to to normalize raw data in both
            Ref and Dev datasets by grid area.
            Default value: False

        enforce_units : boolean
            Set this flag to True force an error if the variables in
            the Ref and Dev datasets have different units.
            Default value: True

        flip_ref : boolean
            Set this flag to True to flip the vertical dimension of
            3D variables in the Ref dataset.
            Default value: False

        flip_dev : boolean
            Set this flag to True to flip the vertical dimension of
            3D variables in the Dev dataset.
            Default value: False

        use_cmap_RdBu : boolean
            Set this flag to True to use a blue-white-red colormap for
            plotting raw reference and development datasets.
            Default value: False

        verbose : logical
            Set this flag to True to enable informative printout.
            Default value: False

        log_color_scale: boolean         
            Set this flag to True to enable plotting data (not diffs)
            on a log color scale.
            Default value: False

        log_yaxis : boolean
            Set this flag to True if you wish to create zonal mean
            plots with a log-pressure Y-axis.
            Default value: False

        extra_title_txt : str
            Specifies extra text (e.g. a date string such as "Jan2016")
            for the top-of-plot title.
            Default value: None

        sigdiff_list: list of str
            Returns a list of all quantities having significant 
            differences (where |max(fractional difference)| > 0.1).
            Default value: []

    Returns:
    --------
        Nothing

    Example:
    --------
        >>> import matplotlib.pyplot as plt
        >>> import xarray as xr
        >>> from gcpy import benchmark
        >>> refds = xr.open_dataset('path/to/ref.nc4')
        >>> devds = xr.open_dataset('path/to/dev.nc4')
        >>> varlist = ['SpeciesConc_O3', 'SpeciesConc_NO']
        >>> benchmark.compare_zonal_mean( refds, '12.3.2', devds, 'bug fix', varlist=varlist )
        >>> plt.show()
    """

    # TODO: refactor this function and single level plot function. There is a lot of overlap and
    # repeated code that could be abstracted.
    warnings.showwarning = warning_format
    if not isinstance(refdata, xr.Dataset):
        raise TypeError("The refdata argument must be an xarray Dataset!")

    if not isinstance(devdata, xr.Dataset):
        raise TypeError("The devdata argument must be an xarray Dataset!")

    # If no varlist is passed, plot all 3D variables in the dataset
    if varlist == None:
        quiet = not verbose
        vardict = core.compare_varnames(refdata, devdata, quiet=quiet)
        varlist = vardict["commonvars3D"]
        print("Plotting all 3D variables")
    n_var = len(varlist)

    # Exit out if there are no 3D variables
    if not n_var:
        print("WARNING: no 3D variables to plot zonal mean for!")
        return

    # If no weightsdir is passed, set to current directory in case it is needed
    if weightsdir == None:
        weightsdir = "."

    # If no pdf name passed, then do not save to PDF
    savepdf = True
    if pdfname == "":
        savepdf = False
    # Cleanup previous temp PDFs
    for i in range(n_var):
        try:
            os.remove(pdfname + "BENCHMARKFIGCREATION.pdf" + str(i))
        except:
            continue
    # Get mid-point pressure and edge pressures for this grid (assume 72-level)
    pmid = GEOS_72L_grid.p_mid()
    pedge = GEOS_72L_grid.p_edge()

    # Get indexes of pressure subrange (full range is default)
    pedge_ind = np.where((pedge <= np.max(pres_range)) & (pedge >= np.min(pres_range)))
    pedge_ind = pedge_ind[0]
    # Pad edges if subset does not include surface or TOA so data spans entire subrange
    if min(pedge_ind) != 0:
        pedge_ind = np.append(min(pedge_ind) - 1, pedge_ind)
    if max(pedge_ind) != 72:
        pedge_ind = np.append(pedge_ind, max(pedge_ind) + 1)
    # pmid indexes do not include last pedge index
    pmid_ind = pedge_ind[:-1]
    nlev = len(pmid_ind)

    # Convert levels to pressures in ref and dev data
    if refdata.sizes["lev"] == 72:
        refdata["lev"] = pmid
    elif refdata.sizes["lev"] == 73:
        refdata["lev"] = pedge
    else:
        print(
            "ERROR: compare_zonal_mean implemented for 72 or 73 levels only. Other values found in ref."
        )
        return
    refdata["lev"].attrs["units"] = "hPa"
    refdata["lev"].attrs["long_name"] = "level pressure"

    if devdata.sizes["lev"] == 72:
        devdata["lev"] = pmid
    elif devdata.sizes["lev"] == 73:
        devdata["lev"] = pedge
    else:
        print(
            "ERROR: compare_zonal_mean implemented for 72 or 73 levels only. Other value found in dev."
        )
        return
    devdata["lev"].attrs["units"] = "hPa"
    devdata["lev"].attrs["long_name"] = "level pressure"

    # ==================================================================
    # Reduce pressure range if reduced range passed as input. Indices
    # must be flipped if flipping vertical axis.
    # ==================================================================

    pmid_ind_ref = pmid_ind
    pmid_ind_dev = pmid_ind
    pmid_ind_flipped = 72 - pmid_ind[::-1] - 1
    if flip_ref:
        pmid_ind_ref = pmid_ind_flipped
    if flip_dev:
        pmid_ind_dev = pmid_ind_flipped

    refdata = refdata.isel(lev=pmid_ind_ref)
    devdata = devdata.isel(lev=pmid_ind_dev)

    # ==================================================================
    # Determine input grid resolutions and types
    #
    # GCC output and GCHP output using pre-v1.0.0 MAPL have lat and
    # lon dims.  GCHP output in v1.0.0 MAPL has XDim and YDim instead.
    # ==================================================================

    refres, refgridtype = get_input_res(refdata)
    devres, devgridtype = get_input_res(devdata)

    # ==================================================================
    # Determine comparison grid resolution (if not passed)
    # ==================================================================

    # If no cmpres is passed then choose highest resolution between ref and dev.
    # If both datasets are cubed sphere then default to 1x1.25 for comparison.
    # If cmpres pass as cubed-sphere, over-ride to be 1x1.25 lat-lon with
    # a warning.
    cmpgridtype = "ll"
    if cmpres == None:
        if refres == devres and refgridtype == "ll":
            cmpres = refres
        elif refgridtype == "ll" and devgridtype == "ll":
            cmpres = min([refres, devres])
        else:
            cmpres = "1x1.25"
    elif "x" not in cmpres:
        print(
            "WARNING: zonal mean comparison grid must be lat-lon. Defaulting to 1x1.25"
        )
        cmpres = "1x1.25"

    # Determine what, if any, need regridding.
    regridref = refres != cmpres
    regriddev = devres != cmpres
    regridany = regridref or regriddev

    # ==================================================================
    # Make grids (ref, dev, and comparison)
    # ==================================================================

    [refgrid, regrid_list] = call_make_grid(refres, refgridtype, True, False)
    [devgrid, devgrid_list] = call_make_grid(devres, devgridtype, True, False)
    [cmpgrid, cmpgrid_list] = call_make_grid(cmpres, cmpgridtype, True, True)

    # ==================================================================
    # Make regridders, if applicable
    # TODO: Add CS to CS regridders
    # ==================================================================

    if regridref:
        if refgridtype == "ll":
            refregridder = make_regridder_L2L(
                refres, cmpres, weightsdir=weightsdir, reuse_weights=True
            )
        else:
            refregridder_list = make_regridder_C2L(
                refres, cmpres, weightsdir=weightsdir, reuse_weights=True
            )
    if regriddev:
        if devgridtype == "ll":
            devregridder = make_regridder_L2L(
                devres, cmpres, weightsdir=weightsdir, reuse_weights=True
            )
        else:
            devregridder_list = make_regridder_C2L(
                devres, cmpres, weightsdir=weightsdir, reuse_weights=True
            )

    ds_refs = [None] * n_var
    ds_devs = [None] * n_var
    for i in range(n_var):
        varname = varlist[i]
        units_ref, units_dev = check_units(refdata, devdata, varname)        

        # ==============================================================
        # Slice the data, allowing for the
        # possibility of no time dimension (bpch)
        # ==============================================================

        # Ref
        vdims = refdata[varname].dims
        if "time" in vdims:
            ds_refs[i] = refdata[varname].isel(time=itime)
        else:
            ds_refs[i] = refdata[varname]

        # Dev
        vdims = devdata[varname].dims
        if "time" in vdims:
            ds_devs[i] = devdata[varname].isel(time=itime)
        else:
            ds_devs[i] = devdata[varname]

        # ==============================================================
        # Reshape cubed sphere data if using MAPL v1.0.0+
        # TODO: update function to expect data in this format
        # ==============================================================
        ds_refs[i] = reshape_MAPL_CS(ds_refs[i], refdata[varname].dims)
        ds_devs[i] = reshape_MAPL_CS(ds_devs[i], devdata[varname].dims)

        # Flip in the vertical if applicable
        if flip_ref:
            ds_refs[i].data = ds_refs[i].data[::-1, :, :]
        if flip_dev:
            ds_devs[i].data = ds_devs[i].data[::-1, :, :]

    ds_ref_cmps = [None] * n_var
    ds_dev_cmps = [None] * n_var

    for i in range(n_var):
        # ==============================================================
        # Get comparison data sets, regridding input slices if needed
        # ==============================================================

        ds_ref = ds_refs[i]
        ds_dev = ds_devs[i]

        # Ref
        if regridref:
            if refgridtype == "ll":
                # regrid ll to ll
                ds_ref_cmps[i] = refregridder(ds_ref)
            else:
                # regrid cs to ll
                ds_ref_reshaped = ds_ref.data.reshape(nlev, 6, refres, refres).swapaxes(
                    0, 1
                )
                ds_ref_cmps[i] = np.zeros(
                    [nlev, cmpgrid["lat"].size, cmpgrid["lon"].size]
                )
                for j in range(6):
                    regridder = refregridder_list[j]
                    ds_ref_cmps[i] += regridder(ds_ref_reshaped[j])
        else:
            ds_ref_cmps[i] = ds_ref

        # Dev
        if regriddev:
            if devgridtype == "ll":
                # regrid ll to ll
                ds_dev_cmps[i] = devregridder(ds_dev)
            else:
                # regrid cs to ll
                ds_dev_reshaped = ds_dev.data.reshape(nlev, 6, devres, devres).swapaxes(
                    0, 1
                )
                ds_dev_cmps[i] = np.zeros(
                    [nlev, cmpgrid["lat"].size, cmpgrid["lon"].size]
                )
                for j in range(6):
                    regridder = devregridder_list[j]
                    ds_dev_cmps[i] += regridder(ds_dev_reshaped[j])
        else:
            ds_dev_cmps[i] = ds_dev

    # ==================================================================
    # Create pdf, if savepdf is passed as True
    # ==================================================================

    # Universal plot setup
    xtick_positions = np.arange(-90, 91, 30)
    xticklabels = ["{}$\degree$".format(x) for x in xtick_positions]

    if savepdf:
        print("Creating {} for {} variables".format(pdfname, n_var))
        pdf = PdfPages(pdfname)
        pdf.close()

    # ==================================================================
    # Loop over variables
    # ==================================================================

    # Loop over variables
    print_units_warning = True

    # This loop is written as a function so it can be called in parallel
    def createfig(ivar):

        # Suppress harmless run-time warnings (mostly about underflow)
        warnings.filterwarnings('ignore', category=RuntimeWarning)
        warnings.filterwarnings('ignore', category=UserWarning)

        if savepdf:
            print("{} ".format(ivar), end="")
        varname = varlist[ivar]
        varndim_ref = refdata[varname].ndim
        varndim_dev = devdata[varname].ndim

        # Convert mol/mol units to ppb and ensure ref and dev units match
        units_ref, units_dev = check_units(refdata, devdata, varname)

        # ==============================================================
        # Area normalization, if any
        # ==============================================================

        # if normalizing by area, adjust units to be per m2, and
        # adjust title string
        units = units_ref
        varndim = varndim_ref
        subtitle_extra = ""

        # if normalizing by area, transform on the native grid and adjust units and subtitle string
        exclude_list = ["WetLossConvFrac", "Prod_", "Loss_"]
        if normalize_by_area and not any(s in varname for s in exclude_list):
            ds_ref.values = ds_ref.values / refdata["AREAM2"].values[np.newaxis, :, :]
            ds_dev.values = ds_dev.values / devdata["AREAM2"].values[np.newaxis, :, :]
            units = "{} m-2".format(units)
            subtitle_extra = ", Normalized by Area"

        # ==============================================================
        # Assign data variables
        # ==============================================================
        ds_ref = ds_refs[ivar]
        ds_dev = ds_devs[ivar]
        ds_ref_cmp = ds_ref_cmps[ivar]
        ds_dev_cmp = ds_dev_cmps[ivar]

        # ==============================================================
        # Calculate zonal mean
        # ==============================================================

        # Ref
        if refgridtype == "ll":
            zm_ref = ds_ref.mean(dim="lon")
        else:
            zm_ref = ds_ref_cmp.mean(axis=2)

        # Dev
        if devgridtype == "ll":
            zm_dev = ds_dev.mean(dim="lon")
        else:
            zm_dev = ds_dev_cmp.mean(axis=2)

        # Comparison
        zm_dev_cmp = ds_dev_cmp.mean(axis=2)
        zm_ref_cmp = ds_ref_cmp.mean(axis=2)

        # ==============================================================
        # Get min and max values for use in the colorbars
        # and also flag if Ref and/or Dev are all zero or all NaN
        # ==============================================================

        # Ref
        vmin_ref = float(zm_ref.min())
        vmax_ref = float(zm_ref.max())

        # Dev
        vmin_dev = float(zm_dev.min())
        vmax_dev = float(zm_dev.max())

        # Comparison
        vmin_cmp = np.min([zm_ref_cmp.min(), zm_dev_cmp.min()])
        vmax_cmp = np.max([zm_ref_cmp.max(), zm_dev_cmp.max()])

        # Take min/max across all grids
        vmin_abs = np.min([vmin_ref, vmin_dev, vmin_cmp])
        vmax_abs = np.max([vmax_ref, vmax_dev, vmax_cmp])

        if verbose:
            print("vmin_ref: {}".format(vmin_ref))
            print("vmin_dev: {}".format(vmin_dev))
            print("vmin_cmp: {}".format(vmin_cmp))
            print("vmin_abs: {}".format(vmin_abs))
            print("vmax_ref: {}".format(vmax_ref))
            print("vmax_dev: {}".format(vmax_dev))
            print("vmax_cmp: {}".format(vmax_cmp))
            print("vmax_abs: {}".format(vmax_abs))

        # ==============================================================
        # Test if Ref and/or Dev contain all zeroes or all NaNs.
        # This will have implications as to how we set min and max
        # values for the color ranges below.
        # ==============================================================

        ref_is_all_zero, ref_is_all_nan = all_zero_or_nan(ds_ref.values)
        dev_is_all_zero, dev_is_all_nan = all_zero_or_nan(ds_dev.values)

        # ==============================================================
        # Calculate zonal mean difference
        # ==============================================================

        zm_diff = np.array(zm_dev_cmp) - np.array(zm_ref_cmp)

        # Test if abs. diff is zero everywhere or NaN everywhere
        absdiff_is_all_zero, absdiff_is_all_nan = all_zero_or_nan(zm_diff)

        # Absolute maximum difference value
        diffabsmax = max([np.abs(zm_diff.min()), np.abs(zm_diff.max())])

        # ==============================================================
        # Calculate fractional difference, set divides by zero to Nan
        # ==============================================================

        zm_fracdiff = (np.array(zm_dev_cmp) - np.array(zm_ref_cmp)) / np.array(
            zm_ref_cmp
        )
        zm_fracdiff = np.where(zm_fracdiff == np.inf, np.nan, zm_fracdiff)

        # Test if the frac. diff is zero everywhere or NaN everywhere
        fracdiff_is_all_zero = not np.any(zm_fracdiff)
        fracdiff_is_all_nan = np.isnan(zm_fracdiff).all()

        # ==============================================================
        # Create 3x2 figure
        # ==============================================================

        # Create figs and axes objects
        figs, ((ax0, ax1), (ax2, ax3), (ax4, ax5)) = plt.subplots(
            3, 2, figsize=[12, 15.3]
        )
        # Ensure subplots don't overlap when invoking plt.show()
        if not savepdf:
            plt.subplots_adjust(hspace=0.4)
        # Give the plot a title
        offset = 0.96
        fontsize = 25
        if extra_title_txt is not None:
            figs.suptitle(
                "{}, Zonal Mean ({})".format(varname, extra_title_txt),
                fontsize=fontsize,
                y=offset,
            )
        else:
            figs.suptitle("{}, Zonal Mean".format(varname), fontsize=fontsize, y=offset)

        # ==============================================================
        # Set color map objects.  Use gray for NaNs (no worries,
        # because zonal means are always plotted on lat-alt grids).
        #
        # Use shallow copy (copy.copy() to create color map objects,
        # in order to avoid set_bad() from being applied to the base
        # color table. See: https://docs.python.org/3/library/copy.html
        # ==============================================================

        if use_cmap_RdBu:
            cmap1 = copy.copy(mpl.cm.RdBu_r)
        else:
            cmap1 = copy.copy(WhGrYlRd)
        cmap1.set_bad("gray")

        cmap_plot = copy.copy(mpl.cm.RdBu_r)
        cmap_plot.set_bad(color="gray")

        # ==============================================================
        # Set titles for plots
        # ==============================================================

        if refgridtype == "ll":
            ref_title = "{} (Ref){}\n{}".format(refstr, subtitle_extra, refres)
            dev_title = "{} (Dev){}\n{}".format(devstr, subtitle_extra, devres)
        else:
            ref_title = "{} (Ref){}\nc{} regridded from c{}".format(
                refstr, subtitle_extra, cmpres, refres
            )
            dev_title = "{} (Dev){}\nc{} regridded from c{}".format(
                devstr, subtitle_extra, cmpres, devres
            )

        if regridany:
            absdiff_dynam_title = "Difference ({})\nDev - Ref, Dynamic Range".format(
                cmpres
            )
            absdiff_fixed_title = "Difference ({})\nDev - Ref, Restricted Range [5%,95%]".format(
                cmpres
            )
            fracdiff_dynam_title = "Fractional Difference ({})\n(Dev-Ref)/Ref, Dynamic Range".format(
                cmpres
            )
            fracdiff_fixed_title = "Fractional Difference ({})\n(Dev-Ref)/Ref, Fixed Range".format(
                cmpres
            )
        else:
            absdiff_dynam_title = "Difference\nDev - Ref, Dynamic Range"
            absdiff_fixed_title = "Difference\nDev - Ref, Restricted Range [5%,95%]"
            fracdiff_dynam_title = "Fractional Difference\n(Dev-Ref)/Ref, Dynamic Range"
            fracdiff_fixed_title = "Fractional Difference\n(Dev-Ref)/Ref, Fixed Range"

        # ==============================================================
        # Bundle variables for 6 parallel plotting calls
        # 0 = Ref                 1 = Dev
        # 2 = Dynamic abs diff    3 = Restricted abs diff
        # 4 = Dynamic frac diff   5 = Restricted frac diff
        # ==============================================================

        subplots = [
            "ref",
            "dev",
            "dyn_abs_diff",
            "res_abs_diff",
            "dyn_frac_diff",
            "res_frac_diff",
        ]
        
        all_zeros = [
            ref_is_all_zero,
            dev_is_all_zero,
            absdiff_is_all_zero,
            absdiff_is_all_zero,
            fracdiff_is_all_zero,
            fracdiff_is_all_zero,
        ]

        all_nans = [
            ref_is_all_nan,
            dev_is_all_nan,
            absdiff_is_all_nan,
            absdiff_is_all_nan,
            fracdiff_is_all_nan,
            fracdiff_is_all_nan,
        ]

        plot_vals = [zm_ref, zm_dev, zm_diff, zm_diff, zm_fracdiff, zm_fracdiff]

        axs = [ax0, ax1, ax2, ax3, ax4, ax5]

        cmaps = [cmap1, cmap1, cmap_plot, cmap_plot, cmap_plot, cmap_plot]

        rowcols = [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1)]

        titles = [
            ref_title,
            dev_title,
            absdiff_dynam_title,
            absdiff_fixed_title,
            fracdiff_dynam_title,
            fracdiff_fixed_title,
        ]

        grids = [refgrid, devgrid, cmpgrid, cmpgrid, cmpgrid, cmpgrid]

        if refgridtype != "ll":
            grids[0] = cmpgrid
        if devgridtype != "ll":
            grids[1] = cmpgrid

        extents = [None, None, None, None, None, None]

        masked = ["ZM", "ZM", "ZM", "ZM", "ZM", "ZM"]

        unit_list = [units, units, "unitless", "unitless", "unitless", "unitless"]

        other_all_nans = [dev_is_all_nan, ref_is_all_nan, False, False, False, False]

        gridtypes = [
            cmpgridtype,
            cmpgridtype,
            cmpgridtype,
            cmpgridtype,
            cmpgridtype,
            cmpgridtype,
        ]

        mins = [vmin_ref, vmin_dev, vmin_abs]
        maxs = [vmax_ref, vmax_dev, vmax_abs]

        # Plot
        for i in range(6):
            sixplot(
                subplots[i],
                all_zeros[i],
                all_nans[i],
                plot_vals[i],
                grids[i],
                axs[i],
                rowcols[i],
                titles[i],
                cmaps[i],
                unit_list[i],
                extents[i],
                masked[i],
                other_all_nans[i],
                gridtypes[i],
                mins,
                maxs,
                use_cmap_RdBu,
                match_cbar,
                verbose,
                log_color_scale,
                pedge,
                pedge_ind,
                log_yaxis,
                plot_type="zonal_mean"
            )

        # ==============================================================
        # Update the list of variables with significant differences.
        # Criterion: abs(max(fracdiff)) > 0.1
        # Do not include NaNs in the criterion, because these indicate
        # places where fracdiff could not be computed (div-by-zero).
        # ==============================================================
        if np.abs(np.nanmax(zm_fracdiff)) > 0.1:
            sigdiff_list.append(varname)

        # ==============================================================
        # Add this page of 6-panel plots to the PDF file
        # ==============================================================
        if savepdf:
            pdf = PdfPages(pdfname + "BENCHMARKFIGCREATION.pdf" + str(ivar))
            pdf.savefig(figs)
            pdf.close()
            plt.close(figs)

    # for i in range(n_var):
    #    createfig(i)

    # do not attempt nested thread parallelization due to issues with matplotlib
    if current_process().name != "MainProcess":
        n_job = 1

    # This seems like it shouldn't be here (bmy, 2/21/20)
    #Parallel(n_jobs = n_job) (delayed(createfig)(i) for i in range(n_var))
        
    if savepdf:
        Parallel(n_jobs = n_job) (delayed(createfig)(i) for i in range(n_var))
    else:
        #disable parallel plotting to allow interactive figure plotting
        for i in range(n_var):
            createfig(i)

    # ==================================================================
    # Finish
    # ==================================================================
    if savepdf:
        print("Closed PDF")
        merge = PdfFileMerger()
        for i in range(n_var):
            merge.append(pdfname + "BENCHMARKFIGCREATION.pdf" + str(i))
            os.remove(pdfname + "BENCHMARKFIGCREATION.pdf" + str(i))
        merge.write(pdfname)
        merge.close()
        warnings.showwarning = warning_format


def get_emissions_varnames(commonvars, template=None):
    """
    Will return a list of emissions diagnostic variable names that
    contain a particular search string.

    Args:
    -----
        commonvars : list of strs
            A list of commmon variable names from two data sets.
            (This can be obtained with method gcpy.core.compare_varnames)

        template : str
            String template for matching variable names corresponding
            to emission diagnostics by sector.

    Returns:
    --------
        varnames : list of strs
            A list of variable names corresponding to emission
            diagnostics for a given species and sector.

    Example:
    --------
        >>> import gcpy
        >>> commonvars = ['EmisCO_Anthro', 'EmisNO_Anthro', 'AREA']
        >>> varnames = gcpy.get_emissions_varnames(commonvars, "Emis")
        >>> print(varnames)
        ['EmisCO_Anthro', 'EmisNO_Anthro']
    """

    # Make sure the commonvars list has at least one element
    if len(commonvars) == 0:
        raise ValueError("No valid variable names were passed!")

    # Define template for emission diagnostics by sector
    if template is None:
        raise ValueError("The template argument was not passed!")

    # Find all emission diagnostics for the given species
    varnames = core.filter_names(commonvars, template)

    return varnames


def create_display_name(diagnostic_name):
    """
    Converts a diagnostic name to a more easily digestible name
    that can be used as a plot title or in a table of totals.

    Args:
    -----
        diagnostic_name : str
            Name of the diagnostic to be formatted

    Returns:
    --------
        display_name : str
            Formatted name that can be used as plot titles or in tables
            of emissions totals.

    Remarks:
    --------
        Assumes that diagnostic names will start with either "Emis"
        (for emissions by category) or "Inv" (for emissions by inventory).
        This should be an OK assumption to make since this routine is
        specifically geared towards model benchmarking.

    Example:
    --------
        >>> import gcpy
        >>> diag_name = "EmisCO_Anthro"
        >>> display_name = gcpy.create_display_name(diag_name)
        >>> print(display_name)
        CO Anthro
    """

    # Initialize
    display_name = diagnostic_name

    if "SpeciesRst" in display_name:
        display_name = display_name.split("_")[1]

    # Special handling for Inventory totals
    if "INV" in display_name.upper():
        display_name = display_name.replace("_", " ")

    # Replace text
    for v in ["Emis", "EMIS", "emis", "Inv", "INV", "inv"]:
        display_name = display_name.replace(v, "")

    # Replace underscores
    display_name = display_name.replace("_", " ")

    return display_name


def print_totals(ref, refstr, dev, devstr, f, mass_tables=False, masks=None):
    """
    Computes and prints Ref and Dev totals (as well as the difference
    Dev - Ref) for two xarray DataArray objects.

    Args:
    -----
        ref : xarray DataArray
            The first DataArray to be compared (aka "Reference")

        refstr : str
            A string that can be used to identify refdata.
            (e.g. a model version number or other identifier).

        dev : xarray DataArray
            The second DataArray to be compared (aka "Development")

        devstr : str
            A string that can be used to identify devdata
            (e.g. a model version number or other identifier).

        f : file
            File object denoting a text file where output will be directed.

    Keyword Arguments (optional):
    -----------------------------
        mass_tables : bool
            Set this switch to True if you would like to print out
            totals of global mass.
            Default value: False (i.e. print out emissions totals)

        masks : dict of xarray DataArray
            Dictionary containing the tropospheric mask arrays
            for Ref and Dev.  If this keyword argument is passed,
            then print_totals will print tropospheric totals
            NOTE: This option is only used if mass_tables=True.
            Default value: None (i.e. print whole-atmosphere totals)

    Remarks:
    --------
        This is an internal method.  It is meant to be called from method
        create_total_emissions_table or create_global_mass_table instead of
        being called directly.
    """

    # ==================================================================
    # Initialization and error checks
    # ==================================================================

    # Make sure that both Ref and Dev are xarray DataArray objects
    if not isinstance(ref, xr.DataArray):
        raise TypeError("The ref argument must be an xarray DataArray!")
    if not isinstance(dev, xr.DataArray):
        raise TypeError("The dev argument must be an xarray DataArray!")

    # Determine if either Ref or Dev have all NaN values:
    ref_is_all_nan = np.isnan(ref.values).all()
    dev_is_all_nan = np.isnan(dev.values).all()

    # If Ref and Dev do not contain all NaNs, then make sure
    # that Ref and Dev have the same units before proceeding.
    if (not ref_is_all_nan) and (not dev_is_all_nan):
        if ref.units != dev.units:
            msg = 'Ref has units "{}", but Dev array has units "{}"'.format(
                ref.units, dev.units
            )
            raise ValueError(msg)

    # ==================================================================
    # Get the diagnostic name and units
    # ==================================================================
    if dev_is_all_nan:
        diagnostic_name = ref.name
        units = ref.units
    else:
        diagnostic_name = dev.name
        units = dev.units

    # Create the display name by editing the diagnostic name
    display_name = create_display_name(diagnostic_name)

    # Special handling for totals
    if "_TOTAL" in diagnostic_name.upper():
        print("-" * 79, file=f)

    # ==================================================================
    # Sum the Ref array (or set to NaN if missing)
    # ==================================================================
    if ref_is_all_nan:
        total_ref = np.nan
    else:
        if masks is None:
            total_ref = np.sum(ref.values)
        else:
            arr = np.ma.masked_array(ref.values, masks["Ref_TropMask"])
            total_ref = np.sum(arr)

    # ==================================================================
    # Sum the Dev array (or set to NaN if missing)
    # ==================================================================
    if dev_is_all_nan:
        total_dev = np.nan
    else:
        if masks is None:
            total_dev = np.sum(dev.values)
        else:
            arr = np.ma.masked_array(dev.values, masks["Dev_TropMask"])
            total_dev = np.sum(arr)

    # ==================================================================
    # Compute differences (or set to NaN if missing)
    # ==================================================================
    if ref_is_all_nan or dev_is_all_nan:
        diff = np.nan
    else:
        diff = total_dev - total_ref

    # ==================================================================
    # Compute % differences (or set to NaN if missing)
    # If ref is very small, near zero, also set the % diff to NaN
    # ==================================================================
    if mass_tables:
        if np.isnan(total_ref) or np.isnan(total_dev):
            pctdiff = np.nan
        else:
            pctdiff = ((total_dev - total_ref) / total_ref) * 100.0
            if total_ref < 1.0e-15:
                pctdiff = np.nan

    # ==================================================================
    # Write output to file
    # ==================================================================
    if mass_tables:
        print(
            "{} : {:18.6f}  {:18.6f}  {:13.6f}  {:8.3f}".format(
                display_name.ljust(12), total_ref, total_dev, diff, pctdiff
            ),
            file=f,
        )
    else:
        print(
            "{} : {:18.6f}  {:18.6f}  {:13.6f} {}".format(
                display_name.ljust(21), total_ref, total_dev, diff, units
            ),
            file=f,
        )


def create_total_emissions_table(
    refdata,
    refstr,
    devdata,
    devstr,
    species,
    outfilename,
    interval=[2678400.0],
    template="Emis{}_",
    ref_area_varname="AREA",
    dev_area_varname="AREA",
):
    """
    Creates a table of emissions totals (by sector and by inventory)
    for a list of species in contained in two data sets.  The data sets,
    which typically represent output from two differnet model versions,
    are usually contained in netCDF data files.

    Args:
    -----
        refdata : xarray Dataset
            The first data set to be compared (aka "Reference" or "Ref").

        refstr : str
            A string that can be used to identify refdata
            (e.g. a model version number or other identifier).

        devdata : xarray Dataset
            The second data set to be compared (aka "Development" or "Dev").

        devstr: str
            A string that can be used to identify the data set specified
            by devfile (e.g. a model version number or other identifier).

        species : dict
            Dictionary containing the name of each species and the target
            unit that emissions will be converted to. The format of
            species is as follows:

                { species_name : target_unit", etc. }

            where "species_name" and "target_unit" are strs.

        outfilename : str
            Name of the text file which will contain the table of
            emissions totals.

    Keyword Args (optional):
    ------------------------
        interval : float
            The length of the data interval in seconds. By default, interval
            is set to the number of seconds in a 31-day month (86400 * 31),
            which corresponds to typical benchmark simulation output.

        template : str
            Template for the diagnostic names that are contained both
            "Reference" and "Development" data sets.  If not specified,
            template will be set to "Emis{}", where {} will be replaced
            by the species name.

        ref_area_varname : str
            Name of the variable containing the grid box surface areas
            (in m2) in the ref dataset.
            Default value: 'AREA'

        dev_area_varname : str
            Name of the variable containing the grid box surface areas
            (in m2) in the dev dataset.
            Default value: 'AREA'

    Remarks:
    --------
        This method is mainly intended for model benchmarking purposes,
        rather than as a general-purpose tool.

        Species properties (such as molecular weights) are read from a
        YAML file called "species_database.yml".

    Example:
    --------
        Print the total of CO and ACET emissions in two different
        data sets, which represent different model versions:

        >>> include gcpy
        >>> include xarray as xr
        >>> reffile = '~/output/12.1.1/HEMCO_diagnostics.201607010000.nc'
        >>> refstr = '12.1.1'
        >>> refdata = xr.open_dataset(reffile)
        >>> devfile = '~/output/12.2.0/HEMCO_sa.diagnostics.201607010000.nc'
        >>> devstr = '12.2.0'
        >>> devdata = xr.open_dataset(devfile)
        >>> outfilename = '12.2.0_emission_totals.txt'
        >>> species = { 'CO' : 'Tg', 'ACET', 'Tg C', 'ALK4', 'Gg C' }
        >>> create_total_emissions_table(refdata, refstr, devdata, devstr,
            species, outfilename)
    """

    # ==================================================================
    # Initialization
    # ==================================================================

    # Make sure refdata and devdata are both xarray Dataset objects
    if not isinstance(refdata, xr.Dataset):
        raise TypeError("The refdata argument must be an xarray Dataset!")
    if not isinstance(devdata, xr.Dataset):
        raise TypeError("The devdata argument must be an xarray Dataset!")

    # Make sure that the area variable is present in both refdata and devdata
    if ref_area_varname not in refdata.data_vars.keys():
        raise ValueError(
            "Area variable {} is not in the ref Dataset!".format(ref_area_varname)
        )
    if dev_area_varname not in devdata.data_vars.keys():
        raise ValueError(
            "Area variable {} is not in the dev Dataset!".format(dev_area_varname)
        )

    # Load a YAML file containing species properties (such as
    # molecular weights), which we will need for unit conversions.
    # This is located in the "data" subfolder of this folder where
    # this benchmark.py file is found.
    properties_path = os.path.join(os.path.dirname(__file__), "species_database.yml")
    properties = yaml_load_file(open(properties_path))

    # ==================================================================
    # Get the list of emission variables for which we will print totals
    # ==================================================================

    # Make sure that Ref and Dev datasets have the same variables.
    # Variables that are in Ref but not in Dev will be added to Dev
    # with all missing values (NaNs). And vice-versa.
    [refdata, devdata] = add_missing_variables(refdata, devdata)

    # Find all common variables between the two datasets
    # and get the lists of variables only in Ref and only in Dev,
    vardict = core.compare_varnames(refdata, devdata, quiet=True)
    cvars = vardict["commonvars"]
    refonly = vardict["refonly"]
    devonly = vardict["devonly"]

    # =================================================================
    # Open the file for output
    # =================================================================
    try:
        f = open(outfilename, "w")
    except FileNotFoundError:
        msg = "Could not open {} for writing!".format(outfilename)
        raise FileNotFoundError(msg)

    # =================================================================
    # Loop through all of the species are in species_dict
    # =================================================================
    for species_name, target_units in species.items():

        # Get a list of emission variable names for each species
        diagnostic_template = template.format(species_name)
        varnames = get_emissions_varnames(cvars, diagnostic_template)

        # Also add variables that might be in either Ref or Dev
        # but not the other.  This will allow us to print totals
        # for all species (and print NaN for the missing ones).
        if len(refonly) > 0:
            matching = [v for v in refonly if diagnostic_template in v]
            varnames = varnames + matching
        if len(devonly) > 0:
            matching = [v for v in devonly if diagnostic_template in v]
            varnames = varnames + matching

        # Sort the list again  to account for new variables added above
        varnames.sort()

        # If no emissions are found, then skip to next species
        if len(varnames) == 0:
            print("No emissions found for {} ... skippping".format(species_name))
            continue

        # Check if there is a total emissions variable in the list
        vartot = [v for v in varnames if "_TOTAL" in v.upper()]

        # Push the total variable to the last list element
        # so that it will be printed last of all
        if len(vartot) == 1:
            varnames.append(varnames.pop(varnames.index(vartot[0])))

        # Title strings
        if "Inv" in template:
            print("Computing inventory totals for {}".format(species_name))
            title1 = "### Emissions totals for inventory {}".format(species_name)
        else:
            print("Computing emissions totals for {}".format(species_name))
            title1 = "### Emissions totals for species {}".format(species_name)

        title2 = "### Ref = {}; Dev = {}".format(refstr, devstr)

        # Print header to file
        print("#" * 79, file=f)
        print("{}{}".format(title1.ljust(76), "###"), file=f)
        print("{}{}".format(title2.ljust(76), "###"), file=f)
        print("#" * 79, file=f)
        print(
            "{}{}{}{}".format(
                " ".ljust(22), "Ref".rjust(20), "Dev".rjust(20), "Dev - Ref".rjust(15)
            ),
            file=f,
        )

        # =============================================================
        # Loop over all emissions variables corresponding to this
        # species and print their totals in Ref and Dev to the file.
        # =============================================================
        for v in varnames:

            if "Inv" in template:
                spc_name = v.split("_")[1]
            else:
                spc_name = species_name

            # Get a list of properties for the given species
            species_properties = properties.get(spc_name)

            # If no properties are found, then skip to next species
            if species_properties is None:
                print("No properties found for {} ... skippping".format(spc_name))
                continue

            # Convert units of Ref and Dev and save to numpy ndarray objects
            # (or set to NaN if the variable is not found in Ref or Dev)
            if v in refonly and v not in devonly:

                # Convert units of Ref
                refarray = convert_units(
                    refdata[v],
                    spc_name,
                    species_properties,
                    target_units,
                    interval,
                    refdata[ref_area_varname],
                )

                # Set Dev to NaN (missing values) everywhere
                devarray = core.create_dataarray_of_nan(
                    name=refdata[v].name,
                    sizes=devdata.sizes,
                    coords=devdata.coords,
                    attrs=refdata[v].attrs,
                )

            elif v in devonly and v not in refonly:

                # Convert units of Dev
                devarray = convert_units(
                    devdata[v],
                    spc_name,
                    species_properties,
                    target_units,
                    interval,
                    devdata[dev_area_varname],
                )

                # Set Ref to NaN (missing values) everywhere
                refarray = core.create_dataarray_of_nan(
                    name=devdata[v].name,
                    sizes=refdata.sizes,
                    coords=refdata.coords,
                    attrs=devdata[v].attrs,
                )

            else:

                # Convert units of both Ref and Dev
                refarray = convert_units(
                    refdata[v],
                    spc_name,
                    species_properties,
                    target_units,
                    interval,
                    refdata[ref_area_varname],
                )
                devarray = convert_units(
                    devdata[v],
                    spc_name,
                    species_properties,
                    target_units,
                    interval,
                    devdata[dev_area_varname],
                )

            # ==========================================================
            # Print emission totals for Ref and Dev
            # ==========================================================
            print_totals(refarray, refstr, devarray, devstr, f)

        # Add newlines before going to the next species
        print(file=f)
        print(file=f)

    # =================================================================
    # Close file
    # =================================================================
    f.close()


def create_global_mass_table(
    refdata,
    refstr,
    devdata,
    devstr,
    varlist,
    met_and_masks,
    label,
    trop_only=False,
    outfilename="GlobalMass_TropStrat.txt",
    verbose=False,
):
    """
    Creates a table of global masses for a list of species in contained in
    two data sets.  The data sets,  which typically represent output from two
    differnet model versions, are usually contained in netCDF data files.

    Args:
    -----
        refdata : xarray Dataset
            The first data set to be compared (aka "Reference").

        refstr : str
            A string that can be used to identify refdata
            (e.g. a model version number or other identifier).

        devdata : xarray Dataset
            The second data set to be compared (aka "Development").

        devstr: str
            A string that can be used to identify the data set specified
            by devfile (e.g. a model version number or other identifier).

        varlist : list of strings
            List of species concentation variable names to include
            in the list of global totals.

        met_and_masks : dict of xarray DataArray
            Dictionary containing the meterological variables and
            masks for the Ref and Dev datasets.
        
        label : str
            Label to go in the header string.  Can be used to
            pass the month & year.

    Keyword Args (optional):
    ------------------------
        trop_only : book
            Set this switch to True if you wish to print totals
            only for the troposphere.
            Default value: False (i.e. print whole-atmosphere totals).

        outfilename : str
            Name of the text file which will contain the table of
            emissions totals.
            Default value: "GlobalMass_TropStrat.txt"

        verbose : bool
            Set this switch to True if you wish to print out extra
            informational messages.
            Default value: False

    Remarks:
    --------
        This method is mainly intended for model benchmarking purposes,
        rather than as a general-purpose tool.

        Species properties (such as molecular weights) are read from a
        YAML file called "species_database.yml".

        The area variable for GEOS-Chem "Classic" will be "AREA",
        but for GCHP it will be "Met_AREAM2".
    """

    # ==================================================================
    # Initialization
    # ==================================================================

    # Make sure refdata and devdata are xarray Dataset objects
    if not isinstance(refdata, xr.Dataset):
        raise TypeError("The refdata argument must be an xarray Dataset!")
    if not isinstance(devdata, xr.Dataset):
        raise TypeError("The devdata argument must be an xarray Dataset!")

    # Make sure required arguments are passed
    if varlist is None:
        raise ValueError('The "varlist" argument was not passed!')
    if met_and_masks is None:
        raise ValueError('The "met_and_masks" argument was not passed!')

    # Load a YAML file containing species properties (such as
    # molecular weights), which we will need for unit conversions.
    # This is located in the "data" subfolder of this current directory.2
    properties_path = os.path.join(os.path.dirname(__file__), "species_database.yml")
    properties = yaml_load_file(open(properties_path))

    # ==================================================================
    # Open file for output
    # ==================================================================

    # Create file
    try:
        f = open(outfilename, "w")
    except FileNotFoundError:
        msg = "Could not open {} for writing!".format(outfilename)
        raise FileNotFoundError(msg)

    # Title strings
    if trop_only:
        title1 = "### Global mass (Gg) {} (Trop only)".format(label)
    else:
        title1 = "### Global mass (Gg) {} (Trop + Strat)".format(label)
    title2 = "### Ref = {}; Dev = {}".format(refstr, devstr)

    # Print header to file
    print("#" * 79, file=f)
    print("{}{}".format(title1.ljust(76), "###"), file=f)
    print("{}{}".format(title2.ljust(76), "###"), file=f)
    print("#" * 79, file=f)
    print(
        "{}{}{}{}{}".format(
            " ".ljust(13),
            "Ref".rjust(20),
            "Dev".rjust(20),
            "Dev - Ref".rjust(15),
            "% diff".rjust(10),
        ),
        file=f,
    )

    # ==================================================================
    # Print global masses for all species
    #
    # NOTE: By this point, all species will be in both Ref and Dev'
    # because we have added them in the calling routine
    # ==================================================================
    for v in varlist:

        # Get the species name
        spc_name = v.split("_")[1]

        # Get a list of properties for the given species
        species_properties = properties.get(spc_name)

        # If no properties are found, then skip to next species
        if species_properties is None:
            print("No properties found for {} ... skippping".format(spc_name))
            continue

        # Specify target units
        target_units = "Gg"
        mol_wt_g = species_properties.get("MW_g")
        if mol_wt_g is None:
            #            print('No molecular weight found for {} ... skippping'.format(
            #                  spc_name))
            continue

        # ==============================================================
        # Convert units of Ref and save to a DataArray
        # (or skip if Ref contains NaNs everywhere)
        # ==============================================================
        if not np.isnan(refdata[v].values).all():
            refarray = convert_units(
                refdata[v],
                spc_name,
                species_properties,
                target_units,
                area_m2=met_and_masks["Ref_Area"],
                delta_p=met_and_masks["Ref_Delta_P"],
                box_height=met_and_masks["Ref_BxHeight"],
            )

        # ==============================================================
        # Convert units of Dev and save to a DataArray
        # (or skip if Dev contains NaNs everywhere)
        # ==============================================================
        if not np.isnan(devdata[v].values).all():
            devarray = convert_units(
                devdata[v],
                spc_name,
                species_properties,
                target_units,
                area_m2=met_and_masks["Dev_Area"],
                delta_p=met_and_masks["Dev_Delta_P"],
                box_height=met_and_masks["Dev_BxHeight"],
            )

        # ==============================================================
        # Print global masses for Ref and Dev
        # (we will mask out tropospheric boxes in print_totals)
        # ==============================================================
        if trop_only:
            print_totals(
                refarray,
                refstr,
                devarray,
                devstr,
                f,
                mass_tables=True,
                masks=met_and_masks,
            )
        else:
            print_totals(refarray, refstr, devarray, devstr, f, mass_tables=True)

    # ==================================================================
    # Close files
    # ==================================================================
    f.close()


def create_budget_table(
    devdata,
    devstr,
    region,
    species,
    varnames,
    outfilename,
    interval=[2678400.0],
    template="Budget_{}",
):
    """
    Creates a table of budgets by species and component for a data set.

    Args:
    -----
        devdata : xarray Dataset
            The second data set to be compared (aka "Development").

        devstr: str
            A string that can be used to identify the data set specified
            by devfile (e.g. a model version number or other identifier).

        region : str
            Name of region for which budget will be computed.

        species : List of strings
            List of species  to include in budget tables.

        varnames : List of strings
            List of variable names in the budget diagnostics.

        outfilename : str
            Name of the text file which will contain the table of
            emissions totals.

    Keyword Args (optional):
    ------------------------
        interval : list of float
            The length of the data interval in seconds. By default, interval
            is set to [2678400.0], which is the number of seconds in July
            (our 1-month benchmarking month).

        template : str
            Template for the diagnostic names that are contained in the
            data set. If not specified, template will be set to "Budget_{}",
            where {} will be replaced by the species name.     
            
    Remarks:
    --------
        This method is mainly intended for model benchmarking purposes,
        rather than as a general-purpose tool.
    """

    # ==================================================================
    # Initialization
    # ==================================================================

    # Error check arguments
    if not isinstance(devdata, xr.Dataset):
        raise TypeError("The devdata argument must be an xarray Dataset!")

    # Open file for output
    try:
        f = open(outfilename, "w")
    except FileNotFoundError:
        msg = "Could not open {} for writing!".format(outfilename)
        raise FileNotFoundError(msg)

    # ==================================================================
    # Loop over species
    # ==================================================================
    for spc_name in species:

        # Title string
        title = "### {} budget totals for species {}".format(devstr, spc_name)

        # Write header to file
        print("#" * 79, file=f)
        print("{}{}".format(title.ljust(76), "###"), file=f)
        print("#" * 79, file=f)

        # Get variable names for this species
        spc_vars = [v for v in varnames if v.endswith("_" + spc_name)]

        for v in spc_vars:

            # Component name
            comp_name = v.replace("Budget", "")
            comp_name = comp_name.replace("_" + spc_name, "")
            comp_name = comp_name.replace(region, "")

            # Convert from kg/s to Tg
            devarray = devdata[v] * interval * 1e-9
            units = "Tg"

            # Compute sum
            total_dev = np.sum(devarray.values)

            # Write output
            print(
                "{} : {:13.6e} {}".format(comp_name.ljust(12), total_dev, units), file=f
            )

        # Add new lines before going to the next species
        print(file=f)
        print(file=f)

    # Close file
    f.close()


def get_species_categories(benchmark_type="FullChemBenchmark"):
    """
    Returns the list of benchmark categories that each species
    belongs to.  This determines which PDF files will contain the
    plots for the various species.

    Args:
    -----
        benchmark_type : str
            Specifies the type of the benchmark (either 
            FullChemBenchmark (default) or TransportTracersBenchmark).

    Returns:
        spc_cat_dict : dict
            A nested dictionary of categories (and sub-categories)
            and the species belonging to each.

    NOTE: The benchmark categories are specified in YAML file
    benchmark_species.yml.
    """
    yamlfile = os.path.join(os.path.dirname(__file__), spc_categories)
    with open(yamlfile, "r") as f:
        spc_cat_dict = yaml.load(f.read())
    return spc_cat_dict[benchmark_type]


def archive_species_categories(dst):
    """
    Writes the list of benchmark categories to a YAML file
    named "benchmark_species.yml".

    Args:
    -----
        dst : str
            Name of the folder where the YAML file containing
            benchmark categories ("benchmark_species.yml")
            will be written.
    """

    src = os.path.join(os.path.dirname(__file__), spc_categories)
    print("Archiving {} in {}".format(spc_categories, dst))
    shutil.copyfile(src, os.path.join(dst, spc_categories))


def make_benchmark_plots(
    ref,
    refstr,
    dev,
    devstr,
    dst="./1mo_benchmark",
    subdst=None,
    overwrite=False,
    verbose=False,
    collection="SpeciesConc",
    benchmark_type="FullChemBenchmark",
    plot_by_spc_cat=True,
    restrict_cats=[],
    plots=["sfc", "500hpa", "zonalmean"],
    use_cmap_RdBu=False,
    log_color_scale=False,
    sigdiff_files=None,
    n_job=-1,
):
    """
    Creates PDF files containing plots of species concentration
    for model benchmarking purposes.

    Args:
    -----
        ref: str
            Path name for the "Ref" (aka "Reference") data set.

        refstr : str
            A string to describe ref (e.g. version number)

        dev : str
            Path name for the "Dev" (aka "Development") data set.
            This data set will be compared against the "Reference"
            data set.

        devstr : str
            A string to describe dev (e.g. version number)

    Keyword Args (optional):
    ------------------------
        dst : str
            A string denoting the destination folder where a PDF
            file containing plots will be written.
            Default value: ./1mo_benchmark

        subdst : str
            A string denoting the sub-directory of dst where PDF
            files containing plots will be written.  In practice,
            subdst is only needed for the 1-year benchmark output,
            and denotes a date string (such as "Jan2016") that
            corresponds to the month that is being plotted.
            Default value: None

        overwrite : boolean
            Set this flag to True to overwrite files in the
            destination folder (specified by the dst argument).
            Default value: False.

        verbose : boolean
            Set this flag to True to print extra informational output.
            Default value: False.

        plot_by_spc_cat: logical
            Set this flag to False to send plots to one file rather
            than separate file per category.
            Default value: True

        restrict_cats : list of strings
            List of benchmark categories in benchmark_categories.yml to make
            plots for. If empty, plots are made for all categories.
            Default value: empty

        plots : list of strings
            List of plot types to create.
            Default value: ['sfc', '500hpa', 'zonalmean']

        log_color_scale: boolean         
            Set this flag to True to enable plotting data (not diffs)
            on a log color scale.
            Default value: False

        sigdiff_files : list of str
            Filenames that will contain the lists of species having
            significant differences in the 'sfc', '500hpa', and
            'zonalmean' plots.  These lists are needed in order to
            fill out the benchmark approval forms.
            Default value: None
    """

    # NOTE: this function could use some refactoring;
    # abstract processing per category?

    # ==================================================================
    # Initialization and data read
    # ==================================================================
    if os.path.isdir(dst) and not overwrite:
        print(
            "Directory {} exists. Pass overwrite=True to overwrite files in that directory, if any.".format(
                dst
            )
        )
        return
    elif not os.path.isdir(dst):
        os.mkdir(dst)

    # Define extra title text (usually a date string)
    # for the top-title of the plot
    if subdst is not None:
        extra_title_txt = subdst
    else:
        extra_title_txt = None

    # Ref dataset
    try:
        refds = xr.open_dataset(ref, drop_variables=skip_these_vars)
    except FileNotFoundError:
        raise FileNotFoundError("Could not find Ref file: {}".format(ref))

    # Dev dataset
    try:
        devds = xr.open_dataset(dev, drop_variables=skip_these_vars)
    except FileNotFoundError:
        raise FileNotFoundError("Could not find Dev file: {}!".format(dev))

    # If sending plots to one file then do all plots here and return
    if not plot_by_spc_cat:
        [refds, devds] = add_missing_variables(refds, devds)
        var_prefix = 'SpeciesConc_'
        varlist = [k for k in refds.data_vars.keys() if var_prefix in k]
        varlist.sort()
        # Surface
        pdfname = os.path.join(dst,'SpeciesConc_Sfc.pdf')
        compare_single_level(refds, refstr, devds, devstr, 
                             varlist=varlist,
                             pdfname=pdfname,
                             use_cmap_RdBu=use_cmap_RdBu,
                             log_color_scale=log_color_scale,
                             extra_title_txt=extra_title_txt)
        add_bookmarks_to_pdf(pdfname, varlist, remove_prefix=var_prefix,
                             verbose=verbose)
        # 500 hPa
        pdfname = os.path.join(dst,'SpeciesConc_500hPa.pdf')
        compare_single_level(refds, refstr, devds, devstr,
                             ilev=22,
                             varlist=varlist,
                             pdfname=pdfname,
                             use_cmap_RdBu=use_cmap_RdBu,
                             log_color_scale=log_color_scale,
                             extra_title_txt=extra_title_txt)

        add_bookmarks_to_pdf(pdfname, varlist, remove_prefix=var_prefix,
                             verbose=verbose)
        # Zonal mean
        pdfname = os.path.join(dst,'SpeciesConc_ZnlMn.pdf')
        compare_zonal_mean(refds, refstr, devds, devstr,
                           varlist=varlist,
                           pdfname=pdfname,
                           use_cmap_RdBu=use_cmap_RdBu,
                           log_color_scale=log_color_scale,
                           extra_title_txt=extra_title_txt)

        add_bookmarks_to_pdf(pdfname, varlist, remove_prefix=var_prefix,
                             verbose=verbose)
        return

    # FullChemBenchmark has lumped species (TransportTracers does not)
    if "FullChem" in benchmark_type:
        refds = core.add_lumped_species_to_dataset(refds, verbose=verbose)
        devds = core.add_lumped_species_to_dataset(devds, verbose=verbose)
        core.archive_lumped_species_definitions(dst)

    # Get the list of species categories
    catdict = get_species_categories(benchmark_type)
    archive_species_categories(dst)

    # Make sure that Ref and Dev datasets have the same variables.
    # Variables that are in Ref but not in Dev will be added to Dev
    # with all missing values (NaNs). And vice-versa.
    [refds, devds] = add_missing_variables(refds, devds)

    # Collection prefix
    coll_prefix = collection.strip() + "_"

    # ==================================================================
    # Create the plots!
    # ==================================================================

    # Use dictionaries to maintain order of significant difference categories
    dict_sfc = {}
    dict_500 = {}
    dict_zm = {}

    def createplots(i, filecat):
        cat_diff_dict = {'sfc' : [], '500' : [], 'zm' : []}
        # Suppress harmless run-time warnings from all threads
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        warnings.filterwarnings("ignore", category=UserWarning)

        # If restrict_cats list is passed,
        # skip all categories except those in the list
        if restrict_cats and filecat not in restrict_cats:
            return

        # Create a directory for each category.
        # If subdst is passed, then create a subdirectory in each
        # category directory (e.g. as for the 1-year benchmark).
        catdir = os.path.join(dst, filecat)
        if not os.path.isdir(catdir):
            os.mkdir(catdir)
        if subdst is not None:
            catdir = os.path.join(catdir, subdst)
            if not os.path.isdir(catdir):
                os.mkdir(catdir)

        varlist = []
        warninglist = []
        for subcat in catdict[filecat]:
            for spc in catdict[filecat][subcat]:
                varname = coll_prefix + spc
                if varname not in refds.data_vars or \
                   varname not in devds.data_vars:
                    warninglist.append(varname)
                    continue
                varlist.append(varname)
        if warninglist != []:
            print(
                "\n\nWarning: variables in {} category not in dataset: {}".format(
                    filecat, warninglist
                )
            )

        # -----------------------
        # Surface plots
        # -----------------------
        if "sfc" in plots:

            if subdst is not None:
                pdfname = os.path.join(
                    catdir, "{}_Surface_{}.pdf".format(filecat, subdst)
                )
            else:
                pdfname = os.path.join(catdir, "{}_Surface.pdf".format(filecat))

            diff_sfc = []
            compare_single_level(
                refds,
                refstr,
                devds,
                devstr,
                varlist=varlist,
                ilev=0,
                pdfname=pdfname,
                use_cmap_RdBu=use_cmap_RdBu,
                log_color_scale=log_color_scale,
                extra_title_txt=extra_title_txt,
                sigdiff_list=diff_sfc,
            )
            diff_sfc[:] = [v.replace(coll_prefix, "") for v in diff_sfc]
            cat_diff_dict['sfc'] = diff_sfc
            add_nested_bookmarks_to_pdf(
                pdfname, filecat, catdict,
                warninglist, remove_prefix=coll_prefix
            )

        # -----------------------
        # 500 hPa plots
        # -----------------------
        if "500hpa" in plots:

            if subdst is not None:
                pdfname = os.path.join(
                    catdir, "{}_500hPa_{}.pdf".format(filecat, subdst)
                )
            else:
                pdfname = os.path.join(catdir, "{}_500hPa.pdf".format(filecat))

            diff_500 = []
            compare_single_level(
                refds,
                refstr,
                devds,
                devstr,
                varlist=varlist,
                ilev=22,
                pdfname=pdfname,
                use_cmap_RdBu=use_cmap_RdBu,
                log_color_scale=log_color_scale,
                extra_title_txt=extra_title_txt,
                sigdiff_list=diff_500,
            )
            diff_500[:] = [v.replace(coll_prefix, "") for v in diff_500]
            #dict_500[filecat] = diff_500
            cat_diff_dict['500'] = diff_500
            add_nested_bookmarks_to_pdf(
                pdfname, filecat, catdict,
                warninglist, remove_prefix=coll_prefix
            )

        # -----------------------
        # Zonal mean plots
        # -----------------------
        if "zonalmean" in plots or "zm" in plots:

            if subdst is not None:
                pdfname = os.path.join(
                    catdir, "{}_FullColumn_ZonalMean_{}.pdf".format(filecat, subdst)
                )
            else:
                pdfname = os.path.join(
                    catdir, "{}_FullColumn_ZonalMean.pdf".format(filecat)
                )

            diff_zm = []
            compare_zonal_mean(
                refds,
                refstr,
                devds,
                devstr,
                varlist=varlist,
                pdfname=pdfname,
                use_cmap_RdBu=use_cmap_RdBu,
                log_color_scale=log_color_scale,
                extra_title_txt=extra_title_txt,
                sigdiff_list=diff_zm,
            )
            diff_zm[:] = [v.replace(coll_prefix, "") for v in diff_zm]
            #dict_zm = diff_zm
            cat_diff_dict['zm'] = diff_zm
            add_nested_bookmarks_to_pdf(
                pdfname, filecat, catdict,
                warninglist, remove_prefix=coll_prefix
            )

            # Strat_ZonalMean plots will use a log-pressure Y-axis, with
            # a range of 1..100 hPa, as per GCSC request. (bmy, 8/13/19)
            if subdst is not None:
                pdfname = os.path.join(
                    catdir, "{}_Strat_ZonalMean_{}.pdf".format(filecat, subdst)
                )
            else:
                pdfname = os.path.join(catdir, "{}_Strat_ZonalMean.pdf".format(filecat))

            compare_zonal_mean(
                refds,
                refstr,
                devds,
                devstr,
                varlist=varlist,
                pdfname=pdfname,
                use_cmap_RdBu=use_cmap_RdBu,
                pres_range=[1, 100],
                log_yaxis=True,
                extra_title_txt=extra_title_txt,
                log_color_scale=log_color_scale,
            )
            add_nested_bookmarks_to_pdf(
                pdfname, filecat, catdict,
                warninglist, remove_prefix=coll_prefix
            )
        return {filecat : cat_diff_dict}
    # Create the plots in parallel
    results = Parallel(n_jobs=n_job)(
        delayed(createplots)(i, filecat) for i, filecat in enumerate(catdict)
    )

    dict_sfc = {list(result.keys())[0] : result[list(result.keys())[0]]['sfc'] for result in results}
    dict_500 = {list(result.keys())[0] : result[list(result.keys())[0]]['500'] for result in results}
    dict_zm  = {list(result.keys())[0] : result[list(result.keys())[0]]['zm']  for result in results}    
    
    # ==============================================================
    # Write the list of species having significant differences,
    # which we need to fill out the benchmark approval forms.
    # ==============================================================
    if sigdiff_files != None:
        for filename in sigdiff_files:
            if "sfc" in plots:
                if "sfc" in filename:
                    with open(filename, "a+") as f:
                        for c, diff_list in dict_sfc.items():
                            print("* {}: ".format(c), file=f, end="")
                            for v in diff_list:
                                print("{} ".format(v), file=f, end="")
                            print(file=f)
                        f.close()

            if "500hpa" in plots:
                if "500hpa" in filename:
                    with open(filename, "a+") as f:
                        for c, diff_list in dict_500.items():
                            print("* {}: ".format(c), file=f, end="")
                            for v in diff_list:
                                print("{} ".format(v), file=f, end="")
                            print(file=f)
                        f.close()

            if "zonalmean" in plots or "zm" in plots:
                if "zonalmean" in filename or "zm" in filename:
                    with open(filename, "a+") as f:
                        for c, diff_list in dict_zm.items():
                            print("* {}: ".format(c), file=f, end="")
                            for v in diff_list:
                                print("{} ".format(v), file=f, end="")
                            print(file=f)
                        f.close()


def make_benchmark_emis_plots(
    ref,
    refstr,
    dev,
    devstr,
    dst="./1mo_benchmark",
    subdst=None,
    plot_by_spc_cat=False,
    plot_by_hco_cat=False,
    overwrite=False,
    verbose=False,
    flip_ref=False,
    flip_dev=False,
    log_color_scale=False,
    sigdiff_files=None,
    n_job=-1,
):
    """
    Creates PDF files containing plots of emissions for model
    benchmarking purposes. This function is compatiblity with benchmark 
    simulation output only. It is not compatible with transport tracers
    emissions diagnostics.

    Args:
    -----
        ref: str
            Path name for the "Ref" (aka "Reference") data set.

        refstr : str
            A string to describe ref (e.g. version number)

        dev : str
            Path name for the "Dev" (aka "Development") data set.
            This data set will be compared against the "Reference"
            data set.

        devstr : str
            A string to describe dev (e.g. version number)

    Keyword Args (optional):
    ------------------------
        dst : str
            A string denoting the destination folder where
            PDF files containing plots will be written.
            Default value: './1mo_benchmark

        subdst : str
            A string denoting the sub-directory of dst where PDF
            files containing plots will be written.  In practice,
            and denotes a date string (such as "Jan2016") that
            corresponds to the month that is being plotted.
            Default value: None

        plot_by_spc_cat : boolean
            Set this flag to True to separate plots into PDF files
            according to the benchmark species categories (e.g. Oxidants,
            Aerosols, Nitrogen, etc.)  These categories are specified
            in the YAML file benchmark_species.yml.
            Default value: False

        plot_by_hco_cat : boolean
            Set this flag to True to separate plots into PDF files
            according to HEMCO emissions categories (e.g. Anthro,
            Aircraft, Bioburn, etc.)
            Default value: False

        overwrite : boolean
            Set this flag to True to overwrite files in the
            destination folder (specified by the dst argument).
            Default value: False

        verbose : boolean
            Set this flag to True to print extra informational output.
            Default value: False

        flip_ref : boolean
            Set this flag to True to reverse the vertical level
            ordering in the "Ref" dataset (in case "Ref" starts
            from the top of atmosphere instead of the surface).
            Default value: False

        flip_dev : boolean
            Set this flag to True to reverse the vertical level
            ordering in the "Dev" dataset (in case "Dev" starts
            from the top of atmosphere instead of the surface).
            Default value: False

        log_color_scale: boolean         
            Set this flag to True to enable plotting data (not diffs)
            on a log color scale.
            Default value: False

         sigdiff_files : list of str
            Filenames that will contain the lists of species having
            significant differences in the 'sfc', '500hpa', and
            'zonalmean' plots.  These lists are needed in order to
            fill out the benchmark approval forms.
            Default value: None

    Remarks:
    --------
        (1) If both plot_by_spc_cat and plot_by_hco_cat are
            False, then all emission plots will be placed into the
            same PDF file.

        (2) Emissions that are 3-dimensional will be plotted as
            column sums.
    """
    # =================================================================
    # Initialization and data read
    # =================================================================

    # Create destination folder if it does not exist
    if os.path.isdir(dst) and not overwrite:
        print(
            "Directory {} exists. Pass overwrite=True to overwrite files in that directory, if any.".format(
                dst
            )
        )
        return
    elif not os.path.isdir(dst):
        os.mkdir(dst)

    # Create the "Emissions" category folder.  If subdst is passed,
    # then create a sub-folder (needed for the 1-year benchmarks).
    emisdir = os.path.join(dst, "Emissions")
    if not os.path.isdir(emisdir):
        os.mkdir(emisdir)
    if subdst is not None:
        emisdir = os.path.join(emisdir, subdst)
        if not os.path.isdir(emisdir):
            os.mkdir(emisdir)
        extra_title_txt = subdst
    else:
        extra_title_txt = None

    # Ref dataset
    try:
        refds = xr.open_dataset(ref, drop_variables=skip_these_vars)
    except FileNotFoundError:
        raise FileNotFoundError("Could not find Ref file: {}".format(ref))

    # Dev dataset
    try:
        devds = xr.open_dataset(dev, drop_variables=skip_these_vars)
    except FileNotFoundError:
        raise FileNotFoundError("Could not find Dev file: {}".format(dev))

    # Make sure that Ref and Dev datasets have the same variables.
    # Variables that are in Ref but not in Dev will be added to Dev
    # with all missing values (NaNs). And vice-versa.
    [refds, devds] = add_missing_variables(refds, devds)

    # Combine 2D and 3D variables into an overall list
    quiet = not verbose
    vardict = core.compare_varnames(refds, devds, quiet=quiet)
    vars2D = vardict["commonvars2D"]
    vars3D = vardict["commonvars3D"]
    varlist = vars2D + vars3D

    # ==================================================================
    # Compute column sums for 3D emissions
    # Make sure not to clobber the DataArray attributes
    # ==================================================================
    with xr.set_options(keep_attrs=True):
        for v in vars3D:
            if "lev" in refds[v].dims:
                refds[v] = refds[v].sum(dim="lev")
            if "lev" in devds[v].dims:
                devds[v] = devds[v].sum(dim="lev")

    # ==================================================================
    # If inputs plot_by* are both false, plot all emissions in same file
    # ==================================================================
    if not plot_by_spc_cat and not plot_by_hco_cat:
        if subdst is not None:
            pdfname = os.path.join(emisdir, "Emissions_{}.pdf".format(subdst))
        else:
            pdfname = os.path.join(emisdir, "Emissions.pdf")

        compare_single_level(
            refds,
            refstr,
            devds,
            devstr,
            varlist=varlist,
            pdfname=pdfname,
            log_color_scale=log_color_scale,
            extra_title_txt=extra_title_txt,
        )
        add_bookmarks_to_pdf(pdfname, varlist, remove_prefix="Emis", verbose=verbose)
        return

    # Get emissions variables (non-inventory), categories, and species
    emis_vars = [v for v in varlist if v[:4] == "Emis"]
    emis_cats = sorted(set([v.split("_")[1] for v in emis_vars]))
    emis_spc = sorted(set([v.split("_")[0][4:] for v in emis_vars]))

    # This is fixed in 12.3.2, comment out for now (bmy, 5/1/19)
    #    # Handle Bioburn and BioBurn as same categories (temporary until 12.3.1)
    #    emis_cats.remove('BioBurn')

    # Sort alphabetically (assume English characters)
    emis_vars.sort(key=str.lower)

    # ==================================================================
    # if plot_by_hco_cat is true, make a file for each HEMCO emissions
    # category that is in the diagnostics file
    #
    # Also write the list of emission quantities that have significant
    # diffs.  We'll need that to fill out the benchmark forms.
    # ==================================================================

    if plot_by_hco_cat:
        emisspcdir = os.path.join(dst, "Emissions")
        if not os.path.isdir(emisspcdir):
            os.mkdir(emisspcdir)
        if subdst is not None:
            emisspcdir = os.path.join(emisspcdir, subdst)
            if not os.path.isdir(emisspcdir):
                os.mkdir(emisspcdir)

        # for c in emis_cats:
        def createfile_hco_cat(c):
            # Handle cases of bioburn and bioBurn (temporary until 12.3.1)
            if c == "Bioburn":
                varnames = [
                    k for k in emis_vars if any(b in k for b in ["Bioburn", "BioBurn"])
                ]
            else:
                varnames = [k for k in emis_vars if c in k]

            # Create the PDF name.  If subdst is passed, then also add
            # subdst to the file name (e.g. as for 1-year benchmarks).
            if subdst is not None:
                pdfname = os.path.join(
                    emisspcdir, "{}_Emissions_{}.pdf".format(c, subdst)
                )
            else:
                pdfname = os.path.join(emisspcdir, "{}_Emissions.pdf".format(c))
            diff_dict = {}
            diff_emis = []
            compare_single_level(
                refds,
                refstr,
                devds,
                devstr,
                varlist=varnames,
                ilev=0,
                pdfname=pdfname,
                log_color_scale=log_color_scale,
                extra_title_txt=extra_title_txt,
                sigdiff_list=diff_emis,
            )

            add_bookmarks_to_pdf(
                pdfname, varnames, remove_prefix="Emis", verbose=verbose
            )
            # Save the list of quantities with significant differences for
            # this category into the diff_dict dictionary for use below
            diff_emis[:] = [v.replace("Emis", "") for v in diff_emis]
            diff_emis[:] = [v.replace("_" + c, "") for v in diff_emis]
            diff_dict[c] = diff_emis
            return diff_dict

        results = Parallel(n_jobs=n_job)(delayed(createfile_hco_cat)(c) for c in emis_cats)

        dict_emis = {list(result.keys())[0] : result[list(result.keys())[0]] for result in results}
        
        # =============================================================
        # Write the list of species having significant differences,
        # which we need to fill out the benchmark approval forms.
        # =============================================================
        if sigdiff_files != None:
            for filename in sigdiff_files:
                if "emis" in filename:
                    with open(filename, "w+") as f:
                        for c, diff_list in dict_emis.items():
                            print("* {}: ".format(c), file=f, end="")
                            for v in diff_list:
                                print("{} ".format(v), file=f, end="")
                            print(file=f)
                        f.close()

    # ==================================================================
    # if plot_by_spc_cat is true, make a file for each benchmark
    # species category with emissions in the diagnostics file
    # ==================================================================
    if plot_by_spc_cat:

        catdict = get_species_categories()
        warninglist = (
            []
        )  # in case any emissions are skipped (for use in nested pdf bookmarks)
        allcatspc = (
            []
        )  # for checking if emissions species not defined in benchmark category file
        emisdict = {}  # used for nested pdf bookmarks
        # for i, filecat in enumerate(catdict):
        def createfile_bench_cat(filecat):
            # Get emissions for species in this benchmark category
            varlist = []
            emisdict[filecat] = {}
            for subcat in catdict[filecat]:
                for spc in catdict[filecat][subcat]:
                    allcatspc.append(spc)
                    if spc in emis_spc:
                        emisdict[filecat][spc] = []
                        emisvars = [v for v in emis_vars if spc == v.split("_")[0][4:]]
                        for var in emisvars:
                            emisdict[filecat][spc].append(var.replace("Emis", ""))
                            varlist.append(var)
            if not varlist:
                print(
                    "\nWarning: no emissions species in benchmark species category {}".format(
                        filecat
                    )
                )
                return

            # Use same directory structure as for concentration plots
            catdir = os.path.join(dst, filecat)
            if not os.path.isdir(catdir):
                os.mkdir(catdir)
            if subdst is not None:
                catdir = os.path.join(catdir, subdst)
                if not os.path.isdir(catdir):
                    os.mkdir(catdir)

            # Create emissions file for this benchmark species category
            # If subdst is passed, add it to the pdf name (e.g. as
            # is needed for the 1-year benchmarks).
            if subdst is not None:
                pdfname = os.path.join(
                    catdir, "{}_Emissions_{}.pdf".format(filecat, subdst)
                )
            else:
                pdfname = os.path.join(catdir, "{}_Emissions.pdf".format(filecat))
            print(pdfname)
            # Create the PDF
            compare_single_level(
                refds,
                refstr,
                devds,
                devstr,
                varlist=varlist,
                ilev=0,
                pdfname=pdfname,
                flip_ref=flip_ref,
                flip_dev=flip_dev,
                log_color_scale=log_color_scale,
                extra_title_txt=extra_title_txt,
            )
            add_nested_bookmarks_to_pdf(pdfname, filecat, emisdict, warninglist)

        Parallel(n_jobs=n_job)(
            delayed(createfile_bench_cat)(filecat) for i, filecat in enumerate(catdict)
        )

        # Give warning if emissions species is not assigned a benchmark category
        for spc in emis_spc:
            if spc not in allcatspc:
                print(
                    "Warning: species {} has emissions diagnostics but is not in benchmark_categories.yml".format(
                        spc
                    )
                )


def make_benchmark_emis_tables(
    reflist,
    refstr,
    devlist,
    devstr,
    dst="./1mo_benchmark",
    overwrite=False,
    interval=[2678400.0],
):
    """
    Creates a text file containing emission totals by species and
    category for benchmarking purposes.

    Args:
    -----
        reflist: list of str
             List with the path names of the emissions file, or emissions
             and met field files, that will constitute the "Ref" (aka 
             "Reference") data set. If two files are passed in the list,
             the met field file must be second.

        refstr : str
            A string to describe ref (e.g. version number)

        devlist : list of str
             List with the path names of the emissions file, or emissions
             and met field files, that will constitute the "Dev" (aka
             "Development") data set. If two files are passed in the list,
             the met field file must be second. The "Dev" data set will
             be compared against the "Ref" data set.

        devstr : str
            A string to describe dev (e.g. version number)

    Keyword Args (optional):
    ------------------------
        dst : str
            A string denoting the destination folder where the file
            containing emissions totals will be written.
            Default value: ./1mo_benchmark

        overwrite : boolean
            Set this flag to True to overwrite files in the
            destination folder (specified by the dst argument).
            Default value : False

        interval : list of float
            The length of the data interval in seconds. By default, interval
            is set to [2678400.0], which is the number of seconds in July
            (our 1-month benchmarking month).
    """

    # ==================================================================
    # Initialization
    # ==================================================================

    # Create destination folder
    if os.path.isdir(dst) and not overwrite:
        print(
            "Directory {} exists. Pass overwrite=True to overwrite files in that directory, if any.".format(
                dst
            )
        )
        return
    elif not os.path.isdir(dst):
        os.mkdir(dst)

    # Create the "Emissions" category folder if it does not exist
    emisdir = os.path.join(dst, "Emissions")
    if not os.path.isdir(emisdir):
        os.mkdir(emisdir)

    # ==================================================================
    # Read data from netCDF into Dataset objects
    # ==================================================================

    # Read the Ref dataset and make sure that the area variables are present
    if len(reflist) == 1:
        reflist = [reflist]
    refds = xr.open_mfdataset(reflist, drop_variables=skip_these_vars)
    refds = core.check_for_area(refds)

    # Read the Dev dataset and make sure that area variables are present
    if len(devlist) == 1:
        devlist = [devlist]
    devds = xr.open_mfdataset(devlist, drop_variables=skip_these_vars)
    devds = core.check_for_area(devds)

    # ==================================================================
    # Create table of emissions
    # ==================================================================

    # Emissions species dictionary
    species = yaml_load_file(
        open(os.path.join(os.path.dirname(__file__), emission_spc))
    )
    inventories = yaml_load_file(
        open(os.path.join(os.path.dirname(__file__), emission_inv))
    )

    # Destination files
    file_emis_totals = os.path.join(emisdir, "Emission_totals.txt")
    file_inv_totals = os.path.join(emisdir, "Inventory_totals.txt")

    # Create table of emissions by species
    create_total_emissions_table(
        refds,
        refstr,
        devds,
        devstr,
        species,
        file_emis_totals,
        interval,
        template="Emis{}_",
    )

    # Create table of emissions by inventory
    create_total_emissions_table(
        refds,
        refstr,
        devds,
        devstr,
        inventories,
        file_inv_totals,
        interval,
        template="Inv{}_",
    )


def make_benchmark_jvalue_plots(
    ref,
    refstr,
    dev,
    devstr,
    varlist=None,
    dst="./1mo_benchmark",
    subdst=None,
    local_noon_jvalues=False,
    plots=["sfc", "500hpa", "zonalmean"],
    overwrite=False,
    verbose=False,
    flip_ref=False,
    flip_dev=False,
    log_color_scale=False,
    sigdiff_files=None,
):
    """
    Creates PDF files containing plots of J-values for model
    benchmarking purposes.

    Args:
    -----
        ref: str
            Path name for the "Ref" (aka "Reference") data set.

        refstr : str
            A string to describe ref (e.g. version number)
    
        dev : str
            Path name for the "Dev" (aka "Development") data set.
            This data set will be compared against the "Reference"
            data set.

        devstr : str
            A string to describe dev (e.g. version number)
    
    Keyword Args (optional):
    ------------------------
        varlist : list of str
            List of J-value variables to plot.  If not passed, 
            then all J-value variables common to both dev 
            and ref will be plotted.  The varlist argument can be
            a useful way of restricting the number of variables
            plotted to the pdf file when debugging.
            Default value: None

        dst : str
            A string denoting the destination folder where a
            PDF file  containing plots will be written.
            Default value: ./1mo_benchmark.

        subdst : str
            A string denoting the sub-directory of dst where PDF
            files containing plots will be written.  In practice,
            subdst is only needed for the 1-year benchmark output,
            and denotes a date string (such as "Jan2016") that
            corresponds to the month that is being plotted.
            Default value: None

        local_noon_jvalues : boolean
            Set this flag to plot local noon J-values.  This will
            divide all J-value variables by the JNoonFrac counter,
            which is the fraction of the time that it was local noon
            at each location.
            Default value : False

        plots : list of strings
            List of plot types to create.
            Default value: ['sfc', '500hpa', 'zonalmean']

        overwrite : boolean
            Set this flag to True to overwrite files in the
            destination folder (specified by the dst argument).
            Default value: False.

        verbose : boolean
            Set this flag to True to print extra informational output.
            Default value: False

        flip_ref : boolean
            Set this flag to True to reverse the vertical level
            ordering in the "Ref" dataset (in case "Ref" starts
            from the top of atmosphere instead of the surface).
            Default value: False

        flip_dev : boolean
            Set this flag to True to reverse the vertical level
            ordering in the "Dev" dataset (in case "Dev" starts
            from the top of atmosphere instead of the surface).
            Default value: False

        log_color_scale: boolean         
            Set this flag to True if you wish to enable plotting data 
            (not diffs) on a log color scale.
            Default value: False

        sigdiff_files : list of str
            Filenames that will contain the lists of J-values having
            significant differences in the 'sfc', '500hpa', and
            'zonalmean' plots.  These lists are needed in order to
            fill out the benchmark approval forms.
            Default value: None

    Remarks:
    --------
         Will create 4 files containing J-value plots:
            (1 ) Surface values
            (2 ) 500 hPa values
            (3a) Full-column zonal mean values.
            (3b) Stratospheric zonal mean values
         These can be toggled on/off with the plots keyword argument.

         At present, we do not yet have the capability to split the
         plots up into separate files per category (e.g. Oxidants,
         Aerosols, etc.).  This is primarily due to the fact that 
         we archive J-values from GEOS-Chem for individual species
         but not family species.  We could attempt to add this 
         functionality later if there is sufficient demand. 
    """

    # ==================================================================
    # Initialization
    # ==================================================================

    # Create the destination folder if it does not exist
    if os.path.isdir(dst) and not overwrite:
        print(
            "Directory {} exists. Pass overwrite=True to overwrite files in tht directory, if any.".format(
                dst
            )
        )
        return
    elif not os.path.isdir(dst):
        os.mkdir(dst)

    # Ref dataset
    try:
        refds = xr.open_dataset(ref, drop_variables=skip_these_vars)
    except FileNotFoundError:
        raise FileNotFoundError("Could not find Ref file: {}".format(ref))

    # Dev dataset
    try:
        devds = xr.open_dataset(dev, drop_variables=skip_these_vars)
    except FileNotFoundError:
        raise FileNotFoundError("Could not find Dev file: {}".format(dev))

    # Make sure that Ref and Dev datasets have the same variables.
    # Variables that are in Ref but not in Dev will be added to Dev
    # with all missing values (NaNs). And vice-versa.
    [refds, devds] = add_missing_variables(refds, devds)

    # Get a list of the 3D variables in both datasets
    if varlist == None:
        quiet = not verbose
        vardict = core.compare_varnames(refds, devds, quiet=quiet)
        cmn = vardict["commonvars3D"]

    # ==================================================================
    # Local noon or continuously-averaged J-values?
    # ==================================================================
    if local_noon_jvalues:

        # Get a list of local noon J-value variables
        # (or use the varlist passed via tha argument list)
        prefix = "JNoon_"
        if varlist == None:
            varlist = [v for v in cmn if prefix in v]

        # Make sure JNoonFrac (fraction of times it was local noon
        # in each column) is present in both Ref and Dev datasets
        if not "JNoonFrac" in cmn:
            msg = "JNoonFrac is not common to Ref and Dev datasets!"
            raise ValueError(msg)

        # JNoon_* are cumulative sums of local noon J-values; we need
        # to divide these by JNoonFrac to get the average value
        refds = core.divide_dataset_by_dataarray(refds, refds["JNoonFrac"], varlist)
        devds = core.divide_dataset_by_dataarray(devds, devds["JNoonFrac"], varlist)

        # Subfolder of dst where PDF files will be printed
        catdir = "JValuesLocalNoon"

    else:

        # Get a list of continuously averaged J-value variables
        # (or use the varlist passed via tha argument list)
        prefix = "Jval_"
        if varlist == None:
            varlist = [v for v in cmn if prefix in v]

        # Subfolder of dst where PDF files will be printed
        catdir = "JValues"

    # ==================================================================
    # Create the plots
    # ==================================================================

    # Make the output folder if it doesn't exist.  If subdst is passed,
    # then create a sub-folder of this directory (e.g. which is needed
    # for the 1-year benchmarks)
    jvdir = os.path.join(dst, catdir)
    if not os.path.isdir(jvdir):
        os.mkdir(jvdir)
    if subdst is not None:
        jvdir = os.path.join(jvdir, subdst)
        if not os.path.isdir(jvdir):
            os.mkdir(jvdir)
        extra_title_txt = subdst
    else:
        extra_title_txt = None

    # Surface plots
    if "sfc" in plots:
        if subdst is not None:
            pdfname = os.path.join(jvdir, "{}Surface_{}.pdf".format(prefix, subdst))
        else:
            pdfname = os.path.join(jvdir, "{}Surface.pdf".format(prefix))

        diff_sfc = []
        compare_single_level(
            refds,
            refstr,
            devds,
            devstr,
            varlist=varlist,
            ilev=0,
            pdfname=pdfname,
            flip_ref=flip_ref,
            flip_dev=flip_dev,
            log_color_scale=log_color_scale,
            extra_title_txt=extra_title_txt,
            sigdiff_list=diff_sfc,
        )
        diff_sfc[:] = [v.replace(prefix, "") for v in diff_sfc]
        add_bookmarks_to_pdf(pdfname, varlist, remove_prefix=prefix, verbose=verbose)

    # 500hPa plots
    if "500hpa" in plots:
        if subdst is not None:
            pdfname = os.path.join(jvdir, "{}500hPa_{}.pdf".format(prefix, subdst))
        else:
            pdfname = os.path.join(jvdir, "{}500hPa.pdf".format(prefix))

        diff_500 = []
        compare_single_level(
            refds,
            refstr,
            devds,
            devstr,
            varlist=varlist,
            ilev=22,
            pdfname=pdfname,
            flip_ref=flip_ref,
            flip_dev=flip_dev,
            log_color_scale=log_color_scale,
            extra_title_txt=extra_title_txt,
            sigdiff_list=diff_500,
        )
        diff_500[:] = [v.replace(prefix, "") for v in diff_500]
        add_bookmarks_to_pdf(pdfname, varlist, remove_prefix=prefix, verbose=verbose)

    # Full-column zonal mean plots
    if "zonalmean" in plots:
        if subdst is not None:
            pdfname = os.path.join(
                jvdir, "{}FullColumn_ZonalMean_{}.pdf".format(prefix, subdst)
            )
        else:
            pdfname = os.path.join(jvdir, "{}FullColumn_ZonalMean.pdf".format(prefix))

        diff_zm = []
        compare_zonal_mean(
            refds,
            refstr,
            devds,
            devstr,
            varlist=varlist,
            pdfname=pdfname,
            flip_ref=flip_ref,
            flip_dev=flip_dev,
            log_color_scale=log_color_scale,
            extra_title_txt=extra_title_txt,
            sigdiff_list=diff_zm,
        )
        diff_zm[:] = [v.replace(prefix, "") for v in diff_zm]
        add_bookmarks_to_pdf(pdfname, varlist, remove_prefix=prefix, verbose=verbose)

        # Strat_ZonalMean plots will use a log-pressure Y-axis, with
        # a range of 1..100 hPa, as per GCSC request. (bmy, 8/13/19)
        if subdst is not None:
            pdfname = os.path.join(
                jvdir, "{}Strat_ZonalMean_{}.pdf".format(prefix, subdst)
            )
        else:
            pdfname = os.path.join(jvdir, "{}Strat_ZonalMean.pdf".format(prefix))

        compare_zonal_mean(
            refds,
            refstr,
            devds,
            devstr,
            varlist=varlist,
            pdfname=pdfname,
            pres_range=[0, 100],
            log_yaxis=True,
            flip_ref=flip_ref,
            flip_dev=flip_dev,
            extra_title_txt=extra_title_txt,
            log_color_scale=log_color_scale,
        )
        add_bookmarks_to_pdf(pdfname, varlist, remove_prefix=prefix, verbose=verbose)

        # ==============================================================
        # Write the lists of J-values that have significant differences,
        # which we need to fill out the benchmark approval forms.
        # ==============================================================
        if sigdiff_files != None:
            for filename in sigdiff_files:
                if "sfc" in plots:
                    if "sfc" in filename:
                        with open(filename, "a+") as f:
                            print("* J-Values: ", file=f, end="")
                            for v in diff_sfc:
                                print("{} ".format(v), file=f, end="")
                            print(file=f)
                            f.close()

                if "500" in plots:
                    if "500" in filename:
                        with open(filename, "a+") as f:
                            print("* J-Values: ", file=f, end="")
                            for v in diff_500:
                                print("{} ".format(v), file=f, end="")
                            print(file=f)
                            f.close()

                if "zonalmean" in plots or "zm" in plots:
                    if "zonalmean" in filename or "zm" in filename:
                        with open(filename, "a+") as f:
                            print("* J-Values: ", file=f, end="")
                            for v in diff_zm:
                                print("{} ".format(v), file=f, end="")
                            print(file=f)
                            f.close()


def make_benchmark_aod_plots(
    ref,
    refstr,
    dev,
    devstr,
    varlist=None,
    dst="./1mo_benchmark",
    subdst=None,
    overwrite=False,
    verbose=False,
    log_color_scale=False,
    sigdiff_files=None,
):
    """
    Creates PDF files containing plots of column aerosol optical
    depths (AODs) for model benchmarking purposes.

    Args:
    -----
        ref: str
            Path name for the "Ref" (aka "Reference") data set.

        refstr : str
            A string to describe ref (e.g. version number)

        dev : str
            Path name for the "Dev" (aka "Development") data set.
            This data set will be compared against the "Reference"
            data set.

        devstr : str
            A string to describe dev (e.g. version number)

    Keyword Args (optional):
    ------------------------
        varlist : list of str
            List of AOD variables to plot.  If not passed, then all
            AOD variables common to both Dev and Ref will be plotted.
            Use the varlist argument to restrict the number of
            variables plotted to the pdf file when debugging.
            Default value: None

        dst : str
            A string denoting the destination folder where a
            PDF file  containing plots will be written.
            Default value: ./1mo_benchmark.

        subdst : str
            A string denoting the sub-directory of dst where PDF
            files containing plots will be written.  In practice,
            subdst is only needed for the 1-year benchmark output,
            and denotes a date string (such as "Jan2016") that
            corresponds to the month that is being plotted.
            Default value: None

        overwrite : boolean
            Set this flag to True to overwrite files in the
            destination folder (specified by the dst argument).
            Default value: False.

        verbose : boolean
            Set this flag to True to print extra informational output.
            Default value: False

        log_color_scale: boolean         
            Set this flag to True to enable plotting data (not diffs)
            on a log color scale.
            Default value: False

        sigdiff_files : list of str
            Filenames that will contain the list of quantities having
            having significant differences in the column AOD plots.
            These lists are needed in order to fill out the benchmark
            approval forms.
            Default value: None
    """
    # ==================================================================
    # Initialization and also read data
    # ==================================================================

    # Create the destination directory if it does not exist
    if os.path.isdir(dst) and not overwrite:
        print(
            "Directory {} exists. Pass overwrite=True to overwrite files in tht directory, if any.".format(
                dst
            )
        )
        return
    elif not os.path.isdir(dst):
        os.mkdir(dst)

    # Create the "Aerosols" directory as a subfolder of dst.
    # If subdst is passed, then create a subdirectory of the "Aerosols"
    # directory (e.g. which is needed for the 1-year benchmarks).
    aoddir = os.path.join(dst, "Aerosols")
    if not os.path.isdir(aoddir):
        os.mkdir(aoddir)
    if subdst is not None:
        aoddir = os.path.join(aoddir, subdst)
        if not os.path.isdir(aoddir):
            os.mkdir(aoddir)
        extra_title_txt = subdst
    else:
        extra_title_txt = None

    # Read the Ref dataset
    try:
        refds = xr.open_dataset(ref, drop_variables=skip_these_vars)
    except FileNotFoundError:
        raise FileNotFoundError("Could not find Ref file: {}".format(ref))

    # Read the Dev dataset
    try:
        devds = xr.open_dataset(dev, drop_variables=skip_these_vars)
    except FileNotFoundError:
        raise FileNotFoundError("Could not find Dev file: {}".format(dev))

    # NOTE: GCHP diagnostic variable exports are defined before the
    # input.geos file is read.  This means "WL1" will not have been
    # replaced with "550nm" in the variable names.  Do this string
    # replace operation here, so that we can compare GCC and GCHP
    # data directly. (bmy, 4/29/19)
    with xr.set_options(keep_attrs=True):

        # Rename variables in the Ref dataset
        old2new = {}
        for v in refds.data_vars.keys():
            if "WL1" in v:
                newname = v.replace("WL1", "550nm")
                old2new[v] = newname
        refds = refds.rename(old2new)

        # Rename variables in the Dev dataset
        old2new = {}
        for v in devds.data_vars.keys():
            if "WL1" in v:
                newname = v.replace("WL1", "550nm")
                old2new[v] = newname
        devds = devds.rename(old2new)

    # Make sure that Ref and Dev datasets have the same variables.
    # Variables that are in Ref but not in Dev will be added to Dev
    # with all missing values (NaNs). And vice-versa.
    [refds, devds] = add_missing_variables(refds, devds)

    # Find common AOD variables in both datasets
    # (or use the varlist passed via keyword argument)
    if varlist == None:
        quiet = not verbose
        vardict = core.compare_varnames(refds, devds, quiet)
        cmn3D = vardict["commonvars3D"]
        varlist = [v for v in cmn3D if "AOD" in v and "_bin" not in v]

    # Dictionary and list for new display names
    newvars = yaml_load_file(open(os.path.join(os.path.dirname(__file__), aod_spc)))
    newvarlist = []

    # ==================================================================
    # Compute the total AOD by summing over the constituent members
    # ==================================================================

    # Take one of the variables so we can use its dims, coords,
    # attrs to create the DataArray object for total AOD
    v = varlist[0]

    # Create a DataArray to hold total column AOD
    # This is the same shape as the DataArray objects in refds
    reftot = xr.DataArray(
        np.zeros(refds[v].values.shape),
        name="AODTotal",
        dims=refds[v].dims,
        coords=refds[v].coords,
        attrs=refds[v].attrs,
    )

    # Create a DataArray to hold total column AOD
    # This is the same shape as the DataArray objects in devds
    devtot = xr.DataArray(
        np.zeros(devds[v].values.shape),
        name="AODTotal",
        dims=devds[v].dims,
        coords=devds[v].coords,
        attrs=devds[v].attrs,
    )

    # Save the variable attributes so that we can reattach them
    refattrs = reftot.attrs
    devattrs = devtot.attrs

    # Compute the sum of all AOD variables
    for v in varlist:
        reftot = reftot + refds[v]
        devtot = devtot + devds[v]

    # Reattach the variable attributes
    reftot.name = "AODTotal"
    reftot.attrs = refattrs
    reftot.attrs["long_name"] = "Total aerosol optical depth"
    devtot.name = "AODTotal"
    devtot.attrs = devattrs
    devtot.attrs["long_name"] = "Total aerosol optical depth"

    # Merge these variables back into the dataset
    refds = xr.merge([refds, reftot])
    devds = xr.merge([devds, devtot])

    # Also add AODTotal to the list
    varlist.append("AODTotal")

    # ==================================================================
    # Compute column AODs
    # Create a new DataArray for each column AOD variable,
    # using the new display name, and preserving attributes.
    # Merge the new DataArrays back into the DataSets.
    # ==================================================================
    for v in varlist:

        # Get the new name for each AOD variable (it's easier to display)
        if v in newvars:
            newname = newvars[v]
            newvarlist.append(newname)
        else:
            raise ValueError("Could not find a display name for {}".format(v))

        # Don't clobber existing DataArray and Dataset attributes
        with xr.set_options(keep_attrs=True):

            # Add column AOD of newname to Ref
            array = refds[v].sum(dim="lev")
            array.name = newname
            array.attrs["units"] = "1"
            refds = xr.merge([refds, array])

            # Add column AOD of newname to Dev
            array = devds[v].sum(dim="lev")
            array.name = newname
            array.attrs["units"] = "1"
            devds = xr.merge([devds, array])

    # ==================================================================
    # Create the plots
    # ==================================================================
    if subdst is not None:
        pdfname = os.path.join(aoddir, "Aerosols_ColumnOptDepth_{}.pdf".format(subdst))
    else:
        pdfname = os.path.join(aoddir, "Aerosols_ColumnOptDepth.pdf")

    diff_aod = []
    compare_single_level(
        refds,
        refstr,
        devds,
        devstr,
        varlist=newvarlist,
        ilev=0,
        pdfname=pdfname,
        log_color_scale=log_color_scale,
        extra_title_txt=extra_title_txt,
        sigdiff_list=diff_aod,
    )
    diff_aod[:] = [v.replace("Column_AOD_", "") for v in diff_aod]
    add_bookmarks_to_pdf(
        pdfname, newvarlist, remove_prefix="Column_AOD_", verbose=verbose
    )

    # ==================================================================
    # Write the list of AOD quantities having significant differences,
    # which we will need to fill out the benchmark forms.
    # ==================================================================
    if sigdiff_files != None:
        for filename in sigdiff_files:
            if "sfc" in filename:
                with open(filename, "a+") as f:
                    print("* Column AOD: ", file=f, end="")
                    for v in diff_aod:
                        print("{} ".format(v), file=f, end="")
                    print(file=f)
                    f.close()


def make_benchmark_mass_tables(
    reflist,
    refstr,
    devlist,
    devstr,
    varlist=None,
    dst="./1mo_benchmark",
    subdst=None,
    overwrite=False,
    verbose=False,
    label="at end of simulation",
):
    """
    Creates a text file containing global mass totals by species and
    category for benchmarking purposes.

    Args:
    -----
        reflist : list of str
            List of files (i.e. pathnames) that will constitute
            the "Ref" (aka "Reference") data set.

        refstr : str
            A string to describe ref (e.g. version number)

        dev : list of str
            List of files (i.e. pathnames) that will constitute
            the "Dev" (aka "Development") data set.  The "Dev"
            data set will be compared against the "Ref" data set.

        devstr : str
            A string to describe dev (e.g. version number)

    Keyword Args (optional):
    ------------------------
        varlist : list of str
            List of variables to include in the list of totals.
            If omitted, then all variables that are found in either
            "Ref" or "Dev" will be included.  The varlist argument
            can be a useful way of reducing the number of
            variables during debugging and testing.
            Default value: None

        dst : str
            A string denoting the destination folder where the file
            containing emissions totals will be written.
            Default value: ./1mo_benchmark

        subdst : str
            A string denoting the sub-directory of dst where PDF
            files containing plots will be written.  In practice,
            subdst is only needed for the 1-year benchmark output,
            and denotes a date string (such as "Jan2016") that
            corresponds to the month that is being plotted.
            Default value: None

        overwrite : boolean
            Set this flag to True to overwrite files in the
            destination folder (specified by the dst argument).
            Default value : False

        verbose : boolean
            Set this flag to True to print extra informational output.
            Default value: False.
    """

    # ==================================================================
    # Define destination directory
    # ==================================================================
    if os.path.isdir(dst) and not overwrite:
        print("Directory {} exists. Pass overwrite=True to overwrite files in that directory, if any.".format(dst))
        return
    elif not os.path.isdir(dst):
        os.makedirs(dst)

    # ==================================================================
    # Read data from netCDF into Dataset objects
    # ==================================================================

    # Ref
    try:
        refds = xr.open_mfdataset(reflist, drop_variables=skip_these_vars)
    except FileNotFoundError:
        raise FileNotFoundError("Error opening Ref files: {}".format(reflist))

    # Dev dataset
    try:
        devds = xr.open_mfdataset(devlist, drop_variables=skip_these_vars)
    except FileNotFoundError:
        raise FileNotFoundError("Error opening Dev files: {}!".format(devlist))

    # ==================================================================
    # Make sure that all necessary meteorological variables are found
    # ==================================================================

    # Find the area variables in Ref and Dev
    ref_area = core.get_area_from_dataset(refds)
    dev_area = core.get_area_from_dataset(devds)

    # Find required meteorological variables in Ref
    # (or exit with an error if we can't find them)
    metvar_list = ["Met_DELPDRY", "Met_BXHEIGHT", "Met_TropLev"]
    refmet = core.get_variables_from_dataset(refds, metvar_list)
    devmet = core.get_variables_from_dataset(devds, metvar_list)

    # ==================================================================
    # Make sure that all necessary species are found
    # ==================================================================

    # If varlist has not been passed as an argument, then use all
    # species concentration variables that are present in Ref or Dev.
    # For species that are in Ref but not in Dev, create a
    # variable of missing values (NaNs). in Dev.  Ditto for Ref.
    if varlist is None:
        [refds, devds] = add_missing_variables(refds, devds)
        quiet = not verbose
        vardict = core.compare_varnames(refds, devds, quiet=quiet)
        varlist = vardict["commonvars3D"]

    # Only use the species concentration variables in the restart
    # file, as we will pass the meteorology variables separately
    # to routine create_global_mass_table.
    varlist = [v for v in varlist if "SpeciesRst_" in v]
    varlist.sort()

    # ==================================================================
    # Create the mask arrays for the troposphere for Ref and Dev
    # ==================================================================
    ref_tropmask = get_troposphere_mask(refmet)
    dev_tropmask = get_troposphere_mask(devmet)

    # ==================================================================
    # Create a dictionary to hold all of the meterological
    # variables and mask variables that we need to pass down
    # ==================================================================
    met_and_masks = {
        "Ref_Area": ref_area,
        "Dev_Area": dev_area,
        "Ref_Delta_P": refmet["Met_DELPDRY"],
        "Dev_Delta_P": devmet["Met_DELPDRY"],
        "Ref_BxHeight": refmet["Met_BXHEIGHT"],
        "Dev_BxHeight": devmet["Met_BXHEIGHT"],
        #                     'Ref_AirVol'   : refmet['Met_AIRVOL'],
        #                     'Dev_AirVol'   : devmet['Met_AIRVOL'],
        "Ref_TropMask": ref_tropmask,
        "Dev_TropMask": dev_tropmask,
    }

    # ==================================================================
    # Create global mass table
    # ==================================================================
    if subdst is not None:
        mass_file = os.path.join(
            dst, "{}_GlobalMass_TropStrat_{}.txt".format(devstr, subdst)
        )
    else:
        mass_file = os.path.join(dst, "{}_GlobalMass_TropStrat.txt".format(devstr))

    create_global_mass_table(
        refds,
        refstr,
        devds,
        devstr,
        varlist,
        met_and_masks,
        label,
        outfilename=mass_file,
        verbose=verbose,
    )

    # ==================================================================
    # Create tropospheric mass table
    # ==================================================================
    if subdst is not None:
        mass_file = os.path.join(
            dst, "{}_GlobalMass_Trop_{}.txt".format(devstr, subdst)
        )
    else:
        mass_file = os.path.join(dst, "{}_GlobalMass_Trop.txt".format(devstr))

    create_global_mass_table(
        refds,
        refstr,
        devds,
        devstr,
        varlist,
        met_and_masks,
        label,
        outfilename=mass_file,
        trop_only=True,
        verbose=verbose,
    )


def make_benchmark_budget_tables(
        dev, devstr, dst="./1mo_benchmark", overwrite=False, subdst=None, interval=[2678400.0], n_job=-1
):
    """
    Creates a text file containing budgets by species for benchmarking
    purposes.

    Args:
    -----
        dev : str
            Path names for the "Dev" (aka "Development") data set.

        devstr : str
            A string to describe dev (e.g. version number)

    Keyword Args (optional):
    ------------------------
        dst : str
            A string denoting the destination folder where the file
            containing emissions totals will be written.
            Default value: ./1mo_benchmark

        overwrite : boolean
            Set this flag to True to overwrite files in the
            destination folder (specified by the dst argument).
            Default value : False

        interval : float
            Specifies the averaging period in seconds, which is used
            to convert fluxes (e.g. kg/m2/s) to masses (e.g kg).
            Default value : None

        subdst : str
            A string denoting the sub-directory of dst where PDF
            files containing plots will be written.  In practice,
            subdst is only needed for the 1-year benchmark output,
            and denotes a date string (such as "Jan2016") that
            corresponds to the month that is being plotted.
            Default value: None

    """

    # ==================================================================
    # Define destination directory
    # ==================================================================
    if os.path.isdir(dst) and not overwrite:
        print(
            "Directory {} exists. Pass overwrite=True to overwrite files in that directory, if any.".format(
                dst
            )
        )
        return
    elif not os.path.isdir(dst):
        os.mkdir(dst)
    budgetdir = os.path.join(dst, "Budget")
    if not os.path.isdir(budgetdir):
        os.mkdir(budgetdir)
    if subdst is not None:
        budgetdir = os.path.join(budgetdir, subdst)
        if not os.path.isdir(budgetdir):
            os.mkdir(budgetdir)

    # ==================================================================
    # Read data from netCDF into Dataset objects
    # ==================================================================

    # Dev
    try:
        devds = xr.open_dataset(dev, drop_variables=skip_these_vars)
    except FileNotFoundError:
        raise FileNotFoundError("Could not find Dev file: {}".format(dev))

    # ==================================================================
    # Create budget table
    # ==================================================================

    # If the averaging interval (in seconds) is not specified,
    # then assume July 2016 = 86400 seconds * 31 days
    #    if interval == None:
    #        interval = 86400.0 * 31.0

    # Get budget variable and regions
    budget_vars = [k for k in devds.data_vars.keys() if k[:6] == "Budget"]
    budget_regions = sorted(set([v.split("_")[0][-4:] for v in budget_vars]))

    # for region in budget_regions:
    def createfile(region):

        # Suppress harmless run-time warnings from all threads
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        warnings.filterwarnings("ignore", category=UserWarning)

        # Destination file
        file_budget = os.path.join(budgetdir, "Budget_" + region + ".txt")

        # Get variable names and species for this region
        region_vars = [k for k in budget_vars if region in k]
        region_spc = sorted(set([v.split("_")[1] for v in region_vars]))

        # Write to file
        create_budget_table(
            devds,
            devstr,
            region,
            region_spc,
            region_vars,
            file_budget,
            interval,
            template="Budget_{}",
        )

    # Create budget tables in parallel
    Parallel(n_jobs=n_job)(delayed(createfile)(region) for region in budget_regions)


def make_benchmark_oh_metrics(
    reflist,
    refstr,
    devlist,
    devstr,
    dst="./1mo_benchmark",
    overwrite=False,
):
    """
    Creates a text file containing metrics of global mean OH, MCF lifetime,
    and CH4 lifetime for benchmarking purposes.

    Args:
    -----
        reflist: list of str
            List with the path names of files that will constitute the
            "Ref" (aka "Reference") data set.

        refstr : str
            A string to describe ref (e.g. version number)

        devlist : list of str
            List with the path names of files that will constitute the
            "Dev" (aka "Development") data set.  The "Dev" data set will be
            compared against the "Ref" data set.

        devstr : str
            A string to describe dev (e.g. version number)

    Keyword Args (optional):
    ------------------------
        dst : str
            A string denoting the destination folder where the file
            containing emissions totals will be written.
            Default value: ./1mo_benchmark

        overwrite : boolean
            Set this flag to True to overwrite files in the
            destination folder (specified by the dst argument).
            Default value : False
    """

    # ==================================================================
    # Define destination directory
    # ==================================================================
    if os.path.isdir(dst) and not overwrite:
        print(
            "Directory {} exists. Pass overwrite=True to overwrite files in that directory, if any.".format(dst))
        return
    elif not os.path.isdir(dst):
        os.makedirs(dst)

    # ==================================================================
    # Read data from netCDF into Dataset objects
    # ==================================================================

    # Ref
    try:
        refds = xr.open_mfdataset(reflist, drop_variables=skip_these_vars)
    except FileNotFoundError:
        raise FileNotFoundError(
            "Could not find one of the Ref files: {}".format(reflist)
        )

    # Dev
    try:
        devds = xr.open_mfdataset(devlist, drop_variables=skip_these_vars)
    except FileNotFoundError:
        raise FileNotFoundError(
            "Could not find one of the Dev files: {}".format(devlist)
        )

    # Make sure that required variables are found
    if "OHconcAfterChem" not in refds.data_vars.keys():
        raise ValueError('Could not find "OHconcAfterChem" in Ref!')
    if "OHconcAfterChem" not in devds.data_vars.keys():
        raise ValueError('Could not find "OHconcAfterChem" in Dev!')

    # ==================================================================
    # Make sure that all necessary variables are found
    # ==================================================================

    # Find the area variables in Ref and Dev
    ref_area = core.get_area_from_dataset(refds)
    dev_area = core.get_area_from_dataset(devds)

    # Find required meteorological variables in Ref
    # (or exit with an error if we can't find them)
    metvar_list = [
        "Met_AD",
        "Met_AIRDEN",
        "Met_BXHEIGHT",
        "Met_T",
        "Met_TropLev",
        "FracOfTimeInTrop",
    ]
    refmet = core.get_variables_from_dataset(refds, metvar_list)
    devmet = core.get_variables_from_dataset(devds, metvar_list)

    # Create the mask arrays for the troposphere for Ref and Dev
    ref_tropmask = get_troposphere_mask(refmet)
    dev_tropmask = get_troposphere_mask(devmet)

    # Get the OH concentration
    ref_oh = refds["OHconcAfterChem"]
    dev_oh = devds["OHconcAfterChem"]

    # ==================================================================
    # Open file for output
    # ==================================================================

    # Create file
    outfilename = os.path.join(dst, "{}_OH_metrics.txt".format(devstr))
    try:
        f = open(outfilename, "w")
    except FileNotFoundError:
        raise FileNotFoundError("Could not open {} for writing!".format(
            outfilename))

    # ==================================================================
    # Compute mass-weighted OH in the troposphere
    # ==================================================================

    # Physical constants
    Avo = 6.022140857e23  # molec/mol
    mw_air = 28.97  # g/mole dry air
    g0 = 9.80665  # m/s2

    # Ref
    ref_oh_trop = np.ma.masked_array(ref_oh.values, ref_tropmask)
    ref_airmass_trop = np.ma.masked_array(refmet["Met_AD"].values, ref_tropmask)
    ref_oh_mass = ref_oh_trop * ref_airmass_trop
    ref_total_ohmass = np.sum(ref_oh_mass)
    ref_total_airmass = np.sum(ref_airmass_trop)
    ref_mean_oh = (ref_total_ohmass / ref_total_airmass) / 1e5

    # Dev
    dev_oh_trop = np.ma.masked_array(dev_oh.values, dev_tropmask)
    dev_airmass_trop = np.ma.masked_array(devmet["Met_AD"].values, dev_tropmask)
    dev_oh_mass = dev_oh_trop * dev_airmass_trop
    dev_total_ohmass = np.sum(dev_oh_mass)
    dev_total_airmass = np.sum(dev_airmass_trop)
    dev_mean_oh = (dev_total_ohmass / dev_total_airmass) / 1e5

    oh_diff = dev_mean_oh - ref_mean_oh
    oh_pctdiff = ((dev_mean_oh - ref_mean_oh) / ref_mean_oh) * 100.0

    # Title strings
    title1 = "### Global mass-weighted OH concentration [1e5 molec/cm3]"
    title2 = "### Ref = {}; Dev = {}".format(refstr, devstr)

    # Print header to file
    print("#" * 79, file=f)
    print("{}{}".format(title1.ljust(76), "###"), file=f)
    print("{}{}".format(title2.ljust(76), "###"), file=f)
    print("#" * 79, file=f)

    # Write results to file
    print(
        "{}{}{}{}".format(
            "  Ref".ljust(15),
            "Dev".ljust(13),
            "Dev - Ref".ljust(13),
            "% diff".ljust(11),
        ),
        file=f,
    )
    print(
        "{:11.6f}  {:11.6f}  {:11.6f}  {:9.4f}".format(
            ref_mean_oh, dev_mean_oh, oh_diff, oh_pctdiff
        ),
        file=f,
    )

    # ==================================================================
    # Compute MCF and CH4 lifetimes
    # ==================================================================

    # Select only boxes that are purely tropospheric
    # This excludes influence from the stratosphere
    ref_timetrop_mask = refmet["FracOfTimeInTrop"].values != 1.0
    dev_timetrop_mask = devmet["FracOfTimeInTrop"].values != 1.0

    # Get grid box volumes [cm3] (trop + strat)
    ref_vol = (refmet["Met_BXHEIGHT"] * ref_area) * 1e6
    dev_vol = (devmet["Met_BXHEIGHT"] * dev_area) * 1e6

    # Get grid box volumes [cm3] (trop only)
    ref_vol_trop = np.ma.masked_array(ref_vol.values, ref_timetrop_mask)
    dev_vol_trop = np.ma.masked_array(dev_vol.values, dev_timetrop_mask)

    # Get MCF and CH4 density [molec/cm3] (trop + strat)
    # Assume that species is evenly distributed in air, with
    # a mixing ratio of 1. Thus species density = air density.
    ref_dens = refmet["Met_AIRDEN"] / 1e6
    dev_dens = devmet["Met_AIRDEN"] / 1e6

    # Get MCF and CH4 density [molec/cm3] (trop only)
    ref_dens_trop = np.ma.masked_array(ref_dens.values, ref_timetrop_mask)
    dev_dens_trop = np.ma.masked_array(dev_dens.values, dev_timetrop_mask)

    # Get temperature [K] (trop only)
    ref_temp = np.ma.masked_array(refmet["Met_T"].values, ref_timetrop_mask)
    dev_temp = np.ma.masked_array(devmet["Met_T"].values, dev_timetrop_mask)

    # Compute Arrhenius parameter K [cm3/molec/s]
    ref_mcf_k = 1.64e-12 * np.exp(-1520e0 / ref_temp)
    dev_mcf_k = 1.64e-12 * np.exp(-1520e0 / dev_temp)
    ref_ch4_k = 2.45e-12 * np.exp(-1775e0 / ref_temp)
    dev_ch4_k = 2.45e-12 * np.exp(-1775e0 / dev_temp)

    # Numerator: Total atmospheric (trop+strat) burden
    ref_num = np.sum(ref_dens.values * ref_vol.values)
    dev_num = np.sum(dev_dens.values * dev_vol.values)

    # Denominator: Loss rate in troposphere
    ref_mcf_denom = np.sum(ref_mcf_k * ref_oh_trop * ref_dens_trop * ref_vol_trop)

    dev_mcf_denom = np.sum(dev_mcf_k * dev_oh_trop * dev_dens_trop * dev_vol_trop)
    ref_ch4_denom = np.sum(ref_ch4_k * ref_oh_trop * ref_dens_trop * ref_vol_trop)
    dev_ch4_denom = np.sum(dev_ch4_k * dev_oh_trop * dev_dens_trop * dev_vol_trop)

    # Compute lifetimes [years]
    sec_to_year = 365.25 * 86400.0
    ref_mcf_lifetime = (ref_num / ref_mcf_denom) / sec_to_year
    dev_mcf_lifetime = (dev_num / dev_mcf_denom) / sec_to_year
    ref_ch4_lifetime = (ref_num / ref_ch4_denom) / sec_to_year
    dev_ch4_lifetime = (dev_num / dev_ch4_denom) / sec_to_year

    # Compute differences
    mcf_diff = dev_mcf_lifetime - ref_mcf_lifetime
    ch4_diff = dev_ch4_lifetime - ref_ch4_lifetime

    mcf_pctdiff = ((dev_mcf_lifetime - ref_mcf_lifetime) / ref_mcf_lifetime) * 100.0
    ch4_pctdiff = ((dev_ch4_lifetime - ref_ch4_lifetime) / ref_ch4_lifetime) * 100.0

    # Title strings
    title1 = "### MCF lifetime w/r/t tropospheric OH [years]"
    title2 = "### Ref = {}; Dev = {}".format(refstr, devstr)

    # Print header to file
    print("", file=f)
    print("#" * 79, file=f)
    print("{}{}".format(title1.ljust(76), "###"), file=f)
    print("{}{}".format(title2.ljust(76), "###"), file=f)
    print("#" * 79, file=f)

    # Write results to file
    print(
        "{}{}{}{}".format(
            "  Ref".ljust(15),
            "Dev".ljust(13),
            "Dev - Ref".ljust(13),
            "% diff".ljust(11),
        ),
        file=f,
    )
    print(
        "{:11.6f}  {:11.6f}  {:11.6f}  {:9.4f}".format(
            ref_mcf_lifetime, dev_mcf_lifetime, mcf_diff, mcf_pctdiff
        ),
        file=f,
    )

    # Title strings
    title1 = "### CH4 lifetime w/r/t tropospheric OH [years]"
    title2 = "### Ref = {}; Dev = {}".format(refstr, devstr)

    # Print header to file
    print("", file=f)
    print("#" * 79, file=f)
    print("{}{}".format(title1.ljust(76), "###"), file=f)
    print("{}{}".format(title2.ljust(76), "###"), file=f)
    print("#" * 79, file=f)

    # Write results to file
    print(
        "{}{}{}{}".format(
            "  Ref".ljust(15),
            "Dev".ljust(13),
            "Dev - Ref".ljust(13),
            "% diff".ljust(11),
        ),
        file=f,
    )
    print(
        "{:11.6f}  {:11.6f}  {:11.6f}  {:9.4f}".format(
            ref_ch4_lifetime, dev_ch4_lifetime, ch4_diff, ch4_pctdiff
        ),
        file=f,
    )


def add_bookmarks_to_pdf(pdfname, varlist, remove_prefix="", verbose=False):
    """
    Adds bookmarks to an existing PDF file.

    Args:
    -----
        pdfname : str
            Name of an existing PDF file of species or emission plots
            to which bookmarks will be attached.

        varlist : list
            List of variables, which will be used to create the
            PDF bookmark names.

    Keyword Args (optional):
    ------------------------
        remove_prefix : str0
            Specifies a prefix to remove from each entry in varlist
            when creating bookmarks.  For example, if varlist has
            a variable name "SpeciesConc_NO", and you specify
            remove_prefix="SpeciesConc_", then the bookmark for
            that variable will be just "NO", etc.

         verbose : boolean
            Set this flag to True to print extra informational output.
            Default value: False
    """

    # Setup
    pdfobj = open(pdfname, "rb")
    input = PdfFileReader(pdfobj, overwriteWarnings=False)
    output = PdfFileWriter()

    for i, varname in enumerate(varlist):
        bookmarkname = varname.replace(remove_prefix, "")
        if verbose:
            print("Adding bookmark for {} with name {}".format(varname, bookmarkname))
        output.addPage(input.getPage(i))
        output.addBookmark(bookmarkname, i)
        output.setPageMode("/UseOutlines")

    # Write to temp file
    pdfname_tmp = pdfname + "_with_bookmarks.pdf"
    outputstream = open(pdfname_tmp, "wb")
    output.write(outputstream)
    outputstream.close()

    # Rename temp file with the target name
    os.rename(pdfname_tmp, pdfname)
    pdfobj.close()


def add_nested_bookmarks_to_pdf(
    pdfname, category, catdict, warninglist, remove_prefix=""
):

    """
    Add nested bookmarks to PDF.

    Args:
    -----
        pdfname : str
            Path of PDF to add bookmarks to

        category : str
            Top-level key name in catdict that maps to contents of PDF

        catdict : dictionary
            Dictionary containing key-value pairs where one top-level
            key matches category and has value fully describing pages
            in PDF. The value is a dictionary where keys are level 1
            bookmark names, and values are lists of level 2 bookmark
            names, with one level 2 name per PDF page.  Level 2 names
            must appear in catdict in the same order as in the PDF.

        warninglist : list of strings
            Level 2 bookmark names to skip since not present in PDF.

    Keyword Args (optional):
    ------------------------
        remove_prefix : str
            Prefix to be remove from warninglist names before comparing with 
            level 2 bookmark names in catdict.
            Default value: empty string (warninglist names match names 
            in catdict) 
    """

    # ==================================================================
    # Setup
    # ==================================================================
    pdfobj = open(pdfname, "rb")
    input = PdfFileReader(pdfobj, overwriteWarnings=False)
    output = PdfFileWriter()
    warninglist = [k.replace(remove_prefix, "") for k in warninglist]

    # ==================================================================
    # Loop over the subcategories in this category; make parent bookmark
    # ==================================================================
    i = -1
    for subcat in catdict[category]:

        # First check that there are actual variables for
        # this subcategory; otherwise skip
        numvars = 0
        if catdict[category][subcat]:
            for varname in catdict[category][subcat]:
                if varname in warninglist:
                    continue
                else:
                    numvars += 1
        else:
            continue
        if numvars == 0:
            continue

        # There are non-zero variables to plot in this subcategory
        i = i + 1
        output.addPage(input.getPage(i))
        parent = output.addBookmark(subcat, i)
        output.setPageMode("/UseOutlines")
        first = True

        # Loop over variables in this subcategory; make children bookmarks
        for varname in catdict[category][subcat]:
            if varname in warninglist:
                print("Warning: skipping {}".format(varname))
                continue
            if first:
                output.addBookmark(varname, i, parent)
                first = False
            else:
                i = i + 1
                output.addPage(input.getPage(i))
                output.addBookmark(varname, i, parent)
                output.setPageMode("/UseOutlines")

    # ==================================================================
    # Write to temp file
    # ==================================================================
    pdfname_tmp = pdfname + "_with_bookmarks.pdf"
    outputstream = open(pdfname_tmp, "wb")
    output.write(outputstream)
    outputstream.close()

    # Rename temp file with the target name
    os.rename(pdfname_tmp, pdfname)
    pdfobj.close()


def add_missing_variables(refdata, devdata, **kwargs):
    """
    Compares two xarray Datasets, "Ref", and "Dev".  For each variable
    that is present  in "Ref" but not in "Dev", a DataArray of missing
    values (i.e. NaN) will be added to "Dev".  Similarly, for each
    variable that is present in "Dev" but not in "Ref", a DataArray
    of missing values will be added to "Ref".

    This routine is mostly intended for benchmark purposes, so that we
    can represent variables that were removed from a new GEOS-Chem
    version by missing values in the benchmark plots.

    Args:
    -----
        refdata : xarray Dataset
            The "Reference" (aka "Ref") dataset.

        devdata : xarray Dataset
            The "Development" (aka "Dev") dataset

    Returns:
    --------
        refdata, devdata : xarray Datasets
            The returned "Ref" and "Dev" datasets, with
            placeholder missing value variables added.
    """
    # ==================================================================
    # Initialize
    # ==================================================================

    # Make sure that refdata and devdata are both xarray Dataset objects
    if not isinstance(refdata, xr.Dataset):
        raise TypeError("The refdata object must be an xarray Dataset!")
    if not isinstance(devdata, xr.Dataset):
        raise TypeError("The refdata object must be an xarray Dataset!")

    # Find common variables as well as variables only in one or the other
    vardict = core.compare_varnames(refdata, devdata, quiet=True)
    refonly = vardict["refonly"]
    devonly = vardict["devonly"]

    # Don't clobber any DataArray attributes
    with xr.set_options(keep_attrs=True):

        # ==============================================================
        # For each variable that is in refdata but not in devdata,
        # add a new DataArray to devdata with the same sizes but
        # containing all NaN's.  This will allow us to represent those
        # variables as missing values # when we plot against refdata.
        # ==============================================================
        for v in refonly:
            dr = core.create_dataarray_of_nan(
                name=refdata[v].name,
                sizes=devdata.sizes,
                coords=devdata.coords,
                attrs=refdata[v].attrs,
                **kwargs
            )
            devdata = xr.merge([devdata, dr])

        # ==============================================================
        # For each variable that is in devdata but not in refdata,
        # add a new DataArray to refdata with the same sizes but
        # containing all NaN's.  This will allow us to represent those
        # variables as missing values # when we plot against devdata.
        # ==================================================================
        for v in devonly:
            dr = core.create_dataarray_of_nan(
                name=devdata[v].name,
                sizes=refdata.sizes,
                coords=refdata.coords,
                attrs=devdata[v].attrs,
                **kwargs
            )
            refdata = xr.merge([refdata, dr])

    return refdata, devdata


def get_troposphere_mask(ds):
    """
    Returns a mask array for picking out the tropospheric grid boxes.
    
    Args:
    -----
        ds : xarray Dataset
            Dataset containing certain met field variables (i.e.
            Met_TropLev, Met_BXHEIGHT).

    Returns:
    --------
        tropmask : numpy ndarray
            Tropospheric mask.  False denotes grid boxes that are
            in the troposphere and True in the stratosphere
            (as per Python masking logic).
    """

    # ==================================================================
    # Initialization
    # ==================================================================

    # Make sure ds is an xarray Dataset object
    if not isinstance(ds, xr.Dataset):
        raise TypeError("The ds argument must be an xarray Dataset!")

    # Make sure certain variables are found
    if "Met_BXHEIGHT" not in ds.data_vars.keys():
        raise ValueError("Met_BXHEIGHT could not be found!")
    if "Met_TropLev" not in ds.data_vars.keys():
        raise ValueError("Met_TropLev could not be found!")

    # Mask of tropospheric grid boxes in the Ref dataset
    shape = core.get_shape_of_data(np.squeeze(ds["Met_BXHEIGHT"]))

    # ==================================================================
    # Create the mask arrays for the troposphere
    #
    # Convert the Met_TropLev DataArray objects to numpy ndarrays of
    # integer.  Also subtract 1 to convert from Fortran to Python
    # array index notation.
    # ==================================================================

    if len(shape) == 4:
        # --------------------------------------------------------------
        # There are multiple time slices
        # --------------------------------------------------------------

        # Create the tropmask array with dims (time, lev, lat*lon)
        tropmask = np.ones((shape[0], shape[1], np.prod(np.array(shape[2:]))), bool)

        # Loop over each time
        for t in range(tropmask.shape[0]):

            # Pick the tropopause level and make a 1-D array
            values = ds["Met_TropLev"].isel(time=t).values
            lev = np.int_(np.squeeze(values) - 1)
            lev_1d = lev.flatten()

            # Create the tropospheric mask array
            for x in range(tropmask.shape[2]):
                tropmask[t, 0 : lev_1d[x], x] = False

    elif len(shape) == 3:
        # --------------------------------------------------------------
        # There is only one time slice
        # --------------------------------------------------------------

        # Create the tropmask array with dims (lev, lat*lon)
        tropmask = np.ones((shape[0], np.prod(np.array(shape[1:]))), bool)

        # Pick the tropopause level and make a 1-D array
        values = ds["Met_TropLev"].values
        lev = np.int_(np.squeeze(values) - 1)
        lev_1d = lev.flatten()

        # Create the tropospheric mask array
        for x in range(tropmask.shape[1]):
            tropmask[0 : lev_1d[x], x] = False

    # Reshape into the same shape as Met_BxHeight
    return tropmask.reshape(shape)


def check_units(refdata, devdata, varname):
    # If units are mol/mol then convert to ppb
    conc_units = ["mol mol-1 dry", "mol/mol", "mol mol-1"]
    if refdata[varname].units.strip() in conc_units:
        refdata[varname].attrs["units"] = "ppbv"
        refdata[varname].values = refdata[varname].values * 1e9
    if devdata[varname].units.strip() in conc_units:
        devdata[varname].attrs["units"] = "ppbv"
        devdata[varname].values = devdata[varname].values * 1e9

    # Binary diagnostic concentrations have units ppbv. Change to ppb.
    if refdata[varname].units.strip() == "ppbv":
        refdata[varname].attrs["units"] = "ppb"
    if devdata[varname].units.strip() == "ppbv":
        devdata[varname].attrs["units"] = "ppb"

    # Check that units match
    units_ref = refdata[varname].units.strip()
    units_dev = devdata[varname].units.strip()
    if units_ref != units_dev:
        print_units_warning = True
        if print_units_warning:
            print("WARNING: ref and dev concentration units do not match!")
            print("Ref units: {}".format(units_ref))
            print("Dev units: {}".format(units_dev))
        if enforce_units:
            # if enforcing units, stop the program if
            # units do not match
            assert units_ref == units_dev, "Units do not match for {}!".format(varname)
        else:
            # if not enforcing units, just keep going after
            # only printing warning once
            print_units_warning = False

    return units_ref, units_dev


def reshape_MAPL_CS(ds, vdims):
    # Reshape cubed sphere data if using MAPL v1.0.0+
    if "nf" in vdims and "Xdim" in vdims and "Ydim" in vdims:
        ds = ds.stack(lat=("nf", "Ydim"))
        ds = ds.rename({"Xdim": "lon"})
        if "lev" in ds.dims:
            ds = ds.transpose("lev", "lat", "lon")
        else:
            ds = ds.transpose("lat", "lon")
    return ds

