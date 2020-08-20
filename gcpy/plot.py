import os
import copy
import yaml
import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
from matplotlib.backends.backend_pdf import PdfPages
from PyPDF2 import PdfFileWriter, PdfFileReader, PdfFileMerger
from .grid import GEOS_72L_grid, GEOS_47L_grid, get_vert_grid, get_pressure_indices, \
     pad_pressure_edges, convert_lev_to_pres, get_grid_extents, call_make_grid, get_input_res
from .regrid import regrid_comparison_data, create_regridders, dict_72_to_47, reduce_72_to_47
from .util import reshape_MAPL_CS, get_diff_of_diffs, get_nan_mask, all_zero_or_nan, slice_by_lev_and_time, compare_varnames
from .units import check_units, data_unit_is_mol_per_mol
from .constants import MW_AIR_g
from joblib import Parallel, delayed, cpu_count, parallel_backend
from multiprocessing import current_process
from shutil import rmtree
from tempfile import TemporaryDirectory
import warnings

# Save warnings format to undo overwriting built into PyPDF2
warning_format = warnings.showwarning

# Suppress numpy divide by zero warnings to prevent output spam
np.seterr(divide="ignore", invalid="ignore")

current_dir = os.path.dirname(__file__)

rgb_WhGrYlRd = np.genfromtxt(current_dir+'/colormaps/WhGrYlRd.txt',
                             delimiter=' ')
WhGrYlRd = mcolors.ListedColormap(rgb_WhGrYlRd/255.0)


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
    plot_type="single_level",
    ratio_log=False
):

    """
    Plotting function to be called from compare_single_level or
    compare_zonal_mean.
    Can also be called on its own
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
    xticklabels: list of str
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
            [vmin, vmax] = [-absmax, absmax]
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
                [vmin, vmax] = [1/fracdiffabsmax, fracdiffabsmax]
                #if vmin > 0.5:
                #    vmin = 0.5
                #if vmax < 2:
                #    vmax = 2
            else:
                [vmin, vmax] = [0.5, 2]
    if verbose:
        print("Subplot ({}) vmin, vmax: {}, {}".format(rowcol, vmin, vmax))

    # Normalize colors (put into range [0..1] for matplotlib methods)
    if subplot in ("ref", "dev"):
        norm = normalize_colors(
            vmin, vmax, is_difference=use_cmap_RdBu,
            log_color_scale=log_color_scale, ratio_log=ratio_log
        )
    elif subplot in ("dyn_abs_diff", "res_abs_diff"):
        norm = normalize_colors(vmin, vmax, is_difference=True)
    else:
        #remove NaNs for compatibility with color normalization
        plot_val = get_nan_mask(plot_val)
        norm = normalize_colors(vmin, vmax, is_difference=True, log_color_scale=True,ratio_log=ratio_log)
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
        elif subplot in ("dyn_frac_diff", "res_frac_diff") and np.all(np.isin(plot_val, [1])):
            cb.set_ticklabels(["Ref and Dev equal throughout domain"])
        elif subplot in ("dyn_frac_diff", "res_frac_diff"):
            if subplot is "dyn_frac_diff" and vmin != 0.5 and vmax != 2.0:
                if vmin > 0.1 and vmax < 10:
                    cb.locator = mticker.MaxNLocator(nbins=4)
                    cb.update_ticks()
                    cb.formatter = mticker.ScalarFormatter()
                else:
                    cb.formatter = mticker.LogFormatter(base=10)
                    cb.update_ticks()
                    cb.locator = mticker.LogLocator(base=10, subs='all')
                cb.update_ticks()
            else:
                cb.formatter = mticker.ScalarFormatter()
                cb.set_ticks([0.5,0.75,1,1.5,2.0])
        else:
            if (vmax - vmin) < 0.1 or (vmax - vmin) > 100:
                cb.locator = mticker.MaxNLocator(nbins=4)

    try:
        cb.formatter.set_useOffset(False)
    except:
        #not all automatically chosen colorbar formatters properly handle the above method
        pass
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
    refmet=None,
    devmet=None,
    weightsdir='.',
    pdfname="",
    cmpres=None,
    match_cbar=True,
    normalize_by_area=False,
    enforce_units=True,
    convert_to_ugm3=False,
    flip_ref=False,
    flip_dev=False,
    use_cmap_RdBu=False,
    verbose=False,
    log_color_scale=False,
    extra_title_txt=None,
    plot_extent = [-1000, -1000, -1000, -1000],
    n_job=-1,
    sigdiff_list=[],
    second_ref=None,
    second_dev=None,
    spcdb_dir=os.path.dirname(__file__)
):
    """
    Create single-level 3x2 comparison map plots for variables common
    in two xarray Datasets. Optionally save to PDF.
    Args:
    -----
        refdata : xarray dataset
            Dataset used as reference in comparison
        refstr  : str OR list of str
            String description for reference data to be used in plots
            OR list containing [ref1str, ref2str] for diff-of-diffs plots
        devdata : xarray dataset
            Dataset used as development in comparison
        devstr  : str OR list of str
            String description for development data to be used in plots
            OR list containing [dev1str, dev2str] for diff-of-diffs plots
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
        refmet : xarray dataset
            Dataset containing ref meteorology
            Default value: None
        devmet : xarray dataset
            Dataset containing dev meteorology  
            Default value: None
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
            raw data by grid area. Input ref and dev datasets must include
            AREA variable in m2 if normalizing by area.
            Default value: False
        enforce_units : boolean
            Set this flag to True to force an error if Ref and Dev
            variables have different units.
            Default value: True
        convert_to_ugm3 : boolean
            Whether to convert data units to ug/m3 for plotting.
            Default value: False
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
        plot_extent : list
            Defines the extent of the region to be plotted in form 
            [minlon, maxlon, minlat, maxlat]. Default value plots extent of input grids.
            Default value: [-1000, -1000, -1000, -1000]            
        n_job : int
            Defines the number of simultaneous workers for parallel plotting.
            Set to 1 to disable parallel plotting. Value of -1 allows the application to decide.
            Default value: -1
        sigdiff_list: list of str
            Returns a list of all quantities having significant
            differences (where |max(fractional difference)| > 0.1).
            Default value: []
        second_ref : xarray Dataset
            A dataset of the same model type / grid as refdata, to be used in diff-of-diffs plotting.
            Default value: None
        second_dev : xarray Dataset
            A dataset of the same model type / grid as devdata, to be used in diff-of-diffs plotting.
            Default value: None
        spcdb_dir  : str 
            Directory containing species_database.yml file.
            Default value: Path of GCPy code repository
    """
    warnings.showwarning = warning_format

    # Error check arguments
    if not isinstance(refdata, xr.Dataset):
        raise TypeError("The refdata argument must be an xarray Dataset!")

    if not isinstance(devdata, xr.Dataset):
        raise TypeError("The devdata argument must be an xarray Dataset!")
                
    # Determine if doing diff-of-diffs
    if second_ref is not None and second_dev is not None:
        diff_of_diffs = True
    else:
        diff_of_diffs = False

    # Prepare diff-of-diffs datasets if needed
    if diff_of_diffs:
        refdata, fracrefdata = get_diff_of_diffs(refdata, second_ref)
        devdata, fracdevdata = get_diff_of_diffs(devdata, second_dev)
        if type(refstr) is list:
            frac_refstr = '{} / {}'.format(refstr[1], refstr[0])
            refstr = '{} - {}'.format(refstr[1], refstr[0])
        else:
            frac_refstr = 'Ref2 / Ref1'
        if type(devstr) is list:
            frac_devstr = '{} / {}'.format(devstr[1], devstr[0])
            devstr = '{} - {}'.format(devstr[1], devstr[0])
        else:
            frac_devstr = 'Dev2 / Dev1'
    # If no varlist is passed, plot all (surface only for 3D)
    if varlist is None:
        quiet = not verbose
        vardict = compare_varnames(refdata, devdata, quiet=quiet)
        varlist = vardict["commonvars3D"] + vardict["commonvars2D"]
        print("Plotting all common variables")
    n_var = len(varlist)

    # If no PDF name passed, then do not save to PDF
    savepdf = True
    if pdfname == "":
        savepdf = False
    # Cleanup previous temporary PDFs produced during parallelization
    #for i in range(n_var):
    #    try:
    #        os.remove(os.path.join(temp_dir, pdfname + "BENCHMARKFIGCREATION.pdf" + str(i)))
    #    except:
    #        continue

    # If converting to ug/m3, load the species database
    if convert_to_ugm3:
        properties_path = os.path.join(spcdb_dir, "species_database.yml")
        properties = yaml.load(open(properties_path), Loader=yaml.FullLoader)

    # Get grid info and regrid if necessary
    [refres, refgridtype, devres, devgridtype, cmpres, cmpgridtype, regridref, 
     regriddev, regridany, refgrid, devgrid, cmpgrid, refregridder,
     devregridder, refregridder_list, devregridder_list] = create_regridders(
         refdata,
         devdata,
         weightsdir,
         cmpres=cmpres
     )

    # ==============================================================
    # Handle grid extents for lat-lon grids
    # ==============================================================

    # Get lat/lon extents, if applicable
    refminlon, refmaxlon, refminlat, refmaxlat = get_grid_extents(refgrid)
    devminlon, devmaxlon, devminlat, devmaxlat = get_grid_extents(devgrid)
    if cmpgridtype is not "cs":
        cmpminlon, cmpmaxlon, cmpminlat, cmpmaxlat = get_grid_extents(cmpgrid)
    else:
        cmpminlon, cmpmaxlon, cmpminlat, cmpmaxlat = [-180, 180, -90, 90]

    # Set plot bounds for non cubed-sphere regridding and plotting
    #if cmpgridtype == "ll":
    ref_extent = (refminlon, refmaxlon, refminlat, refmaxlat)
    dev_extent = (devminlon, devmaxlon, devminlat, devmaxlat)
    cmp_extent = (cmpminlon, cmpmaxlon, cmpminlat, cmpmaxlat)
    # Trim data to extent of comparison grid if using lat-lon grids.
    # Get comparison date extents in same midpoint format as lat-lon grid.
    cmp_mid_minlon, cmp_mid_maxlon, cmp_mid_minlat, cmp_mid_maxlat = \
                                    get_grid_extents(cmpgrid,edges=False)
    if refgridtype == "ll" and ref_extent != cmp_extent:
        refdata = refdata.\
                      where(refdata.lon >= cmp_mid_minlon, drop=True).\
                      where(refdata.lon <= cmp_mid_maxlon, drop=True).\
                      where(refdata.lat >= cmp_mid_minlat, drop=True).\
                      where(refdata.lat <= cmp_mid_maxlat, drop=True)
        if diff_of_diffs:
            fracrefdata = fracrefdata.\
                            where(fracrefdata.lon >= cmp_mid_minlon,drop=True).\
                            where(fracrefdata.lon <= cmp_mid_maxlon,drop=True).\
                            where(fracrefdata.lat >= cmp_mid_minlat,drop=True).\
                            where(fracrefdata.lat <= cmp_mid_maxlat,drop=True)
    if devgridtype == "ll" and dev_extent != cmp_extent:
        devdata = devdata.\
                      where(devdata.lon >= cmp_mid_minlon, drop=True).\
                      where(devdata.lon <= cmp_mid_maxlon, drop=True).\
                      where(devdata.lat >= cmp_mid_minlat, drop=True).\
                      where(devdata.lat <= cmp_mid_maxlat, drop=True)
        if diff_of_diffs:
            fracdevdata = fracdevdata.\
                            where(fracdevdata.lon >= cmp_mid_minlon,drop=True).\
                            where(fracdevdata.lon <= cmp_mid_maxlon,drop=True).\
                            where(fracdevdata.lat >= cmp_mid_minlat,drop=True).\
                            where(fracdevdata.lat <= cmp_mid_maxlat,drop=True)

    # ==============================================================
    # Loop over all variables
    # ==============================================================
    ds_refs = [None] * n_var
    frac_ds_refs = [None] * n_var
    ds_devs = [None] * n_var
    frac_ds_devs = [None] * n_var
    for i in range(n_var):

        varname = varlist[i]
        ref_dims = refdata[varname].dims
        dev_dims = devdata[varname].dims

        # ==============================================================
        # Slice the data, allowing for no time dimension (bpch)
        # ==============================================================
        # Ref
        ds_refs[i] = slice_by_lev_and_time(
            refdata,
            varname,
            itime,
            ilev,
            flip_ref
        )
        if diff_of_diffs:
            frac_ds_refs[i] = slice_by_lev_and_time(
                fracrefdata,
                varname,
                itime,
                ilev,
                flip_ref
            )
            fracref_dims = fracrefdata[varname].dims
        # Dev
        ds_devs[i] = slice_by_lev_and_time(
            devdata,
            varname,
            itime,
            ilev,
            flip_dev
        )
        if diff_of_diffs:
            frac_ds_devs[i] = slice_by_lev_and_time(
                fracdevdata,
                varname,
                itime,
                ilev,
                flip_dev
            )
            fracdev_dims = fracdevdata[varname].dims

        # ==================================================================
        #  Handle units as needed
        # ==================================================================
        
        # Convert to ppb if units string is variation of mol/mol
        if data_unit_is_mol_per_mol(ds_refs[i]):
            ds_refs[i].values = ds_refs[i].values * 1e9
            ds_refs[i].attrs["units"] = "ppb"
        if data_unit_is_mol_per_mol(ds_devs[i]):
            ds_devs[i].values = ds_devs[i].values * 1e9
            ds_devs[i].attrs["units"] = "ppb"

        # If units string is ppbv (true for bpch data) then rename units
        if ds_refs[i].units.strip() == "ppbv":
            ds_refs[i].attrs["units"] = "ppb"
        if ds_devs[i].units.strip() == "ppbv":
            ds_devs[i].attrs["units"] = "ppb"

        # If units string is W/m2 (may be true for bpch data) then rename units
        if ds_refs[i].units.strip() == "W/m2":
            ds_refs[i].attrs["units"] = "W m-2"
        if ds_devs[i].units.strip() == "W/m2":
            ds_devs[i].attrs["units"] = "W m-2"

        # If units string is UNITLESS (may be true for bpch data) then rename units
        if ds_refs[i].units.strip() == "UNITLESS":
            ds_refs[i].attrs["units"] = "1"
        if ds_devs[i].units.strip() == "UNITLESS":
            ds_devs[i].attrs["units"] = "1"

        # Check that units are the same in ref and dev. Will exit with
        # an error if do not match and enforce_units is true (default).
        check_units(ds_refs[i], ds_devs[i])

        # Convert from ppb to ug/m3 if convert_to_ugm3 is passed as true
        if convert_to_ugm3:

            # Error checks: must pass met, not normalize by area, and be in ppb
            if refmet is None or devmet is None:
                msg = "Met mata ust be passed to convert units to ug/m3."
                raise ValueError(msg)
            elif normalize_by_area:
                msg = "Normalizing by area is now allowed if plotting ug/m3"
                raise ValueError(msg)
            elif ds_refs[i].units != "ppb" or ds_devs[i].units != "ppb":
                msg = "Units must be mol/mol if converting to ug/m3."
                raise ValueError(msg)

            # Slice air density data by lev and time
            # (assume same format and dimensions as refdata and devdata)
            ref_airden = slice_by_lev_and_time(
                refmet,
                "Met_AIRDEN",
                itime,
                ilev,
                False
            )
            dev_airden = slice_by_lev_and_time(
                devmet,
                "Met_AIRDEN",
                itime,
                ilev,
                False
            )

            # Get a list of properties for the given species
            spc_name = varname.replace(varname.split("_")[0]+"_", "")
            species_properties = properties.get(spc_name)

            # If no properties are found, then exit with an error.
            # Otherwise, get the molecular weight in g/mol.
            if species_properties is None:
                # Hack lumped species until we implement a solution
                if spc_name in ["Simple_SOA","Complex_SOA"]:
                    spc_mw_g = 150.0
                else:
                    msg = "No properties found for {}. Cannot convert" \
                          + " to ug/m3."
                    raise ValueError(msg.format(spc_name))
            else:
                spc_mw_g = species_properties.get("MW_g")
                if spc_mw_g is None:
                    msg = "Molecular weight not found for for species {}!" \
                          + " Cannot convert to ug/m3.".format(spc_name)
                    raise ValueError(msg)

            # Convert values from ppb to ug/m3:
            # ug/m3 = mol/mol * mol/g air * kg/m3 air * 1e3g/kg
            #         * g/mol spc * 1e6ug/g
            #       = ppb * air density * (spc MW / air MW)
            ds_refs[i].values = ds_refs[i].values * ref_airden.values \
                                * ( spc_mw_g / MW_AIR_g )
            ds_devs[i].values = ds_devs[i].values * dev_airden.values \
                                * ( spc_mw_g / MW_AIR_g )

            # Update units string
            ds_refs[i].attrs["units"] = "\u03BCg/m3" # ug/m3 using mu
            ds_devs[i].attrs["units"] = "\u03BCg/m3"

        # ==============================================================
        # Reshape cubed sphere data if using MAPL v1.0.0+
        # TODO: update function to expect data in this format
        # ==============================================================
        ds_refs[i] = reshape_MAPL_CS(ds_refs[i])
        ds_devs[i] = reshape_MAPL_CS(ds_devs[i])
        if diff_of_diffs:
            frac_ds_refs[i] = reshape_MAPL_CS(frac_ds_refs[i])
            frac_ds_devs[i] = reshape_MAPL_CS(frac_ds_devs[i])

    # ==================================================================
    # Get the area variables if normalize_by_area=True. They can be
    # either in the main datasets as variable AREA or in the optionally
    # passed meteorology datasets as Met_AREAM2.
    # ==================================================================
    if normalize_by_area:
        # ref
        if "AREA" in refdata.data_vars.keys():
            ref_area = refdata["AREA"]
        elif refmet is not None:
            if "Met_AREAM2" in refmet.data_vars.keys():
                ref_area = refmet["Met_AREAM2"]
        else:
            msg = "normalize_by_area = True but AREA not " \
                + "present in the Ref dataset and ref met with Met_AREAM2" \
                + " not passed!"
            raise ValueError(msg)
        if "time" in ref_area.dims:
            ref_area = ref_area.isel(time=0)
        if refgridtype == 'cs':
            ref_area = reshape_MAPL_CS(ref_area)

        # dev
        if "AREA" in devdata.data_vars.keys():
            dev_area = devdata["AREA"]
        elif devmet is not None:
            if "Met_AREAM2" in devmet.data_vars.keys():
                dev_area = devmet["Met_AREAM2"]
        else:
            msg = "normalize_by_area = True but AREA not " \
                + "present in the Dev dataset and dev met with Met_AREAM2" \
                | " not passed!"
            raise ValueError(msg)
        if "time" in dev_area.dims:
            dev_area = dev_area.isel(time=0)
        if devgridtype == 'cs':
            dev_area = reshape_MAPL_CS(dev_area)

        # Make sure the areas do not have a lev dimension
        if "lev" in ref_area.dims:
            ref_area = ref_area.isel(lev=0)
        if "lev" in dev_area.dims:
            dev_area = dev_area.isel(lev=0)

    # ==================================================================
    # Create arrays for each variable in Ref and Dev datasets
    # and do any necessary regridding. 'cmp' stands for comparison
    # and represents ref and dev data regridded as needed to a common
    # grid type and resolution for use in difference and ratio plots.
    # ==================================================================
    ds_ref_cmps = [None] * n_var
    ds_dev_cmps = [None] * n_var
    frac_ds_ref_cmps = [None] * n_var
    frac_ds_dev_cmps = [None] * n_var

    global_cmp_grid = call_make_grid(cmpres, cmpgridtype, False, True)[0]
    cmpminlon_ind = np.where(global_cmp_grid["lon"] >= cmpminlon)[0][0]
    cmpmaxlon_ind = np.where(global_cmp_grid["lon"] <= cmpmaxlon)[0][-1]
    cmpminlat_ind = np.where(global_cmp_grid["lat"] >= cmpminlat)[0][0]
    cmpmaxlat_ind = np.where(global_cmp_grid["lat"] <= cmpmaxlat)[0][-1]

    for i in range(n_var):
        ds_ref = ds_refs[i]
        ds_dev = ds_devs[i]

        # Do area normalization before regridding if normalize_by_area is True.
        # Assumes units are the same in ref and dev. If enforce_units is passed
        # as false then normalization may not be correct.
        if normalize_by_area:
            exclude_list = ["WetLossConvFrac", "Prod_", "Loss_"]
            if not any(s in varname for s in exclude_list):
                ds_ref.values = ds_ref.values / ref_area.values
                ds_dev.values = ds_dev.values / dev_area.values
                ds_refs[i] = ds_ref
                ds_devs[i] = ds_dev
                if diff_of_diffs:
                    frac_ds_refs[i] = frac_ds_refs[i].values / ref_area.values
                    frac_ds_devs[i] = frac_ds_devs[i].values / dev_area.values

        # Ref
        ds_ref_cmps[i] = regrid_comparison_data(
            ds_ref,
            refres,
            regridref,
            refregridder,
            refregridder_list,
            global_cmp_grid,
            refgridtype,
            cmpminlat_ind,
            cmpmaxlat_ind,
            cmpminlon_ind,
            cmpmaxlon_ind
        )

        # Dev
        ds_dev_cmps[i] = regrid_comparison_data(
            ds_dev,
            devres,
            regriddev,
            devregridder,
            devregridder_list,
            global_cmp_grid,
            devgridtype,
            cmpminlat_ind,
            cmpmaxlat_ind,
            cmpminlon_ind,
            cmpmaxlon_ind
        )

        # Diff of diffs
        if diff_of_diffs:
            frac_ds_ref_cmps[i] = regrid_comparison_data(
                frac_ds_refs[i],
                refres,
                regridref,
                refregridder,
                refregridder_list,
                global_cmp_grid,
                refgridtype,
                cmpminlat_ind,
                cmpmaxlat_ind,
                cmpminlon_ind, 
                cmpmaxlon_ind
            )
            frac_ds_dev_cmps[i] = regrid_comparison_data(
                frac_ds_devs[i],
                devres,
                regriddev,
                devregridder,
                devregridder_list,
                global_cmp_grid,
                devgridtype,
                cmpminlat_ind,
                cmpmaxlat_ind,
                cmpminlon_ind,
                cmpmaxlon_ind
            )


    # =================================================================
    # Define function to create a single page figure to be called
    # in a parallel loop
    # =================================================================
    def createfig(ivar, temp_dir):

        # Suppress harmless run-time warnings (mostly about underflow)
        warnings.filterwarnings('ignore', category=RuntimeWarning)
        warnings.filterwarnings('ignore', category=UserWarning)

        if savepdf and verbose:
            print("{} ".format(ivar), end="")
        varname = varlist[ivar]
        varndim_ref = refdata[varname].ndim
        varndim_dev = devdata[varname].ndim

        ds_ref = ds_refs[ivar]
        ds_dev = ds_devs[ivar]

        # ==============================================================
        # Set units and subtitle, including modification if normalizing
        # area. Note if enforce_units is False (non-default) then
        # units on difference plots will be wrong.
        # ==============================================================
        cmn_units = ds_ref.attrs["units"]
        subtitle_extra = ""
        varndim = varndim_ref
        if normalize_by_area:
            exclude_list = ["WetLossConvFrac", "Prod_", "Loss_"]
            if not any(s in varname for s in exclude_list):
                if "/" in cmn_units:
                    cmn_units = "{}/m2".format(cmn_units)
                else:
                    cmn_units = "{} m-2".format(cmn_units)
                ds_ref.attrs["units"] = cmn_units
                ds_dev.attrs["units"] = cmn_units
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
        frac_ds_ref_cmp = frac_ds_ref_cmps[ivar]
        frac_ds_dev_cmp = frac_ds_dev_cmps[ivar]

        # Reshape comparison cubed sphere data, if any
        if cmpgridtype == "cs":
            ds_ref_cmp_reshaped = ds_ref_cmp.data.reshape(6, cmpres, cmpres)
            ds_dev_cmp_reshaped = ds_dev_cmp.data.reshape(6, cmpres, cmpres)
            if frac_ds_ref_cmp is not None:
                frac_ds_ref_cmp_reshaped = \
                                frac_ds_ref_cmp.data.reshape(6, cmpres, cmpres)
            if frac_ds_dev_cmp is not None:
                frac_ds_dev_cmp_reshaped = \
                                frac_ds_dev_cmp.data.reshape(6, cmpres, cmpres)

        # ==============================================================
        # Get min and max values for use in the colorbars
        # ==============================================================

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
        # boundary line, as described here: https://stackoverflow.com/
        # questions/46527456/preventing-spurious-horizontal-lines-for-
        # ungridded-pcolormesh-data
        if cmpgridtype == "cs":
            absdiff = np.ma.masked_where(np.abs(cmpgrid["lon"] - 180) < 2,
                                         absdiff)

        # ==============================================================
        # Calculate fractional difference, set divides by zero to NaN
        # ==============================================================
        if cmpgridtype == "ll":
            # Replace fractional difference plots with absolute difference
            # of fractional datasets if necessary
            if frac_ds_dev_cmp is not None and frac_ds_ref_cmp is not None:
                fracdiff = np.array(frac_ds_dev_cmp) -       \
                           np.array(frac_ds_ref_cmp)
            else:
                fracdiff = np.abs(np.array(ds_dev_cmp)) /    \
                           np.abs(np.array(ds_ref_cmp))
        else:
            if frac_ds_dev_cmp is not None and frac_ds_ref_cmp is not None:
                fracdiff = frac_ds_dev_cmp_reshaped -        \
                           frac_ds_ref_cmp_reshaped
            else:
                fracdiff = np.abs(ds_dev_cmp_reshaped) /     \
            np.abs(ds_ref_cmp_reshaped)

        # Replace Infinity values with NaN
        fracdiff = np.where(np.abs(fracdiff) == np.inf, np.nan, fracdiff)
        fracdiff[np.abs(fracdiff>1e308)] = np.nan

        # Test if the frac. diff. is zero everywhere or NaN everywhere
        fracdiff_is_all_zero = not np.any(fracdiff) or       \
                               (np.nanmin(fracdiff) == 0 and \
                                np.nanmax(fracdiff) ==0)
        fracdiff_is_all_nan = np.isnan(fracdiff).all() or ref_is_all_zero

        # Absolute max value of fracdiff, excluding NaNs
        fracdiffabsmax = max([np.abs(np.nanmin(fracdiff)),
                              np.abs(np.nanmax(fracdiff))])

        # For cubed-sphere, take special care to avoid a spurious
        # boundary line, as described here: https://stackoverflow.com/
        # questions/46527456/preventing-spurious-horizontal-lines-for-
        # ungridded-pcolormesh-data
        if cmpgridtype == "cs":
            fracdiff = np.ma.masked_where(np.abs(cmpgrid["lon"] - 180) < 2, \
                                          fracdiff)

        # ==============================================================
        # Create 3x2 figure
        # ==============================================================

        # Create figures and axes objects
        # Also define the map projection that will be shown
        figs, ((ax0, ax1), (ax2, ax3), (ax4, ax5)) = plt.subplots(
            3, 2, figsize=[12, 14],
            subplot_kw={"projection": ccrs.PlateCarree()}
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
                    "{}, {}".format(varname, levstr),
                    fontsize=fontsize, y=offset
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
        # Set titles for plots
        # ==============================================================

        if refgridtype == "ll":
            ref_title = "{} (Ref){}\n{}".format(refstr, subtitle_extra, refres)
        else:
            ref_title = "{} (Ref){}\nc{}".format(refstr, subtitle_extra, refres)

        if devgridtype == "ll":
            dev_title = "{} (Dev){}\n{}".format(devstr, subtitle_extra, devres)
        else:
            dev_title = "{} (Dev){}\nc{}".format(devstr, subtitle_extra, devres)

        if regridany:
            absdiff_dynam_title = \
                "Difference ({})\nDev - Ref, Dynamic Range".format(cmpres)
            absdiff_fixed_title = \
                "Difference ({})\nDev - Ref, Restricted Range [5%,95%]".\
                format(cmpres)
            if diff_of_diffs:
                fracdiff_dynam_title = \
                        "Difference ({}), Dynamic Range\n{} - {}".\
                        format(cmpres,frac_devstr,frac_refstr)
                fracdiff_fixed_title = \
                        "Difference ({}), Restricted Range [5%,95%]\n{} - {}".\
                        format(cmpres, frac_devstr, frac_refstr)
            else:
                fracdiff_dynam_title = \
                        "Ratio ({})\nDev/Ref, Dynamic Range".format(cmpres)
                fracdiff_fixed_title = \
                        "Ratio ({})\nDev/Ref, Fixed Range".format(cmpres)
        else:
            absdiff_dynam_title = "Difference\nDev - Ref, Dynamic Range"
            absdiff_fixed_title = \
                            "Difference\nDev - Ref, Restricted Range [5%,95%]"
            if diff_of_diffs:
                fracdiff_dynam_title = \
                        "Difference, Dynamic Range\n{} - {}".\
                        format(frac_devstr,frac_refstr)
                fracdiff_fixed_title = \
                        "Difference, Restricted Range [5%,95%]\n{} - {}".\
                        format(frac_devstr, frac_refstr)
            else:
                fracdiff_dynam_title = "Ratio \nDev/Ref, Dynamic Range"
                fracdiff_fixed_title = "Ratio \nDev/Ref, Fixed Range"

        # ==============================================================
        # Bundle variables for 6 parallel plotting calls
        # 0 = Ref                 1 = Dev
        # 2 = Dynamic abs diff    3 = Restricted abs diff
        # 4 = Dynamic frac diff   5 = Restricted frac diff
        # ==============================================================

        subplots = [
            "ref", "dev",
            "dyn_abs_diff", "res_abs_diff",
            "dyn_frac_diff", "res_frac_diff",
        ]
        if diff_of_diffs:
            subplots =  ["ref", "dev",
                         "dyn_abs_diff", "res_abs_diff",
                         "dyn_abs_diff", "res_abs_diff"]

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
            extents = [cmp_extent, cmp_extent,
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
            cmaps = [ref_cmap, dev_cmap, cmap_gray, \
                     cmap_gray, cmap_gray, cmap_gray]
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

        unit_list = [ds_ref.units, ds_dev.units, cmn_units,
                     cmn_units, "unitless", "unitless"]

        other_all_nans = [dev_is_all_nan, ref_is_all_nan, \
                          False, False, False, False]

        mins = [vmin_ref, vmin_dev, vmin_abs]
        maxs = [vmax_ref, vmax_dev, vmax_abs]

        ratio_logs = [False, False, False, False, True, True]

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
                plot_type="single_level",
                ratio_log=ratio_logs[i]
            )

        # ==============================================================
        # Add this page of 6-panel plots to a PDF file
        # ==============================================================
        if savepdf:

            folders = pdfname.split('/')
            full_path = ''
            for folder in folders[:-1]:
                full_path = os.path.join(temp_dir, full_path, folder)
                if not os.path.isdir(full_path):
                    os.mkdir(full_path)
            pdf = PdfPages(os.path.join(temp_dir, pdfname + "BENCHMARKFIGCREATION.pdf" + str(ivar)))
            pdf.savefig(figs)
            pdf.close()
            plt.close(figs)
        # ==============================================================
        # Update the list of variables with significant differences.
        # Criterion: abs(1 - max(fracdiff)) > 0.1
        # Do not include NaNs in the criterion, because these indicate
        # places where fracdiff could not be computed (div-by-zero).
        # ==============================================================
        if np.abs(1 - np.nanmax(fracdiff)) > 0.1:
            sigdiff_list.append(varname)
            return varname
        else:
            return

    # ==================================================================
    # Call figure generation function in a parallel loop over variables
    # ==================================================================
    #do not attempt nested thread parallelization due to issues with matplotlib
    if current_process().name != "MainProcess":
        n_job = 1

    if not savepdf:
        #disable parallel plotting to allow interactive figure plotting
        for i in range(n_var):
            createfig(i)

    else:
        with TemporaryDirectory() as temp_dir:
            results = Parallel(n_jobs = n_job) (delayed(createfig)(i, temp_dir) \
                                                for i in range(n_var))
            #update sig diffs after parallel calls
            if current_process().name == "MainProcess":
                for varname in results:
                    if type(varname) is str:
                        sigdiff_list.append(varname)

            # ==================================================================
            # Finish
            # ==================================================================
            if verbose:
                print("Closed PDF")
            merge = PdfFileMerger()
            #print("Creating {} for {} variables".format(pdfname, n_var))
            pdf = PdfPages(pdfname)
            pdf.close()
            for i in range(n_var):
                merge.append(os.path.join(temp_dir, pdfname + "BENCHMARKFIGCREATION.pdf" + str(i)))
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
    refmet=None,
    devmet=None,
    weightsdir='.',
    pdfname="",
    cmpres=None,
    match_cbar=True,
    pres_range=[0, 2000],
    normalize_by_area=False,
    enforce_units=True,
    convert_to_ugm3=False,
    flip_ref=False,
    flip_dev=False,
    use_cmap_RdBu=False,
    verbose=False,
    log_color_scale=False,
    log_yaxis=False,
    extra_title_txt=None,
    n_job=-1,
    sigdiff_list=[],
    second_ref=None,
    second_dev=None,
    spcdb_dir=os.path.dirname(__file__)
):

    """
    Create single-level 3x2 comparison zonal-mean plots for variables
    common in two xarray Daatasets. Optionally save to PDF.
    Args:
    -----
        refdata : xarray dataset
            Dataset used as reference in comparison
        refstr  : str OR list of str
            String description for reference data to be used in plots
            OR list containing [ref1str, ref2str] for diff-of-diffs plots
        devdata : xarray dataset
            Dataset used as development in comparison
        devstr  : str OR list of str
            String description for development data to be used in plots
            OR list containing [dev1str, dev2str] for diff-of-diffs plots
    Keyword Args (optional):
    ------------------------
        varlist : list of strings
            List of xarray dataset variable names to make plots for
            Default value: None (will compare all common 3D variables)
        itime : integer
            Dataset time dimension index using 0-based system
            Default value: 0
        refmet : xarray dataset
            Dataset containing ref meteorology
            Default value: None
        devmet : xarray dataset
            Dataset containing dev meteorology  
            Default value: None
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
            Ref and Dev datasets by grid area. Input ref and dev datasets
            must include AREA variable in m2 if normalizing by area.
            Default value: False
        enforce_units : boolean
            Set this flag to True force an error if the variables in
            the Ref and Dev datasets have different units.
            Default value: True
        convert_to_ugm3 : str
            Whether to convert data units to ug/m3 for plotting.
            Default value: False
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
        n_job : int
            Defines the number of simultaneous workers for parallel plotting.
            Set to 1 to disable parallel plotting. Value of -1 allows the application to decide.
            Default value: -1
        sigdiff_list: list of str
            Returns a list of all quantities having significant
            differences (where |max(fractional difference)| > 0.1).
            Default value: []
        second_ref : xarray Dataset
            A dataset of the same model type / grid as refdata, to be used in diff-of-diffs plotting.
            Default value: None    
        second_dev : xarray Dataset
            A dataset of the same model type / grid as devdata, to be used in diff-of-diffs plotting.
            Default value: None
        spcdb_dir  : str 
            Directory containing species_database.yml file.
            Default value: Path of GCPy code repository
    """
    warnings.showwarning = warning_format
    if not isinstance(refdata, xr.Dataset):
        raise TypeError("The refdata argument must be an xarray Dataset!")

    if not isinstance(devdata, xr.Dataset):
        raise TypeError("The devdata argument must be an xarray Dataset!")

    # Determine if doing diff-of-diffs
    if second_ref is not None and second_dev is not None:
        diff_of_diffs = True
    else:
        diff_of_diffs = False

    # Prepare diff-of-diffs datasets if needed
    if diff_of_diffs:
        refdata, fracrefdata = get_diff_of_diffs(refdata, second_ref)
        devdata, fracdevdata = get_diff_of_diffs(devdata, second_dev)
        if type(refstr) is list:
            frac_refstr = '{} / {}'.format(refstr[1], refstr[0])
            refstr = '{} - {}'.format(refstr[1], refstr[0])
        else:
            frac_refstr = 'Ref2 / Ref1'
        if type(devstr) is list:
            frac_devstr = '{} / {}'.format(devstr[1], devstr[0])
            devstr = '{} - {}'.format(devstr[1], devstr[0])
        else:
            frac_devstr = 'Dev2 / Dev1'
        
    # If no varlist is passed, plot all 3D variables in the dataset
    if varlist is None:
        quiet = not verbose
        vardict = compare_varnames(refdata, devdata, quiet=quiet)
        varlist = vardict["commonvars3D"]
        print("Plotting all 3D variables")
    n_var = len(varlist)

    # Exit out if there are no 3D variables
    if not n_var:
        print("WARNING: no 3D variables to plot zonal mean for!")
        return

    # If no PDF name passed, then do not save to PDF
    savepdf = True
    if pdfname == "":
        savepdf = False
    # Cleanup previous temporary PDFs produced during parallelization
    #for i in range(n_var):
    #    try:
    #        os.remove(os.path.join(temp_dir, pdfname + "BENCHMARKFIGCREATION.pdf" + str(i)))
    #    except:
    #        continue

    # If converting to ug/m3, load the species database
    if convert_to_ugm3:
        properties_path = os.path.join(spcdb_dir, "species_database.yml")
        properties = yaml.load(open(properties_path), Loader=yaml.FullLoader)

    # Get mid-point pressure and edge pressures for this grid (assume 72-level)
    ref_pedge, ref_pmid, ref_grid_cat = get_vert_grid(refdata)
    dev_pedge, dev_pmid, dev_grid_cat = get_vert_grid(devdata)

    # Get indexes of pressure subrange (full range is default)
    ref_pedge_ind = get_pressure_indices(ref_pedge, pres_range)
    dev_pedge_ind = get_pressure_indices(dev_pedge, pres_range)


    # Pad edges if subset does not include surface or TOA so data spans
    # entire subrange
    ref_pedge_ind = pad_pressure_edges(ref_pedge_ind, refdata.sizes["lev"])
    dev_pedge_ind = pad_pressure_edges(dev_pedge_ind, devdata.sizes["lev"])

    # pmid indexes do not include last pedge index
    ref_pmid_ind = ref_pedge_ind[:-1]
    dev_pmid_ind = dev_pedge_ind[:-1]
    # Convert levels to pressures in ref and dev data
    refdata = convert_lev_to_pres(refdata, ref_pmid, ref_pedge)
    devdata = convert_lev_to_pres(devdata, dev_pmid, dev_pedge)
    
    if diff_of_diffs:
        fracrefdata = convert_lev_to_pres(fracrefdata, ref_pmid, ref_pedge)
        fracdevdata = convert_lev_to_pres(fracdevdata, dev_pmid, dev_pedge)

    # ==================================================================
    # Reduce pressure range if reduced range passed as input. Indices
    # must be flipped if flipping vertical axis.
    # ==================================================================
    #this may require checking for 48 / 73 levels
    ref_pmid_ind_flipped = refdata.sizes["lev"] - ref_pmid_ind[::-1] - 1
    dev_pmid_ind_flipped = devdata.sizes["lev"] - dev_pmid_ind[::-1] - 1
    if flip_ref:
        ref_pmid_ind = pmid_ind_flipped
    if flip_dev:
        dev_pmid_ind = pmid_ind_flipped

    refdata = refdata.isel(lev=ref_pmid_ind)
    devdata = devdata.isel(lev=dev_pmid_ind)
    if diff_of_diffs:
        fracrefdata = fracrefdata.isel(lev=ref_pmid_ind)
        fracdevdata = fracdevdata.isel(lev=dev_pmid_ind)

    [refres, refgridtype, devres, devgridtype, cmpres, cmpgridtype,
     regridref, regriddev, regridany, refgrid, devgrid, cmpgrid,
     refregridder, devregridder,refregridder_list, devregridder_list] = \
            create_regridders(
                refdata, 
                devdata, 
                weightsdir=weightsdir,
                cmpres=cmpres,
                zm=True
            )

    ref_to_47 = False
    dev_to_47 = False

    # Determine whether one grid must be reduced to 47 levels and assign
    # plotting variables
    if ref_grid_cat is 72 and dev_grid_cat is 47:
        ref_to_47 = True
        pedge = dev_pedge
        pedge_ind = dev_pedge_ind
    elif ref_grid_cat is 47 and dev_grid_cat is 72:
        dev_to_47 = True
        pedge = ref_pedge
        pedge_ind = ref_pedge_ind
    else:
        pedge = ref_pedge
        pedge_ind = ref_pedge_ind

    if ref_to_47 or dev_to_47:
        conv_dict = dict_72_to_47()

    # ==================================================================
    # Loop over all variables
    # ==================================================================
    ds_refs = [None] * n_var
    frac_ds_refs = [None] * n_var
    ds_devs = [None] * n_var
    frac_ds_devs = [None] * n_var
    for i in range(n_var):

        varname = varlist[i]

        # ==================================================================
        # Slice the data, allowing for no time dimension (bpch)
        # ==================================================================

        # Ref
        if "time" in refdata[varname].dims:
            ds_refs[i] = refdata[varname].isel(time=itime)
            if diff_of_diffs:
                frac_ds_refs[i] = fracrefdata[varname].isel(time=itime)
        else:
            ds_refs[i] = refdata[varname]
            if diff_of_diffs:
                frac_ds_refs[i] = fracrefdata[varname]
        
        # Dev
        if "time" in devdata[varname].dims:
            ds_devs[i] = devdata[varname].isel(time=itime)
            if diff_of_diffs:
                frac_ds_devs[i] = fracdevdata[varname].isel(time=itime)

        else:
            ds_devs[i] = devdata[varname]
            if diff_of_diffs:
                frac_ds_devs[i] = fracdevdata[varname]

        # ==================================================================
        #  Handle units as needed
        # ==================================================================
        
        # Convert to ppb if units string is variation of mol/mol
        if data_unit_is_mol_per_mol(ds_refs[i]):
            ds_refs[i].values = ds_refs[i].values * 1e9
            ds_refs[i].attrs["units"] = "ppb"
        if data_unit_is_mol_per_mol(ds_devs[i]):
            ds_devs[i].values = ds_devs[i].values * 1e9
            ds_devs[i].attrs["units"] = "ppb"

        # If units string is ppbv (true for bpch data) then rename units
        if ds_refs[i].units.strip() == "ppbv":
            ds_refs[i].attrs["units"] = "ppb"
        if ds_devs[i].units.strip() == "ppbv":
            ds_devs[i].attrs["units"] = "ppb"

        # If units string is W/m2 (may be true for bpch data) then rename units
        if ds_refs[i].units.strip() == "W/m2":
            ds_refs[i].attrs["units"] = "W m-2"
        if ds_devs[i].units.strip() == "W/m2":
            ds_devs[i].attrs["units"] = "W m-2"

        # If units string is UNITLESS (may be true for bpch data) then rename units
        if ds_refs[i].units.strip() == "UNITLESS":
            ds_refs[i].attrs["units"] = "1"
        if ds_devs[i].units.strip() == "UNITLESS":
            ds_devs[i].attrs["units"] = "1"

        # Check that units are the same in ref and dev. Will exit with
        # an error if do not match and enforce_units is true (default).
        check_units(ds_refs[i], ds_devs[i])
            
        # Convert from ppb to ug/m3 if convert_to_ugm3 is passed as true
        if convert_to_ugm3:

            # Error checks: must pass met, not normalize by area, and be in ppb
            if refmet is None or devmet is None:
                msg = "Met mata ust be passed to convert units to ug/m3."
                raise ValueError(msg)
            elif normalize_by_area:
                msg = "Normalizing by area is now allowed if plotting ug/m3"
                raise ValueError(msg)
            elif ds_refs[i].units != "ppb" or ds_devs[i].units != "ppb":
                msg = "Units must be mol/mol if converting to ug/m3."
                raise ValueError(msg)

            # Slice air density data by time and lev
            # (assume same format and dimensions as refdata and devdata)
            if "time" in refmet["Met_AIRDEN"].dims:
                ref_airden = refmet["Met_AIRDEN"].isel(time=itime,
                                                       lev=ref_pmid_ind)
            else:
                ref_airden = refmet["Met_AIRDEN"].isel(lev=ref_pmid_ind)
            if "time" in devmet["Met_AIRDEN"].dims:
                dev_airden = devmet["Met_AIRDEN"].isel(time=itime,
                                                       lev=dev_pmid_ind)
            else:
                dev_airden = devmet["Met_AIRDEN"].isel(lev=dev_pmid_ind)

            # Get a list of properties for the given species
            spc_name = varname.replace(varname.split("_")[0]+"_", "")
            species_properties = properties.get(spc_name)

            # If no properties are found, then exit with an error.
            # Otherwise, get the molecular weight in g/mol.
            if species_properties is None:
                # Hack lumped species until we implement a solution
                if spc_name in ["Simple_SOA","Complex_SOA"]:
                    spc_mw_g = 150.0
                else:
                    msg = "No properties found for {}. Cannot convert" \
                          + " to ug/m3."
                    raise ValueError(msg.format(spc_name))
            else:
                # Get the species molecular weight in g/mol
                spc_mw_g = species_properties.get("MW_g")
                if spc_mw_g is None:
                    msg = "Molecular weight not found for for species {}!" \
                          + " Cannot convert to ug/m3.".format(spc_name)
                    raise ValueError(msg)

            # Convert values from ppb to ug/m3:
            # ug/m3 = 1e-9ppb * mol/g air * kg/m3 air * 1e3g/kg
            #         * g/mol spc * 1e6ug/g
            #       = ppb * air density * (spc MW / air MW)
            ds_refs[i].values = ds_refs[i].values * ref_airden.values \
                                * ( spc_mw_g / MW_AIR_g )
            ds_devs[i].values = ds_devs[i].values * dev_airden.values \
                                * ( spc_mw_g / MW_AIR_g )

            # Update units string
            ds_refs[i].attrs["units"] = "\u03BCg/m3" # ug/m3 using mu
            ds_devs[i].attrs["units"] = "\u03BCg/m3"

        # ==============================================================
        # Reshape cubed sphere data if using MAPL v1.0.0+
        # TODO: update function to expect data in this format
        # ==============================================================
        ds_refs[i] = reshape_MAPL_CS(ds_refs[i])
        ds_devs[i] = reshape_MAPL_CS(ds_devs[i])
        if diff_of_diffs:
            frac_ds_refs[i] = reshape_MAPL_CS(frac_ds_refs[i])
            frac_ds_devs[i] = reshape_MAPL_CS(frac_ds_devs[i])        

        # Flip in the vertical if applicable
        if flip_ref:
            ds_refs[i].data = ds_refs[i].data[::-1, :, :]
            if diff_of_diffs:
                frac_ds_refs[i].data = frac_ds_refs[i].data[::-1, :, :]
        if flip_dev:
            ds_devs[i].data = ds_devs[i].data[::-1, :, :]
            if diff_of_diffs:
                frac_ds_devs[i].data = frac_ds_devs[i].data[::-1, :, :]

    # ==================================================================
    # Get the area variables if normalize_by_area=True. They can be
    # either in the main datasets as variable AREA or in the optionally
    # passed meteorology datasets as Met_AREAM2.
    # ==================================================================
    if normalize_by_area:
        if "AREA" in refdata.data_vars.keys():
            ref_area = refdata["AREA"]
        elif refmet is not None:
            if "Met_AREAM2" in refmet.data_vars.keys():
                ref_area = refmet["Met_AREAM2"]
        else:
            msg = "normalize_by_area = True but AREA not " \
                + "present in the Ref dataset and ref met with Met_AREAM2" \
                + " not passed!"
            raise ValueError(msg)
        if "time" in ref_area.dims:
            ref_area = ref_area.isel(time=0)
        if refgridtype == 'cs':
            ref_area = reshape_MAPL_CS(ref_area)

        if "AREA" in devdata.data_vars.keys():
            dev_area = devdata["AREA"]
        elif devmet is not None:
            if "Met_AREAM2" in devmet.data_vars.keys():
                dev_area = devmet["Met_AREAM2"]
        else:
            msg = "normalize_by_area = True but AREA not " \
                + "present in the Dev dataset and dev met with Met_AREAM2" \
                | " not passed!"
            raise ValueError(msg)
        if "time" in dev_area.dims:
            dev_area = dev_area.isel(time=0)
        if devgridtype == 'cs':
            dev_area = reshape_MAPL_CS(dev_area)

        # Make sure the areas do not have a lev dimension
        if "lev" in ref_area.dims:
            ref_area = ref_area.isel(lev=0)
        if "lev" in dev_area.dims:
            dev_area = dev_area.isel(lev=0)

    # ==================================================================
    # Create arrays for each variable in the Ref and Dev dataset
    # and regrid to the comparison grid.
    # ==================================================================
    ds_ref_cmps = [None] * n_var
    ds_dev_cmps = [None] * n_var
    frac_ds_ref_cmps = [None] * n_var
    frac_ds_dev_cmps = [None] * n_var

    for i in range(n_var):
        ds_ref = ds_refs[i]
        ds_dev = ds_devs[i]
        frac_ds_ref = frac_ds_refs[i]
        frac_ds_dev = frac_ds_devs[i]
        # ==============================================================
        # Reduce variables to 47L from 72L if necessary for comparison
        # ==============================================================
        if ref_to_47:
            ds_ref = reduce_72_to_47(
                ds_ref, 
                conv_dict, 
                ref_pmid_ind, 
                dev_pmid_ind
            )
            if diff_of_diffs:
                frac_ds_ref = reduce_72_to_47(
                    frac_ds_ref,
                    ref_pmid_ind,
                    dev_pmid_ind
                )
        if dev_to_47:
            ds_dev = reduce_72_to_47(
                ds_dev,
                conv_dict,
                dev_pmid_ind,
                ref_pmid_ind
            )
            if diff_of_diffs:
                frac_ds_dev = reduce_72_to_47(
                    frac_ds_dev,
                    ref_pmid_ind,
                    dev_pmid_ind
                )

        ref_nlev = len(ds_ref['lev'])
        dev_nlev = len(ds_dev['lev'])

        # Do area normalization before regridding if normalize_by_area=True
        if normalize_by_area:
            exclude_list = ["WetLossConvFrac", "Prod_", "Loss_"]
            if not any(s in varname for s in exclude_list):
                ds_ref.values = ds_ref.values / ref_area.values
                ds_dev.values = ds_dev.values / dev_area.values
                ds_refs[i] = ds_ref
                ds_devs[i] = ds_dev
                if diff_of_diffs:
                    frac_ds_ref.values = frac_ds_ref.values / ref_area.values
                    frac_ds_refs[i] = frac_ds_ref
                    frac_ds_dev.values = frac_ds_dev.values / dev_area.values
                    frac_ds_devs[i] = frac_ds_dev
                
        # Ref
        ds_ref_cmps[i] = regrid_comparison_data(
            ds_ref,
            refres,
            regridref,
            refregridder,
            refregridder_list,
            cmpgrid,
            refgridtype,
            nlev=ref_nlev
        )
        if diff_of_diffs:
            frac_ds_ref_cmps[i] = regrid_comparison_data(
                frac_ds_refs[i],
                refres,
                regridref,
                refregridder,
                refregridder_list,
                cmpgrid,
                refgridtype,
                nlev=ref_nlev
            )
        # Dev
        ds_dev_cmps[i] = regrid_comparison_data(
            ds_dev,
            devres,
            regriddev,
            devregridder,
            devregridder_list,
            cmpgrid,
            devgridtype,
            nlev=dev_nlev
        )
        if diff_of_diffs:
            frac_ds_dev_cmps[i] = regrid_comparison_data(
                frac_ds_devs[i],
                devres,
                regriddev,
                devregridder,
                devregridder_list,
                cmpgrid,
                devgridtype,
                nlev=dev_nlev
            )

    # Universal plot setup
    xtick_positions = np.arange(-90, 91, 30)
    xticklabels = ["{}$\degree$".format(x) for x in xtick_positions]

    # ==================================================================
    # Define function to create a single page figure to be called
    # in a parallel loop
    # ==================================================================
    def createfig(ivar, temp_dir):

        # Suppress harmless run-time warnings (mostly about underflow)
        warnings.filterwarnings('ignore', category=RuntimeWarning)
        warnings.filterwarnings('ignore', category=UserWarning)

        if savepdf and verbose:
            print("{} ".format(ivar), end="")
        varname = varlist[ivar]
        varndim_ref = refdata[varname].ndim
        varndim_dev = devdata[varname].ndim

        # ==============================================================
        # Assign data variables
        # ==============================================================
        ds_ref = ds_refs[ivar]
        ds_dev = ds_devs[ivar]
        ds_ref_cmp = ds_ref_cmps[ivar]
        ds_dev_cmp = ds_dev_cmps[ivar]
        frac_ds_ref_cmp = frac_ds_ref_cmps[ivar]
        frac_ds_dev_cmp = frac_ds_dev_cmps[ivar]

        # ==============================================================
        # Area normalization units and subtitle
        # Set units and subtitle, including modification if normalizing
        # area. Note if enforce_units is False (non-default) then
        # units on difference plots will be wrong.
        # ==============================================================
        cmn_units = ds_ref.attrs["units"]
        varndim = varndim_ref
        subtitle_extra = ""
        if normalize_by_area:
            exclude_list = ["WetLossConvFrac", "Prod_", "Loss_"]
            if not any(s in varname for s in exclude_list):
                if "/" in cmn_units:
                    cmn_units = "{}/m2".format(cmn_units)
                else:
                    cmn_units = "{} m-2".format(cmn_units)
                ds_ref.attrs["units"] = cmn_units
                ds_dev.attrs["units"] = cmn_units
                subtitle_extra = ", Normalized by Area"

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
        if diff_of_diffs:
            frac_zm_dev_cmp = frac_ds_dev_cmp.mean(axis=2)
            frac_zm_ref_cmp = frac_ds_ref_cmp.mean(axis=2)
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
        if diff_of_diffs:
            zm_fracdiff = np.array(frac_zm_dev_cmp) -       \
                          np.array(frac_zm_ref_cmp)
        else:
            zm_fracdiff = np.abs(np.array(zm_dev_cmp)) /    \
                          np.abs(np.array(zm_ref_cmp))
        zm_fracdiff = np.where(np.abs(zm_fracdiff) ==       \
                               np.inf, np.nan, zm_fracdiff)
        zm_fracdiff[zm_fracdiff>1e308] = np.nan
        # Test if the frac. diff is zero everywhere or NaN everywhere
        fracdiff_is_all_zero = not np.any(zm_fracdiff) or       \
                               (np.nanmin(zm_fracdiff) == 0 and \
                                np.nanmax(zm_fracdiff) ==0)
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
            figs.suptitle("{}, Zonal Mean".format(varname),
                          fontsize=fontsize, y=offset)

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
        else:
            ref_title = "{} (Ref){}\n{} regridded from c{}".format(
                refstr, subtitle_extra, cmpres, refres
            )

        if devgridtype == "ll":
            dev_title = "{} (Dev){}\n{}".format(devstr, subtitle_extra, devres)
        else:
            dev_title = "{} (Dev){}\n{} regridded from c{}".format(
                devstr, subtitle_extra, cmpres, devres)

        if regridany:
            absdiff_dynam_title = \
                    "Difference ({})\nDev - Ref, Dynamic Range".format(cmpres)
            absdiff_fixed_title = \
                    "Difference ({})\nDev - Ref, Restricted Range [5%,95%]".\
                    format(cmpres)
            if diff_of_diffs:
                fracdiff_dynam_title = \
                    "Difference ({}), Dynamic Range\n{} - {}".\
                    format(cmpres,frac_devstr,frac_refstr)
                fracdiff_fixed_title = \
                    "Difference ({}), Restricted Range [5%,95%]\n{} - {}".\
                    format(cmpres, frac_devstr, frac_refstr)
            else:
                fracdiff_dynam_title = \
                    "Ratio ({})\nDev/Ref, Dynamic Range".format(cmpres)
                fracdiff_fixed_title = \
                    "Ratio ({})\nDev/Ref, Fixed Range".format(cmpres)
        else:
            absdiff_dynam_title = "Difference\nDev - Ref, Dynamic Range"
            absdiff_fixed_title = \
                    "Difference\nDev - Ref, Restricted Range [5%,95%]"
            if diff_of_diffs:
                fracdiff_dynam_title = \
                    "Difference, Dynamic Range\n{} - {}".\
                    format(frac_devstr,frac_refstr)
                fracdiff_fixed_title = \
                    "Difference, Restricted Range [5%,95%]\n{} - {}".\
                    format(frac_devstr, frac_refstr)
            else:
                fracdiff_dynam_title = "Ratio \nDev/Ref, Dynamic Range"
                fracdiff_fixed_title = "Ratio \nDev/Ref, Fixed Range"

        # ==============================================================
        # Bundle variables for 6 parallel plotting calls
        # 0 = Ref                 1 = Dev
        # 2 = Dynamic abs diff    3 = Restricted abs diff
        # 4 = Dynamic frac diff   5 = Restricted frac diff
        # ==============================================================

        subplots = [
            "ref", "dev",
            "dyn_abs_diff", "res_abs_diff",
            "dyn_frac_diff", "res_frac_diff",
        ]
        if diff_of_diffs:
            subplots =  ["ref", "dev",
                         "dyn_abs_diff", "res_abs_diff",
                         "dyn_abs_diff", "res_abs_diff"]

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

        plot_vals = [zm_ref, zm_dev, zm_diff, zm_diff,
                     zm_fracdiff, zm_fracdiff]

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

        unit_list = [ds_ref.units, ds_dev.units, cmn_units, cmn_units,
                     "unitless", "unitless"]

        other_all_nans = [dev_is_all_nan, ref_is_all_nan, \
                          False, False, False, False]

        gridtypes = [
            cmpgridtype,
            cmpgridtype,
            cmpgridtype,
            cmpgridtype,
            cmpgridtype,
            cmpgridtype,
        ]

        pedges = [ref_pedge, dev_pedge, pedge, pedge, pedge, pedge]

        pedge_inds = [ref_pedge_ind, dev_pedge_ind, pedge_ind, \
                      pedge_ind, pedge_ind, pedge_ind]

        mins = [vmin_ref, vmin_dev, vmin_abs]
        maxs = [vmax_ref, vmax_dev, vmax_abs]

        ratio_logs = [False, False, False, False, True, True]
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
                pedges[i],
                pedge_inds[i],
                log_yaxis,
                plot_type="zonal_mean",
                xtick_positions=xtick_positions,
                xticklabels=xticklabels,
                ratio_log=ratio_logs[i]
            )

        # ==============================================================
        # Add this page of 6-panel plots to the PDF file
        # ==============================================================
        if savepdf:
            folders = pdfname.split('/')
            full_path = ''
            for folder in folders[:-1]:
                full_path = os.path.join(temp_dir, full_path, folder)
                if not os.path.isdir(full_path):
                    os.mkdir(full_path)

            pdf = PdfPages(os.path.join(temp_dir, pdfname + "BENCHMARKFIGCREATION.pdf" + str(ivar)))
            pdf.savefig(figs)
            pdf.close()
            plt.close(figs)

        # ==============================================================
        # Update the list of variables with significant differences.
        # Criterion: abs(1 - max(fracdiff)) > 0.1
        # Do not include NaNs in the criterion, because these indicate
        # places where fracdiff could not be computed (div-by-zero).
        # ==============================================================
        if np.abs(1 - np.nanmax(zm_fracdiff)) > 0.1:
            sigdiff_list.append(varname)
            return varname
        else:
            return

    # ==================================================================
    # Call figure generation function in a parallel loop over variables
    # ==================================================================
    # do not attempt nested thread parallelization due to issues with
    # matplotlib
    if current_process().name != "MainProcess":
        n_job = 1

    if not savepdf:
        #disable parallel plotting to allow interactive figure plotting
        for i in range(n_var):
            createfig(i)

    else:
        with TemporaryDirectory() as temp_dir:
            results = Parallel(n_jobs = n_job) (delayed(createfig)(i, temp_dir) \
                                                for i in range(n_var))
            #update sig diffs after parallel calls
            if current_process().name == "MainProcess":
                for varname in results:
                    if type(varname) is str:
                        sigdiff_list.append(varname)

            # ==================================================================
            # Finish
            # ==================================================================
            if verbose:
                print("Closed PDF")
            merge = PdfFileMerger()
            #print("Creating {} for {} variables".format(pdfname, n_var))
            pdf = PdfPages(pdfname)
            pdf.close()
            for i in range(n_var):
                merge.append(os.path.join(temp_dir, pdfname + "BENCHMARKFIGCREATION.pdf" + str(i)))
            merge.write(pdfname)
            merge.close()
            warnings.showwarning = warning_format

def normalize_colors(vmin, vmax, is_difference=False, log_color_scale=False, ratio_log=False):
    """
    Normalizes a data range to the colormap range used by matplotlib
    functions.
    For log-color scales, special handling is done to prevent
    taking the log of data that is all zeroes.
    Args:
    -----
        vmin : float
            Minimum value of the data range.
        vmax : float
            Maximum value of the data range.
    Keyword Args:
    -------------
        is_difference : boolean
            Set this switch to denote that we are using a difference
            color scale (i.e. with zero in the middle of the range).
            Default value: False
        log_color_scale : boolean
            Logical flag to denote that we are using a logarithmic
            color scale instead of a linear color scale.
            Default value: False
    Returns:
    --------
        norm : matplotlib Norm
            The normalized matplotlib color range, stored in
            a matplotlib Norm object.
    Remarks:
    --------
         For log color scales, we will use a range of 3 orders of
         magnitude (i.e. from vmax/1e3 to vmax).
    """

    #Define class for logarithmic non-symmetric color scheme
    class MidpointLogNorm(mcolors.LogNorm):
        def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
            mcolors.LogNorm.__init__(self,vmin=vmin,vmax=vmax,clip=clip)
            self.midpoint=midpoint
        def __call__(self,value,clip=None):
            result, is_scalar = self.process_value(value)
            x = [np.log(self.vmin), np.log(self.midpoint), np.log(self.vmax)]
            y = [0, 0.5, 1]
            return np.ma.array(np.interp(np.log(value), x, y), mask=result.mask, copy = False)
            
    if (abs(vmin) == 0 and abs(vmax) == 0) or (np.isnan(vmin) and np.isnan(vmax)):
        # If the data is zero everywhere (vmin=vmax=0) or undefined
        # everywhere (vmin=vmax=NaN), then normalize the data range
        # so that the color corresponding to zero (white) will be
        # placed in the middle of the colorbar, where we will
        # add a single tick.
        if is_difference:
            return mcolors.Normalize(vmin=-1.0, vmax=1.0)
        else:
            return mcolors.Normalize(vmin=0.0, vmax=1.0)

    else:
        
        # For log color scales, assume a range 3 orders of magnitude
        # below the maximum value.  Otherwise use a linear scale.
        if log_color_scale and not ratio_log:
            return mcolors.LogNorm(vmin=vmax / 1e3, vmax=vmax)
        elif log_color_scale:
            return MidpointLogNorm(vmin=vmin, vmax=vmax, midpoint=1)
        else:
            return mcolors.Normalize(vmin=vmin, vmax=vmax)

def gcplot(plot_vals,
           ax=None,
           plot_type="single_level",
           grid={},
           gridtype="",
           title="fill",
           comap=WhGrYlRd,
           norm=[],
           unit="",
           extent=(None, None, None, None),
           masked_data=None,
           use_cmap_RdBu=False,
           log_color_scale=False,
           add_cb=True,
           pres_range=[0, 2000],
           pedge=np.full((1, 1), -1),
           pedge_ind=np.full((1,1), -1),
           log_yaxis=False,
           xtick_positions=[],
           xticklabels=[],
           proj=ccrs.PlateCarree()
):
    """
    Core plotting routine -- creates a single plot panel.
    Args:
    -----
        plot_vals : xarray DataArray or numpy array
            Single data variable GEOS-Chem output to plot
    Keyword Args (Optional):
    -----
    
        ax : matplotlib axes
            Axes object to plot information. Will create a new axes if None is passed.
        plot_type : str
            Either "single_level" or "zonal_mean"
        grid : dict
            Dictionary mapping plot_vals to plottable coordinates
        gridtype : str
            "ll" for lat/lon or "cs" for cubed-sphere
        title : str
            Title to put at top of plot    
        comap : matplotlib Colormap
            Colormap for plotting data values
        norm : list
            List with range [0..1] normalizing color range for matplotlib methods
        unit : ""
            Units of plotted data
        extent : tuple (minlon, maxlon, minlat, maxlat)
            Describes minimum and maximum latitude and longitude of input data
        masked_data : numpy array
            Masked area for cubed-sphere plotting
        use_cmap_RdBu : boolean
            Set this flag to True to use a blue-white-red colormap
        log_color_scale : boolean
            Set this flag to True to use a log-scale colormap
        add_cb : boolean
            Set this flag to True to add a colorbar to the plot 
        pres_range : list(int)
            Range from minimum to maximum pressure for zonal mean plotting
        pedge : numpy array
            Edge pressures of vertical grid cells in plot_vals
        pedge_ind : numpy array
            Index of edge pressure values within pressure range in plot_vals
        log_yaxis : boolean
            Set this flag to True to enable log scaling of pressure in zonal mean plots
        xtick_positions : list(float)
            Locations of lat/lon or lon ticks on plot
        xticklabels: list(str)
            Labels for lat/lon ticks
    Returns:
    -----
    
    plot : matplotlib plot
        Plot object created from input
    """

    #Eliminate 1D level or time dimensions
    plot_vals=plot_vals.squeeze()
    data_is_xr = type(plot_vals) is xr.DataArray
    if extent == (None, None, None, None) or extent == None:
        extent = get_grid_extents(grid)
        #convert to -180 to 180 grid if needed (necessary if going cross-dateline later)
        if extent[0] > 180 or extent[1] > 180:
            extent = [((extent[0]+180)%360)-180, ((extent[1]+180)%360)-180, extent[2], extent[3]]
    #Account for cross-dateline extent
    if extent[0] > extent[1]:
        if gridtype == "ll":
            proj = ccrs.PlateCarree(central_longitude=180)
            extent[0] = extent[0]%360-180        
            extent[1] = extent[1]%360-180
            plot_vals = plot_vals.assign_coords(lon=plot_vals.lon%360-180)
            plot_vals = plot_vals.sortby(plot_vals.lon)
        else:
            proj = ccrs.PlateCarree(central_longitude=180)
            extent[0] = extent[0]%360-180        
            extent[1] = extent[1]%360-180
            grid["lon_b"]=grid["lon_b"]-180
            
    if xtick_positions == []:
        #if plot_type == "single_level":
        #    xtick_positions = np.arange(extent[0], extent[1], (extent[1]-extent[0])/12)
        if plot_type == "zonal_mean":
            xtick_positions = np.arange(-90, 90, 30)

    if xticklabels == []:
        xticklabels = ["{}$\degree$".format(x) for x in xtick_positions]

    if unit == "" and data_is_xr:
        unit = plot_vals.units.strip()

    if ax == None:
        if plot_type == "zonal_mean":
            ax = plt.axes()
        if plot_type == "single_level":
            ax = plt.axes(projection = proj)

    if title == "fill" and data_is_xr:
        title = plot_vals.name

    #Generate grid if not passed
    if grid == {}:
        res, gridtype = get_input_res(plot_vals)

        if plot_type == 'single_level':
            [grid, _] = call_make_grid(res, gridtype, False, False)

        else: #zonal mean
            if np.all(pedge_ind == -1) or np.all(pedge == -1):
                
                # Get mid-point pressure and edge pressures for this grid (assume 72-level)
                pedge, pmid, grid_cat = get_vert_grid(plot_vals)
                
                # Get indexes of pressure subrange (full range is default)
                pedge_ind = get_pressure_indices(pedge, pres_range)
                
                # Pad edges if subset does not include surface or TOA so data spans
                # entire subrange
                pedge_ind = pad_pressure_edges(pedge_ind, plot_vals.sizes["lev"])
                
                # pmid indexes do not include last pedge index
                pmid_ind = pedge_ind[:-1]
                
                # Convert levels to pressures in ref and dev data
                plot_vals = convert_lev_to_pres(plot_vals, pmid, pedge)
                
                #get proper levels
                plot_vals = plot_vals.isel(lev=pmid_ind)            
                
            [input_res, input_gridtype, _, _,
             new_res, new_gridtype, regrid, _, _, _, _, 
             grid, regridder, _, regridder_list, _] = create_regridders(
                                                      plot_vals, 
                                                      plot_vals, 
                                                      weightsdir='.',
                                                      cmpres=None,
                                                      zm=True
                                                   )
            if gridtype == 'cs':
                plot_vals = reshape_MAPL_CS(plot_vals)
                nlev = len(plot_vals['lev'])
            
                # Ref
                plot_vals = regrid_comparison_data(
                    plot_vals,
                    input_res,
                    regrid,
                    regridder,
                    regridder_list,
                    grid,
                    input_gridtype,
                    nlev=nlev
                )
            #calculate zonal means
            plot_vals = plot_vals.mean(axis=2)


    data_is_xr = type(plot_vals) is xr.DataArray
    # Normalize colors (put into range [0..1] for matplotlib methods)
    if norm == []:
        if data_is_xr:
            vmin = plot_vals.data.min()
            vmax = plot_vals.data.max()
        elif type(plot_vals) is np.ndarray:
            vmin = np.min(plot_vals)
            vmax = np.max(plot_vals)
        norm = normalize_colors(
            vmin, vmax, is_difference=use_cmap_RdBu, log_color_scale=log_color_scale
        )


    # Create plot
    ax.set_title(title)
    if plot_type == "zonal_mean":
        #Zonal mean plot
        plot = ax.pcolormesh(
            grid["lat_b"], pedge[pedge_ind], plot_vals, cmap=comap, norm=norm
        )
        ax.set_aspect("auto")
        ax.set_ylabel("Pressure (hPa)")
        if log_yaxis:
            ax.set_yscale("log")
            ax.yaxis.set_major_formatter(
                mticker.FuncFormatter(lambda y, _: "{:g}".format(y))
            )
        ax.invert_yaxis()
        ax.set_xticks(xtick_positions)
        ax.set_xticklabels(xticklabels)

    elif gridtype == "ll":
        #Lat/Lon single level
        if type(plot_vals) is xr.DataArray:
            [minlon, maxlon, minlat, maxlat] = extent
            #filter data by bounds of extent
            plot_vals = plot_vals.where(plot_vals.lon>=minlon, 
                                        drop=True).where(plot_vals.lon<=maxlon,
                                                         drop=True).where(plot_vals.lat>=minlat, 
                                                                          drop=True).where(plot_vals.lat<=maxlat, drop=True)
        else:
            #for numpy arrays
            [minlon, maxlon, minlat, maxlat] = extent
            minlon_ind = np.where(grid["lon"] >= minlon)[0][0]
            maxlon_ind = np.where(grid["lon"] <= maxlon)[0][-1]
            minlat_ind = np.where(grid["lat"] >= minlat)[0][0]
            maxlat_ind = np.where(grid["lat"] <= maxlat)[0][-1]
            #assume lat comes first in indexing
            plot_vals = plot_vals[minlat_ind:maxlat_ind+1,minlon_ind:maxlon_ind+1].squeeze()
        # Create a lon/lat plot
        plot = ax.imshow(
            plot_vals, extent=extent, transform=proj, cmap=comap, norm=norm
        )
        ax.coastlines()
        ax.set_xticks(xtick_positions)
        ax.set_xticklabels(xticklabels)
        
    else:
        #Cubed-sphere single level
        ax.coastlines()
        try:
            if masked_data == None:
                masked_data = np.ma.masked_where(np.abs(grid["lon"] - 180) < 2, plot_vals.data.reshape(6, res, res))
        except ValueError:
            #Comparison of numpy arrays throws errors
            pass
        [minlon,maxlon,minlat,maxlat] = extent
        #Catch issue with plots extending into both the western and eastern hemisphere
        if np.max(grid["lon_b"] > 180):
            grid["lon_b"] = (((grid["lon_b"]+180)%360)-180)
        for j in range(6):
            plot = ax.pcolormesh(
                grid["lon_b"][j, :, :],
                grid["lat_b"][j, :, :],
                masked_data[j, :, :],
                transform=proj,
                cmap=comap,
                norm=norm
            )
        ax.set_xlim(minlon, maxlon)
        ax.set_ylim(minlat, maxlat)
        ax.set_xticks(xtick_positions)
        ax.set_xticklabels(xticklabels)

    if add_cb == True:
        cb = plt.colorbar(plot, ax=ax, orientation="horizontal", pad=0.10)
        cb.mappable.set_norm(norm)
        if data_is_xr:
            all_zero, all_nan = all_zero_or_nan(plot_vals.values)
        else:
            all_zero, all_nan = all_zero_or_nan(plot_vals)
        if all_zero or all_nan:
            if use_cmap_RdBu:
                cb.set_ticks([0.0])
            else:
                cb.set_ticks([0.5])
            if all_nan:
                cb.set_ticklabels(["Undefined throughout domain"])
            else:
                cb.set_ticklabels(["Zero throughout domain"])
        else:
            if log_color_scale:
                cb.formatter = mticker.LogFormatter(base=10)
            else:
                if (vmax - vmin) < 0.1 or (vmax - vmin) > 100:
                    cb.locator = mticker.MaxNLocator(nbins=4)

        try:
            cb.formatter.set_useOffset(False)
        except:
            #not all automatically chosen colorbar formatters properly handle the above method
            pass
        cb.update_ticks()
        cb.set_label(unit)

    return plot

