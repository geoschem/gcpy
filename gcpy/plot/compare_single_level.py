"""
compare_single_level.py: Function to create a six-panel plot comparing
quantities at a single model level for two different model versions.
Called from the GEOS-Chem benchmarking scripts and from the
compare_diags.py example script.
"""
import os
import copy
import warnings
from multiprocessing import current_process
from tempfile import TemporaryDirectory
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
from joblib import Parallel, delayed
from pypdf import PdfMerger
from gcpy.grid import get_grid_extents, call_make_grid
from gcpy.regrid import regrid_comparison_data, create_regridders
from gcpy.util import reshape_MAPL_CS, get_diff_of_diffs, \
    all_zero_or_nan, slice_by_lev_and_time, compare_varnames, \
    read_config_file, verify_variable_type
from gcpy.units import check_units, data_unit_is_mol_per_mol
from gcpy.constants import MW_AIR_g
from gcpy.plot.core import gcpy_style, six_panel_subplot_names, \
    _warning_format, WhGrYlRd
from gcpy.plot.six_plot import six_plot

# Suppress numpy divide by zero warnings to prevent output spam
np.seterr(divide="ignore", invalid="ignore")

# Use a style sheet to control plot attributes
plt.style.use(gcpy_style)


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
        extent=None,
        n_job=-1,
        sigdiff_list=None,
        second_ref=None,
        second_dev=None,
        spcdb_dir=os.path.dirname(__file__),
        sg_ref_path='',
        sg_dev_path='',
        ll_plot_func='imshow',
        **extra_plot_args
):
    """
    Create single-level 3x2 comparison map plots for variables common
    in two xarray Datasets. Optionally save to PDF.

    Args:
        refdata: xarray dataset
            Dataset used as reference in comparison
        refstr: str
            String description for reference data to be used in plots
        devdata: xarray dataset
            Dataset used as development in comparison
        devstr: str
            String description for development data to be used in plots

    Keyword Args (optional):
        varlist: list of strings
            List of xarray dataset variable names to make plots for
            Default value: None (will compare all common variables)
        ilev: integer
            Dataset level dimension index using 0-based system.
            Indexing is ambiguous when plotting differing vertical grids
            Default value: 0
        itime: integer
            Dataset time dimension index using 0-based system
            Default value: 0
        refmet: xarray dataset
            Dataset containing ref meteorology
            Default value: None
        devmet: xarray dataset
            Dataset containing dev meteorology
            Default value: None
        weightsdir: str
            Directory path for storing regridding weights
            Default value: None (will create/store weights in
            current directory)
        pdfname: str
            File path to save plots as PDF
            Default value: Empty string (will not create PDF)
        cmpres: str
            String description of grid resolution at which
            to compare datasets
            Default value: None (will compare at highest resolution
            of ref and dev)
        match_cbar: bool
            Set this flag to True if you wish to use the same colorbar
            bounds for the Ref and Dev plots.
            Default value: True
        normalize_by_area: bool
            Set this flag to True if you wish to normalize the Ref
            and Dev raw data by grid area. Input ref and dev datasets
            must include AREA variable in m2 if normalizing by area.
            Default value: False
        enforce_units: bool
            Set this flag to True to force an error if Ref and Dev
            variables have different units.
            Default value: True
        convert_to_ugm3: bool
            Whether to convert data units to ug/m3 for plotting.
            Default value: False
        flip_ref: bool
            Set this flag to True to flip the vertical dimension of
            3D variables in the Ref dataset.
            Default value: False
        flip_dev: bool
            Set this flag to True to flip the vertical dimension of
            3D variables in the Dev dataset.
            Default value: False
        use_cmap_RdBu: bool
            Set this flag to True to use a blue-white-red colormap
            for plotting the raw data in both the Ref and Dev datasets.
            Default value: False
        verbose: bool
            Set this flag to True to enable informative printout.
            Default value: False
        log_color_scale: bool
            Set this flag to True to plot data (not diffs)
            on a log color scale.
            Default value: False
        extra_title_txt: str
            Specifies extra text (e.g. a date string such as "Jan2016")
            for the top-of-plot title.
            Default value: None
        extent: list
            Defines the extent of the region to be plotted in form
            [minlon, maxlon, minlat, maxlat].
            Default value plots extent of input grids.
            Default value: [-1000, -1000, -1000, -1000]
        n_job: int
            Defines the number of simultaneous workers for parallel
            plotting.  Set to 1 to disable parallel plotting.
            Value of -1 allows the application to decide.
            Default value: -1
        sigdiff_list: list of str
            Returns a list of all quantities having significant
            differences (where |max(fractional difference)| > 0.1).
            Default value: None
        second_ref: xarray Dataset
            A dataset of the same model type / grid as refdata,
            to be used in diff-of-diffs plotting.
            Default value: None
        second_dev: xarray Dataset
            A dataset of the same model type / grid as devdata,
            to be used in diff-of-diffs plotting.
            Default value: None
        spcdb_dir: str
            Directory containing species_database.yml file.
            Default value: Path of GCPy code repository
        sg_ref_path: str
            Path to NetCDF file containing stretched-grid info
            (in attributes) for the ref dataset
            Default value: '' (will not be read in)
        sg_dev_path: str
            Path to NetCDF file containing stretched-grid info
            (in attributes) for the dev dataset
            Default value: '' (will not be read in)
        ll_plot_func: str
            Function to use for lat/lon single level plotting with
            possible values 'imshow' and 'pcolormesh'. imshow is much
            faster but is slightly displaced when plotting from
            dateline to dateline and/or pole to pole.
            Default value: 'imshow'
        extra_plot_args: various
            Any extra keyword arguments are passed through the
            plotting functions to be used in calls to pcolormesh() (CS)
            or imshow() (Lat/Lon).
    """
    warnings.showwarning = _warning_format
    # Error check arguments
    verify_variable_type(refdata, xr.Dataset)
    verify_variable_type(devdata, xr.Dataset)

    # Create empty lists for keyword arguments
    if extent is None:
        extent = [-1000, -1000, -1000, -1000]
    if sigdiff_list is None:
        sigdiff_list = []

    # Determine if doing diff-of-diffs
    diff_of_diffs = second_ref is not None and second_dev is not None

    # Prepare diff-of-diffs datasets if needed
    if diff_of_diffs:
        refdata, devdata = refdata.load(), devdata.load()
        second_ref, second_dev = second_ref.load(), second_dev.load()

#        # If needed, use fake time dim in case dates are different
#        # in datasets.  This needs more work for case of single versus
#        # multiple times.
#        aligned_time = [np.datetime64('2000-01-01')] * refdata.dims['time']
#        refdata = refdata.assign_coords({'time': aligned_time})
#        devdata = devdata.assign_coords({'time': aligned_time})
#        second_ref = second_ref.assign_coords({'time': aligned_time})
#        second_dev = second_dev.assign_coords({'time': aligned_time})

        refdata, fracrefdata = get_diff_of_diffs(refdata, second_ref)
        devdata, fracdevdata = get_diff_of_diffs(devdata, second_dev)
        frac_refstr = 'GCC_dev / GCC_ref'
        frac_devstr = 'GCHP_dev / GCHP_ref'
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
    if convert_to_ugm3:
        properties = read_config_file(
            os.path.join(
                spcdb_dir,
                "species_database.yml"
            ),
            quiet=True
        )

    sg_ref_params = [1, 170, -90]
    sg_dev_params = [1, 170, -90]
    # Get stretched-grid info if passed
    if sg_ref_path != '':
        sg_ref_attrs = xr.open_dataset(sg_ref_path).attrs
        sg_ref_params = [
            sg_ref_attrs['stretch_factor'],
            sg_ref_attrs['target_longitude'],
            sg_ref_attrs['target_latitude']]

    if sg_dev_path != '':
        sg_dev_attrs = xr.open_dataset(sg_dev_path).attrs
        sg_dev_params = [
            sg_dev_attrs['stretch_factor'],
            sg_dev_attrs['target_longitude'],
            sg_dev_attrs['target_latitude']]

    # Get grid info and regrid if necessary
    [refres, refgridtype, devres, devgridtype, cmpres, cmpgridtype, regridref,
     regriddev, regridany, refgrid, devgrid, cmpgrid, refregridder,
     devregridder, refregridder_list, devregridder_list] = create_regridders(
         refdata,
         devdata,
         weightsdir,
         cmpres=cmpres,
         sg_ref_params=sg_ref_params,
         sg_dev_params=sg_dev_params
    )

    # ==============================================================
    # Handle grid extents for lat-lon grids
    # ==============================================================

    # Get lat/lon extents, if applicable
    refminlon, refmaxlon, refminlat, refmaxlat = get_grid_extents(refgrid)
    devminlon, devmaxlon, devminlat, devmaxlat = get_grid_extents(devgrid)

    if -1000 not in extent:
        cmpminlon, cmpmaxlon, cmpminlat, cmpmaxlat = extent
    else:
        # Account for 0-360 coordinate scale
        uniform_refminlon, uniform_refmaxlon = refminlon, refmaxlon
        uniform_devminlon, uniform_devmaxlon = devminlon, devmaxlon
        if uniform_refmaxlon > 185:
            uniform_refminlon, uniform_refmaxlon = -180, 180
        if uniform_devmaxlon > 185:
            uniform_devminlon, uniform_devmaxlon = -180, 180

        cmpminlon, cmpmaxlon, cmpminlat, cmpmaxlat = \
            [np.max([(uniform_refminlon+180%360)-180, uniform_devminlon]),
             np.min([uniform_refmaxlon, uniform_devmaxlon]),
             np.max([refminlat, devminlat]),
             np.min([refmaxlat, devmaxlat])]

    # Set plot bounds for non cubed-sphere regridding and plotting
    # Pylint says ref_extent and dev_extent are not used
    #  -- Bob Yantosca (15 Aug 2023)
    #ref_extent = (refminlon, refmaxlon, refminlat, refmaxlat)
    #dev_extent = (devminlon, devmaxlon, devminlat, devmaxlat)
    cmp_extent = (cmpminlon, cmpmaxlon, cmpminlat, cmpmaxlat)
    # ==============================================================
    # Loop over all variables
    # ==============================================================
    ds_refs = [None] * n_var
    frac_ds_refs = [None] * n_var
    ds_devs = [None] * n_var
    frac_ds_devs = [None] * n_var
    for i in range(n_var):
        varname = varlist[i]
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

        # If units string is UNITLESS (may be true for bpch data) then rename
        # units
        if ds_refs[i].units.strip() == "UNITLESS":
            ds_refs[i].attrs["units"] = "1"
        if ds_devs[i].units.strip() == "UNITLESS":
            ds_devs[i].attrs["units"] = "1"

        # Compare units of ref and dev. The check_units function will throw an
        # error if units do not match and enforce_units is True.
        check_units(ds_refs[i], ds_devs[i], enforce_units)

        # Convert from ppb to ug/m3 if convert_to_ugm3 is passed as true
        if convert_to_ugm3:

            # Error checks: must pass met, not normalize by area, and be in ppb
            if refmet is None or devmet is None:
                msg = "Met mata ust be passed to convert units to ug/m3."
                raise ValueError(msg)
            if normalize_by_area:
                msg = "Normalizing by area is not allowed if plotting ug/m3"
                raise ValueError(msg)
            if ds_refs[i].units != "ppb" or ds_devs[i].units != "ppb":
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
            spc_name = varname.replace(varname.split("_")[0] + "_", "")
            species_properties = properties.get(spc_name)

            # If no properties are found, then exit with an error.
            # Otherwise, get the molecular weight in g/mol.
            if species_properties is None:
                # Hack lumped species until we implement a solution
                if spc_name in ["Simple_SOA", "Complex_SOA"]:
                    spc_mw_g = 150.0
                else:
                    msg = f"No properties found for {spc_name}. Cannot convert" \
                          + " to ug/m3."
                    raise ValueError(msg)
            else:
                spc_mw_g = species_properties.get("MW_g")
                if spc_mw_g is None:
                    msg = f"Molecular weight not found for species {spc_name}!" \
                          + " Cannot convert to ug/m3."
                    raise ValueError(msg)

            # Convert values from ppb to ug/m3:
            # ug/m3 = mol/mol * mol/g air * kg/m3 air * 1e3g/kg
            #         * g/mol spc * 1e6ug/g
            #       = ppb * air density * (spc MW / air MW)
            ds_refs[i].values = ds_refs[i].values * ref_airden.values \
                * (spc_mw_g / MW_AIR_g)
            ds_devs[i].values = ds_devs[i].values * dev_airden.values \
                * (spc_mw_g / MW_AIR_g)

            # Update units string
            ds_refs[i].attrs["units"] = "\u03BCg/m3"  # ug/m3 using mu
            ds_devs[i].attrs["units"] = "\u03BCg/m3"

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

    # ==============================================================
    # Reshape cubed sphere data if using MAPL v1.0.0+
    # TODO: update function to expect data in this format
    # ==============================================================

    for i in range(n_var):
        ds_refs[i] = reshape_MAPL_CS(ds_refs[i])
        ds_devs[i] = reshape_MAPL_CS(ds_devs[i])
        #ds_ref_cmps[i] = reshape_MAPL_CS(ds_ref_cmps[i])
        #ds_dev_cmps[i] = reshape_MAPL_CS(ds_dev_cmps[i])
        if diff_of_diffs:
            frac_ds_refs[i] = reshape_MAPL_CS(frac_ds_refs[i])
            frac_ds_devs[i] = reshape_MAPL_CS(frac_ds_devs[i])
            #frac_ds_ref_cmps[i] = reshape_MAPL_CS(frac_ds_ref_cmps[i])
            #frac_ds_dev_cmps[i] = reshape_MAPL_CS(frac_ds_dev_cmps[i])


    # ==================================================================
    # Create arrays for each variable in Ref and Dev datasets
    # and do any necessary horizontal regridding. 'cmp' stands for comparison
    # and represents ref and dev data regridded as needed to a common
    # grid type and resolution for use in difference and ratio plots.
    # ==================================================================
    ds_ref_cmps = [None] * n_var
    ds_dev_cmps = [None] * n_var
    frac_ds_ref_cmps = [None] * n_var
    frac_ds_dev_cmps = [None] * n_var

    global_cmp_grid = call_make_grid(cmpres, cmpgridtype)[0]
    # Get grid limited to cmp_extent for comparison datasets
    # Do not do this for cross-dateline plotting
    if cmp_extent[0] < cmp_extent[1]:
        regional_cmp_extent = cmp_extent
    else:
        regional_cmp_extent = [-180, 180, -90, 90]

    regional_cmp_grid = call_make_grid(cmpres, cmpgridtype,
                                       in_extent=[-180,180,-90,90],
                                       out_extent=regional_cmp_extent)[0]

    # Get comparison data extents in same midpoint format as lat-lon grid.
    cmp_mid_minlon, cmp_mid_maxlon, cmp_mid_minlat, cmp_mid_maxlat = \
        get_grid_extents(regional_cmp_grid, edges=False)

    cmpminlon_ind = np.where(global_cmp_grid["lon"] >= cmp_mid_minlon)[0][0]
    cmpmaxlon_ind = np.where(global_cmp_grid["lon"] <= cmp_mid_maxlon)[0][-1]
    cmpminlat_ind = np.where(global_cmp_grid["lat"] >= cmp_mid_minlat)[0][0]
    cmpmaxlat_ind = np.where(global_cmp_grid["lat"] <= cmp_mid_maxlat)[0][-1]

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
        ref_cs_res = refres
        dev_cs_res = devres
        if cmpgridtype == "cs":
            ref_cs_res = cmpres
            dev_cs_res = cmpres
        # Ref
        ds_ref_cmps[i] = regrid_comparison_data(
            ds_ref,
            ref_cs_res,
            regridref,
            refregridder,
            refregridder_list,
            global_cmp_grid,
            refgridtype,
            cmpgridtype,
            cmpminlat_ind,
            cmpmaxlat_ind,
            cmpminlon_ind,
            cmpmaxlon_ind
        )
        # Dev
        ds_dev_cmps[i] = regrid_comparison_data(
            ds_dev,
            dev_cs_res,
            regriddev,
            devregridder,
            devregridder_list,
            global_cmp_grid,
            devgridtype,
            cmpgridtype,
            cmpminlat_ind,
            cmpmaxlat_ind,
            cmpminlon_ind,
            cmpmaxlon_ind
        )
        # Diff of diffs
        if diff_of_diffs:
            frac_ds_ref_cmps[i] = regrid_comparison_data(
                frac_ds_refs[i],
                ref_cs_res,
                regridref,
                refregridder,
                refregridder_list,
                global_cmp_grid,
                refgridtype,
                cmpgridtype,
                cmpminlat_ind,
                cmpmaxlat_ind,
                cmpminlon_ind,
                cmpmaxlon_ind
            )
            frac_ds_dev_cmps[i] = regrid_comparison_data(
                frac_ds_devs[i],
                dev_cs_res,
                regriddev,
                devregridder,
                devregridder_list,
                global_cmp_grid,
                devgridtype,
                cmpgridtype,
                cmpminlat_ind,
                cmpmaxlat_ind,
                cmpminlon_ind,
                cmpmaxlon_ind
            )
    # =================================================================
    # Define function to create a single page figure to be called
    # in a parallel loop
    # =================================================================
    def createfig(ivar, temp_dir=''):

        # Suppress harmless run-time warnings (mostly about underflow)
        warnings.filterwarnings('ignore', category=RuntimeWarning)
        warnings.filterwarnings('ignore', category=UserWarning)

        if savepdf and verbose:
            print(f"{ivar} ", end="")
        varname = varlist[ivar]

        ds_ref = ds_refs[ivar]
        ds_dev = ds_devs[ivar]

        # ==============================================================
        # Set units and subtitle, including modification if normalizing
        # area. Note if enforce_units is False (non-default) then
        # units on difference plots will be wrong.
        # ==============================================================
        cmn_units = ds_ref.attrs["units"]
        subtitle_extra = ""
        if normalize_by_area:
            exclude_list = ["WetLossConvFrac", "Prod_", "Loss_"]
            if not any(s in varname for s in exclude_list):
                if "/" in cmn_units:
                    cmn_units = f"{cmn_units}/m2"
                else:
                    cmn_units = f"{cmn_units} m-2"
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
            def call_reshape(cmp_data):
                new_data = None
                if isinstance(cmp_data, xr.DataArray):
                    new_data = cmp_data.data.reshape(6, cmpres, cmpres)
                elif isinstance(cmp_data, np.ndarray):
                    new_data = cmp_data.reshape(6, cmpres, cmpres)
                return new_data

            ds_ref_cmp_reshaped = call_reshape(ds_ref_cmp)
            ds_dev_cmp_reshaped = call_reshape(ds_dev_cmp)
            frac_ds_ref_cmp_reshaped = call_reshape(frac_ds_ref_cmp)
            frac_ds_dev_cmp_reshaped = call_reshape(frac_ds_dev_cmp)

        # ==============================================================
        # Get min and max values for use in the colorbars
        # ==============================================================

        # Choose from values within plot extent
        if -1000 not in extent:
            min_max_extent = extent
        else:
            min_max_extent = cmp_extent
        # Find min and max lon
        min_max_minlon = np.min([min_max_extent[0], min_max_extent[1]])
        min_max_maxlon = np.max([min_max_extent[0], min_max_extent[1]])
        min_max_minlat = min_max_extent[2]
        min_max_maxlat = min_max_extent[3]

        def get_extent_for_colors(dset, minlon, maxlon, minlat, maxlat):
            ds_new = dset.copy()
            lat_var='lat'
            lon_var='lon'
            # Account for cubed-sphere data
            if 'lons' in ds_new.coords:
                lat_var='lats'
                lon_var='lons'
            if ds_new['lon'].max() > 190:
                minlon=minlon%360
                maxlon=maxlon%360
                # account for global plot
                if minlon == maxlon and maxlon == 180:
                    minlon = 0
                    maxlon = 360
            # account for cross dateline
            if minlon > maxlon:
                minlon, maxlon = maxlon, minlon

            # Add .compute() to force evaluation of ds_new[lon_var]
            # See https://github.com/geoschem/gcpy/issues/254
            # Also note: This may return as a dask.array.Array object
            return ds_new.where(\
                ds_new[lon_var].compute() >= minlon, drop=True).\
                where(ds_new[lon_var].compute() <= maxlon, drop=True).\
                where(ds_new[lat_var].compute() >= minlat, drop=True).\
                where(ds_new[lat_var].compute() <= maxlat, drop=True)

        ds_ref_reg = get_extent_for_colors(
            ds_ref,
            min_max_minlon,
            min_max_maxlon,
            min_max_minlat,
            min_max_maxlat
        )
        ds_dev_reg = get_extent_for_colors(
            ds_dev,
            min_max_minlon,
            min_max_maxlon,
            min_max_minlat,
            min_max_maxlat
        )

        # Ref
        vmin_ref = float(np.nanmin(ds_ref_reg.data))
        vmax_ref = float(np.nanmax(ds_ref_reg.data))

        # Dev
        vmin_dev = float(np.nanmin(ds_dev_reg.data))
        vmax_dev = float(np.nanmax(ds_dev_reg.data))

# Pylint says that these are unused variables, so comment out
#  -- Bob Yantosca (15 Aug 2023)
#        # Comparison
#        if cmpgridtype == "cs":
#            vmin_ref_cmp = float(np.nanmin(ds_ref_cmp))
#            vmax_ref_cmp = float(np.nanmax(ds_ref_cmp))
#            vmin_dev_cmp = float(np.nanmin(ds_dev_cmp))
#            vmax_dev_cmp = float(np.nanmax(ds_dev_cmp))
#            vmin_cmp = np.nanmin([vmin_ref_cmp, vmin_dev_cmp])
#            vmax_cmp = np.nanmax([vmax_ref_cmp, vmax_dev_cmp])
#        else:
#            vmin_cmp = np.nanmin([np.nanmin(ds_ref_cmp), np.nanmin(ds_dev_cmp)])
#            vmax_cmp = np.nanmax([np.nanmax(ds_ref_cmp), np.nanmax(ds_dev_cmp)])

        # Get overall min & max
        vmin_abs = np.nanmin([vmin_ref, vmin_dev])#, vmin_cmp])
        vmax_abs = np.nanmax([vmax_ref, vmax_dev])#, vmax_cmp])
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
        fracdiff[np.abs(fracdiff > 1e308)] = np.nan

        # Test if the frac. diff. is zero everywhere or NaN everywhere
        fracdiff_is_all_zero = not np.any(fracdiff) or       \
            (np.nanmin(fracdiff) == 0 and
             np.nanmax(fracdiff) == 0)
        fracdiff_is_all_nan = np.isnan(fracdiff).all() or ref_is_all_zero

        # For cubed-sphere, take special care to avoid a spurious
        # boundary line, as described here: https://stackoverflow.com/
        # questions/46527456/preventing-spurious-horizontal-lines-for-
        # ungridded-pcolormesh-data
        if cmpgridtype == "cs":
            fracdiff = np.ma.masked_where(np.abs(cmpgrid["lon"] - 180) < 2,
                                          fracdiff)

        # ==============================================================
        # Create 3x2 figure
        # ==============================================================

        # Create figures and axes objects
        # Also define the map projection that will be shown
        if extent[0] > extent[1]:
            proj = ccrs.PlateCarree(central_longitude=180)
        else:
            proj = ccrs.PlateCarree()
        figs, ((ax0, ax1), (ax2, ax3), (ax4, ax5)) = plt.subplots(
            3, 2, figsize=[12, 14],
            subplot_kw={"projection": proj}
        )
        # Ensure subplots don't overlap when invoking plt.show()
        if not savepdf:
            plt.subplots_adjust(hspace=0.4)
        # Give the figure a title
        offset = 0.96
        if "lev" in ds_ref.dims and "lev" in ds_dev.dims:
            if ilev == 0:
                levstr = "Surface"
            elif ilev == 22:
                levstr = "500 hPa"
            else:
                levstr = "Level " + str(ilev - 1)
            if extra_title_txt is not None:
                figs.suptitle(
                    f"{varname}, {levstr} ({extra_title_txt})",
                    y=offset,
                )
            else:
                figs.suptitle(
                    f"{varname}, {levstr}",
                    y=offset
                )
        elif (
            "lat" in ds_ref.dims
            and "lat" in ds_dev.dims
            and "lon" in ds_ref.dims
            and "lon" in ds_dev.dims
        ):
            if extra_title_txt is not None:
                figs.suptitle(
                    f"{varname} ({extra_title_txt})",
                    y=offset,
                )
            else:
                figs.suptitle(
                    f"{varname}",
                    y=offset)
        else:
            print(f"Incorrect dimensions for {varname}!")

        # ==============================================================
        # Set colormaps for data plots
        #
        # Use shallow copy (copy.copy() to create color map objects,
        # in order to avoid set_bad() from being applied to the base
        # color table. See: https://docs.python.org/3/library/copy.html
        # ==============================================================

        # Colormaps for 1st row (Ref and Dev)
        if use_cmap_RdBu:
            cmap_toprow_nongray = copy.copy(mpl.colormaps["RdBu_r"])
            cmap_toprow_gray = copy.copy(mpl.colormaps["RdBu_r"])
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
        cmap_nongray = copy.copy(mpl.colormaps["RdBu_r"])
        cmap_gray = copy.copy(mpl.colormaps["RdBu_r"])
        cmap_gray.set_bad(color="gray")

        # ==============================================================
        # Set titles for plots
        # ==============================================================

        if refgridtype == "ll":
            ref_title = f"{refstr} (Ref){subtitle_extra}\n{refres}"
        else:
            ref_title = f"{refstr} (Ref){subtitle_extra}\nc{refres}"

        if devgridtype == "ll":
            dev_title = f"{devstr} (Dev){subtitle_extra}\n{devres}"
        else:
            dev_title = f"{devstr} (Dev){subtitle_extra}\nc{devres}"
        if regridany:
            absdiff_dynam_title = \
                f"Difference ({cmpres})\nDev - Ref, Dynamic Range"
            absdiff_fixed_title = \
                f"Difference ({cmpres})\nDev - Ref, Restricted Range [5%,95%]"
            if diff_of_diffs:
                fracdiff_dynam_title = \
                    f"Difference ({cmpres}), " + \
                    f"Dynamic Range\n{frac_devstr} - {frac_refstr}"
                fracdiff_fixed_title = \
                    f"Difference ({cmpres}), " + \
                    f"Restricted Range [5%,95%]\n{frac_devstr} - {frac_refstr}"
            else:
                fracdiff_dynam_title = \
                    f"Ratio ({cmpres})\nDev/Ref, Dynamic Range"
                fracdiff_fixed_title = \
                    f"Ratio ({cmpres})\nDev/Ref, Fixed Range"
        else:
            absdiff_dynam_title = "Difference\nDev - Ref, Dynamic Range"
            absdiff_fixed_title = \
                "Difference\nDev - Ref, Restricted Range [5%,95%]"
            if diff_of_diffs:
                fracdiff_dynam_title = \
                    f"Difference, Dynamic Range\n{frac_devstr} - {frac_refstr}"
                fracdiff_fixed_title = \
                    "Difference, Restricted Range " + \
                    f"[5%,95%]\n{frac_devstr} - {frac_refstr}"
            else:
                fracdiff_dynam_title = "Ratio \nDev/Ref, Dynamic Range"
                fracdiff_fixed_title = "Ratio \nDev/Ref, Fixed Range"

        # ==============================================================
        # Bundle variables for 6 parallel plotting calls
        # 0 = Ref                 1 = Dev
        # 2 = Dynamic abs diff    3 = Restricted abs diff
        # 4 = Dynamic frac diff   5 = Restricted frac diff
        # ==============================================================

        subplots = six_panel_subplot_names(diff_of_diffs)

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
        if -1000 not in extent:
            extents = [extent[:], extent[:],
                       extent[:], extent[:],
                       extent[:], extent[:]]
        else:
            plot_extent = [np.max([cmp_extent[0], -180]),
                           np.min([cmp_extent[1], 180]),
                           cmp_extent[2], cmp_extent[3]]
            extents = [plot_extent[:], plot_extent[:],
                       plot_extent[:], plot_extent[:],
                       plot_extent[:], plot_extent[:]]
        plot_vals = [ds_ref, ds_dev, absdiff, absdiff, fracdiff, fracdiff]
        grids = [refgrid, devgrid, regional_cmp_grid.copy(), regional_cmp_grid.copy(),
                 regional_cmp_grid.copy(), regional_cmp_grid.copy()]
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
            cmaps = [ref_cmap, dev_cmap, cmap_gray,
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

        other_all_nans = [dev_is_all_nan, ref_is_all_nan,
                          False, False, False, False]

        mins = [vmin_ref, vmin_dev, vmin_abs]
        maxs = [vmax_ref, vmax_dev, vmax_abs]

        ratio_logs = [False, False, False, False, True, True]

        # Plot
        for i in range(6):
            six_plot(
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
                ratio_log=ratio_logs[i],
                proj=proj,
                ll_plot_func=ll_plot_func,
                **extra_plot_args
            )


        # ==============================================================
        # Add this page of 6-panel plots to a PDF file
        # ==============================================================
        if savepdf:
            folders = pdfname.split('/')
            pdfname_temp = folders[-1] + "BENCHMARKFIGCREATION.pdf" + str(ivar)
            full_path = temp_dir
            for folder in folders[:-1]:
                full_path = os.path.join(full_path, folder)
                if not os.path.isdir(full_path):
                    try:
                        os.mkdir(full_path)
                    except FileExistsError:
                        pass
            pdf = PdfPages(os.path.join(full_path, pdfname_temp))
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
        return ""

    # ==================================================================
    # Call figure generation function in a parallel loop over variables
    # ==================================================================
    # do not attempt nested thread parallelization due to issues with
    # matplotlib
    if current_process().name != "MainProcess":
        n_job = 1

    if not savepdf:
        # disable parallel plotting to allow interactive figure plotting
        for i in range(n_var):
            createfig(i)

    else:
        with TemporaryDirectory() as temp_dir:
            # ---------------------------------------
            # Turn off parallelization if n_job=1
            if n_job != 1:
                results = Parallel(n_jobs=n_job)(
                    delayed(createfig)(i, temp_dir)
                    for i in range(n_var)
                )
            else:
                results = []
                for i in range(n_var):
                    results.append(createfig(i, temp_dir))
            # ---------------------------------------

            # update sig diffs after parallel calls
            if current_process().name == "MainProcess":
                for varname in results:
                    if isinstance(varname, str):
                        sigdiff_list.append(varname)

            # ==========================================================
            # Finish
            # ==========================================================
            if verbose:
                print("Closed PDF")
            merge = PdfMerger()
            #print(f"Creating {pdfname} for {n_var} variables")
            pdf = PdfPages(pdfname)
            pdf.close()
            for i in range(n_var):
                temp_pdfname = pdfname
                if pdfname[0] == '/':
                    temp_pdfname = temp_pdfname[1:]
                merge.append(
                    os.path.join(
                        str(temp_dir),
                        temp_pdfname +
                        "BENCHMARKFIGCREATION.pdf" +
                        str(i)))
            merge.write(pdfname)
            merge.close()
            warnings.showwarning = _warning_format
