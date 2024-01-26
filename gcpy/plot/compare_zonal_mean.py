"""
Creates a six-panel comparison plot of zonal means from two different
GEOS-Chem model versions.  Called from the GEOS-Chem benchmarking scripts
and from the compare_diags.py example script.
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
from joblib import Parallel, delayed
from pypdf import PdfMerger
from gcpy.grid import get_vert_grid, get_pressure_indices, \
    pad_pressure_edges, convert_lev_to_pres
from gcpy.regrid import regrid_comparison_data, create_regridders, gen_xmat, \
    regrid_vertical
from gcpy.util import reshape_MAPL_CS, get_diff_of_diffs, \
    all_zero_or_nan, compare_varnames, \
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
        pres_range=None,
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
        sigdiff_list=None,
        second_ref=None,
        second_dev=None,
        spcdb_dir=os.path.dirname(__file__),
        sg_ref_path='',
        sg_dev_path='',
        ref_vert_params=None,
        dev_vert_params=None,
        **extra_plot_args
):
    """
    Creates 3x2 comparison zonal-mean plots for variables
    common in two xarray Datasets. Optionally save to PDF.

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
            Default value: None (will compare all common 3D variables)
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
            of Ref and Dev)
        match_cbar: bool
            Set this flag to True to use same the colorbar bounds
            for both Ref and Dev plots.
            Default value: True
        pres_range: list of two integers
            Pressure range of levels to plot [hPa]. The vertical axis
            will span the outer pressure edges of levels that contain
            pres_range endpoints.
            Default value: [0, 2000]
        normalize_by_area: bool
            Set this flag to True to to normalize raw data in both
            Ref and Dev datasets by grid area. Input ref and dev
            datasets must include AREA variable in m2 if normalizing
            by area.
            Default value: False
        enforce_units: bool
            Set this flag to True force an error if the variables in
            the Ref and Dev datasets have different units.
            Default value: True
        convert_to_ugm3: str
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
            Set this flag to True to use a blue-white-red colormap for
            plotting raw reference and development datasets.
            Default value: False
        verbose: logical
            Set this flag to True to enable informative printout.
            Default value: False
        log_color_scale: bool
            Set this flag to True to enable plotting data (not diffs)
            on a log color scale.
            Default value: False
        log_yaxis: bool
            Set this flag to True if you wish to create zonal mean
            plots with a log-pressure Y-axis.
            Default value: False
        extra_title_txt: str
            Specifies extra text (e.g. a date string such as "Jan2016")
            for the top-of-plot title.
            Default value: None
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
        ref_vert_params: list(AP, BP) of list-like types
            Hybrid grid parameter A in hPa and B (unitless).
            Needed if ref grid is not 47 or 72 levels.
            Default value: None
        dev_vert_params: list(AP, BP) of list-like types
            Hybrid grid parameter A in hPa and B (unitless).
            Needed if dev grid is not 47 or 72 levels.
            Default value: None
        extra_plot_args: various
            Any extra keyword arguments are passed through the
            plotting functions to be used in calls to pcolormesh()
            (CS) or imshow() (Lat/Lon).
    """
    warnings.showwarning = _warning_format
    verify_variable_type(refdata, xr.Dataset)
    verify_variable_type(devdata, xr.Dataset)

    # Create empty lists for keyword arguments
    if sigdiff_list is None:
        sigdiff_list = []
    if ref_vert_params is None:
        ref_vert_params = [[], []]
    if dev_vert_params is None:
        dev_vert_params = [[], []]
    if pres_range is None:
        pres_range = [0, 2000]

    # Determine if doing diff-of-diffs
    diff_of_diffs = second_ref is not None and second_dev is not None

    # Prepare diff-of-diffs datasets if needed
    if diff_of_diffs:
        refdata, devdata = refdata.load(), devdata.load()
        second_ref, second_dev = second_ref.load(), second_dev.load()

#        # If needed, use fake time dim in case dates are different in datasets.
#        # This needs more work for case of single versus multiple times.
#        aligned_time = np.datetime64('2000-01-01')
#        refdata = refdata.assign_coords({'time' : [aligned_time]})
#        devdata = devdata.assign_coords({'time' : [aligned_time]})
#        second_ref = second_ref.assign_coords({'time' : [aligned_time]})
#        second_dev = second_dev.assign_coords({'time' : [aligned_time]})

        refdata, fracrefdata = get_diff_of_diffs(refdata, second_ref)
        devdata, fracdevdata = get_diff_of_diffs(devdata, second_dev)

        frac_refstr = 'GCC_dev / GCC_ref'
        frac_devstr = 'GCHP_dev / GCHP_ref'

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
    # If converting to ug/m3, load the species database
    if convert_to_ugm3:
        properties = read_config_file(
            os.path.join(
                spcdb_dir,
                "species_database.yml"
            ),
            quiet=True
        )

    # Get mid-point pressure and edge pressures for this grid
    ref_pedge, ref_pmid, _ = get_vert_grid(refdata, *ref_vert_params)
    dev_pedge, dev_pmid, _ = get_vert_grid(devdata, *dev_vert_params)

    # Get indexes of pressure subrange (full range is default)
    ref_pedge_ind = get_pressure_indices(ref_pedge, pres_range)
    dev_pedge_ind = get_pressure_indices(dev_pedge, pres_range)

    # Pad edges if subset does not include surface or TOA so data spans
    # entire subrange
    ref_pedge_ind = pad_pressure_edges(
        ref_pedge_ind,
        refdata.sizes["lev"],
        np.size(ref_pmid))
    dev_pedge_ind = pad_pressure_edges(
        dev_pedge_ind,
        devdata.sizes["lev"],
        np.size(dev_pmid))

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
    # this may require checking for 48 / 73 levels
    ref_pmid_ind_flipped = refdata.sizes["lev"] - ref_pmid_ind[::-1] - 1
    dev_pmid_ind_flipped = devdata.sizes["lev"] - dev_pmid_ind[::-1] - 1
    if flip_ref:
        ref_pmid_ind = ref_pmid_ind_flipped
    if flip_dev:
        dev_pmid_ind = dev_pmid_ind_flipped

    refdata = refdata.isel(lev=ref_pmid_ind)
    devdata = devdata.isel(lev=dev_pmid_ind)
    if diff_of_diffs:
        fracrefdata = fracrefdata.isel(lev=ref_pmid_ind)
        fracdevdata = fracdevdata.isel(lev=dev_pmid_ind)

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

    [refres, refgridtype, devres, devgridtype, cmpres, cmpgridtype,
     regridref, regriddev, regridany, refgrid, devgrid, cmpgrid,
     refregridder, devregridder, refregridder_list, devregridder_list] = \
        create_regridders(
        refdata,
        devdata,
        weightsdir=weightsdir,
        cmpres=cmpres,
        zm=True,
        sg_ref_params=sg_ref_params,
        sg_dev_params=sg_dev_params
    )

    # Use smaller vertical grid as target for vertical regridding
    # NOTE: Convert target_index from numpy.int64 to int to conform
    # to the Python style guide (as per Pylint).
    #  -- Bob Yantosca (21 Sep 2023)
    target_index = int(np.array([len(ref_pedge), len(dev_pedge)]).argmin())
    pedge = [ref_pedge, dev_pedge][target_index]
    pedge_ind = [ref_pedge_ind, dev_pedge_ind][target_index]

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

        # If units string is UNITLESS (may be true for bpch data) then rename
        # units
        if ds_refs[i].units.strip() == "UNITLESS":
            ds_refs[i].attrs["units"] = "1"
        if ds_devs[i].units.strip() == "UNITLESS":
            ds_devs[i].attrs["units"] = "1"

        # Compare units of ref and dev. The check_units function will throw an error
        # if the units do not match and enforce_units is True.
        check_units(ds_refs[i], ds_devs[i], enforce_units)

        # Convert from ppb to ug/m3 if convert_to_ugm3 is passed as true
        if convert_to_ugm3:

            # Error checks: must pass met, not normalize by area, and be in ppb
            if refmet is None or devmet is None:
                msg = "Met mata ust be passed to convert units to ug/m3."
                raise ValueError(msg)
            if normalize_by_area:
                msg = "Normalizing by area is now allowed if plotting ug/m3"
                raise ValueError(msg)
            if ds_refs[i].units != "ppb" or ds_devs[i].units != "ppb":
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
                # Get the species molecular weight in g/mol
                spc_mw_g = species_properties.get("MW_g")
                if spc_mw_g is None:
                    msg = f"Molecular weight not found for for species {spc_name}!" \
                          + " Cannot convert to ug/m3."
                    raise ValueError(msg)

            # Convert values from ppb to ug/m3:
            # ug/m3 = 1e-9ppb * mol/g air * kg/m3 air * 1e3g/kg
            #         * g/mol spc * 1e6ug/g
            #       = ppb * air density * (spc MW / air MW)
            ds_refs[i].values = ds_refs[i].values * ref_airden.values \
                * (spc_mw_g / MW_AIR_g)
            ds_devs[i].values = ds_devs[i].values * dev_airden.values \
                * (spc_mw_g / MW_AIR_g)

            # Update units string
            ds_refs[i].attrs["units"] = "\u03BCg/m3"  # ug/m3 using mu
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
    # store units in case data changes from DataArray to numpy array
    ref_units = [None] * n_var
    dev_units = [None] * n_var

    # regrid vertically if necessary
    if len(ref_pedge) != len(pedge):
        xmat = gen_xmat(ref_pedge[ref_pedge_ind], pedge[pedge_ind])
    elif len(dev_pedge) != len(pedge):
        xmat = gen_xmat(dev_pedge[dev_pedge_ind], pedge[pedge_ind])

    for i in range(n_var):

        ds_ref = ds_refs[i]
        ds_dev = ds_devs[i]
        frac_ds_ref = frac_ds_refs[i]
        frac_ds_dev = frac_ds_devs[i]
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

        # save units for later use
        ref_units[i] = ds_ref.attrs["units"]
        dev_units[i] = ds_dev.attrs["units"]

        ref_nlev = len(ds_ref['lev'])
        dev_nlev = len(ds_dev['lev'])

        # Regrid variables horizontally
        # Ref
        ds_ref = regrid_comparison_data(
            ds_ref,
            refres,
            regridref,
            refregridder,
            refregridder_list,
            cmpgrid,
            refgridtype,
            cmpgridtype,
            nlev=ref_nlev
        )
        if diff_of_diffs:
            frac_ds_ref = regrid_comparison_data(
                frac_ds_ref,
                refres,
                regridref,
                refregridder,
                refregridder_list,
                cmpgrid,
                cmpgridtype,
                refgridtype,
                nlev=ref_nlev
            )
        # Dev
        ds_dev = regrid_comparison_data(
            ds_dev,
            devres,
            regriddev,
            devregridder,
            devregridder_list,
            cmpgrid,
            devgridtype,
            cmpgridtype,
            nlev=dev_nlev
        )
        if diff_of_diffs:
            frac_ds_dev = regrid_comparison_data(
                frac_ds_dev,
                devres,
                regriddev,
                devregridder,
                devregridder_list,
                cmpgrid,
                devgridtype,
                cmpgridtype,
                nlev=dev_nlev
            )

        # store regridded CS data before dealing with vertical regridding
        if refgridtype == "cs":
            ds_refs[i] = ds_ref
            frac_ds_refs[i] = frac_ds_ref
        if devgridtype == "cs":
            ds_devs[i] = ds_dev
            frac_ds_devs[i] = frac_ds_dev

        # Reduce variables to smaller vert grid if necessary for comparison
        if len(ref_pedge) != len(pedge):
            ds_ref = regrid_vertical(ds_ref, xmat, dev_pmid[dev_pmid_ind])
            if diff_of_diffs:
                frac_ds_ref = regrid_vertical(frac_ds_ref, xmat, dev_pmid[dev_pmid_ind])

        if len(dev_pedge) != len(pedge):
            ds_dev = regrid_vertical(ds_dev, xmat, ref_pmid[ref_pmid_ind])
            if diff_of_diffs:
                frac_ds_dev = regrid_vertical(frac_ds_dev, xmat, ref_pmid[ref_pmid_ind])
        ds_ref_cmps[i] = ds_ref
        ds_dev_cmps[i] = ds_dev
        if diff_of_diffs:
            frac_ds_ref_cmps[i] = frac_ds_ref
            frac_ds_dev_cmps[i] = frac_ds_dev
    # Universal plot setup
    xtick_positions = np.arange(-90, 91, 30)
    xticklabels = [rf"{x}$\degree$" for x in xtick_positions]

    # ==================================================================
    # Define function to create a single page figure to be called
    # in a parallel loop
    # ==================================================================
    def createfig(ivar, temp_dir=''):

        # Suppress harmless run-time warnings (mostly about underflow)
        warnings.filterwarnings('ignore', category=RuntimeWarning)
        warnings.filterwarnings('ignore', category=UserWarning)

        if savepdf and verbose:
            print(f"{ivar} ", end="")
        varname = varlist[ivar]

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
        cmn_units = ref_units[ivar]
        subtitle_extra = ""
        if normalize_by_area:
            exclude_list = ["WetLossConvFrac", "Prod_", "Loss_"]
            if not any(s in varname for s in exclude_list):
                if "/" in cmn_units:
                    cmn_units = f"{cmn_units}/m2"
                else:
                    cmn_units = f"{cmn_units} m-2"
                ref_units[ivar] = cmn_units
                dev_units[ivar] = cmn_units
                subtitle_extra = ", Normalized by Area"

        # ==============================================================
        # Calculate zonal mean
        # ==============================================================
        # Ref
        if refgridtype == "ll":
            zm_ref = ds_ref.mean(dim="lon")
        else:
            zm_ref = ds_ref.mean(axis=2)

        # Dev
        if devgridtype == "ll":
            zm_dev = ds_dev.mean(dim="lon")
        else:
            zm_dev = ds_dev.mean(axis=2)
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

        # ==============================================================
        # Test if Ref and/or Dev contain all zeroes or all NaNs.
        # This will have implications as to how we set min and max
        # values for the color ranges below.
        # ==============================================================
        ref_values = ds_ref.values if isinstance(ds_ref, xr.DataArray) else ds_ref
        dev_values = ds_dev.values if isinstance(ds_dev, xr.DataArray) else ds_dev
        ref_is_all_zero, ref_is_all_nan = all_zero_or_nan(ref_values)
        dev_is_all_zero, dev_is_all_nan = all_zero_or_nan(dev_values)

        # ==============================================================
        # Calculate zonal mean difference
        # ==============================================================

        zm_diff = np.array(zm_dev_cmp) - np.array(zm_ref_cmp)

        # Test if abs. diff is zero everywhere or NaN everywhere
        absdiff_is_all_zero, absdiff_is_all_nan = all_zero_or_nan(zm_diff)

        # ==============================================================
        # Calculate fractional difference, set divides by zero to Nan
        # ==============================================================
        if diff_of_diffs:
            zm_fracdiff = np.array(frac_zm_dev_cmp) -       \
                np.array(frac_zm_ref_cmp)
        else:
            zm_fracdiff = np.abs(np.array(zm_dev_cmp)) /    \
                np.abs(np.array(zm_ref_cmp))
        zm_fracdiff = np.where(np.abs(zm_fracdiff) ==
                               np.inf, np.nan, zm_fracdiff)
        zm_fracdiff[zm_fracdiff > 1e308] = np.nan
        # Test if the frac. diff is zero everywhere or NaN everywhere
        fracdiff_is_all_zero = not np.any(zm_fracdiff) or       \
            (np.nanmin(zm_fracdiff) == 0 and
             np.nanmax(zm_fracdiff) == 0)
        fracdiff_is_all_nan = np.isnan(zm_fracdiff).all()

        # ==============================================================
        # Create 3x2 figure
        # ==============================================================

        # Create figs and axes objects
        figs, ((ax0, ax1), (ax2, ax3), (ax4, ax5)) = plt.subplots(
            3, 2, figsize=[12, 15.3]
        )
        # Add extra adding so that plots don't bump into each other.
        # For zonal mean plots, we need to leave extra padding at the
        # left (for the Y-axis label) and at the bottom (for the colrobar).
        plt.subplots_adjust(
            left=0.10,    # Fraction of page width, from left edge
            right=0.925,  # Fraction of page width, from left edge
            bottom=0.05,  # Fraction of page height, from bottom edge
            wspace=0.25,  # Horizontal spacing btw subplots (frac of width)
            hspace=0.35   # Vertical spacing btw subplots (fract of height)
        )
        # Give the plot a title
        offset = 0.96
        if extra_title_txt is not None:
            figs.suptitle(
                f"{varname}, Zonal Mean ({extra_title_txt})",
                y=offset,
            )
        else:
            figs.suptitle(
                f"{varname}, Zonal Mean",
                y=offset
            )

        # ==============================================================
        # Set color map objects.  Use gray for NaNs (no worries,
        # because zonal means are always plotted on lat-alt grids).
        #
        # Use shallow copy (copy.copy() to create color map objects,
        # in order to avoid set_bad() from being applied to the base
        # color table. See: https://docs.python.org/3/library/copy.html
        # ==============================================================

        if use_cmap_RdBu:
            cmap1 = copy.copy(mpl.colormaps["RdBu_r"])
        else:
            cmap1 = copy.copy(WhGrYlRd)
        cmap1.set_bad("gray")

        cmap_plot = copy.copy(mpl.colormaps["RdBu_r"])
        cmap_plot.set_bad(color="gray")

        # ==============================================================
        # Set titles for plots
        # ==============================================================

        if refgridtype == "ll":
            ref_title = f"{refstr} (Ref){subtitle_extra}\n{refres}"
        else:
            ref_title = f"{refstr} (Ref){subtitle_extra}\n{cmpres} regridded from c{refres}"

        if devgridtype == "ll":
            dev_title = f"{devstr} (Dev){subtitle_extra}\n{devres}"
        else:
            dev_title = f"{devstr} (Dev){subtitle_extra}\n{cmpres} regridded from c{devres}"

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

        unit_list = [ref_units[ivar], dev_units[ivar], cmn_units, cmn_units,
                     "unitless", "unitless"]

        other_all_nans = [dev_is_all_nan, ref_is_all_nan,
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

        pedge_inds = [ref_pedge_ind, dev_pedge_ind, pedge_ind,
                      pedge_ind, pedge_ind, pedge_ind]

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
                pedges[i],
                pedge_inds[i],
                log_yaxis,
                plot_type="zonal_mean",
                xtick_positions=xtick_positions,
                xticklabels=xticklabels,
                ratio_log=ratio_logs[i],
                **extra_plot_args
            )

        # ==============================================================
        # Add this page of 6-panel plots to the PDF file
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
        if np.abs(1 - np.nanmax(zm_fracdiff)) > 0.1:
            sigdiff_list.append(varname)
            return varname
        return ""

    # ==================================================================
    # Call figure generation function in a parallel loop over variables
    #
    # ==================================================================

    # Disable parallelization if this routine is already being
    # called in parallel.  This is due to issues with matplotlib.
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
            #print("Creating {} for {} variables".format(pdfname, n_var))
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
