#!/usr/bin/env python3
"""
Specific utilities for creating plots from GEOS-Chem benchmark simulations.
"""
import os
import warnings
import itertools
import gc
import numpy as np
import pandas as pd
import xarray as xr
import sparselt.esmf
import sparselt.xr
from joblib import Parallel, delayed
from tabulate import tabulate
from gcpy.regrid import create_regridders
from gcpy.grid import get_troposphere_mask
from gcpy.util import \
    add_bookmarks_to_pdf, add_missing_variables, add_nested_bookmarks_to_pdf, \
    array_equals, compare_varnames, create_blank_dataarray, dataset_reader, \
    dataset_mean, divide_dataset_by_dataarray, get_area_from_dataset, \
    get_emissions_varnames, get_filepath, get_variables_from_dataset, \
    insert_text_into_file, make_directory, print_totals, read_config_file, \
    read_species_metadata, rename_and_flip_gchp_rst_vars, \
    replace_whitespace, unique_values, verify_variable_type, wrap_text
from gcpy.units import convert_units
from gcpy.constants import \
    COL_WIDTH, ENCODING, MW_AIR_g, SKIP_THESE_VARS, TABLE_WIDTH
from gcpy.plot.compare_single_level import compare_single_level
from gcpy.plot.compare_zonal_mean import compare_zonal_mean
from gcpy.benchmark.modules.benchmark_utils import \
    AOD_SPC, EMISSION_SPC, EMISSION_INV, add_lumped_species_to_dataset, \
    archive_lumped_species_definitions, get_species_categories, \
    archive_species_categories, rename_speciesconc_to_speciesconcvv

# Suppress numpy divide by zero warnings to prevent output spam
np.seterr(divide="ignore", invalid="ignore")


def create_global_mass_table(
        refdata,
        refstr,
        devdata,
        devstr,
        varlist,
        met_and_masks,
        spcdb_files,
        ref_hdr_label="",
        dev_hdr_label="",
        trop_only=False,
        outfilename="GlobalMass_TropStrat.txt",
        verbose=False,
):
    """
    Creates a table of global masses for a list of species in contained in
    two data sets.  The data sets,  which typically represent output from two
    different model versions, are usually contained in netCDF data files.

    Args:
        refdata: xarray Dataset
            The first data set to be compared (aka "Reference").
        refstr: str
            A string that can be used to identify refdata
            (e.g. a model version number or other identifier).
        devdata: xarray Dataset
            The second data set to be compared (aka "Development").
        devstr: str
            A string that can be used to identify the data set specified
            by devfile (e.g. a model version number or other identifier).
        varlist: list of strings
            List of species concentation variable names to include
            in the list of global totals.
        met_and_masks: dict of xarray DataArray
            Dictionary containing the meterological variables and
            masks for the Ref and Dev datasets.
        spcdb_files : list
            Paths to species_database.yml files in Ref & Dev rundirs

    Keyword Args (optional):
        ref_hdr_label : str
            Label for Ref, placed after refstr in the file header
            Default value: ""
        dev_hdr_label : str
            Label for Dev, placed after devstr in the file header
            Default value: ""
        trop_only: bool
            Set this switch to True if you wish to print totals
            only for the troposphere.
            Default value: False (i.e. print whole-atmosphere totals).
        outfilename: str
            Name of the text file which will contain the table of
            emissions totals.
            Default value: "GlobalMass_TropStrat.txt"
        verbose: bool
            Set this switch to True if you wish to print out extra
            informational messages.
            Default value: False

    Remarks:
        This method is mainly intended for model benchmarking purposes,
        rather than as a general-purpose tool.

        Species metadata (such as molecular weights) are read from a
        YAML file called "species_database.yml".
    """

    # ==================================================================
    # Initialization
    # ==================================================================
    verify_variable_type(refdata, xr.Dataset)
    verify_variable_type(devdata, xr.Dataset)

    # Make sure required arguments are passed
    if varlist is None:
        raise ValueError('The "varlist" argument was not passed!')
    if met_and_masks is None:
        raise ValueError('The "met_and_masks" argument was not passed!')

    # Read the species database files in the Ref & Dev rundirs, and
    # return a dict containing metadata for each.
    ref_metadata, dev_metadata = read_species_metadata(
        spcdb_files,
        quiet=True
    )

    # Replace whitespace in the ref and dev labels
    refstr = replace_whitespace(refstr)
    devstr = replace_whitespace(devstr)

    # ==================================================================
    # Open file for output
    # ==================================================================

    # Create file
    try:
        f = open(outfilename, "w", encoding=ENCODING)
    except (IOError, OSError, FileNotFoundError) as e:
        raise e(f"Could not open {outfilename} for writing!") from e

    # Define a list for differences
    diff_list = []

    # Title strings
    title1 = f"### Global mass (Gg) (Trop + Strat)"
    if trop_only:
        title1 = f"### Global mass (Gg) (Trop only)"
    title2 = f"### Ref = {refstr} {ref_hdr_label}"
    title3 = f"### Dev = {devstr} {dev_hdr_label}"

    # Write a placeholder to the file that denotes where
    # the list of species with differences will be written
    placeholder = "@%% insert diff status here %%@"

    # Print header to file
    print("#" * TABLE_WIDTH, file=f)
    print(f"{title1 : <{TABLE_WIDTH-3}}{'###'}", file=f)
    print(f"{'###'  : <{TABLE_WIDTH-3}}{'###'}", file=f)
    print(f"{title2 : <{TABLE_WIDTH-3}}{'###'}", file=f)
    print(f"{title3 : <{TABLE_WIDTH-3}}{'###'}", file=f)
    print(f"{'###'  : <{TABLE_WIDTH-3}}{'###'}", file=f)
    print(f"{placeholder}", file=f)
    print("#" * TABLE_WIDTH, file=f)

    # Column headers
    print(f"{'' : <{COL_WIDTH-1}}{'Ref' : >{COL_WIDTH}}{'Dev' : >{COL_WIDTH}}{'Dev - Ref' : >{COL_WIDTH}}{'% diff' : >{COL_WIDTH}} {'diffs'}", file=f)

    # ==================================================================
    # Print global masses for all species
    #
    # NOTE: By this point, all secies will be in both Ref and Dev'
    # because we have added them in the calling routine
    # ==================================================================
    for v in varlist:

        # Get the species name
        spc_name = v.split("_")[1]

        # Get metadta for the given species
        ref_species_metadata = ref_metadata.get(spc_name)
        dev_species_metadata = dev_metadata.get(spc_name)
        if ref_species_metadata is None and dev_species_metadata is None:
            if verbose:
                msg = f"No metadata found for {spc_name} ... skippping"
                print(msg)
            continue

        # Specify target units
        target_units = "Gg"

        # ==============================================================
        # Convert units of Ref and save to a DataArray
        # (or skip if Ref contains NaNs everywhere)
        # ==============================================================
        refarray = refdata[v]
        if not np.isnan(refdata[v].values).all():
            refarray = convert_units(
                refarray,
                spc_name,
                ref_species_metadata,
                target_units,
                area_m2=met_and_masks["Ref_Area"],
                delta_p=met_and_masks["Ref_Delta_P"],
                box_height=met_and_masks["Ref_BxHeight"],
            )

        # ==============================================================
        # Convert units of Dev and save to a DataArray
        # (or skip if Dev contains NaNs everywhere)
        # ==============================================================
        devarray = devdata[v]
        if not np.isnan(devdata[v].values).all():
            devarray = convert_units(
                devarray,
                spc_name,
                dev_species_metadata,
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
                devarray,
                f,
                diff_list,
                masks=met_and_masks
            )
        else:
            print_totals(
                refarray,
                devarray,
                f,
                diff_list
            )

    # ==================================================================
    # Cleanup and quit
    # ==================================================================

    # Close file
    f.close()

    # Reopen file and replace placeholder text by diff_text
    insert_text_into_file(
        filename=outfilename,
        search_text=placeholder,
        replace_text=diff_list_to_text(
            refstr,
            devstr,
            diff_list,
            fancy_format=True
        ),
        width=TABLE_WIDTH
    )


def create_mass_accumulation_table(
        refdatastart,
        refdataend,
        refstr,
        refperiodstr,
        devdatastart,
        devdataend,
        devstr,
        devperiodstr,
        varlist,
        met_and_masks,
        label,
        spcdb_files,
        trop_only=False,
        outfilename="GlobalMassAccum_TropStrat.txt",
        verbose=False,
):
    """
    Creates a table of global mass accumulation for a list of species in
    two data sets.  The data sets, which typically represent output from two
    different model versions, are usually contained in netCDF data files.

    Args:
        refdatastart: xarray Dataset
            The first data set to be compared (aka "Reference").
        refdataend: xarray Dataset
            The first data set to be compared (aka "Reference").
        refstr: str
            A string that can be used to identify refdata
            (e.g. a model version number or other identifier).
        refperiodstr: str
            Ref simulation period start and end
        devdatastart: xarray Dataset
            The second data set to be compared (aka "Development").
        devdataend: xarray Dataset
            The second data set to be compared (aka "Development").
        devstr: str
            A string that can be used to identify the data set specified
            by devfile (e.g. a model version number or other identifier).
        devperiodstr: str
            Ref simulation period start and end
        varlist: list of strings
            List of species concentation variable names to include
            in the list of global totals.
        met_and_masks: dict of xarray DataArray
            Dictionary containing the meterological variables and
            masks for the Ref and Dev datasets.
        label: str
            Label to go in the header string.  Can be used to
            pass the month & year.
        spcdb_files : list
            Paths to species_database.yml files in Ref & Dev rundirs

    Keyword Args (optional):
        trop_only: bool
            Set this switch to True if you wish to print totals
            only for the troposphere.
            Default value: False (i.e. print whole-atmosphere totals).
        outfilename: str
            Name of the text file which will contain the table of
            emissions totals.
            Default value: "GlobalMass_TropStrat.txt"
        verbose: bool
            Set this switch to True if you wish to print out extra
            informational messages.
            Default value: False

    Remarks:
        This method is mainly intended for model benchmarking purposes,
        rather than as a general-purpose tool.

        Species metadata (such as molecular weights) are read from a
        YAML file called "species_database.yml".
    """

    # ==================================================================
    # Initialization
    # ==================================================================
    verify_variable_type(refdatastart, xr.Dataset)
    verify_variable_type(refdataend, xr.Dataset)
    verify_variable_type(devdatastart, xr.Dataset)
    verify_variable_type(devdataend, xr.Dataset)

    # Make sure required arguments are passed
    if varlist is None:
        raise ValueError('The "varlist" argument was not passed!')
    if met_and_masks is None:
        raise ValueError('The "met_and_masks" argument was not passed!')

    # Read the species database files in the Ref & Dev rundirs, and
    # return a dict containing metadata for each.
    ref_metadata, dev_metadata = read_species_metadata(
        spcdb_files,
        quiet=True
    )

    # Replace whitespace in the ref and dev labels
    refstr = replace_whitespace(refstr)
    devstr = replace_whitespace(devstr)

    # ==================================================================
    # Open file for output
    # ==================================================================

    # Create file
    try:
        f = open(outfilename, "w", encoding=ENCODING)
    except (IOError, OSError, FileNotFoundError) as e:
        raise e(f"Could not open {outfilename} for writing!") from e

    # Define a list for differences
    diff_list = []

    # Title strings
    title1 = f"### Global mass accumulation (Gg) {label} (Trop + Strat)"
    if trop_only:
        title1 = f"### Global mass accumulation (Gg) {label} (Trop only)"
    title2 = "### Computed as change in instantaneous mass across period"
    title3 = f"### Ref = {refstr}"
    title4 = f"### Dev = {devstr}"
    title5 = f"### Ref period: {refperiodstr}"
    title6 = f"### Dev period: {devperiodstr}"

    # Write a placeholder to the file that denotes where
    # the list of species with differences will be written
    placeholder = "@%% insert diff status here %%@"

    # Print header to file
    print("#" * TABLE_WIDTH, file=f)
    print(f"{title1 : <{TABLE_WIDTH-3}}{'###'}", file=f)
    print(f"{'###'  : <{TABLE_WIDTH-3}}{'###'}", file=f)
    print(f"{title2 : <{TABLE_WIDTH-3}}{'###'}", file=f)
    print(f"{'###'  : <{TABLE_WIDTH-3}}{'###'}", file=f)
    print(f"{title3 : <{TABLE_WIDTH-3}}{'###'}", file=f)
    print(f"{title4 : <{TABLE_WIDTH-3}}{'###'}", file=f)
    print(f"{'###'  : <{TABLE_WIDTH-3}}{'###'}", file=f)
    print(f"{title5 : <{TABLE_WIDTH-3}}{'###'}", file=f)
    print(f"{title6 : <{TABLE_WIDTH-3}}{'###'}", file=f)
    print(f"{'###'  : <{TABLE_WIDTH-3}}{'###'}", file=f)
    print(f"{placeholder}", file=f)
    print("#" * TABLE_WIDTH, file=f)

    # Column headers
    print(f"{'' : <{COL_WIDTH-1}}{'Ref' : >{COL_WIDTH}}{'Dev' : >{COL_WIDTH}}{'Dev - Ref' : >{COL_WIDTH}}{'% diff' : >{COL_WIDTH}} {'diffs'}", file=f)

    # ==================================================================
    # Print global masses for all species
    #
    # NOTE: By this point, all secies will be in both Ref and Dev'
    # because we have added them in the calling routine
    # ==================================================================
    for v in varlist:

        # Get the species name
        spc_name = v.split("_")[1]

        # Get a list of metadata for the given species
        ref_species_metadata = ref_metadata.get(spc_name)
        dev_species_metadata = dev_metadata.get(spc_name)
        if ref_species_metadata is None and dev_species_metadata is None:
            if verbose:
                msg = f"No properties found for {spc_name} ... skippping"
                print(msg)
            continue

        # Specify target units
        target_units = "Gg"

        # ==============================================================
        # Convert units of Ref and save to a DataArray
        # (or skip if Ref contains NaNs everywhere)
        # ==============================================================
        refarrays = refdatastart[v]
        if not np.isnan(refdatastart[v].values).all():
            refarrays = convert_units(
                refarrays,
                spc_name,
                ref_species_metadata,
                target_units,
                area_m2=met_and_masks["Refs_Area"],
                delta_p=met_and_masks["Refs_Delta_P"],
                box_height=met_and_masks["Refs_BxHeight"],
            )

        refarraye = refdataend[v]
        if not np.isnan(refdataend[v].values).all():
            refarraye = convert_units(
                refarraye,
                spc_name,
                ref_species_metadata,
                target_units,
                area_m2=met_and_masks["Refe_Area"],
                delta_p=met_and_masks["Refe_Delta_P"],
                box_height=met_and_masks["Refe_BxHeight"],
            )

        refarray = refarrays
        refarray.values = refarraye.values - refarrays.values

        # ==============================================================
        # Convert units of Dev and save to a DataArray
        # (or skip if Dev contains NaNs everywhere)
        # ==============================================================
        devarrays = devdatastart[v]
        if not np.isnan(devdatastart[v].values).all():
            devarrays = convert_units(
                devarrays,
                spc_name,
                dev_species_metadata,
                target_units,
                area_m2=met_and_masks["Devs_Area"],
                delta_p=met_and_masks["Devs_Delta_P"],
                box_height=met_and_masks["Devs_BxHeight"],
            )
        #print('devarrays: {}'.format(devarrays.values))

        devarraye = devdataend[v]
        if not np.isnan(devdataend[v].values).all():
            devarraye = convert_units(
                devarraye,
                spc_name,
                dev_species_metadata,
                target_units,
                area_m2=met_and_masks["Deve_Area"],
                delta_p=met_and_masks["Deve_Delta_P"],
                box_height=met_and_masks["Deve_BxHeight"],
            )

        devarray = devarrays
        devarray.values = devarraye.values - devarrays.values

        # ==============================================================
        # Print global masses for Ref and Dev
        # (we will mask out tropospheric boxes in print_totals)
        # ==============================================================
        # ewl: for now trop_only is always false for accumulation table
        if trop_only:
            print_totals(
                refarray,
                devarray,
                f,
                diff_list,
                masks=met_and_masks
            )
        else:
            print_totals(
                refarray,
                devarray,
                f,
                diff_list
            )

    # ==================================================================
    # Cleanup and quit
    # ==================================================================

    # Close file
    f.close()

    # Reopen file and replace placeholder text by diff_text
    insert_text_into_file(
        filename=outfilename,
        search_text=placeholder,
        replace_text=diff_list_to_text(
            refstr,
            devstr,
            diff_list,
            fancy_format=True
        ),
        width=TABLE_WIDTH
    )


def make_benchmark_jvalue_plots(
        ref,
        refstr,
        dev,
        devstr,
        spcdb_files,
        varlist=None,
        dst="./benchmark",
        subdst=None,
        local_noon_jvalues=False,
        cmpres=None,
        plots=["sfc", "500hpa", "zonalmean"],
        overwrite=False,
        verbose=False,
        flip_ref=False,
        flip_dev=False,
        log_color_scale=False,
        sigdiff_files=None,
        weightsdir='.',
        n_job=-1,
        time_mean=False,
):
    """
    Creates PDF files containing plots of J-values for model
    benchmarking purposes.

    Args:
        ref: str
            Path name for the "Ref" (aka "Reference") data set.
        refstr: str
            A string to describe ref (e.g. version number)
        dev: str
            Path name for the "Dev" (aka "Development") data set.
            This data set will be compared against the "Reference"
            data set.
        devstr: str
            A string to describe dev (e.g. version number)
        spcdb_files : list
            Paths to species_database.yml files in Ref & Dev rundirs

    Keyword Args (optional):
        varlist: list of str
            List of J-value variables to plot.  If not passed,
            then all J-value variables common to both dev
            and ref will be plotted.  The varlist argument can be
            a useful way of restricting the number of variables
            plotted to the pdf file when debugging.
            Default value: None
        dst: str
            A string denoting the destination folder where a
            PDF file  containing plots will be written.
            Default value: ./benchmark.
        subdst: str
            A string denoting the sub-directory of dst where PDF
            files containing plots will be written.  In practice,
            subdst is only needed for the 1-year benchmark output,
            and denotes a date string (such as "Jan2016") that
            corresponds to the month that is being plotted.
            Default value: None
        local_noon_jvalues: bool
            Set this flag to plot local noon J-values.  This will
            divide all J-value variables by the JNoonFrac counter,
            which is the fraction of the time that it was local noon
            at each location.
            Default value: False
        cmpres: string
            Grid resolution at which to compare ref and dev data, e.g. '1x1.25'
        plots: list of strings
            List of plot types to create.
            Default value: ['sfc', '500hpa', 'zonalmean']
        overwrite: bool
            Set this flag to True to overwrite files in the
            destination folder (specified by the dst argument).
            Default value: False.
        verbose: bool
            Set this flag to True to print extra informational output.
            Default value: False
        flip_ref: bool
            Set this flag to True to reverse the vertical level
            ordering in the "Ref" dataset (in case "Ref" starts
            from the top of atmosphere instead of the surface).
            Default value: False
        flip_dev: bool
            Set this flag to True to reverse the vertical level
            ordering in the "Dev" dataset (in case "Dev" starts
            from the top of atmosphere instead of the surface).
            Default value: False
        log_color_scale: bool
            Set this flag to True if you wish to enable plotting data
            (not diffs) on a log color scale.
            Default value: False
        sigdiff_files: list of str
            Filenames that will contain the lists of J-values having
            significant differences in the 'sfc', '500hpa', and
            'zonalmean' plots.  These lists are needed in order to
            fill out the benchmark approval forms.
            Default value: None
        weightsdir: str
            Directory in which to place (and possibly reuse) xESMF regridder
            netCDF files.
            Default value: '.'
        n_job: int
            Defines the number of simultaneous workers for parallel plotting.
            Set to 1 to disable parallel plotting. Value of -1 allows the
            application to decide.
            Default value: -1
        time_mean : bool
            Determines if we should average the datasets over time
            Default value: False

    Remarks:
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

    # Create the directory for output
    make_directory(dst, overwrite)

    # Replace whitespace in the ref and dev labels
    refstr = replace_whitespace(refstr)
    devstr = replace_whitespace(devstr)

    # Get the function that will read file(s) into a Dataset
    reader = dataset_reader(time_mean, verbose=verbose)

    # Ref dataset
    try:
        refds = reader(ref, drop_variables=SKIP_THESE_VARS)
    except (OSError, IOError, FileNotFoundError) as e:
        raise e(f"Could not find Ref file: {ref}") from e

    # Dev dataset
    try:
        devds = reader(dev, drop_variables=SKIP_THESE_VARS)
    except (OSError, IOError, FileNotFoundError) as e:
        raise e(f"Could not find Ref file: {dev}") from e

    # Compute mean of data over the time dimension (if time_mean=True)
    if time_mean:
        refds = dataset_mean(refds)
        devds = dataset_mean(devds)

    # Make sure that Ref and Dev datasets have the same variables.
    # Variables that are in Ref but not in Dev will be added to Dev
    # with all missing values (NaNs). And vice-versa.
    [refds, devds] = add_missing_variables(refds, devds)

    # Create regridding files if necessary
    [_ for _ in create_regridders(refds, devds, weightsdir=weightsdir)]

    # Get a list of the 3D variables in both datasets
    if varlist is None:
        vardict = compare_varnames(refds, devds, quiet=not verbose)
        cmn = vardict["commonvars3D"]

    # ==================================================================
    # Local noon or continuously-averaged J-values?
    # ==================================================================
    if local_noon_jvalues:

        # Get a list of local noon J-value variables
        # (or use the varlist passed via tha argument list)
        prefix = "JNoon_"
        if varlist is None:
            varlist = [v for v in cmn if prefix in v]

        # Make sure JNoonFrac (fraction of times it was local noon
        # in each column) is present in both Ref and Dev datasets
        if "JNoonFrac" not in cmn:
            msg = "JNoonFrac is not common to Ref and Dev datasets!"
            raise ValueError(msg)

        # JNoon_* are cumulative sums of local noon J-values; we need
        # to divide these by JNoonFrac to get the average value
        refds = divide_dataset_by_dataarray(refds,
                                                 refds["JNoonFrac"], varlist)
        devds = divide_dataset_by_dataarray(devds,
                                                 devds["JNoonFrac"], varlist)

        # Subfolder of dst where PDF files will be printed
        catdir = "JValuesLocalNoon"

    else:

        # Get a list of continuously averaged J-value variables
        # (or use the varlist passed via tha argument list)
        prefix = "Jval"
        if varlist is None:
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
            pdfname = os.path.join(
                jvdir,
                f"{prefix}_Surface_{subdst}.pdf"
            )
        else:
            pdfname = os.path.join(
                jvdir,
                f"{prefix}_Surface.pdf"
            )

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
            weightsdir=weightsdir,
            n_job=n_job,
            spcdb_files=spcdb_files,
        )
        diff_sfc[:] = [v.replace(prefix, "") for v in diff_sfc]
        add_bookmarks_to_pdf(
            pdfname,
            varlist,
            remove_prefix=prefix,
            verbose=verbose)

    # 500hPa plots
    if "500hpa" in plots:
        if subdst is not None:
            pdfname = os.path.join(
                jvdir,
                f"{prefix}_500hPa_{subdst}.pdf"
            )
        else:
            pdfname = os.path.join(
                jvdir, f"{prefix}_500hPa.pdf"
            )

        diff_500 = []
        compare_single_level(
            refds,
            refstr,
            devds,
            devstr,
            varlist=varlist,
            cmpres=cmpres,
            ilev=22,
            pdfname=pdfname,
            flip_ref=flip_ref,
            flip_dev=flip_dev,
            log_color_scale=log_color_scale,
            extra_title_txt=extra_title_txt,
            sigdiff_list=diff_500,
            weightsdir=weightsdir,
            n_job=n_job,
            spcdb_files=spcdb_files,
        )
        diff_500[:] = [v.replace(prefix, "") for v in diff_500]
        add_bookmarks_to_pdf(
            pdfname,
            varlist,
            remove_prefix=prefix,
            verbose=verbose
        )

    # Full-column zonal mean plots
    if "zonalmean" in plots:
        if subdst is not None:
            pdfname = os.path.join(
                jvdir,
                f"{prefix}_FullColumn_ZonalMean_{subdst}.pdf"
            )
        else:
            pdfname = os.path.join(
                jvdir, f"{prefix}_FullColumn_ZonalMean.pdf"
            )

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
            weightsdir=weightsdir,
            n_job=n_job,
            spcdb_files=spcdb_files,
        )
        diff_zm[:] = [v.replace(prefix, "") for v in diff_zm]
        add_bookmarks_to_pdf(
            pdfname,
            varlist,
            remove_prefix=prefix,
            verbose=verbose)

        # Strat_ZonalMean plots will use a log-pressure Y-axis, with
        # a range of 1..100 hPa, as per GCSC request. (bmy, 8/13/19)
        if subdst is not None:
            pdfname = os.path.join(
                jvdir,
                f"{prefix}_Strat_ZonalMean_{subdst}.pdf"
            )
        else:
            pdfname = os.path.join(
                jvdir,
                f"{prefix}_Strat_ZonalMean.pdf"
            )

        compare_zonal_mean(
            refds,
            refstr,
            devds,
            devstr,
            varlist=varlist,
            pdfname=pdfname,
            pres_range=[1, 100],
            log_yaxis=True,
            flip_ref=flip_ref,
            flip_dev=flip_dev,
            extra_title_txt=extra_title_txt,
            log_color_scale=log_color_scale,
            weightsdir=weightsdir,
            n_job=n_job,
            spcdb_files=spcdb_files,
        )
        add_bookmarks_to_pdf(
            pdfname,
            varlist,
            remove_prefix=prefix,
            verbose=verbose
        )

        # ==============================================================
        # Write the lists of J-values that have significant differences,
        # which we need to fill out the benchmark approval forms.
        # ==============================================================
        if sigdiff_files is not None:
            for filename in sigdiff_files:
                if "sfc" in plots:
                    if "sfc" in filename:
                        with open(filename, "a+", encoding=ENCODING) as f:
                            print("* J-Values: ", file=f, end="")
                            for v in diff_sfc:
                                print(f"{v} ", file=f, end="")
                            print(file=f)
                            f.close()

                if "500" in plots:
                    if "500" in filename:
                        with open(filename, "a+", encoding=ENCODING) as f:
                            print("* J-Values: ", file=f, end="")
                            for v in diff_500:
                                print(f"{v} ", file=f, end="")
                            print(file=f)
                            f.close()

                if "zonalmean" in plots or "zm" in plots:
                    if "zonalmean" in filename or "zm" in filename:
                        with open(filename, "a+", encoding=ENCODING) as f:
                            print("* J-Values: ", file=f, end="")
                            for v in diff_zm:
                                print(f"{v} ", file=f, end="")
                            print(file=f)
                            f.close()

    # -------------------------------------------
    # Clean up
    # -------------------------------------------
    del refds
    del devds
    gc.collect()

def make_benchmark_collection_2d_var_plots(
        ref,
        refstr,
        dev,
        devstr,
        colname,
        spcdb_files,
        var_prefix=None,
        varlist=None,
        dst="./benchmark",
        cmpres=None,
        overwrite=False,
        verbose=False,
        flip_ref=False,
        flip_dev=False,
        log_color_scale=False,
        weightsdir='.',
        n_job=-1,
        time_mean=False,
):
    """
    Creates PDF file containing plots comparing all 2D variables in
    a collection without special handling.

    Args:
        ref: str
            Path name for the "Ref" (aka "Reference") data set.
        refstr: str
            A string to describe ref (e.g. version number)
        dev: str
            Path name for the "Dev" (aka "Development") data set.
            This data set will be compared against the "Reference"
            data set.
        devstr: str
            A string to describe dev (e.g. version number)
        colname: str
            Name of file collection, to be used in PDF name
        spcdb_files : list
            Paths to species_database.yml files in Ref & Dev rundirs

    Keyword Args (optional):
        var_prefix: str
            Variable name prefix in file to search for. If excluded then all 2D variables
            will be plotted regardless of name.
            Default value: None
        varlist: list of str
            List of variables to plot.  If not passed,
            then all variables common to both dev
            and ref will be plotted.  The varlist argument can be
            a useful way of restricting the number of variables
            plotted to the pdf file when debugging.
            Default value: None
        dst: str
            A string denoting the destination folder where a
            PDF file  containing plots will be written.
            Default value: ./benchmark.
        cmpres: string
            Grid resolution at which to compare ref and dev data, e.g. '1x1.25'
        overwrite: bool
            Set this flag to True to overwrite files in the
            destination folder (specified by the dst argument).
            Default value: False.
        verbose: bool
            Set this flag to True to print extra informational output.
            Default value: False
        flip_ref: bool
            Set this flag to True to reverse the vertical level
            ordering in the "Ref" dataset (in case "Ref" starts
            from the top of atmosphere instead of the surface).
            Default value: False
        flip_dev: bool
            Set this flag to True to reverse the vertical level
            ordering in the "Dev" dataset (in case "Dev" starts
            from the top of atmosphere instead of the surface).
            Default value: False
        log_color_scale: bool
            Set this flag to True if you wish to enable plotting data
            (not diffs) on a log color scale.
            Default value: False
        weightsdir: str
            Directory in which to place (and possibly reuse) xESMF regridder
            netCDF files.
            Default value: '.'
        n_job: int
            Defines the number of simultaneous workers for parallel plotting.
            Set to 1 to disable parallel plotting. Value of -1 allows the
            application to decide.
            Default value: -1
        time_mean : bool
            Determines if we should average the datasets over time
            Default value: False

    Remarks:
         Will create 1 file containing Budget plots.

    """

    # ==================================================================
    # Initialization
    # ==================================================================

    # Create the directory for output
    make_directory(dst, overwrite)

    # Replace whitespace in the ref and dev labels
    refstr = replace_whitespace(refstr)
    devstr = replace_whitespace(devstr)

    # Get the function that will read file(s) into a Dataset
    reader = dataset_reader(time_mean, verbose=verbose)

    # Ref dataset
    try:
        refds = reader(ref, drop_variables=SKIP_THESE_VARS)
    except (OSError, IOError, FileNotFoundError) as e:
        raise e(f"Could not find Ref file: {ref}") from e

    # Dev dataset
    try:
        devds = reader(dev, drop_variables=SKIP_THESE_VARS)
    except (OSError, IOError, FileNotFoundError) as e:
        raise e(f"Could not find Ref file: {dev}") from e

    # Compute mean of data over the time dimension (if time_mean=True)
    if time_mean:
        refds = dataset_mean(refds)
        devds = dataset_mean(devds)

    # Make sure that Ref and Dev datasets have the same variables.
    # Variables that are in Ref but not in Dev will be added to Dev
    # with all missing values (NaNs). And vice-versa.
    [refds, devds] = add_missing_variables(refds, devds)

    # Create regridding files if necessary
    [_ for _ in create_regridders(refds, devds, weightsdir=weightsdir)]

    # Get a list of the 2D variables in both datasets
    if varlist is None:
        vardict = compare_varnames(refds, devds, quiet=not verbose)
        cmn = vardict["commonvars2D"]

    # Get a list of variables
    # (or use the varlist passed via tha argument list)
    if varlist is None:
        if var_prefix is None:
            varlist = cmn
        else:
            varlist = [v for v in cmn if var_prefix in v]

    # ==================================================================
    # Create the plots
    # ==================================================================

    # Make the output folder if it doesn't exist
    bdir = os.path.join(dst, colname)
    if not os.path.isdir(bdir):
        os.mkdir(bdir)
    extra_title_txt = None
    pdfname = os.path.join(
        bdir,
        f"{colname}_2D.pdf"
    )

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
        weightsdir=weightsdir,
        n_job=n_job,
        spcdb_files=spcdb_files,
    )
    add_bookmarks_to_pdf(
        pdfname,
        varlist,
        verbose=verbose
    )

    # -------------------------------------------
    # Clean up
    # -------------------------------------------
    del refds
    del devds
    gc.collect()


def make_benchmark_collection_3d_var_plots(
        ref,
        refstr,
        dev,
        devstr,
        colname,
        spcdb_files,
        var_prefix=None,
        varlist=None,
        dst="./benchmark",
        cmpres=None,
        plots=["sfc", "500hpa", "zonalmean"],
        overwrite=False,
        verbose=False,
        flip_ref=False,
        flip_dev=False,
        log_color_scale=False,
        weightsdir='.',
        n_job=-1,
        time_mean=False,
):
    """
    Creates PDF files containing plots comparing all 3D variables in
    a collection without special handling.

    Args:
        ref: str
            Path name for the "Ref" (aka "Reference") data set.
        refstr: str
            A string to describe ref (e.g. version number)
        dev: str
            Path name for the "Dev" (aka "Development") data set.
            This data set will be compared against the "Reference"
            data set.
        devstr: str
            A string to describe dev (e.g. version number)
        colname: str
            Name of file collection, to be used in PDF name
        spcdb_files : list
            Paths to species_database.yml files in Ref & Dev rundirs

    Keyword Args (optional):
        var_prefix: str
            Variable name prefix in file to search for. If excluded then all 2D variables
            will be plotted regardless of name.
            Default value: None
        varlist: list of str
            List of J-value variables to plot.  If not passed,
            then all J-value variables common to both dev
            and ref will be plotted.  The varlist argument can be
            a useful way of restricting the number of variables
            plotted to the pdf file when debugging.
            Default value: None
        dst: str
            A string denoting the destination folder where a
            PDF file  containing plots will be written.
            Default value: ./benchmark.
        cmpres: string
            Grid resolution at which to compare ref and dev data, e.g. '1x1.25'
        plots: list of strings
            List of plot types to create.
            Default value: ['sfc', '500hpa', 'zonalmean']
        overwrite: bool
            Set this flag to True to overwrite files in the
            destination folder (specified by the dst argument).
            Default value: False.
        verbose: bool
            Set this flag to True to print extra informational output.
            Default value: False
        flip_ref: bool
            Set this flag to True to reverse the vertical level
            ordering in the "Ref" dataset (in case "Ref" starts
            from the top of atmosphere instead of the surface).
            Default value: False
        flip_dev: bool
            Set this flag to True to reverse the vertical level
            ordering in the "Dev" dataset (in case "Dev" starts
            from the top of atmosphere instead of the surface).
            Default value: False
        log_color_scale: bool
            Set this flag to True if you wish to enable plotting data
            (not diffs) on a log color scale.
            Default value: False
        weightsdir: str
            Directory in which to place (and possibly reuse) xESMF regridder
            netCDF files.
            Default value: '.'
        n_job: int
            Defines the number of simultaneous workers for parallel plotting.
            Set to 1 to disable parallel plotting. Value of -1 allows the
            application to decide.
            Default value: -1
        time_mean : bool
            Determines if we should average the datasets over time
            Default value: False

    Remarks:
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

    # Create the directory for output
    make_directory(dst, overwrite)

    # Replace whitespace in the ref and dev labels
    refstr = replace_whitespace(refstr)
    devstr = replace_whitespace(devstr)

    # Get the function that will read file(s) into a Dataset
    reader = dataset_reader(time_mean, verbose=verbose)

    # Ref dataset
    try:
        refds = reader(ref, drop_variables=SKIP_THESE_VARS)
    except (OSError, IOError, FileNotFoundError) as e:
        raise e(f"Could not find Ref file: {ref}") from e

    # Dev dataset
    try:
        devds = reader(dev, drop_variables=SKIP_THESE_VARS)
    except (OSError, IOError, FileNotFoundError) as e:
        raise e(f"Could not find Ref file: {dev}") from e

    # Compute mean of data over the time dimension (if time_mean=True)
    if time_mean:
        refds = dataset_mean(refds)
        devds = dataset_mean(devds)

    # Make sure that Ref and Dev datasets have the same variables.
    # Variables that are in Ref but not in Dev will be added to Dev
    # with all missing values (NaNs). And vice-versa.
    [refds, devds] = add_missing_variables(refds, devds)

    # Create regridding files if necessary
    [_ for _ in create_regridders(refds, devds, weightsdir=weightsdir)]

    # Get a list of the 3D variables in both datasets
    if varlist is None:
        vardict = compare_varnames(refds, devds, quiet=not verbose)
        cmn = vardict["commonvars3D"]

    # Get a list of variables
    # (or use the varlist passed via tha argument list)
    if varlist is None:
        if var_prefix is None:
            varlist = cmn
        else:
            varlist = [v for v in cmn if var_prefix in v]

    # ==================================================================
    # Create the plots
    # ==================================================================

    # Make the output folder if it doesn't exist
    bdir = os.path.join(dst, colname)
    if not os.path.isdir(bdir):
        os.mkdir(bdir)
    extra_title_txt = None
    pdfname = os.path.join(
        bdir,
        f"{colname}.pdf"
    )

    # Surface plots
    if "sfc" in plots:
        pdfname = os.path.join(
            bdir,
            f"{colname}_Surface.pdf"
        )

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
            weightsdir=weightsdir,
            n_job=n_job,
            spcdb_files=spcdb_files
        )
        add_bookmarks_to_pdf(
            pdfname,
            varlist,
            verbose=verbose
        )

    # 500hPa plots
    if "500hpa" in plots:
        pdfname = os.path.join(
            bdir,
            f"{colname}_500hPa.pdf"
        )

        compare_single_level(
            refds,
            refstr,
            devds,
            devstr,
            varlist=varlist,
            cmpres=cmpres,
            ilev=22,
            pdfname=pdfname,
            flip_ref=flip_ref,
            flip_dev=flip_dev,
            log_color_scale=log_color_scale,
            extra_title_txt=extra_title_txt,
            weightsdir=weightsdir,
            n_job=n_job,
            spcdb_files=spcdb_files,
        )
        add_bookmarks_to_pdf(
            pdfname,
            varlist,
            verbose=verbose
        )
    # Full-column zonal mean plots
    if "zonalmean" in plots:
        pdfname = os.path.join(
            bdir,
            f"{colname}_ZonalMean.pdf"
        )

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
            weightsdir=weightsdir,
            n_job=n_job,
            spcdb_files=spcdb_files
        )
        add_bookmarks_to_pdf(
            pdfname,
            varlist,
            verbose=verbose
        )

        # Strat_ZonalMean plots will use a log-pressure Y-axis, with
        # a range of 1..100 hPa, as per GCSC request. (bmy, 8/13/19)
        pdfname = os.path.join(
            bdir,
            f"{colname}_Strat_ZonalMean.pdf"
        )

        compare_zonal_mean(
            refds,
            refstr,
            devds,
            devstr,
            varlist=varlist,
            pdfname=pdfname,
            pres_range=[1, 100],
            log_yaxis=True,
            flip_ref=flip_ref,
            flip_dev=flip_dev,
            extra_title_txt=extra_title_txt,
            log_color_scale=log_color_scale,
            weightsdir=weightsdir,
            n_job=n_job,
            spcdb_files=spcdb_files
        )
        add_bookmarks_to_pdf(
            pdfname,
            varlist,
            verbose=verbose
        )

    # -------------------------------------------
    # Clean up
    # -------------------------------------------
    del refds
    del devds
    gc.collect()


def make_benchmark_aod_plots(
        ref,
        refstr,
        dev,
        devstr,
        spcdb_files,
        varlist=None,
        dst="./benchmark",
        subdst=None,
        cmpres=None,
        overwrite=False,
        verbose=False,
        log_color_scale=False,
        sigdiff_files=None,
        weightsdir='.',
        n_job=-1,
        time_mean=False,
):
    """
    Creates PDF files containing plots of column aerosol optical
    depths (AODs) for model benchmarking purposes.

    Args:
        ref: str
            Path name for the "Ref" (aka "Reference") data set.
        refstr: str
            A string to describe ref (e.g. version number)
        dev: str
            Path name for the "Dev" (aka "Development") data set.
            This data set will be compared against the "Reference"
            data set.
        devstr: str
            A string to describe dev (e.g. version number)
        spcdb_files : list
            Paths to species_database.yml files in Ref & Dev rundirs

    Keyword Args (optional):
        varlist: list of str
            List of AOD variables to plot.  If not passed, then all
            AOD variables common to both Dev and Ref will be plotted.
            Use the varlist argument to restrict the number of
            variables plotted to the pdf file when debugging.
            Default value: None
        dst: str
            A string denoting the destination folder where a
            PDF file  containing plots will be written.
            Default value: ./benchmark.
        subdst: str
            A string denoting the sub-directory of dst where PDF
            files containing plots will be written.  In practice,
            subdst is only needed for the 1-year benchmark output,
            and denotes a date string (such as "Jan2016") that
            corresponds to the month that is being plotted.
            Default value: None
        cmpres: string
            Grid resolution at which to compare ref and dev data, e.g. '1x1.25'
        overwrite: bool
            Set this flag to True to overwrite files in the
            destination folder (specified by the dst argument).
            Default value: False.
        verbose: bool
            Set this flag to True to print extra informational output.
            Default value: False
        log_color_scale: bool
            Set this flag to True to enable plotting data (not diffs)
            on a log color scale.
            Default value: False
        sigdiff_files: list of str
            Filenames that will contain the list of quantities having
            having significant differences in the column AOD plots.
            These lists are needed in order to fill out the benchmark
            approval forms.
            Default value: None
        weightsdir: str
            Directory in which to place (and possibly reuse) xESMF regridder
            netCDF files.
            Default value: '.'
        n_job: int
            Defines the number of simultaneous workers for parallel plotting.
            Set to 1 to disable parallel plotting. Value of -1 allows the
            application to decide.
            Default value: -1
        time_mean : bool
            Determines if we should average the datasets over time
            Default value: False
    """
    # ==================================================================
    # Initialization and also read data
    # ==================================================================

    # Create destination plots directory
    make_directory(dst, overwrite)

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

    # Replace whitespace in the ref and dev labels
    refstr = replace_whitespace(refstr)
    devstr = replace_whitespace(devstr)

    # Get the function that will read file(s) into a dataset
    reader = dataset_reader(time_mean, verbose=verbose)

    # Read the Ref dataset
    try:
        refds = reader(ref, drop_variables=SKIP_THESE_VARS)
    except (OSError, IOError, FileNotFoundError) as e:
        raise e(f"Could not find Ref file: {ref}") from e

    # Read the Dev dataset
    try:
        devds = reader(dev, drop_variables=SKIP_THESE_VARS)
    except (OSError, IOError, FileNotFoundError) as e:
        raise e(f"Could not find Ref file: {dev}") from e

    # Compute mean of data over the time dimension (if time_mean=True)
    if time_mean:
        refds = dataset_mean(refds)
        devds = dataset_mean(devds)

    # Create regridding files if necessary
    [_ for _ in create_regridders(refds, devds, weightsdir=weightsdir)]

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
    if varlist is None:
        vardict = compare_varnames(refds, devds, quiet=not verbose)
        cmn3D = vardict["commonvars3D"]
        varlist = [v for v in cmn3D if "AOD" in v]

    # Dictionary and list for new display names
    newvars = read_config_file(
        os.path.join(
            os.path.dirname(__file__),
            AOD_SPC
        ),
        quiet=True
    )
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
    # Avoid double-counting the following:
    # (1) Individual dust AODs, use the column AOD total instead
    # (2) SOA from aqueous isoprene, which is already accounted
    #     for in AODHyg550nm_OCPI.  Also see Github issue:
    #     https://github.com/geoschem/gcpy/issues/65
    for v in varlist:
        if "_bin" in v or "AODSOAfromAqIsoprene550nm" in v:
            continue
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
            raise ValueError(f"Could not find a display name for {v}")

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
        pdfname = os.path.join(
            aoddir,
            f"Aerosols_ColumnOptDepth_{subdst}.pdf"
        )
    else:
        pdfname = os.path.join(
            aoddir,
            "Aerosols_ColumnOptDepth.pdf"
        )

    diff_aod = []
    compare_single_level(
        refds,
        refstr,
        devds,
        devstr,
        varlist=newvarlist,
        cmpres=cmpres,
        ilev=0,
        pdfname=pdfname,
        log_color_scale=log_color_scale,
        extra_title_txt=extra_title_txt,
        sigdiff_list=diff_aod,
        weightsdir=weightsdir,
        n_job=n_job,
        spcdb_files=spcdb_files
    )
    diff_aod[:] = [v.replace("Column_AOD_", "") for v in diff_aod]
    add_bookmarks_to_pdf(
        pdfname,
        newvarlist,
        remove_prefix="Column_AOD_",
        verbose=verbose
    )

    # ==================================================================
    # Write the list of AOD quantities having significant differences,
    # which we will need to fill out the benchmark forms.
    # ==================================================================
    if sigdiff_files is not None:
        for filename in sigdiff_files:
            if "sfc" in filename:
                with open(filename, "a+", encoding=ENCODING) as f:
                    print("* Column AOD: ", file=f, end="")
                    for v in diff_aod:
                        print(f"{v} ", file=f, end="")
                    print(file=f)
                    f.close()

    # -------------------------------------------
    # Clean up
    # -------------------------------------------
    del refds
    del devds
    gc.collect()


def make_benchmark_mass_tables(
        ref,
        refstr,
        dev,
        devstr,
        spcdb_files,
        varlist=None,
        dst="./benchmark",
        subdst=None,
        overwrite=False,
        verbose=False,
        ref_hdr_label="",
        dev_hdr_label="",
        ref_met_extra=None,
        dev_met_extra=None,
):
    """
    Creates a text file containing global mass totals by species and
    category for benchmarking purposes.

    Args:
        reflist: str
            Pathname that will constitute
            the "Ref" (aka "Reference") data set.
        refstr: str
            A string to describe ref (e.g. version number)
        dev: list of str
            Pathname that will constitute
            the "Dev" (aka "Development") data set.  The "Dev"
            data set will be compared against the "Ref" data set.
        devstr: str
            A string to describe dev (e.g. version number)
        spcdb_files : list
            Paths to species_database.yml files in Ref & Dev rundirs

    Keyword Args (optional):
        varlist: list of str
            List of variables to include in the list of totals.
            If omitted, then all variables that are found in either
            "Ref" or "Dev" will be included.  The varlist argument
            can be a useful way of reducing the number of
            variables during debugging and testing.
            Default value: None
        dst: str
            A string denoting the destination folder where the file
            containing emissions totals will be written.
            Default value: ./benchmark
        subdst: str
            A string denoting the sub-directory of dst where PDF
            files containing plots will be written.  In practice,
            subdst is only needed for the 1-year benchmark output,
            and denotes a date string (such as "Jan2016") that
            corresponds to the month that is being plotted.
            Default value: None
        overwrite: bool
            Set this flag to True to overwrite files in the
            destination folder (specified by the dst argument).
            Default value: False
        verbose: bool
            Set this flag to True to print extra informational output.
            Default value: False.
        ref_hdr_label : str
            Label for Ref, placed after refstr in the file header
            Default value: ""
        dev_hdr_label : str
            Label for Dev, placed after devstr in the file header
            Default value: ""
        ref_met_extra: str
            Path to ref Met file containing area data for use with restart files
            which do not contain the Area variable.
            Default value: ''
        dev_met_extra: str
            Path to dev Met file containing area data for use with restart files
            which do not contain the Area variable.
            Default value: ''
    """
    # ==================================================================
    # Initialization
    # ==================================================================

    # Replace whitespace in the ref and dev labels
    refstr = replace_whitespace(refstr)
    devstr = replace_whitespace(devstr)

    # ==================================================================
    # Define destination directory
    # ==================================================================
    if os.path.isdir(dst) and not overwrite:
        msg = "Directory {} exists. Pass overwrite=True to overwrite " \
            + "files in that directory, if any."
        msg = msg.format(dst)
        raise ValueError(msg)
    if not os.path.isdir(dst):
        try:
            os.makedirs(dst)
        except FileExistsError:
            pass

    # ==================================================================
    # Read data from netCDF into Dataset objects
    # ==================================================================

    # Read data
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=xr.SerializationWarning)
        refds = xr.open_dataset(ref, drop_variables=SKIP_THESE_VARS)
        devds = xr.open_dataset(dev, drop_variables=SKIP_THESE_VARS)

    # ==================================================================
    # Update GCHP restart dataset (if any)
    # ==================================================================

    # If the data is from a GCHP restart file, rename variables and
    # flip levels to match the GEOS-Chem Classic naming and level
    # conventions.  Otherwise no changes will be made.
    refds = rename_and_flip_gchp_rst_vars(refds)
    devds = rename_and_flip_gchp_rst_vars(devds)

    # ==================================================================
    # Make sure that all necessary meteorological variables are found
    # ==================================================================
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=xr.SerializationWarning)

        # Find the area variable in Ref
        if ref_met_extra is None:
            ref_area = get_area_from_dataset(refds)
        else:
            ref_area = get_area_from_dataset(
                xr.open_dataset(
                    ref_met_extra,
                    drop_variables=SKIP_THESE_VARS
                )
            )

        # Find the area variable in Dev
        if dev_met_extra is None:
            dev_area = get_area_from_dataset(devds)
        else:
            dev_area = get_area_from_dataset(
                xr.open_dataset(
                    dev_met_extra,
                    drop_variables=SKIP_THESE_VARS
                )
            )

    # Find required meteorological variables in Ref
    # (or exit with an error if we can't find them)
    metvar_list = ["Met_DELPDRY", "Met_BXHEIGHT", "Met_TropLev"]
    refmet = get_variables_from_dataset(refds, metvar_list)
    devmet = get_variables_from_dataset(devds, metvar_list)

    # ==================================================================
    # Make sure that all necessary species are found
    # ==================================================================

    # Get lists of variables names in datasets
    vardict = compare_varnames(refds, devds, quiet=not verbose)
    commonvars = vardict["commonvars3D"]
    refonly = vardict['refonly']
    devonly = vardict['devonly']

    # Narrow down the lists to only include species
    commonspc = [v for v in commonvars if "SpeciesRst_" in v]
    refonlyspc = [v for v in refonly if v.startswith('SpeciesRst_')]
    devonlyspc = [v for v in devonly if v.startswith('SpeciesRst_')]

    # Add ref only species to dev dataset with all nan values
    if refonlyspc:
        for v in refonlyspc:
            devds[v] = devds[commonspc[0]]
            devds[v].data = np.full(devds[v].shape, np.nan)
            devds[v].attrs['units'] = refds[v].units
            commonspc.append(v)

    # Add dev only species to ref dataset with all nan values
    if devonlyspc:
        for v in devonlyspc:
            refds[v] = refds[commonspc[0]]
            refds[v].data = np.full(refds[v].shape, np.nan)
            devds[v].attrs['units'] = refds[v].units
            commonspc.append(v)

    # Set list of variables to print in mass table. If this list was passed
    # as argument, check that all the vars are now in commonspc to ensure
    # in both datasets.
    if varlist:
        for v in varlist:
            if v not in commonspc:
                raise ValueError(
                    f"{dst} folder error: Variable {v} in varlist passed to make_benchmark_mass_tables is not present in Ref and Dev datasets"
                )
    else:
        varlist = commonspc

    # Sort the list of species to be printed alphabetically
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
        "Ref_TropMask": ref_tropmask,
        "Dev_TropMask": dev_tropmask,
    }

    # ==================================================================
    # Create global mass table
    # ==================================================================
    if subdst is not None:
        mass_filename = f"GlobalMass_TropStrat_{subdst}.txt"
    else:
        mass_filename = "GlobalMass_TropStrat.txt"
    mass_file = os.path.join(dst, mass_filename)
    create_global_mass_table(
        refds,
        refstr,
        devds,
        devstr,
        varlist,
        met_and_masks,
        spcdb_files,
        ref_hdr_label=ref_hdr_label,
        dev_hdr_label=dev_hdr_label,
        outfilename=mass_file,
        verbose=verbose,
    )

    # ==================================================================
    # Create tropospheric mass table
    # ==================================================================

    # If a file name has not been specified, then use the "filename"
    # keyword argument.  Otherwise generate a default filename.
    if subdst is not None:
        mass_filename = f"GlobalMass_Trop_{subdst}.txt"
    else:
        mass_filename = 'GlobalMass_Trop.txt'
    mass_file = os.path.join(dst, mass_filename)
    create_global_mass_table(
        refds,
        refstr,
        devds,
        devstr,
        varlist,
        met_and_masks,
        spcdb_files,
        ref_hdr_label=ref_hdr_label,
        dev_hdr_label=dev_hdr_label,
        outfilename=mass_file,
        trop_only=True,
        verbose=verbose,
    )

    # -------------------------------------------
    # Clean up
    # -------------------------------------------
    del refds
    del devds
    gc.collect()


def make_benchmark_mass_accumulation_tables(
        ref_start,
        ref_end,
        refstr,
        refperiodstr,
        dev_start,
        dev_end,
        devstr,
        devperiodstr,
        spcdb_files,
        varlist=None,
        dst="./benchmark",
        subdst=None,
        overwrite=False,
        verbose=False,
        label="at end of simulation",
):
    """
    Creates a text file containing global mass totals by species and
    category for benchmarking purposes.

    Args:
        ref_start: list of str
            Pathname that will constitute
            the "Ref" (aka "Reference") data set.
        ref_end: list of str
            Pathname that will constitute
            the "Ref" (aka "Reference") data set.
        refstr: str
            A string to describe ref (e.g. version number)
        refperiodstr: str
            Ref simulation period start and end
        dev_start: list of str
            Pathname that will constitute
            the "Dev" (aka "Development") data set.  The "Dev"
            data set will be compared against the "Ref" data set.
        dev_end: list of str
            Pathname that will constitute
            the "Dev" (aka "Development") data set.  The "Dev"
            data set will be compared against the "Ref" data set.
        devstr: str
            A string to describe dev (e.g. version number)
        devperiodstr: str
            Dev simulation period start and end
        spcdb_files : list
            Paths to species_database.yml files in Ref & Dev rundirs

    Keyword Args (optional):
        varlist: list of str
            List of variables to include in the list of totals.
            If omitted, then all variables that are found in either
            "Ref" or "Dev" will be included.  The varlist argument
            can be a useful way of reducing the number of
            variables during debugging and testing.
            Default value: None
        dst: str
            A string denoting the destination folder where the file
            containing emissions totals will be written.
            Default value: ./benchmark
        subdst: str
            A string denoting the sub-directory of dst where PDF
            files containing plots will be written.  In practice,
            subdst is only needed for the 1-year benchmark output,
            and denotes a date string (such as "Jan2016") that
            corresponds to the month that is being plotted.
            Default value: None
        overwrite: bool
            Set this flag to True to overwrite files in the
            destination folder (specified by the dst argument).
            Default value: False
        verbose: bool
            Set this flag to True to print extra informational output.
            Default value: False.
    """
    # ==================================================================
    # Initialization
    # ==================================================================

    # Replace whitespace in the ref and dev labels
    refstr = replace_whitespace(refstr)
    devstr = replace_whitespace(devstr)

    # ==================================================================
    # Define destination directory
    # ==================================================================
    if os.path.isdir(dst) and not overwrite:
        msg = "Directory {} exists. Pass overwrite=True to overwrite " \
            + "files in that directory, if any."
        msg = msg.format(dst)
        raise ValueError(msg)
    if not os.path.isdir(dst):
        try:
            os.makedirs(dst)
        except FileExistsError:
            pass

    # ==================================================================
    # Read data from netCDF into Dataset objects
    # ==================================================================

    print('Creating mass accumulation tables from four restart files:')
    print(f'   Ref start: {ref_start}')
    print(f'   Ref end:   {ref_end}')
    print(f'   Dev start: {dev_start}')
    print(f'   Dev end:   {dev_end}')

    # Read data
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=xr.SerializationWarning)
        refSds = xr.open_dataset(ref_start, drop_variables=SKIP_THESE_VARS)
        refEds = xr.open_dataset(ref_end, drop_variables=SKIP_THESE_VARS)
        devSds = xr.open_dataset(dev_start, drop_variables=SKIP_THESE_VARS)
        devEds = xr.open_dataset(dev_end, drop_variables=SKIP_THESE_VARS)

    # ==================================================================
    # Update GCHP restart dataset if needed
    # ==================================================================

    # If the data is from a GCHP restart file, rename variables and
    # flip levels to match the GEOS-Chem Classic naming and level
    # conventions.  Otherwise no changes will be made.
    refSds = rename_and_flip_gchp_rst_vars(refSds)
    refEds = rename_and_flip_gchp_rst_vars(refEds)
    devSds = rename_and_flip_gchp_rst_vars(devSds)
    devEds = rename_and_flip_gchp_rst_vars(devEds)

    # Add area to start restart dataset if area in end but not start
    # Need to consider area variable names used in both GC-Classic and GCHP
    # Should put this in a function (todo)
    refSkeys = refSds.data_vars.keys()
    refEkeys = refEds.data_vars.keys()
    devSkeys = devSds.data_vars.keys()
    devEkeys = devEds.data_vars.keys()
    areaVars = ["Met_AREAM2", "AREA"]
    for areaVar in areaVars:
        if areaVar in refEkeys and areaVar not in refSkeys:
            refSds[areaVar] = refEds[areaVar]
        if areaVar in devEkeys and areaVar not in devSkeys:
            devSds[areaVar] = devEds[areaVar]

    # ==================================================================
    # Make sure that all necessary meteorological variables are found
    # ==================================================================
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=xr.SerializationWarning)

        # Find the area variable in Ref
        refs_area = get_area_from_dataset(refSds)
        refe_area = get_area_from_dataset(refEds)

        # Find the area variable in Dev
        devs_area = get_area_from_dataset(devSds)
        deve_area = get_area_from_dataset(devEds)

    # Find required meteorological variables in Ref
    # (or exit with an error if we can't find them)
    metvar_list = ["Met_DELPDRY", "Met_BXHEIGHT", "Met_TropLev"]
    refsmet = get_variables_from_dataset(refSds, metvar_list)
    refemet = get_variables_from_dataset(refEds, metvar_list)
    devsmet = get_variables_from_dataset(devSds, metvar_list)
    devemet = get_variables_from_dataset(devEds, metvar_list)

    # ==================================================================
    # Make sure that all necessary species are found
    # ==================================================================

    # Get lists of variables names in datasets
    vardict = compare_varnames(refSds, devSds, quiet=not verbose)
    commonvars = vardict["commonvars3D"]
    refonly = vardict['refonly']
    devonly = vardict['devonly']

    # Narrow down the lists to only include species
    commonspc = [v for v in commonvars if "SpeciesRst_" in v]
    refonlyspc = [v for v in refonly if v.startswith('SpeciesRst_')]
    devonlyspc = [v for v in devonly if v.startswith('SpeciesRst_')]

    # Add ref only species to dev dataset with all nan values
    if refonlyspc:
        for v in refonlyspc:
            devSds[v] = devSds[commonspc[0]]
            devSds[v].data = np.full(devSds[v].shape, np.nan)
            devSds[v].attrs['units'] = refSds[v].units
            devEds[v] = devEds[commonspc[0]]
            devEds[v].data = np.full(devEds[v].shape, np.nan)
            devEds[v].attrs['units'] = refEds[v].units
            commonspc.append(v)

    # Add dev only species to ref dataset with all nan values
    if devonlyspc:
        for v in devonlyspc:
            refSds[v] = refSds[commonspc[0]]
            refSds[v].data = np.full(refSds[v].shape, np.nan)
            devSds[v].attrs['units'] = refSds[v].units
            refEds[v] = refEds[commonspc[0]]
            refEds[v].data = np.full(refEds[v].shape, np.nan)
            devEds[v].attrs['units'] = refEds[v].units
            commonspc.append(v)

    # Set list of variables to print in mass table. If this list was passed
    # as argument, check that all the vars are now in commonspc to ensure
    # in both datasets.
    if varlist:
        for v in varlist:
            if v not in commonspc:
                raise ValueError(
                    f"{dst} folder error: Variable {v} in varlist passed to make_benchmark_mass_tables is not present in Ref and Dev datasets"
                )
    else:
        varlist = commonspc

    # Sort the list of species to be printed alphabetically
    varlist.sort()

    # ==================================================================
    # Create the mask arrays for the troposphere for Ref and Dev
    # ==================================================================
    refs_tropmask = get_troposphere_mask(refsmet)
    refe_tropmask = get_troposphere_mask(refemet)
    devs_tropmask = get_troposphere_mask(devsmet)
    deve_tropmask = get_troposphere_mask(devemet)

    # ==================================================================
    # Create a dictionary to hold all of the meterological
    # variables and mask variables that we need to pass down
    # ==================================================================
    met_and_masks = {
        "Refs_Area": refs_area,
        "Refe_Area": refe_area,
        "Devs_Area": devs_area,
        "Deve_Area": deve_area,
        "Refs_Delta_P": refsmet["Met_DELPDRY"],
        "Refe_Delta_P": refemet["Met_DELPDRY"],
        "Devs_Delta_P": devsmet["Met_DELPDRY"],
        "Deve_Delta_P": devemet["Met_DELPDRY"],
        "Refs_BxHeight": refsmet["Met_BXHEIGHT"],
        "Refe_BxHeight": refemet["Met_BXHEIGHT"],
        "Devs_BxHeight": devsmet["Met_BXHEIGHT"],
        "Deve_BxHeight": devemet["Met_BXHEIGHT"],
        "Refs_TropMask": refs_tropmask,
        "Refe_TropMask": refe_tropmask,
        "Devs_TropMask": devs_tropmask,
        "Deve_TropMask": deve_tropmask,
    }

    # ==================================================================
    # Create global mass accumulation table
    # ==================================================================
    if subdst is not None:
        mass_filename = f"GlobalMassAccumulation_TropStrat_{subdst}.txt"
    else:
        mass_filename = "GlobalMassAccumulation_TropStrat.txt"
    mass_file = os.path.join(dst, mass_filename)
    create_mass_accumulation_table(
        refSds,
        refEds,
        refstr,
        refperiodstr,
        devSds,
        devEds,
        devstr,
        devperiodstr,
        varlist,
        met_and_masks,
        label,
        outfilename=mass_file,
        verbose=verbose,
        spcdb_files=spcdb_files,
    )

    ## ==================================================================
    ## Create tropospheric mass table
    ## ==================================================================
    #if subdst is not None:
    #    mass_filename = f"GlobalMassAccumulation_Trop_{subdst}.txt"
    #else:
    #    mass_filename = 'GlobalMassAccumulation_Trop.txt'
    #mass_file = os.path.join(dst, mass_filename)
    #create_mass_accumulation_table(
    #    refSds,
    #    refEds,
    #    refstr,
    #    devSds,
    #    devEds,
    #    devstr,
    #    varlist,
    #    met_and_masks,
    #    label,
    #    outfilename=mass_file,
    #    trop_only=True,
    #    verbose=verbose,
    #    spcdb_files=spcdb_files,
    #)

    # -------------------------------------------
    # Clean up
    # -------------------------------------------
    del refSds
    del refEds
    del devSds
    del devEds
    gc.collect()




def make_benchmark_wetdep_plots(
        ref,
        refstr,
        dev,
        devstr,
        collection,
        spcdb_files,
        dst="./benchmark",
        cmpres=None,
        datestr=None,
        overwrite=False,
        verbose=False,
        benchmark_type="TransportTracersBenchmark",
        plots=["sfc", "500hpa", "zonalmean"],
        log_color_scale=False,
        normalize_by_area=False,
        areas=None,
        refmet=None,
        devmet=None,
        weightsdir='.',
        n_job=-1,
        time_mean=False,
):
    """
    Creates PDF files containing plots of species concentration
    for model benchmarking purposes.

    Args:
        ref: str
            Path name for the "Ref" (aka "Reference") data set.
        refstr: str
            A string to describe ref (e.g. version number)
        dev: str
            Path name for the "Dev" (aka "Development") data set.
            This data set will be compared against the "Reference"
            data set.
        devstr: str
            A string to describe dev (e.g. version number)
        collection: str
            String name of collection to plot comparisons for.
        spcdb_files : list
            Paths to species_database.yml files in Ref & Dev rundirs

    Keyword Args (optional):
        dst: str
            A string denoting the destination folder where a PDF
            file containing plots will be written.
            Default value: ./benchmark
        datestr: str
            A string with date information to be included in both the
            plot pdf filename and as a destination folder subdirectory
            for writing plots
            Default value: None
        benchmark_type: str
            A string denoting the type of benchmark output to plot, options are
            FullChemBenchmark, TransportTracersBenchmark, or CH4Benchmark.
            Default value: "FullChemBenchmark"
        overwrite: bool
            Set this flag to True to overwrite files in the
            destination folder (specified by the dst argument).
            Default value: False.
        verbose: bool
            Set this flag to True to print extra informational output.
            Default value: False.
        plots: list of strings
            List of plot types to create.
            Default value: ['sfc', '500hpa', 'zonalmean']
        normalize_by_area: bool
            Set this flag to true to enable normalization of data
            by surfacea area (i.e. kg s-1 --> kg s-1 m-2).
            Default value: False
        areas: dict of xarray DataArray:
            Grid box surface areas in m2 on Ref and Dev grids.
            Default value: None
        refmet: str
            Path name for ref meteorology
            Default value: None
        devmet: str
            Path name for dev meteorology
            Default value: None
        n_job: int
            Defines the number of simultaneous workers for parallel plotting.
            Set to 1 to disable parallel plotting. Value of -1 allows the
            application to decide.
            Default value: -1
        time_mean : bool
            Determines if we should average the datasets over time
            Default value: False
    """

    # Create destination plot directory
    make_directory(dst, overwrite)

    # Make a collection subdirectory
    targetdst = os.path.join(dst, collection)
    if not os.path.isdir(targetdst):
        os.mkdir(targetdst)

    # If datestr is passed, create a further subdirectory
    if datestr is not None:
        targetdst = os.path.join(targetdst, datestr)
        if not os.path.isdir(targetdst):
            os.mkdir(targetdst)

    # Replace whitespace in the ref and dev labels
    refstr = replace_whitespace(refstr)
    devstr = replace_whitespace(devstr)

    # Get the function that will read file(s) into a dataset
    reader = dataset_reader(time_mean, verbose=verbose)

    # Open datasets
    refds = reader(ref, drop_variables=SKIP_THESE_VARS)
    devds = reader(dev, drop_variables=SKIP_THESE_VARS)

    # Open met datasets if passed as arguments
    refmetds = None
    devmetds = None
    if refmet is not None:
        refmetds = reader(refmet, drop_variables=SKIP_THESE_VARS)
    if devmet is not None:
        devmetds = reader(devmet, drop_variables=SKIP_THESE_VARS)

    # Compute mean of data over the time dimension (if time_mean=True)
    if time_mean:
        refds = dataset_mean(refds)
        devds = dataset_mean(devds)
        if refmet is not None:
            refmetds = dataset_mean(refmetds)
        if devmet is not None:
            devmetds = dataset_mean(devmetds)

    # Make sure that Ref and Dev datasets have the same variables.
    # Variables that are in Ref but not in Dev will be added to Dev
    # with all missing values (NaNs). And vice-versa.
    # Turn this off for now since add_missing_variables inserts GCC area into
    # GCHP files, which causes problems with area normalization (ewl)
    #[refds, devds] = add_missing_variables(refds, devds)

    # Get list of variables in collection
    vardict = compare_varnames(refds, devds, quiet=not verbose)
    varlist = [v for v in vardict["commonvars3D"] if collection + "_" in v]
    varlist.sort()

    # Surface plots
    if "sfc" in plots:
        if datestr is not None:
            plotfilename = f"{collection}_Surface_{datestr}.pdf"
        else:
            plotfilename = f"{collection}_Surface.pdf"
        pdfname = os.path.join(targetdst, plotfilename)
        compare_single_level(
            refds,
            refstr,
            devds,
            devstr,
            varlist=varlist,
            cmpres=cmpres,
            ilev=0,
            refmet=refmetds,
            devmet=devmetds,
            pdfname=pdfname,
            normalize_by_area=normalize_by_area,
            extra_title_txt=datestr,
            weightsdir=weightsdir,
            n_job=n_job,
            spcdb_files=spcdb_files,
        )
        add_bookmarks_to_pdf(
            pdfname,
            varlist,
            remove_prefix=collection + '_',
            verbose=verbose)

    # 500 hPa plots
    if "500hpa" in plots:
        if datestr is not None:
            plotfilename = f"{collection}_500hPa_{datestr}.pdf"
        else:
            plotfilename = f"{collection}_500hPa.pdf"
        pdfname = os.path.join(targetdst, plotfilename)
        compare_single_level(
            refds,
            refstr,
            devds,
            devstr,
            varlist=varlist,
            cmpres=cmpres,
            ilev=22,
            refmet=refmetds,
            devmet=devmetds,
            pdfname=pdfname,
            normalize_by_area=normalize_by_area,
            extra_title_txt=datestr,
            weightsdir=weightsdir,
            n_job=n_job,
            spcdb_files=spcdb_files,
        )
        add_bookmarks_to_pdf(
            pdfname,
            varlist,
            remove_prefix=collection + '_',
            verbose=verbose
        )

    # Zonal mean plots
    if "zonalmean" in plots or "zm" in plots:

        # Full column
        if datestr is not None:
            plotfilename = f"{collection}_FullColumn_ZonalMean_{datestr}.pdf"
        else:
            plotfilename = f"{collection}_FullColumn_ZonalMean.pdf"
        pdfname = os.path.join(targetdst, plotfilename)
        compare_zonal_mean(
            refds,
            refstr,
            devds,
            devstr,
            varlist=varlist,
            refmet=refmetds,
            devmet=devmetds,
            pdfname=pdfname,
            log_color_scale=log_color_scale,
            normalize_by_area=normalize_by_area,
            extra_title_txt=datestr,
            weightsdir=weightsdir,
            n_job=n_job,
            spcdb_files=spcdb_files,
        )
        add_bookmarks_to_pdf(
            pdfname,
            varlist,
            remove_prefix=collection + '_',
            verbose=verbose
        )

        # Stratosphere
        if datestr is not None:
            plotfilename = f"{collection}_Strat_ZonalMean_{datestr}.pdf"
        else:
            plotfilename = f"{collection}_Strat_ZonalMean.pdf"
        pdfname = os.path.join(targetdst, plotfilename)
        compare_zonal_mean(
            refds,
            refstr,
            devds,
            devstr,
            varlist=varlist,
            refmet=refmetds,
            devmet=devmetds,
            pdfname=pdfname,
            pres_range=[1, 100],
            log_yaxis=True,
            extra_title_txt=datestr,
            normalize_by_area=normalize_by_area,
            weightsdir=weightsdir,
            n_job=n_job,
            spcdb_files=spcdb_files
        )
        add_bookmarks_to_pdf(
            pdfname,
            varlist,
            remove_prefix=collection + '_',
            verbose=verbose
        )

    # -------------------------------------------
    # Clean up
    # -------------------------------------------
    del refds
    del devds
    del refmetds
    del devmetds
    gc.collect()


def make_benchmark_aerosol_tables(
        devdir,
        devlist_aero,
        devlist_spc,
        devlist_met,
        devstr,
        year,
        days_per_mon,
        spcdb_files,
        dst='./benchmark',
        overwrite=False,
        is_gchp=False,
):
    """
    Compute FullChemBenchmark aerosol budgets & burdens

    Args:
        devdir: str
            Path to development ("Dev") data directory
        devlist_aero: list of str
            List of Aerosols collection files (different months)
        devlist_spc: list of str
            List of SpeciesConc collection files (different months)
        devlist_met: list of str
            List of meteorology collection files (different months)
        devstr: str
            Descriptive string for datasets (e.g. version number)
        year: str
            The year of the benchmark simulation (e.g. '2016').
        days_per_month: list of int
            List of number of days per month for all months
        spcdb_files : list
            Paths to species_database.yml files in Ref & Dev rundirs

    Keyword Args (optional):
        dst: str
            Directory where budget tables will be created.
            Default value: './benchmark'
        overwrite: bool
            Overwrite burden & budget tables? (default=True)
            Default value: False
        is_gchp: bool
            Whether datasets are for GCHP
            Default value: False

    """
    # Create destination directory
    make_directory(dst, overwrite)

    # List of species (and subsets for the trop & strat)
    species_list = [
        "BCPI", "OCPI", "SO4", "DST1", "DST2", "DST3", "DST4",
        "DSTbin1", "DSTbin2", "DSTbin3", "DSTbin4", "DSTbin5",
        "DSTbin6", "DSTbin7", "SALA", "SALC"
    ]

    # Read the species database files in the Ref & Dev rundirs, and
    # return a dict containing metadata for the union of species.
    # We'll need properties such as mol. wt. for unit conversions, etc.
    _, dev_metadata = read_species_metadata(spcdb_files, quiet=True)

    # Get the list of relevant AOD diagnostics from a YAML file
    ifile= AOD_SPC
    aod = read_config_file(
        os.path.join(
            os.path.dirname(__file__),
            ifile,
        ),
        quiet=True
    )
    ds_aer = xr.open_mfdataset(
        devlist_aero,
        drop_variables=SKIP_THESE_VARS
    )
    ds_spc = xr.open_mfdataset(
        devlist_spc,
        drop_variables=SKIP_THESE_VARS
    )
    ds_met = xr.open_mfdataset(
        devlist_met,
        drop_variables=SKIP_THESE_VARS
    )

    # Rename SpeciesConc_ to SpeciesConcVV_ for consistency with new
    # naming introduced in GEOS-Chem 14.1.0
    ds_spc = rename_speciesconc_to_speciesconcvv(ds_spc)

    # Trim species_list to the variables that are only in the dataset
    varlist = [
        var for var in species_list if f"SpeciesConcVV_{var}" in ds_spc.data_vars
    ]

    # Molecular weights [g mol-1], as taken from the species database
    mw = {}
    full_names = {}
    for var in varlist:
        full_names[var] = dev_metadata[var]["FullName"].strip()
        mw[var] = dev_metadata[var]["MW_g"]
    mw["Air"] = MW_AIR_g

    # Get troposphere mask
    tropmask = get_troposphere_mask(ds_met)

    # Get number of months
    n_mon = len(days_per_mon)

    # ------------------------------------------------------------------
    # Surface area
    # (kludgey but it works - revisit this)
    # ------------------------------------------------------------------

    # Get number of vertical levels
    N_LEVS = ds_spc.dims["lev"]

    if is_gchp:
        area = ds_met["Met_AREAM2"].values
        a = area.shape
        area_m2 = np.zeros([a[0], N_LEVS, a[1], a[2], a[3]])
        for t in range(n_mon):
            for k in range(N_LEVS):
                area_m2[t, k, :, :, :] = area[t, :, :, :]
        total_area_m2 = np.sum(area_m2[0, 0, :, :, :])
    else:
        area = ds_met["AREA"].values
        a = area.shape
        area_m2 = np.zeros([a[0], N_LEVS, a[1], a[2]])
        for t in range(n_mon):
            for k in range(N_LEVS):
                area_m2[t, k, :, :] = area[t, :, :]
        total_area_m2 = np.sum(area_m2[0, 0, :, :])

    # ------------------------------------------------------------------
    # Conversion factors and time increments
    # ------------------------------------------------------------------
    # v/v dry --> Tg
    vv_to_Tg = {}
    for var in varlist:
        if f"SpeciesConcVV_{var}" in ds_spc.data_vars:
            vv_to_Tg[var] = \
                ds_met["Met_AD"].values * (mw[var] / mw["Air"]) * 1e-9

    # Days in the benchmark duration
    days_per_yr = np.sum(days_per_mon)

    # ------------------------------------------------------------------
    # Define function to print tables
    # ------------------------------------------------------------------
    def print_aerosol_metrics(data, varlist, namelist, filename, title, label):

        with open(filename, "w+", encoding=ENCODING) as ofile:

            # Print top header
            print("%" * 79, file=ofile)
            print(f" {title} for {year} in {devstr}", file=ofile)
            print(" (weighted by the number of days per month)", file=ofile)
            print("%" * 79, file=ofile)
            line = "\n" + " " *67 + "Strat         Trop         Strat+Trop\n"
            line += " " * 67 + "-----------   ----------   ----------"
            print(line, file=ofile)

            # Print data
            for var in varlist:
                line = f"{namelist[var] : <41} ({var : <7}) {label} :  {data[var + '_s']:11.9f}   {data[var + '_t']:10.8f}   {data[var + '_f']:10.8f}\n"
                print(line, file=ofile)

    # ------------------------------------------------------------------
    # Compute aerosol burdens [Tg] and print
    # ------------------------------------------------------------------

    # Table info
    filename = f"{dst}/Aerosol_Burdens.txt"
    if n_mon == 12:
        title = "Annual average global aerosol burdens"
    else:
        title = f"Average global aerosol burdens across {n_mon} months"
    label = "burden [Tg]"

    # Initialize
    q = {}
    q_sum_f = np.zeros(n_mon)
    q_sum_t = np.zeros(n_mon)
    q_sum_s = np.zeros(n_mon)
    burdens = {}

    # Define the axes we need to sum over to make monthly sums
    if is_gchp:
        sum_axes = (1, 2, 3, 4)
    else:
        sum_axes = (1, 2, 3)

    # Loop over species
    for var in varlist:

        # Whole-atmosphere and trop-only quantities [g]
        # NOTE: DryDep is by nature trop-only
        varname = f"SpeciesConcVV_{var}"
        q[var + "_f"] = ds_spc[varname].values * vv_to_Tg[var]
        q[var + "_t"] = np.ma.masked_array(q[var + "_f"], tropmask)

        # Compute monthly sums, weighted by the number of days per month
        q_sum_f = np.sum(q[var + "_f"], axis=sum_axes) * days_per_mon
        q_sum_t = np.sum(q[var + "_t"], axis=sum_axes) * days_per_mon
        q_sum_s = q_sum_f - q_sum_t

        # Compute annual averages
        burdens[var + "_f"] = np.sum(q_sum_f) / days_per_yr
        burdens[var + "_t"] = np.sum(q_sum_t) / days_per_yr
        burdens[var + "_s"] = np.sum(q_sum_s) / days_per_yr

    print_aerosol_metrics(burdens, varlist, full_names, filename, title, label)

    # ------------------------------------------------------------------
    # Define function to print tables
    # ------------------------------------------------------------------
    def print_aods(data, varlist, filename, title):

        with open(filename, "w+", encoding=ENCODING) as ofile:

            # Print top header
            print("%" * 79, file=ofile)
            print(f" {title} for {year} in {devstr}", file=ofile)
            print(" (weighted by the number of days per month)", file=ofile)
            print("%" * 79, file=ofile)
            line = "\n" + " " *32 + "Strat         Trop         Strat+Trop\n"
            line += " " * 32 + "-----------   ----------   ----------"
            print(line, file=ofile)

            # Print data
            for var in varlist:
                line = f"{var : <4} column optical depth [1]:  {data[var + '_s']:11.9f}   {data[var + '_t']:10.8f}   {data[var + '_f']:10.8f}\n"
                print(line, file=ofile)

    # ------------------------------------------------------------------
    # Compute average AOD's [Tg] and print
    # ------------------------------------------------------------------

    # Table info
    filename = f"{dst}/Global_Mean_AOD.txt"
    if n_mon == 12:
        title = "Annual average global AODs"
    else:
        title = f"Average global AODs across {n_mon} months"
    label = "mean AOD [1]"

    # Initialize
    q = {}
    q_sum_f = np.zeros(n_mon)
    q_sum_t = np.zeros(n_mon)
    q_sum_s = np.zeros(n_mon)
    aods = {}

    # Define axes to sum over, and total surface area
    if is_gchp:
        sum_axes = (1, 2, 3, 4)
    else:
        sum_axes = (1, 2, 3)

    # List of AOD species to plot
    varlist = [var for var in aod.keys() if "Dust" in var or "Hyg" in var]
    if is_gchp:
        varlist = [var.replace('550nm', 'WL1') for var in varlist]

    # Loop over AOD variables
    aod_list = []
    for varname in varlist:

        # Extract the s
        if "Dust" in varname:
            var = "Dust"
        elif "Total" in varname:
            continue
        else:
            var = varname.split("_")[1]
        aod_list.append(var)

        # Whole-atmosphere AOD [1]
        q[var + "_f"] = ds_aer[varname].values

        # Tropospheric-only AOD [1]
        q[var + "_t"] = np.ma.masked_array(q[var + "_f"], tropmask)

        # Create monthly sums, weighted by the number of days per month
        q_sum_f = np.sum(q[var + "_f"] * area_m2, axis=sum_axes) * days_per_mon
        q_sum_t = np.sum(q[var + "_t"] * area_m2, axis=sum_axes) * days_per_mon
        q_sum_s = q_sum_f - q_sum_t

        # Take annual averages
        aods[var + "_f"] = np.sum(q_sum_f) / total_area_m2 / days_per_yr
        aods[var + "_t"] = np.sum(q_sum_t) / total_area_m2 / days_per_yr
        aods[var + "_s"] = np.sum(q_sum_s) / total_area_m2 / days_per_yr

    print_aods(aods, aod_list, filename, title)

    # ------------------------------------------------------------------
    # Clean up
    # ------------------------------------------------------------------
    del ds_aer
    del ds_spc
    del ds_met
    gc.collect()


def make_benchmark_operations_budget(
        refstr,
        reffiles,
        devstr,
        devfiles,
        ref_interval,
        dev_interval,
        spcdb_files,
        benchmark_type=None,
        label=None,
        col_sections=["Full", "Trop", "PBL", "FixedLevs", "Strat"],
        operations=["Chemistry", "Convection", "EmisDryDep",
                    "Mixing", "Transport", "WetDep"],
        compute_accum=True,
        compute_restart=False,
        require_overlap=False,
        dst='.',
        species=None,
        overwrite=True,
        verbose=False,
):
    """
    Prints the "operations budget" (i.e. change in mass after
    each operation) from a GEOS-Chem benchmark simulation.

    Args:
        refstr: str
            Labels denoting the "Ref" versions
        reffiles: list of str
            Lists of files to read from the "Ref" version.
        devstr: str
            Labels denoting the "Dev" versions
        devfiles: list of str
            Lists of files to read from "Dev" version.
        interval: float
            Number of seconds in the diagnostic interval.
        spcdb_files : list
            Paths to species_database.yml files in Ref & Dev rundirs

    Keyword Args (optional):
        benchmark_type: str
            A string denoting the type of benchmark output to plot, options are
            FullChemBenchmark, TransportTracersBenchmark, or CH4Benchmark.
            Default value: None
        label: str
            Contains the date or date range for each dataframe title.
            Default value: None
        col_sections: list of str
            List of column sections to calculate global budgets for. May
            include Strat eventhough not calculated in GEOS-Chem, but Full
            and Trop must also be present to calculate Strat.
            Default value: ["Full", "Trop", "PBL", "FixedLevs", "Strat"]
        operations: list of str
            List of operations to calculate global budgets for. Accumulation
            should not be included. It will automatically be calculated if
            all GEOS-Chem budget operations are passed and optional arg
            compute_accum is True.
            Default value: ["Chemistry","Convection","EmisDryDep",
                            "Mixing","Transport","WetDep"]
        compute_accum: bool
            Optionally turn on/off accumulation calculation. If True, will
            only compute accumulation if all six GEOS-Chem operations budgets
            are computed. Otherwise a message will be printed warning that
            accumulation will not be calculated.
            Default value: True
        compute_accum: bool
            Optionally turn on/off accumulation calculation. If True, will
            only compute accumulation if all six GEOS-Chem operations budgets
            are computed. Otherwise a message will be printed warning that
            accumulation will not be calculated.
            Default value: True
        compute_restart: bool
            Optionally turn on/off calculation of mass change based on restart
            file. Only functional for "Full" column section.
            Default value: False
        require_overlap: bool
            Whether to calculate budgets for only species that are present in
            both Ref or Dev.
            Default value: False
        dst: str
            Directory where plots & tables will be created.
            Default value: '.' (directory in which function is called)
        species: list of str
            List of species for which budgets will be created.
            Default value: None (all species)
        overwrite: bool
            Denotes whether to overwrite existing budget file.
            Default value: True
        verbose: bool
            Set this switch to True if you wish to print out extra
            informational messages.
            Default value: False
    """
    # Replace whitespace in the ref and dev labels
    refstr = replace_whitespace(refstr)
    devstr = replace_whitespace(devstr)

    # ------------------------------------------
    # Column sections
    # ------------------------------------------

    # Print info. Only allow Strat if Trop and Full are present
    print("Column sections:")
    for col_section in col_sections:
        print(f"  {col_section}")
    n_sections = len(col_sections)
    compute_strat = False
    if "Strat" in col_sections:
        compute_strat = True
        if "Full" not in col_sections or "Trop" not in col_sections:
            msg = "Strat budget table requires Full and Trop column sections"
            raise ValueError(msg)
        print("*** Will compute Strat budget as difference of Full"
              " and Trop ***")

    # Make a subset of column sections not including Strat
    gc_sections = [r for r in col_sections if "Strat" not in r]

    # ------------------------------------------
    # Operations
    # ------------------------------------------

    # Rename the optional argument to be clear it is GEOS-Chem budget
    # operations
    gc_operations = operations

    # Handle whether to compute accumulation.
    all_operations = gc_operations
    if compute_accum and len(gc_operations) == 6:
        all_operations = gc_operations + ["ACCUMULATION"]
    if compute_restart:
        all_operations = gc_operations + ["RESTART"]
    n_ops = len(all_operations)

    # Print info
    print("Operations:")
    for all_operation in all_operations:
        print(f"  {all_operation}")
    if compute_accum:
        if "ACCUMULATION" in all_operations:
            print("*** Will compute ACCUMULATION operation as sum of all "
                  "GEOS-Chem operations ***")
        else:
            print("***Will not compute ACCUMULATION since not all GEOS-Chem"
                  " operation budgets will be computed.")
    if compute_restart:
        print("*** Will compute RESTART operation as mass change "
              "based on simulation start and end restart files ***")

    # ------------------------------------------
    # Read data
    # ------------------------------------------

    # Assume this will be annual budget if interval greater than 3e7 sec
    annual = ref_interval > 3.0e7

    # Read data from disk (either one month or 12 months)
    print('Opening ref and dev data')
    skip_vars = SKIP_THESE_VARS
    if annual:
        ref_ds = xr.open_mfdataset(reffiles, drop_variables=skip_vars)
        dev_ds = xr.open_mfdataset(devfiles, drop_variables=skip_vars)
    else:
        ref_ds = xr.open_dataset(reffiles, drop_variables=skip_vars)
        dev_ds = xr.open_dataset(devfiles, drop_variables=skip_vars)

    # TODO: Add section for reading files for computing mass from restart file

    # ------------------------------------------
    # Species
    # ------------------------------------------

    # Get information about variables in data files
    vardict = compare_varnames(ref_ds, dev_ds, quiet=True)
    refonly = vardict["refonly"]
    devonly = vardict["devonly"]
    cmnvars = vardict["commonvars2D"]

    # Reduce each list to only variables containing "Budget" and not "Strat"
    refonly = [v for v in refonly if "Budget" in v and "Strat" not in v]
    devonly = [v for v in devonly if "Budget" in v and "Strat" not in v]
    cmnvars = [v for v in cmnvars if "Budget" in v and "Strat" not in v]

    # Special handling for fixed level budget diagnostic
    # Get variable name prefix, e.g. Levs1to35. Check that all fixed level
    # vars have the same prefix. Update section names used in table.
    fixedlevvars = [v for v in cmnvars if "Budget" in v and "Levs" in v]
    if fixedlevvars is not None:
        fixedlevnames = [v[v.index('Levs'):].split("_")[0] for v in fixedlevvars]
        if len(set(fixedlevnames)) > 1:
            msg = "Budget fixed level diagnostic name must be constant!"
            raise ValueError(msg)
        col_sections = [v.replace('FixedLevs',fixedlevnames[0]) for v in col_sections]
        gc_sections = [v.replace('FixedLevs',fixedlevnames[0]) for v in gc_sections]

    # Get the species list, depending on if species was passed as argument.
    if species is not None:
        spclist = species
    else:
        # For each column section, get the union or intersection (depending
        # on optional argument require_overlap) of all budget diagnostic
        # species and sort alphabetically. A species is counted even if only
        # present for as few as one operation.
        spclists = {}
        for gc_section in gc_sections:
            refspc = [v.split("_")[1] for v in refonly if gc_section in v]
            devspc = [v.split("_")[1] for v in devonly if gc_section in v]
            cmnspc = [v.split("_")[1] for v in cmnvars if gc_section in v]
            if require_overlap:
                section_spc = list(set(cmnspc))
            else:
                section_spc = list(set(refspc + devspc + cmnspc))
            section_spc.sort()
            spclists[gc_section] = section_spc

        # If calculating Strat, use intersection of Trop and Full
        if "Strat" in col_sections:
            spclists["Strat"] = list(set(spclists["Trop"] + spclists["Full"]))

        # For now, define one species list as species combined across all
        # column sections. Compute budgets of these species for each section.
        # NOTE: eventually would want to be able to compute budget for only
        # species per column section, not all species for all section. If that
        # is implemented then handling of optional species list needs to be
        # revisited.
        spclist = []
        for s in spclists:
            spclist = spclist + spclists[s]

    # Make list items unique and put in alphabetical order
    spclist = list(set(spclist))
    spclist.sort()
    n_spc = len(spclist)

    # ------------------------------------------
    # Concentration units and if a wetdep species
    # ------------------------------------------

    # Load a YAML file containing species properties
    ref_metadata, dev_metadata = read_species_metadata(
        spcdb_files,
        quiet=True
    )

    # Determine what the converted units and conversion factor should be
    # based on benchmark type and species (tracer) name. Assume raw data [kg/s]
    ref_conv_fac = {}
    dev_conv_fac = {}
    units = {}
    is_wetdep = {}
    for spc in spclist:

        # Identify wetdep species
        is_wetdep[spc] = None
        ref_species_metadata = ref_metadata.get(spc)
        dev_species_metadata = dev_metadata.get(spc)
        if ref_species_metadata is not None and \
           dev_species_metadata is not None:
            is_wetdep[spc] = \
                ref_species_metadata.get("Is_WetDep") or \
                dev_species_metadata.get("Is_WetDep")

        # Unit conversion factors and units
        ref_conv_fac[spc] = ref_interval * 1e-6
        dev_conv_fac[spc] = dev_interval * 1e-6
        units[spc] = '[Gg]'
        if benchmark_type is not None:
            if "TransportTracers" in benchmark_type and "Tracer" not in spc:
                ref_conv_fac[spc] = ref_interval
                dev_conv_fac[spc] = dev_interval
                if annual:
                    units[spc] = '[kg/yr]'
                else:
                    units[spc] = '[kg]'
            elif annual:
                ref_conv_fac[spc] = ref_interval * 1e-9
                dev_conv_fac[spc] = dev_interval * 1e-9
                units[spc] = '[Tg/yr]'

    # ------------------------------------------
    # Create dataframe
    # ------------------------------------------
    columns = ["Column_Section", "Species", "Operation", "Ref_raw", "Dev_raw",
               "Units_converted", "Ref", "Dev", "Diff", "Pct_diff"]

    # Make column data to initialize with
    col_section = list(itertools.chain.from_iterable(
        itertools.repeat(x, n_ops * n_spc) for x in col_sections))
    col_spc = list(itertools.chain.from_iterable(
        itertools.repeat(x, n_ops) for x in spclist)) * n_sections
    col_ops = all_operations * n_sections * n_spc

    # Put the column data together into a dictionary
    data = {
        'Species': col_spc,
        'Operation': col_ops,
        'Column_Section': col_section,
    }

    # Create the dataframe from the data dictionary and column names list
    df = pd.DataFrame(data, columns=columns)

    # ------------------------------------------
    # Populate dataframe for GEOS-Chem operations and column sections
    # ------------------------------------------
    print('Calculating budgets for all data operations and column sections...')

    # Loop over sections (only those with data in files)
    for gc_section in gc_sections:

        # Keep track of progress in log
        print(f"  {gc_section}")

        # Loop over species in that section
        for i, spc in enumerate(spclist):

            # Keep track of progress (debugging print)
            #if (i + 1) % 50 == 0:
            #    print(f"  {gc_section}: species {i + 1} of {n_spc}")

            # Loop over operations (only those with data in files)
            for gc_operation in gc_operations:

                # Get the dataframe row to fill. Skip if not found.
                dfrow = (df["Column_Section"] == gc_section) \
                    & (df["Species"] == spc) \
                    & (df["Operation"] == gc_operation)
                if not any(dfrow):
                    continue

                # Get the variable name in the datasets
                varname = "Budget" + gc_operation + gc_section + "_" + spc

                # Calculate Ref and dev raw as global sum
                if varname in ref_ds.data_vars.keys():
                    refraw = ref_ds[varname].values.sum()
                elif gc_operation == "WetDep" and is_wetdep[spc] is None:
                    refraw = 0.0
                else:
                    refraw = np.nan
                if varname in dev_ds.data_vars.keys():
                    devraw = dev_ds[varname].values.sum()
                elif gc_operation == "WetDep" and is_wetdep[spc] is None:
                    devraw = 0.0
                else:
                    devraw = np.nan

                # Convert units
                refconv = refraw * ref_conv_fac[spc]
                devconv = devraw * dev_conv_fac[spc]

                # Calculate diff and % diff from conc with converted units
                if not np.isnan(refconv) and not np.isnan(devconv):
                    diff = devconv - refconv
                    try:
                        pctdiff = diff / refconv * 100
                    except BaseException:
                        pctdiff = np.nan
                else:
                    diff = np.nan
                    pctdiff = np.nan

                # Fill dataframe
                df.loc[dfrow, "Ref_raw"] = refraw
                df.loc[dfrow, "Dev_raw"] = devraw
                df.loc[dfrow, "Units_converted"] = units[spc]
                df.loc[dfrow, "Ref"] = refconv
                df.loc[dfrow, "Dev"] = devconv
                df.loc[dfrow, "Diff"] = diff
                df.loc[dfrow, "Pct_diff"] = pctdiff

    # ------------------------------------------
    # Compute Strat for each data operation (if applicable)
    # ------------------------------------------
    if compute_strat:
        print('Computing Strat budgets from Trop and Full')

        # Loop over species
        for i, spc in enumerate(spclist):

            # Keep track of progress (debugging print)
            #if (i + 1) % 50 == 0:
            #    print(f"  Strat: species {i + 1} of {n_spc}")

            # Loop over operations (only those with data in files)
            for gc_operation in gc_operations:

                # Get the strat dataframe row to fill. Skip if not found.
                dfrow = (df["Column_Section"] == "Strat") \
                    & (df["Species"] == spc) \
                    & (df["Operation"] == gc_operation)
                if not any(dfrow):
                    continue

                # Get the "Full" row to use in the calculation
                dfrow_full = (df["Column_Section"] == "Full") \
                    & (df["Species"] == spc) \
                    & (df["Operation"] == gc_operation)
                if not any(dfrow_full):
                    continue

                # Get the "Trop" row to use in the calculation
                dfrow_trop = (df["Column_Section"] == "Trop") \
                    & (df["Species"] == spc) \
                    & (df["Operation"] == gc_operation)
                if not any(dfrow_trop):
                    continue

                # Calculate strat concentrations as Full minus Trop. Do
                # not bother with raw value computation.
                refstrat = df.loc[dfrow_full, "Ref"].values[0] \
                    - df.loc[dfrow_trop, "Ref"].values[0]
                devstrat = df.loc[dfrow_full, "Dev"].values[0] \
                    - df.loc[dfrow_trop, "Dev"].values[0]

                # Calculate diff and % diff
                if not np.isnan(refstrat) and not np.isnan(devstrat):
                    diff = devstrat - refstrat
                    try:
                        pctdiff = diff / refstrat * 100
                    except BaseException:
                        pctdiff = np.nan
                else:
                    diff = np.nan
                    pctdiff = np.nan

                # Fill dataframe
                df.loc[dfrow, "Units_converted"] = units[spc]
                df.loc[dfrow, "Ref"] = refstrat
                df.loc[dfrow, "Dev"] = devstrat
                df.loc[dfrow, "Diff"] = diff
                df.loc[dfrow, "Pct_diff"] = pctdiff

    # ------------------------------------------
    # Compute Accumulation for each column section (if applicable)
    # ------------------------------------------
    if compute_accum:
        print('Computing ACCUMULATION operation budgets...')

        # Loop over all column sections
        for col_section in col_sections:

            # Keep track of progress in log
            print(f"  {col_section}")

            # Loop over species
            for i, spc in enumerate(spclist):

                # Keep track of progress (debugging print)
                #if (i + 1) % 50 == 0:
                #    print(f"  {col_section}: species {i + 1} of {n_spc}")

                # Get the accumulation dataframe row to fill.Skip if not found.
                dfrow = (df["Column_Section"] == col_section) \
                    & (df["Species"] == spc) \
                    & (df["Operation"] == "ACCUMULATION")
                if not any(dfrow):
                    continue

                # Get the rows to sum
                dfrows_to_sum = (df["Column_Section"] == col_section) \
                    & (df["Species"] == spc) \
                    & (df["Operation"].isin(gc_operations))

                # Sum the concentrations. If there is one or more NaN values
                # then set to NaN since not enough data for accumulation.
                refsum = df.loc[dfrows_to_sum, "Ref"].sum()
                devsum = df.loc[dfrows_to_sum, "Dev"].sum()

                # Calculate diff and % diff
                if not np.isnan(refsum) and not np.isnan(devsum):
                    diff = devsum - refsum
                    try:
                        pctdiff = diff / refsum * 100
                    except BaseException:
                        pctdiff = np.nan
                else:
                    diff = np.nan
                    pctdiff = np.nan

                # Fill dataframe
                df.loc[dfrow, "Units_converted"] = units[spc]
                df.loc[dfrow, "Ref"] = refsum
                df.loc[dfrow, "Dev"] = devsum
                df.loc[dfrow, "Diff"] = diff
                df.loc[dfrow, "Pct_diff"] = pctdiff

    # ------------------------------------------
    # Compute mass change in restarts for each column section (if applicable)
    # ------------------------------------------
    if compute_restart:
        print('Computing RESTART operation budgets...')

        # Read the species database files in the Ref & Dev rundirs,
        # and return a dict containing metadata for each.
        ref_metadata, dev_metadata = read_species_metadata(
            spcdb_files,
            quiet=True
        )

        # Loop over all column sections
        for col_section in col_sections:

            # Loop over species
            for i, spc in enumerate(spclist):

                # Keep track of progress
                if (i + 1) % 50 == 0:
                    print(f"  {col_section}: species {i + 1} of {n_spc}")

                # Get the accumulation dataframe row to fill. Skip if not found
                # of if not Full column section.
                dfrow = (df["Column_Section"] == "Full") \
                    & (df["Species"] == spc) \
                    & (df["Operation"] == "RESTART")
                if not any(dfrow):
                    continue

                # Get ref and dev mass

                # Get species metadata for unit conversion. If none, skip.
                ref_species_metadata = ref_metadata.get(spc)
                dev_species_metadata = dev_metadata.get(spc)
                if ref_species_metadata is None and \
                   dev_species_metadata is None:
                    continue

                # Get molecular weights
                #ref_mol_wt_g = get_molwt_from_metadata(
                #    ref_species_metadata,
                #    spc
                #)
                #dev_mol_wt_g = get_molwt_from_metadata(
                #    ref_species_metadata,
                #    spc
                #)

                # Specify target units
                target_units = "Gg"

                # ==============================================================
                # Convert units of Ref and save to a DataArray
                # (or skip if Ref contains NaNs everywhere)
                # ==============================================================
                refarray = refdata[v]
                if not np.isnan(refdata[v].values).all():
                    refarray = convert_units(
                        refarray,
                        spc,
                        ref_species_metadata,
                        target_units,
                        area_m2=met_and_masks["Ref_Area"],
                        delta_p=met_and_masks["Ref_Delta_P"],
                        box_height=met_and_masks["Ref_BxHeight"],
                    )

                # ==============================================================
                # Convert units of Dev and save to a DataArray
                # (or skip if Dev contains NaNs everywhere)
                # ==============================================================
                devarray = devdata[v]
                if not np.isnan(devdata[v].values).all():
                    devarray = convert_units(
                        devarray,
                        spc_name,
                        dev_species_metadata,
                        target_units,
                        area_m2=met_and_masks["Dev_Area"],
                        delta_p=met_and_masks["Dev_Delta_P"],
                        box_height=met_and_masks["Dev_BxHeight"],
                    )


                # Compute ref mass as end mass minus start mass in ref
                # TODO

                # Compute dev mass as end mass minus start mass in dev
                # TODO - copy above once compete. The rest should just work.

                # Calculate diff and % diff
                if not np.isnan(refmass) and not np.isnan(devmass):
                    diff = devmass - refmass
                    try:
                        pctdiff = diff / refmass * 100
                    except BaseException:
                        pctdiff = np.nan
                else:
                    diff = np.nan
                    pctdiff = np.nan

                # Fill dataframe
                df.loc[dfrow, "Units_converted"] = units[spc]
                df.loc[dfrow, "Ref"] = refmass
                df.loc[dfrow, "Dev"] = devmass
                df.loc[dfrow, "Diff"] = diff
                df.loc[dfrow, "Pct_diff"] = pctdiff

    #  Sanity check write to csv (for testing. Keep commented out otherwise)
    #df.to_csv('df.csv', na_rep='NA')

    # ------------------------------------------
    # Make budget file
    # ------------------------------------------

    # Create the target output directory hierarchy if it doesn't already exist
    make_directory(dst, overwrite)

    # Print budgets to file
    filename = f"{dst}/Budgets_After_Operations.txt"
    with open(filename, "w+", encoding=ENCODING) as f:
        print("#" * 78, file=f)
        if label is not None and benchmark_type is not None:
            print(f"{benchmark_type} budgets for {label}", file=f)
        else:
            print(f"Budgets across {ref_interval}/{dev_interval} sec", file=f)
        print("\n", file=f)
        print("NOTES:", file=f)
        msg = " - When using the non-local mixing scheme (default), "\
              "'Mixing' includes\n   emissions and dry deposition "\
              "applied below the PBL. 'EmisDryDep'\n   therefore only "\
              "captures fluxes above the PBL.\n - When using full mixing, "\
              "'Mixing' and 'EmisDryDep' are fully separated.\n - Budgets "\
              "are calculated as the sum of [kg/s] tendencies\n - Strat "\
              "budget is calculated as Full minus Trop\n - ACCUMULATION "\
              "is calculated as sum of all other operations"
        print(msg, file=f)
        print("#" * 78, file=f)
        print(file=f)

        # Loop over species
        for i, spc in enumerate(spclist):
            print(f"{spc} budgets (Ref={refstr}; Dev={devstr})", file=f)

            # Print a table for each column section
            for col_section in col_sections:

                # Get the dataframe rows. Skip if none found.
                dfrows = (df["Column_Section"] == col_section) \
                    & (df["Species"] == spc) \
                    & (df["Operation"].isin(all_operations))
                if not any(dfrows):
                    continue

                # Print dataframe subset to file
                print(f"{col_section} {units[spc]} : {spc}", file=f)
                print(tabulate(
                    df.loc[dfrows,
                           ["Operation",
                            "Ref",
                            "Dev",
                            "Diff",
                            "Pct_diff"]],
                    headers='keys',
                    tablefmt='psql',
                    showindex=False,
                    floatfmt=(".5f", ".5f", ".5f", ".5f", ".5f"),
                ), file=f
            )
            print("\n", file=f)

    # ------------------------------------------
    # Clean up
    # ------------------------------------------
    del df
    del ref_ds
    del dev_ds
    gc.collect()


def create_benchmark_summary_table(
        refpath,
        refstr,
        refdate,
        devpath,
        devstr,
        devdate,
        collections,
        dst="./benchmark",
        overwrite=False,
        outfilename="Summary.txt",
        verbose=False,
        ref_gchp=False,
        dev_gchp=False
):
    """
    Creates a benchmark summary table that shows which data collections
    have difference.  Useful for scanning the 1-hr and 1-month benchmark
    outputs.

    Args:
        refpath: str
            Path to the first data set to be compared (aka "Ref").
        refstr: str
            A string that can be used to identify refdata
            (e.g. a model version number or other identifier).
        refdate: np.datetime64
            Date/time stamp used by the "Ref" data files.
        ref_gchp: bool
            Set to True if the "Ref" data comes from a GCHP run.
            Default value: False
        devpath: str
            Path to the second data set to be compared (aka "Dev").
        devstr: str
            A string that can be used to identify the data set specified
            by devfile (e.g. a model version number or other identifier).
        dev_gchp: bool
            Set to True if the "Ref" data comes from a GCHP run.
            Default value: False
        devdate: np.datetime64
            Date/time stamp used by the "Dev" data files.
        collections: list of strings
            List of diagnostic collections to examine.

    Keyword Args (optional):
        dst: str
            A string denoting the destination folder where the file
            containing emissions totals will be written.
            Default value: "./benchmark"
        overwrite: bool
            Set this flag to True to overwrite files in the
            destination folder (specified by the dst argument).
            Default value: False
        outfilename: str
            Name of the text file which will contain the table of
            emissions totals.
            Default value: "Summary.txt"
        verbose: bool
            Set this switch to True if you wish to print out extra
            informational messages.
            Default value: False

    Remarks:
        This method is mainly intended for model benchmarking purposes,
        rather than as a general-purpose tool.

        Species properties (such as molecular weights) are read from a
        YAML file called "species_database.yml".
    """

    # ==================================================================
    # Open file for output
    # ==================================================================

    # Replace whitespace in the ref and dev labels
    refstr = replace_whitespace(refstr)
    devstr = replace_whitespace(devstr)

    # Create the directory for output
    make_directory(dst, overwrite)

    # Create file
    try:
        f = open(os.path.join(dst, outfilename), "w", encoding=ENCODING)
    except (IOError, OSError, FileNotFoundError) as e:
        msg = f"Could not open {outfilename} for writing!"
        raise e(msg) from e

    # Title strings
    title1 = "### Benchmark summary table"
    title2 = f"### Ref = {refstr}"
    title3 = f"### Dev = {devstr}"

    # Print header to file
    print("#" * 80, file=f)
    print(f"{title1 : <77}{'###'}", file=f)
    print(f"{'###'  : <77}{'###'}", file=f)
    print(f"{title2 : <77}{'###'}", file=f)
    print(f"{title3 : <77}{'###'}", file=f)
    print("#" * 80, file=f)
    print(file=f)

    # ==================================================================
    # Read data and look differences btw Ref & Dev versions
    # ==================================================================

    # Variables to skip
    skip_vars = SKIP_THESE_VARS
    skip_vars.append("AREA")

    # Pick the proper function to read the data
    reader = dataset_reader(
        multi_files=False,
        verbose=verbose
    )

    # Loop over diagnostic files
    for col in collections:

        # Read Ref data
        refdata = reader(
            get_filepath(
                refpath,
                col,
                refdate,
                is_gchp=ref_gchp
            ),
            drop_variables=skip_vars
        )

        # Get Dev data
        devdata = reader(
            get_filepath(
                devpath,
                col,
                devdate,
                is_gchp=dev_gchp
            ),
            drop_variables=skip_vars
        )

        # Make sure that Ref and Dev datasets have the same variables.
        # Variables that are in Ref but not in Dev will be added to Dev
        # with all missing values (NaNs). And vice-versa.
        [refdata, devdata] = add_missing_variables(
            refdata,
            devdata
        )

        # Find all common variables between the two datasets
        vardict = compare_varnames(
            refdata,
            devdata,
            quiet=True
        )

        # List of differences for this collection
        diff_list = []

        # Keep track of which variables are different
        # NOTE: Use 32-point float for comparisons since this is
        # the precision used for History diagnostics.
        for v in vardict["commonvarsData"]:
            if not array_equals(
                    refdata[v],
                    devdata[v],
                    dtype=np.float32
            ):
                diff_list.append(v)

        # Drop duplicate values from diff_list
        diff_list = unique_values(diff_list, drop=[None])

        if len(diff_list) == 0:
            print("-" *  79, file=f)
            print(f"{col}: {devstr} is identical to {refstr}", file=f)
            print(file=f)
        else:
            print("-" *  79, file=f)
            print(f"{col}: {devstr} differs from {refstr}", file=f)
            print("\n  Diagnostics that differ", file=f)
            for i, v in enumerate(diff_list):
                print(f"    {v}", file=f)
                if i > 10:
                    print(f"    ... and {len(diff_list) - 10} others", file=f)
                    break
            print(file=f)

    # ==================================================================
    # Close files
    # ==================================================================
    f.close()


def diff_list_to_text(
        refstr,
        devstr,
        diff_list,
        fancy_format=False
):
    """
    Converts a list of species/emissions/inventories/diagnostics that
    show differences between GEOS-Chem versions ot a printable text
    string.

    Args:
    -----
    diff_list : list
        List to be converted into text.  "None" values will be dropped.
    fancy_format: bool
        Set to True if you wish output text to be bookended with '###'.

    Returns:
    diff_text : str
        String with concatenated list values.
    """
    verify_variable_type(diff_list, list)

    # Use "Dev" and "Ref" for inserting into a header
    if fancy_format:
        refstr = "Ref"
        devstr = "Dev"

    # Strip out duplicates from diff_list
    # Prepare a message about species differences (or alternate msg)
    diff_list = unique_values(diff_list, drop=[None])

    # Print the text
    n_diff = len(diff_list)
    if n_diff > 0:
        diff_text = f"{devstr} and {refstr} show {n_diff} differences"
    else:
        diff_text = f"{devstr} and {refstr} are identical"

    # If we are placing the text in a header, trim the length of diff_text
    # to fit.  NOTE: TABLE_WIDTH-7 leaves room for the '### ' at the start
    # of the string and the '###' at the end of the string,
    if fancy_format:
        diff_text = f"### {diff_text : <{TABLE_WIDTH-7}}{'###'}"
        diff_text = wrap_text(
            diff_text,
            width=TABLE_WIDTH
        )

    return diff_text.strip()


def diff_of_diffs_toprow_title(config, model):
    """
    Creates the diff-of-diffs plot title for the top row of the
    six-plot output.  If the title string is too long (as empirically
    determined), then a newline will be inserted in order to prevent
    the title strings from overlapping.

    Args:
    -----
    config : dict
       Dictionary containing the benchmark options (as read from a
       YAML file such as 1mo_benchmark.yml, etc.)
    model: str
       The model to plot.  Accepted values are "gcc" or "gchp".

    Returns:
    --------
    title: str
        The plot title string for the diff-of-diff
    """
    verify_variable_type(config, dict)
    verify_variable_type(model, str)
    if not "gcc" in model and not "gchp" in model:
        msg = "The 'model' argument must be either 'gcc' or 'gchp'!"
        raise ValueError(msg)

    title = (
        config["data"]["dev"][model]["version"]
        + " - "
        + config["data"]["ref"][model]["version"]
    )

    if len(title) > 40:
        title = (
            config["data"]["dev"][model]["version"]
            + " -\n"
            + config["data"]["ref"][model]["version"]
        )

    return title


def create_benchmark_sanity_check_table(
        devpath,
        devstr,
        devdate,
        collections,
        dst="./benchmark",
        is_gchp=False,
        overwrite=False,
        outfilename="Diagnostic_Sanity_Check.txt",
        verbose=False,
):
    """
    Creates a diagnostic sanity check table that shows which diagnostic
    variables are zero or NaN everywhere.  This can help to identify
    bugs in diagnostic output.

    Args:
        devpath: str
            Path to the data set to be compared (aka "Dev").
        devstr: str
            A string that can be used to identify the data set specified
            by devfile (e.g. a model version number or other identifier).
        devdate: np.datetime64
            Date/time stamp used by the "Dev" data files.
        collections: list of strings
            List of diagnostic collections to examine.

    Keyword Args (optional):
        dst: str
            A string denoting the destination folder where the file
            containing emissions totals will be written.
            Default value: "./benchmark"
        is_gchp : bool
           Set this flag to true to denote if the data is from GCHP.
        overwrite: bool
            Set this flag to True to overwrite files in the
            destination folder (specified by the dst argument).
            Default value: False
        outfilename: str
            Name of the text file which will contain the table of
            emissions totals.
            Default value: "Summary.txt"
        verbose: bool
            Set this switch to True if you wish to print out extra
            informational messages.
            Default value: False

    Remarks:
        This method is mainly intended for model benchmarking purposes,
        rather than as a general-purpose tool.
    """

    # ==================================================================
    # Initial preparations
    # ==================================================================

    # Replace whitespace in the ref and dev labels
    devstr = replace_whitespace(devstr)

    # Create the directory for output (if necessary)
    make_directory(dst, overwrite)
    outfilename = os.path.join(dst, outfilename)

    # Pick the proper function to read the data
    reader = dataset_reader(
        multi_files=False,
        verbose=verbose
    )

    # Variables to skip
    skip_vars = SKIP_THESE_VARS
    skip_vars.append("AREA")

    # ==================================================================
    # Open output file and write header
    # ==================================================================
    with open(outfilename, "w", encoding=ENCODING) as ofile:

        # Title strings
        title1 = "### Benchmark diagnostic sanity check table"
        title2 = f"### Dev = {devstr}"

        # Print header to file
        print("#" * 80, file=ofile)
        print(f"{title1 : <77}{'###'}", file=ofile)
        print(f"{'###'  : <77}{'###'}", file=ofile)
        print(f"{title2 : <77}{'###'}", file=ofile)
        print("#" * 80, file=ofile)

        # ==============================================================
        # Loop over diagnostic collections and scan files
        # ==============================================================
        for col in collections:

            # Read data into an xr.DataSet object
            file_name = get_filepath(
                devpath,
                col,
                devdate,
                is_gchp=is_gchp,
            )
            dset = reader(
                file_name,
                drop_variables=skip_vars
            )

            # Determine which variables are all zeroes or NaN
            all_zeros_or_nans = []
            for var in dset.data_vars:
                data = dset[var].values
                if np.all(data == 0) or np.all(data == np.nan):
                    all_zeros_or_nans.append(var)

            # ===========================================================
            # Print results for each collection
            # ===========================================================
            print("", file=ofile)
            print("="*80, file=ofile)
            print(f"{os.path.basename(file_name)}", file=ofile)
            print("="*80, file=ofile)
            print("", file=ofile)

            if len(all_zeros_or_nans) == 0:
                print("No variables were all zero or all NaN", file=ofile)
            else:
                print("These variables were all zero or all NaN:", file=ofile)
                for var in all_zeros_or_nans:
                    print(f"   {var}", file=ofile)
