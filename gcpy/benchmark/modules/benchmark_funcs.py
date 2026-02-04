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
from tabulate import tabulate
from gcpy.regrid import create_regridders
from gcpy.grid import get_troposphere_mask
from gcpy.util import \
    add_bookmarks_to_pdf, add_missing_variables, \
    array_equals, compare_varnames, dataset_reader, \
    dataset_mean, get_area_from_dataset, \
    get_filepath, get_variables_from_dataset, \
    insert_text_into_file, make_directory, print_totals, read_config_file, \
    read_species_metadata, rename_and_flip_gchp_rst_vars, \
    replace_whitespace, unique_values, verify_variable_type, wrap_text
from gcpy.units import convert_units
from gcpy.constants import \
    COL_WIDTH, ENCODING, MW_AIR_g, SKIP_THESE_VARS, TABLE_WIDTH
from gcpy.plot.compare_single_level import compare_single_level
from gcpy.plot.compare_zonal_mean import compare_zonal_mean
from gcpy.benchmark.modules.benchmark_utils import \
    AOD_SPC, rename_speciesconc_to_speciesconcvv

# Suppress numpy divide by zero warnings to prevent output spam
np.seterr(divide="ignore", invalid="ignore")


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
