#!/usr/bin/env python3
"""
Generates plots from specified collections in GEOS-Chem benchmark output.
"""
import os
import gc
import numpy as np
from gcpy.regrid import create_regridders
from gcpy.util import \
    add_bookmarks_to_pdf, add_missing_variables, compare_varnames, \
    dataset_reader, dataset_mean, make_directory, replace_whitespace
from gcpy.constants import SKIP_THESE_VARS
from gcpy.plot.compare_single_level import compare_single_level
from gcpy.plot.compare_zonal_mean import compare_zonal_mean

# Suppress numpy divide by zero warnings to prevent output spam
np.seterr(divide="ignore", invalid="ignore")


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
