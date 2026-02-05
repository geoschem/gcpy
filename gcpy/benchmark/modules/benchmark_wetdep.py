#!/usr/bin/env python3
"""
Specific utilities for creating plots from GEOS-Chem benchmark simulations.
"""
import os
import gc
import numpy as np
from gcpy.util import \
    add_bookmarks_to_pdf, compare_varnames, dataset_reader, \
    dataset_mean, make_directory, replace_whitespace
from gcpy.constants import SKIP_THESE_VARS
from gcpy.plot.compare_single_level import compare_single_level
from gcpy.plot.compare_zonal_mean import compare_zonal_mean

# Suppress numpy divide by zero warnings to prevent output spam
np.seterr(divide="ignore", invalid="ignore")


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
