"""
Specific utilities for creating plots from GEOS-Chem benchmark simulations.
"""
import os
import gc
import numpy as np
from gcpy import util
from gcpy.constants import skip_these_vars
from gcpy.plot.compare_single_level import compare_single_level
from gcpy.plot.compare_zonal_mean import compare_zonal_mean
import gcpy.benchmark.modules.benchmark_utils as bmk_util

# Suppress numpy divide by zero warnings to prevent output spam
np.seterr(divide="ignore", invalid="ignore")


def make_benchmark_drydep_plots(
        ref,
        refstr,
        dev,
        devstr,
        collection,
        dst="./benchmark",
        cmpres=None,
        datestr=None,
        overwrite=False,
        verbose=False,
        log_color_scale=False,
        weightsdir=".".
        n_job=-1,
        time_mean=False,
        spcdb_dir=None,
        var_prefix="DryDepVel"
):
    """
    Creates six-panel comparison plots (PDF format) from GEOS-Chem
    benchmark simualtion output.  Can be used with data collections
    that do not require special handling (e.g. concentrations).

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
        n_job: int
            Defines the number of simultaneous workers for parallel plotting.
            Set to 1 to disable parallel plotting. Value of -1 allows the
            application to decide.
            Default value: -1
        spcdb_dir: str
            Directory of species_datbase.yml file
            Default value: Directory of GCPy code repository
        time_mean : bool
            Determines if we should average the datasets over time
            Default value: False
        var_prefix : str
            If 
    """

    # Create destination plot directory
    bmk_util.make_directory(dst, overwrite)
    target_dst = bmk_util.make_collection_subdir(dst, collection, datestr)

    # Defaults for arguments
    if plots is None:
        plots = ["sfc", "500hpa", "zonalmean"]
    if spcdb_dir is None:
        spcdb_dir = os.path.join(os.path.dirname(__file__), "..", "..")

    # Read data
    refds, devds = bmk.util_read_ref_and_dev(ref, dev, time_mean)
    if refmet is not None and devmet is not None:
        refmetds, devmetds = bmk_util.read_ref_and_dev(ref, dev, time_mean)

    # Make sure that Ref and Dev datasets have the same variables.
    # Variables that are in Ref but not in Dev will be added to Dev
    # with all missing values (NaNs). And vice-versa.
    # Turn this off for now since add_missing_variables inserts GCC area into
    # GCHP files, which causes problems with area normalization (ewl)
    #[refds, devds] = add_missing_variables(refds, devds)

    # Get common variables between Ref and Dev
    varlist = bmk_util.get_common_varnames(
        refds, 
        devds, 
        var_prefix=var_prefix
        verbose=verbose
    )
        
    # Surface plots
    plotfilename = f"{collection}_Surface.pdf"
    if datestr is not None:
        plotfilename = f"{collection}_Surface_{datestr}.pdf"
    pdfname = os.path.join(target_dst, plotfilename)
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
        spcdb_dir=spcdb_dir
    )
    util.add_bookmarks_to_pdf(
        pdfname,
        varlist,
        remove_prefix=collection + '_',
        verbose=verbose)

    )


    # -------------------------------------------
    # Clean up
    # -------------------------------------------
    del refds
    del devds
    del refmetds
    del devmetds
    gc.collect()
