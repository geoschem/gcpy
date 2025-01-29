"""
Specific utilities for creating plots from GEOS-Chem benchmark simulations.
"""
import gc
import numpy as np
from gcpy import util
from gcpy.plot.compare_single_level import compare_single_level
from gcpy.benchmark.modules.benchmark_utils import \
    get_common_varnames, make_output_dir, pdf_filename, \
    print_sigdiffs, read_ref_and_dev

# Suppress numpy divide by zero warnings to prevent output spam
np.seterr(divide="ignore", invalid="ignore")


def make_benchmark_drydep_plots(
        ref,
        refstr,
        dev,
        devstr,
        collection="DryDep",
        dst="./benchmark",
        subdst=None,
        cmpres=None,
        overwrite=False,
        verbose=False,
        log_color_scale=False,
        weightsdir=".",
        sigdiff_files=None,
        n_job=-1,
        time_mean=False,
        varlist=None,
        spcdb_dir=None,
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

    Keyword Args (optional):
        collection : str
            Name of the diagnostic collection (e.g. "DryDep")
        dst: str
            A string denoting the destination folder where a PDF
            file containing plots will be written.
            Default value: ./benchmark
        subdst: str
            A string denoting the sub-directory of dst where PDF
            files containing plots will be written.  In practice,
            subdst is only needed for the 1-year benchmark output,
            and denotes a date string (such as "Jan2016") that
            corresponds to the month that is being plotted.
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
        n_job: int
            Defines the number of simultaneous workers for parallel plotting.
            Set to 1 to disable parallel plotting. Value of -1 allows the
            application to decide.
            Default value: -1
        spcdb_dir: str
            Directory of species_datbase.yml file
            Default value: None
        time_mean : bool
            Determines if we should average the datasets over time
            Default value: False
        varlist: list of str
            List of variables to plot.  If varlist is None, then
            all common variables in Ref & Dev will be plotted.
    """
    # Make sure the species database folder is passed
    if spcdb_dir is None:
        msg = "The spcdb_dir argument has not been specified!"
        raise ValueError(msg)

    # Replace whitespace in the ref and dev labels
    refstr = util.replace_whitespace(refstr)
    devstr = util.replace_whitespace(devstr)

    # Create directory for plots (if it doesn't exist)
    dst = make_output_dir(
        dst,
        collection,
        subdst,
        overwrite=overwrite,
    )

    # Read data
    refdata, devdata = read_ref_and_dev(
        ref,
        dev,
        time_mean=time_mean
    )

    # Get common variables between Ref and Dev
    if varlist is None:
        varlist = get_common_varnames(
            refdata,
            devdata,
            prefix="DryDepVel_",
            verbose=verbose
        )

    # Create surface plots
    sigdiff_list = []
    pdfname = pdf_filename(
        dst,
        collection,
        subdst,
        plot_type="Surface"
    )
    compare_single_level(
        refdata,
        refstr,
        devdata,
        devstr,
        varlist=varlist,
        cmpres=cmpres,
        ilev=0,
        pdfname=pdfname,
        log_color_scale=log_color_scale,
        extra_title_txt=subdst,
        sigdiff_list=sigdiff_list,
        weightsdir=weightsdir,
        n_job=n_job,
        spcdb_dir=spcdb_dir
    )
    util.add_bookmarks_to_pdf(
        pdfname,
        varlist,
        remove_prefix=collection + '_',
        verbose=verbose
    )

    # Write significant differences to file (if there are any)
    print_sigdiffs(
        sigdiff_files,
        sigdiff_list,
        sigdiff_type="sfc",
        sigdiff_cat="DryDepVel"
    )

    # -------------------------------------------
    # Clean up
    # -------------------------------------------
    del refdata
    del devdata
    gc.collect()


def drydepvel_species():
    """
    Returns a list of species for the dry deposition velocity
    (DryDepVel) benchmark plots:

    Returns:
    --------
    varnames (list of str): Variable names to plot
    """
    # These are key dry deposition species (as per Mat Evans)
    return ["DryDepVel_ACET", "DryDepVel_HNO3", "DryDepVel_NH3",
            "DryDepVel_NH4", "DryDepVel_NIT", "DryDepVel_NITs",
            "DryDepVel_O3", "DryDepVel_SO4"]
