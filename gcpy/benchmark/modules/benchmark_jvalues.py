#!/usr/bin/env python3
"""
Creates J-Value plots from GEOS-Chem benchmark simulation output.
"""
import os
import gc
import numpy as np
from gcpy.regrid import create_regridders
from gcpy.util import \
    add_bookmarks_to_pdf, add_missing_variables, compare_varnames, \
    dataset_reader, dataset_mean, divide_dataset_by_dataarray, \
    make_directory, replace_whitespace
from gcpy.constants import \
    ENCODING, SKIP_THESE_VARS
from gcpy.plot.compare_single_level import compare_single_level
from gcpy.plot.compare_zonal_mean import compare_zonal_mean

# Suppress numpy divide by zero warnings to prevent output spam
np.seterr(divide="ignore", invalid="ignore")


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
