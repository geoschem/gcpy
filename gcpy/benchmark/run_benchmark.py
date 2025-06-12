#!/usr/bin/env python
"""
run_benchmark.py: Driver script for creating benchmark plots.

Run this script to generate benchmark comparisons between:

    (1) GCC (aka GEOS-Chem "Classic") vs. GCC
    (2) GCHP vs GCC
    (3) GCHP vs GCHP
    (4) GCHP vs GCC diff-of-diffs

You can customize this by editing the settings in the corresponding yaml
configuration file.

Calling sequence:

    $ python -m gcpy.benchmark.run_benchmark <path-to-configuration-file>

Remarks:

    By default, matplotlib will try to open an X window for plotting.
    If you are running this script in an environment where you do not have
    an active X display (such as in a computational queue), then you will
    need to use these commands to disable the X-window functionality.

        import os
        os.environ["QT_QPA_PLATFORM"]="offscreen"

    For more information, please see this issue posted at the ipython site:

        https://github.com/ipython/ipython/issues/10627

    Also, to disable matplotlib from trying to open X windows, you may
    need to set the following environment variable in your shell:

        $ export MPLBACKEND=agg

This script corresponds with GCPy 1.6.2. Edit this version ID if releasing
a new version of GCPy.
"""

# ======================================================================
# Imports and global settings (you should not need to edit these)
# ======================================================================

import os
import sys
import warnings
from datetime import datetime
import numpy as np
from gcpy.util import copy_file_to_dir, get_filepath, read_config_file
from gcpy.date_time import add_months, is_full_year
from gcpy.benchmark.modules.benchmark_funcs import \
    diff_of_diffs_toprow_title, get_species_database_dir, \
    create_benchmark_summary_table, \
    create_benchmark_sanity_check_table, \
    make_benchmark_conc_plots, make_benchmark_emis_plots, \
    make_benchmark_emis_tables, make_benchmark_jvalue_plots, \
    make_benchmark_aod_plots, make_benchmark_mass_tables, \
    make_benchmark_mass_accumulation_tables, \
    make_benchmark_operations_budget, \
    make_benchmark_collection_2d_var_plots, \
    make_benchmark_collection_3d_var_plots, \
    make_benchmark_operations_budget
from gcpy.benchmark.modules.ste_flux import make_benchmark_ste_table
from gcpy.benchmark.modules.oh_metrics import make_benchmark_oh_metrics
from gcpy.benchmark.modules.run_1yr_fullchem_benchmark \
    import run_benchmark as run_1yr_benchmark
from gcpy.benchmark.modules.run_1yr_tt_benchmark \
    import run_benchmark as run_1yr_tt_benchmark
from gcpy.benchmark.modules.benchmark_utils import \
    gcc_vs_gcc_dirs, gchp_vs_gcc_dirs, gchp_vs_gchp_dirs, \
    get_log_filepaths, print_benchmark_info
from gcpy.benchmark.modules.benchmark_drydep \
    import drydepvel_species, make_benchmark_drydep_plots
from gcpy.benchmark.modules.benchmark_scrape_gcclassic_timers import \
    make_benchmark_gcclassic_timing_table
from gcpy.benchmark.modules.benchmark_scrape_gchp_timers import \
    make_benchmark_gchp_timing_table

# Tell matplotlib not to look for an X-window
os.environ["QT_QPA_PLATFORM"] = "offscreen"

# Suppress harmless run-time warnings (mostly about underflow in division)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)


def choose_benchmark_type(config):
    """
    Decides which benchmark to run (default, 1yr, or 1yr_tt)

    Args:
        config : dict
            Contains configuration for 1mon benchmark from yaml file.
    """
    if not (
        config["options"]["bmk_type"] == "FullChemBenchmark"
        or config["options"]["bmk_type"] == "TransportTracersBenchmark"
        or config["options"]["bmk_type"] == "CH4Benchmark"
    ):
        print(
            f"Error: invalid benchmark type {config['options']['bmk_type']}. "
            + "Please enter FullChemBenchmark, TransportTracersBenchmark, or "
            + "CH4Benchmark."
        )
        sys.exit()

    start = np.datetime64(config["data"]["ref"]["gcc"]["bmk_start"])
    end = np.datetime64(config["data"]["ref"]["gcc"]["bmk_end"])
    # determine benchmark type and run relevant script
    if is_full_year(start, end):
        if (
            config["options"]["bmk_type"] == "FullChemBenchmark"
            or config["options"]["bmk_type"] == "CH4Benchmark"
        ):
            run_1yr_benchmark(
                config,
                str(start.astype(datetime).year),
                str(start.astype(datetime).year)
            )
        else:
            run_1yr_tt_benchmark(
                config,
                str(start.astype(datetime).year),
                str(start.astype(datetime).year)
            )
    else:
        run_benchmark_default(config)


def run_benchmark_default(config):
    """
    Runs flexible date benchmark with the given configuration settings.

    Args:
        config : dict
            Contains configuration for 1mon benchmark from yaml file.
    """
    # Folder in which the species_database.yml file is found
    spcdb_dir = get_species_database_dir(config)

    # =====================================================================
    # Data directories
    # For gchp_vs_gcc_refdir use config["data"]["dev"]["gcc"]["version"],
    # not ref (mps, 6/27/19)
    # =====================================================================

    # Diagnostics file directory paths
    s = "outputs_subdir"
    gcc_vs_gcc_refdir, gcc_vs_gcc_devdir = gcc_vs_gcc_dirs(config, s)
    gchp_vs_gcc_refdir, gchp_vs_gcc_devdir = gchp_vs_gcc_dirs(config, s)
    gchp_vs_gchp_refdir, gchp_vs_gchp_devdir = gchp_vs_gchp_dirs(config, s)

    # Restart file directory paths
    s = "restarts_subdir"
    gcc_vs_gcc_refrstdir, gcc_vs_gcc_devrstdir = gcc_vs_gcc_dirs(config, s)
    gchp_vs_gcc_refrstdir, gchp_vs_gcc_devrstdir = gchp_vs_gcc_dirs(config, s)
    gchp_vs_gchp_refrstdir, gchp_vs_gchp_devrstdir = gchp_vs_gchp_dirs(config, s)

    # Log file directory paths
    s = "logs_subdir"
    gcc_vs_gcc_reflogdir, gcc_vs_gcc_devlogdir = gcc_vs_gcc_dirs(config, s)
    gchp_vs_gcc_reflogdir, gchp_vs_gcc_devlogdir = gchp_vs_gcc_dirs(config, s)
    gchp_vs_gchp_reflogdir, gchp_vs_gchp_devlogdir = gchp_vs_gchp_dirs(config, s)

    # =====================================================================
    # Benchmark output directories
    # =====================================================================
    # Results directories
    gcc_vs_gcc_resultsdir = os.path.join(
        config["paths"]["main_dir"],
        config["data"]["dev"]["gcc"]["dir"],
        config["paths"]["results_dir"],
    )
    gchp_vs_gchp_resultsdir = os.path.join(
        config["paths"]["main_dir"],
        config["data"]["dev"]["gchp"]["dir"],
        config["paths"]["results_dir"],
        config["options"]["comparisons"]["gchp_vs_gchp"]["dir"],
    )
    gchp_vs_gcc_resultsdir = os.path.join(
        config["paths"]["main_dir"],
        config["data"]["dev"]["gchp"]["dir"],
        config["paths"]["results_dir"],
        config["options"]["comparisons"]["gchp_vs_gcc"]["dir"],
    )
    diff_of_diffs_resultsdir = os.path.join(
        config["paths"]["main_dir"],
        config["data"]["dev"]["gchp"]["dir"],
        config["paths"]["results_dir"],
        "GCHP_GCC_diff_of_diffs",
    )
    base_gchp_resultsdir = os.path.join(
        config["paths"]["main_dir"],
        config["data"]["dev"]["gchp"]["dir"],
        config["paths"]["results_dir"],
    )

    # make results directories that don't exist
    for (resdir, plotting_type) in zip(
        [
            gcc_vs_gcc_resultsdir,
            base_gchp_resultsdir,
            gchp_vs_gchp_resultsdir,
            gchp_vs_gcc_resultsdir,
            diff_of_diffs_resultsdir,
        ],
        [
            config["options"]["comparisons"]["gcc_vs_gcc"]["run"],
            config["options"]["comparisons"]["gchp_vs_gcc"]["run"]
            or config["options"]["comparisons"]["gchp_vs_gchp"]["run"]
            or config["options"]["comparisons"]["gchp_vs_gcc_diff_of_diffs"]["run"],
            config["options"]["comparisons"]["gchp_vs_gchp"]["run"],
            config["options"]["comparisons"]["gchp_vs_gcc"]["run"],
            config["options"]["comparisons"]["gchp_vs_gcc_diff_of_diffs"]["run"],
        ],
    ):
        if plotting_type and not os.path.exists(resdir):
            os.mkdir(resdir)

            # Copy this script and the config file to each results dir
            if resdir in [
                    gcc_vs_gcc_resultsdir,
                    gchp_vs_gcc_resultsdir,
                    gchp_vs_gchp_resultsdir,
                    diff_of_diffs_resultsdir,
            ]:
                copy_file_to_dir(__file__, resdir)
                copy_file_to_dir(config["configuration_file_name"], resdir)

    gcc_vs_gcc_tablesdir = os.path.join(
        gcc_vs_gcc_resultsdir,
        config["options"]["comparisons"]["gcc_vs_gcc"]["tables_subdir"],
    )
    gchp_vs_gchp_tablesdir = os.path.join(
        gchp_vs_gchp_resultsdir,
        config["options"]["comparisons"]["gchp_vs_gchp"]["tables_subdir"],
    )
    gchp_vs_gcc_tablesdir = os.path.join(
        gchp_vs_gcc_resultsdir,
        config["options"]["comparisons"]["gchp_vs_gcc"]["tables_subdir"],
    )

    # =====================================================================
    # Plot title strings
    # For gchp_vs_gcc_refstr use config["data"]["dev"]["gcc"]["version"], not ref (mps, 6/27/19)
    # =====================================================================
    gcc_vs_gcc_refstr = config["data"]["ref"]["gcc"]["version"]
    gcc_vs_gcc_devstr = config["data"]["dev"]["gcc"]["version"]
    gchp_vs_gcc_refstr = config["data"]["dev"]["gcc"]["version"]
    gchp_vs_gcc_devstr = config["data"]["dev"]["gchp"]["version"]
    gchp_vs_gchp_refstr = config["data"]["ref"]["gchp"]["version"]
    gchp_vs_gchp_devstr = config["data"]["dev"]["gchp"]["version"]
    diff_of_diffs_refstr = diff_of_diffs_toprow_title(config, "gcc")
    diff_of_diffs_devstr = diff_of_diffs_toprow_title(config, "gchp")

    ########################################################################
    ###    THE REST OF THESE SETTINGS SHOULD NOT NEED TO BE CHANGED      ###
    ########################################################################

    # =====================================================================
    # Dates and times
    # =====================================================================

    # Ref start used in diagnostic filename
    gchp_ref_date = np.datetime64(config["data"]["ref"]["gcc"]["bmk_start"])
    gcc_ref_date = np.datetime64(config["data"]["ref"]["gcc"]["bmk_start"])
    # Ref end used in restart filename)
    gcc_end_ref_date = np.datetime64(config["data"]["ref"]["gcc"]["bmk_end"])
    gchp_end_ref_date = np.datetime64(config["data"]["ref"]["gchp"]["bmk_end"])

    # Dev start used in diagnostic filename
    gcc_dev_date = np.datetime64(config["data"]["dev"]["gcc"]["bmk_start"])
    gchp_dev_date = np.datetime64(config["data"]["dev"]["gchp"]["bmk_start"])
    # Dev end used in restart filename
    gcc_end_dev_date = np.datetime64(config["data"]["dev"]["gcc"]["bmk_end"])
    gchp_end_dev_date = np.datetime64(config["data"]["dev"]["gchp"]["bmk_end"])

    # Seconds per month
    gcc_ref_sec_diff = (gcc_end_ref_date - gcc_ref_date).astype("float64")
    gchp_ref_sec_diff = (gchp_end_ref_date - gchp_ref_date).astype("float64")
    gcc_dev_sec_diff = (gcc_end_dev_date - gcc_dev_date).astype("float64")
    gchp_dev_sec_diff = (gchp_end_dev_date - gchp_dev_date).astype("float64")

    # ======================================================================
    # Significant difference filenames
    # ======================================================================

    gcc_vs_gcc_sigdiff = [
        os.path.join(gcc_vs_gcc_resultsdir, "SigDiffs_sfc.txt"),
        os.path.join(gcc_vs_gcc_resultsdir, "SigDiffs_500hpa.txt"),
        os.path.join(gcc_vs_gcc_resultsdir, "SigDiffs_zonalmean.txt"),
        os.path.join(gcc_vs_gcc_resultsdir, "SigDiffs_emissions.txt"),
    ]
    gchp_vs_gcc_sigdiff = [
        os.path.join(gchp_vs_gcc_resultsdir, "SigDiffs_sfc.txt"),
        os.path.join(gchp_vs_gcc_resultsdir, "SigDiffs_500hpa.txt"),
        os.path.join(gchp_vs_gcc_resultsdir, "SigDiffs_zonalmean.txt"),
        os.path.join(gchp_vs_gcc_resultsdir, "SigDiffs_emissions.txt"),
    ]
    gchp_vs_gchp_sigdiff = [
        os.path.join(gchp_vs_gchp_resultsdir, "SigDiffs_sfc.txt"),
        os.path.join(gchp_vs_gchp_resultsdir, "SigDiffs_500hpa.txt"),
        os.path.join(gchp_vs_gchp_resultsdir, "SigDiffs_zonalmean.txt"),
        os.path.join(gchp_vs_gchp_resultsdir, "SigDiffs_emissions.txt"),
    ]

    # ======================================================================
    # Print the list of plots & tables to the screen
    # ======================================================================
    print_benchmark_info(config)

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Create GCC vs GCC benchmark plots and tables
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if config["options"]["comparisons"]["gcc_vs_gcc"]["run"]:

        # ==================================================================
        # GCC vs GCC error checks
        # ==================================================================
        if not np.equal(gcc_ref_sec_diff, gcc_dev_sec_diff):
            print("Skipping GCC vs GCC emissions tables and operations")
            print("budget tables because months are different lengths")
            config["options"]["outputs"]["emis_table"] = False
            config["options"]["outputs"]["ops_budget_table"] = False

        # ==================================================================
        # GCC vs GCC string for month and year
        # (e.g. "2019-07-01T00:00:00 - 2019-08-01T00:00:00")
        # ==================================================================
        if np.equal(gcc_ref_date, gcc_dev_date) and np.equal(
            gcc_end_ref_date, gcc_end_dev_date
        ):
            comparison_str = (
                f"{config['data']['dev']['gcc']['bmk_start']} "
                + f"- {config['data']['dev']['gcc']['bmk_end']}"
            )
        else:
            comparison_str = (
                f"{config['data']['dev']['gcc']['bmk_start']} "
                + f"- {config['data']['dev']['gcc']['bmk_end']}"
                + f" Vs {config['data']['ref']['gcc']['bmk_start']} "
                + f"- {config['data']['ref']['gcc']['bmk_end']}"
            )

        # ==================================================================
        # GCC vs GCC filepaths for StateMet collection data
        # ==================================================================
        refmet = get_filepath(gcc_vs_gcc_refdir, "StateMet", gcc_ref_date)
        devmet = get_filepath(gcc_vs_gcc_devdir, "StateMet", gcc_dev_date)

        # ==================================================================
        # GCC vs GCC species concentration plots
        #
        # Includes lumped species and separates by category if plot_by_spc_cat
        # is true; otherwise excludes lumped species and writes to one file
        # ==================================================================
        if config["options"]["outputs"]["plot_conc"]:
            title = "\n%%% Creating GCC vs. GCC concentration plots %%%"

            # Diagnostic collection files to read
            ref = get_filepath(gcc_vs_gcc_refdir, "SpeciesConc", gcc_ref_date)
            dev = get_filepath(gcc_vs_gcc_devdir, "SpeciesConc", gcc_dev_date)

            # Create plots
            make_benchmark_conc_plots(
                ref,
                gcc_vs_gcc_refstr,
                dev,
                gcc_vs_gcc_devstr,
                refmet=refmet,
                devmet=devmet,
                dst=gcc_vs_gcc_resultsdir,
                weightsdir=config["paths"]["weights_dir"],
                plot_by_spc_cat=config["options"]["outputs"]["plot_options"][
                    "by_spc_cat"],
                benchmark_type=config["options"]["bmk_type"],
                overwrite=True,
                sigdiff_files=gcc_vs_gcc_sigdiff,
                spcdb_dir=spcdb_dir,
                n_job=config["options"]["n_cores"]
            )

        # ==================================================================
        # GCC vs. GCC emissions plots
        # ==================================================================
        if config["options"]["outputs"]["plot_emis"]:
            print("\n%%% Creating GCC vs. GCC emissions plots %%%")

            # Filepaths
            ref = get_filepath(gcc_vs_gcc_refdir, "Emissions", gcc_ref_date)
            dev = get_filepath(gcc_vs_gcc_devdir, "Emissions", gcc_dev_date)

            # Create emissions plots
            make_benchmark_emis_plots(
                ref,
                gcc_vs_gcc_refstr,
                dev,
                gcc_vs_gcc_devstr,
                dst=gcc_vs_gcc_resultsdir,
                weightsdir=config["paths"]["weights_dir"],
                plot_by_spc_cat=config["options"]["outputs"]["plot_options"][
                    "by_spc_cat"],
                plot_by_hco_cat=config["options"]["outputs"]["plot_options"][
                    "by_hco_cat"],
                benchmark_type=config["options"]["bmk_type"],
                overwrite=True,
                sigdiff_files=gcc_vs_gcc_sigdiff,
                spcdb_dir=spcdb_dir,
                n_job=config["options"]["n_cores"]
            )

        # ==================================================================
        # GCC vs. GCC tables of emission and inventory totals
        # ==================================================================
        if config["options"]["outputs"]["emis_table"]:
            print("\n%%% Creating GCC vs. GCC emissions & inventory tables %%%")

            # Filepaths
            ref = get_filepath(gcc_vs_gcc_refdir, "Emissions", gcc_ref_date)
            dev = get_filepath(gcc_vs_gcc_devdir, "Emissions", gcc_dev_date)

            # Create tables
            make_benchmark_emis_tables(
                ref,
                gcc_vs_gcc_refstr,
                dev,
                gcc_vs_gcc_devstr,
                dst=gcc_vs_gcc_resultsdir,
                ref_interval=[gcc_ref_sec_diff],
                dev_interval=[gcc_dev_sec_diff],
                benchmark_type=config["options"]["bmk_type"],
                overwrite=True,
                spcdb_dir=spcdb_dir,
            )

        # ==================================================================
        # GCC vs GCC J-value plots
        # ==================================================================
        if config["options"]["outputs"]["plot_jvalues"]:
            print("\n%%% Creating GCC vs. GCC J-value plots %%%")

            # Diagnostic collection files to read
            ref = get_filepath(gcc_vs_gcc_refdir, "JValues", gcc_ref_date)
            dev = get_filepath(gcc_vs_gcc_devdir, "JValues", gcc_dev_date)

            # Create plots
            make_benchmark_jvalue_plots(
                ref,
                gcc_vs_gcc_refstr,
                dev,
                gcc_vs_gcc_devstr,
                dst=gcc_vs_gcc_resultsdir,
                weightsdir=config["paths"]["weights_dir"],
                overwrite=True,
                sigdiff_files=gcc_vs_gcc_sigdiff,
                spcdb_dir=spcdb_dir,
                n_job=config["options"]["n_cores"]
            )

        # ==================================================================
        # GCC vs GCC budget plots
        # ==================================================================
        if config["options"]["outputs"]["plot_budget"]:
            print("\n%%% Creating GCC vs. GCC Budget plots %%%")

            colname = "Budget"

            # Diagnostic collection files to read
            ref = get_filepath(gcc_vs_gcc_refdir, colname, gcc_ref_date)
            dev = get_filepath(gcc_vs_gcc_devdir, colname, gcc_dev_date)

            # Create plots
            make_benchmark_collection_2d_var_plots(
                ref,
                gcc_vs_gcc_refstr,
                dev,
                gcc_vs_gcc_devstr,
                colname=colname,
                var_prefix=colname,
                dst=gcc_vs_gcc_resultsdir,
                weightsdir=config["paths"]["weights_dir"],
                overwrite=True,
                spcdb_dir=spcdb_dir,
                n_job=config["options"]["n_cores"]
            )

        # ==================================================================
        # GCC vs GCC UVFlux plots
        # ==================================================================
        if config["options"]["outputs"]["plot_uvflux"]:
            print("\n%%% Creating GCC vs. GCC UVFlux plots %%%")

            colname = "UVFlux"

            # Diagnostic collection files to read
            ref = get_filepath(gcc_vs_gcc_refdir, colname, gcc_ref_date)
            dev = get_filepath(gcc_vs_gcc_devdir, colname, gcc_dev_date)

            # Create plots
            make_benchmark_collection_3d_var_plots(
                ref,
                gcc_vs_gcc_refstr,
                dev,
                gcc_vs_gcc_devstr,
                colname=colname,
                dst=gcc_vs_gcc_resultsdir,
                weightsdir=config["paths"]["weights_dir"],
                overwrite=True,
                spcdb_dir=spcdb_dir,
                n_job=config["options"]["n_cores"]
            )

        # ==================================================================
        # GCC vs GCC 2D StateMet plots
        # ==================================================================
        if config["options"]["outputs"]["plot_2d_met"]:
            print("\n%%% Creating GCC vs. GCC 2D StateMet plots %%%")

            # Create plots
            make_benchmark_collection_2d_var_plots(
                refmet,
                gcc_vs_gcc_refstr,
                devmet,
                gcc_vs_gcc_devstr,
                colname="StateMet",
                var_prefix='Met',
                dst=gcc_vs_gcc_resultsdir,
                weightsdir=config["paths"]["weights_dir"],
                overwrite=True,
                spcdb_dir=spcdb_dir,
                n_job=config["options"]["n_cores"]
            )

        # ==================================================================
        # GCC vs GCC 3D StateMet plots
        # ==================================================================
        if config["options"]["outputs"]["plot_3d_met"]:
            print("\n%%% Creating GCC vs. GCC 3D StateMet plots %%%")

            # Create plots
            make_benchmark_collection_3d_var_plots(
                refmet,
                gcc_vs_gcc_refstr,
                devmet,
                gcc_vs_gcc_devstr,
                colname="StateMet",
                var_prefix='Met',
                dst=gcc_vs_gcc_resultsdir,
                weightsdir=config["paths"]["weights_dir"],
                overwrite=True,
                spcdb_dir=spcdb_dir,
                n_job=config["options"]["n_cores"]
            )

        # ==================================================================
        # GCC vs GCC column AOD plots
        # ==================================================================
        if config["options"]["outputs"]["plot_aod"]:
            print("\n%%% Creating GCC vs. GCC column AOD plots %%%")

            # Filepaths
            ref = get_filepath(gcc_vs_gcc_refdir, "Aerosols", gcc_ref_date)
            dev = get_filepath(gcc_vs_gcc_devdir, "Aerosols", gcc_dev_date)

            # Create plots
            make_benchmark_aod_plots(
                ref,
                gcc_vs_gcc_refstr,
                dev,
                gcc_vs_gcc_devstr,
                dst=gcc_vs_gcc_resultsdir,
                weightsdir=config["paths"]["weights_dir"],
                overwrite=True,
                sigdiff_files=gcc_vs_gcc_sigdiff,
                spcdb_dir=spcdb_dir,
                n_job=config["options"]["n_cores"]
            )

        # ==================================================================
        # GCC vs GCC drydep plots
        # ==================================================================
        if config["options"]["outputs"]["plot_drydep"]:
            print("\n%%% Creating GCC vs. GCC drydep plots %%%")

            # Filepaths
            ref = get_filepath(gcc_vs_gcc_refdir, "DryDep", gcc_ref_date)
            dev = get_filepath(gcc_vs_gcc_devdir, "DryDep", gcc_dev_date)

            # Create plots
            make_benchmark_drydep_plots(
                ref,
                gcc_vs_gcc_refstr,
                dev,
                gcc_vs_gcc_devstr,
                dst=gcc_vs_gcc_resultsdir,
                weightsdir=config["paths"]["weights_dir"],
                overwrite=True,
                sigdiff_files=gcc_vs_gcc_sigdiff,
                spcdb_dir=spcdb_dir,
                n_job=config["options"]["n_cores"],
                varlist=drydepvel_species()
            )

        # ==================================================================
        # GCC vs GCC global mass tables
        # ==================================================================
        if config["options"]["outputs"]["mass_table"]:
            print("\n%%% Creating GCC vs. GCC global mass tables %%%")

            # Filepaths
            ref = get_filepath(gcc_vs_gcc_refrstdir, "Restart", gcc_end_ref_date)
            dev = get_filepath(gcc_vs_gcc_devrstdir, "Restart", gcc_end_dev_date)

            # Create tables
            make_benchmark_mass_tables(
                ref,
                config["data"]["ref"]["gcc"]["version"],
                dev,
                config["data"]["dev"]["gcc"]["version"],
                dst=gcc_vs_gcc_tablesdir,
                overwrite=True,
                spcdb_dir=spcdb_dir,
            )

        # ==================================================================
        # GCC vs GCC global mass accumulation tables
        # ==================================================================
        if config["options"]["outputs"]["mass_accum_table"]:
            print("\n%%% Creating GCC vs. GCC mass accumulation tables %%%")

            # Filepaths for start and end restart files
            refs = get_filepath(gcc_vs_gcc_refrstdir, "Restart", gcc_ref_date)
            devs = get_filepath(gcc_vs_gcc_devrstdir, "Restart", gcc_dev_date)
            refe = get_filepath(gcc_vs_gcc_refrstdir, "Restart", gcc_end_ref_date)
            deve = get_filepath(gcc_vs_gcc_devrstdir, "Restart", gcc_end_dev_date)

            # Get period strings
            refs_str = np.datetime_as_string(gcc_ref_date, unit="s")
            devs_str = np.datetime_as_string(gcc_dev_date, unit="s")
            refe_str = np.datetime_as_string(gcc_end_ref_date, unit="s")
            deve_str = np.datetime_as_string(gcc_end_dev_date, unit="s")
            refperiod = refs_str + ' - ' + refe_str
            devperiod = devs_str + ' - ' + deve_str

            # Create tables
            make_benchmark_mass_accumulation_tables(
                refs,
                refe,
                config["data"]["ref"]["gcc"]["version"],
                refperiod,
                devs,
                deve,
                config["data"]["dev"]["gcc"]["version"],
                devperiod,
                overwrite=True,
                dst=gcc_vs_gcc_tablesdir,
                spcdb_dir=spcdb_dir,
            )

        # ==================================================================
        # GCC vs GCC operation budgets tables
        # ==================================================================
        if config["options"]["outputs"]["ops_budget_table"]:
            print("\n%%% Creating GCC vs. GCC operations budget tables %%%")

            # Filepaths
            ref = get_filepath(gcc_vs_gcc_refdir, "Budget", gcc_ref_date)
            dev = get_filepath(gcc_vs_gcc_devdir, "Budget", gcc_dev_date)

            # Create table
            make_benchmark_operations_budget(
                config["data"]["ref"]["gcc"]["version"],
                ref,
                config["data"]["dev"]["gcc"]["version"],
                dev,
                gcc_ref_sec_diff,
                gcc_dev_sec_diff,
                benchmark_type=config["options"]["bmk_type"],
                label=comparison_str,
                dst=gcc_vs_gcc_tablesdir,
            )

        # ==================================================================
        # GCC vs. GCC global mean OH, MCF Lifetime, CH4 Lifetime
        # ==================================================================
        if config["options"]["outputs"]["OH_metrics"]:
            print("\n%%% Creating GCC vs. GCC OH metrics table %%%")

            # Filepaths
            ref = get_filepath(gcc_vs_gcc_refdir, "Metrics", gcc_ref_date)
            dev = get_filepath(gcc_vs_gcc_devdir, "Metrics", gcc_dev_date)

            # Create table
            make_benchmark_oh_metrics(
                ref,
                config["data"]["ref"]["gcc"]["version"],
                dev,
                config["data"]["dev"]["gcc"]["version"],
                dst=gcc_vs_gcc_tablesdir,
                overwrite=True,
                spcdb_dir=spcdb_dir,
            )

        # ==================================================================
        # GCC dev Strat-Trop Exchange
        # ==================================================================
        if config["options"]["outputs"]["ste_table"]:
            print("\n%%% Creating GCC dev Strat-Trop Exchange table %%%")

            # Diagnostic collection files to read
            dev = get_filepath(gcc_vs_gcc_devdir, "AdvFluxVert", gcc_dev_date)

            # Compute monthly and annual average strat-trop exchange of O3
            # only run if comparison is exactly 1 month
            if add_months(gcc_dev_date, 1) == gcc_end_dev_date:
                make_benchmark_ste_table(
                    config["data"]["dev"]["gcc"]["version"],
                    dev,
                    gcc_dev_date.astype(datetime).year,
                    bmk_type=config["options"]["bmk_type"],
                    dst=gcc_vs_gcc_tablesdir,
                    species=["O3"],
                    overwrite=True,
                    month=gcc_dev_date.astype(datetime).month,
                )

        # ==================================================================
        # GCC vs. GCC Benchmark Timing Table
        # ==================================================================
        if config["options"]["outputs"]["timing_table"]:
            print("\n%%% Creating GCC vs. GCC Benchmark Timing table %%%")

            # Filepaths
            ref = get_log_filepaths(
                gcc_vs_gcc_reflogdir,
                config["data"]["ref"]["gcc"]["logs_template"],
                gcc_ref_date,
            )
            dev = get_log_filepaths(
                gcc_vs_gcc_devlogdir,
                config["data"]["dev"]["gcc"]["logs_template"],
                gcc_dev_date,
            )

            # Create the table
            make_benchmark_gcclassic_timing_table(
                ref,
                config["data"]["ref"]["gcc"]["version"],
                dev,
                config["data"]["dev"]["gcc"]["version"],
                dst=gcc_vs_gcc_tablesdir,
                overwrite=True,
            )

        # ==================================================================
        # GCC vs. GCC summary table
        # ==================================================================
        if config["options"]["outputs"]["summary_table"]:
            print("\n%%% Creating GCC vs. GCC summary table %%%")

            # Print summary of which collections are identical
            # between Ref & Dev, and which are not identical.
            create_benchmark_summary_table(
                gcc_vs_gcc_refdir,
                config["data"]["ref"]["gcc"]["version"],
                gcc_ref_date,
                gcc_vs_gcc_devdir,
                config["data"]["dev"]["gcc"]["version"],
                gcc_dev_date,
                collections = [
                    'AerosolMass',
                    'Aerosols',
                    'DryDep',
                    'Emissions',
                    'JValues',
                    'Metrics',
                    'SpeciesConc',
                    'StateMet'
                ],
                dst=gcc_vs_gcc_tablesdir,
                outfilename="Summary.txt",
                overwrite=True,
                verbose=False,
            )

        # ==================================================================
        # GCC vs. GCC diagnostic sanity check table
        # ==================================================================
        if config["options"]["outputs"]["sanity_check_table"]:
            print("\n%%% Creating GCC vs. GCC sanity_check table %%%")

            # Print summary of which collections are identical
            # between Ref & Dev, and which are not identical.
            create_benchmark_sanity_check_table(
                gcc_vs_gcc_devdir,
                config["data"]["dev"]["gcc"]["version"],
                gcc_dev_date,
                collections = [
                    'AerosolMass',
                    'Aerosols',
                    'DryDep',
                    'Emissions',
                    'JValues',
                    'Metrics',
                    'SpeciesConc',
                    'StateMet'
                ],
                dst=gcc_vs_gcc_tablesdir,
                outfilename="Diagnostic_Sanity_Check.txt",
                overwrite=True,
                verbose=False,
            )
            

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Create GCHP vs GCC benchmark plots and tables
    #
    # (1) The GCC version in "GCHP vs GCC" is the Dev of "GCC vs GCC".
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if config["options"]["comparisons"]["gchp_vs_gcc"]["run"]:

        # ==================================================================
        # GCHP vs GCC error checks
        # ==================================================================
        if not np.equal(gcc_dev_sec_diff, gchp_dev_sec_diff):
            print("Skipping GCHP vs GCC emissions tables and operations")
            print("budget tables because months are different lengths")
            config["options"]["outputs"]["emis_table"] = False
            config["options"]["outputs"]["ops_budget_table"] = False

        # ==================================================================
        # GCHP vs GCC string for month and year
        # (e.g.  "2019-07-01T00:00:00 - 2019-08-01T00:00:00")
        # ==================================================================
        if np.equal(gcc_dev_date, gchp_dev_date) and np.equal(
            gcc_end_dev_date, gchp_end_dev_date
        ):
            comparison_str = (
                f"{config['data']['dev']['gcc']['bmk_start']} "
                + f"- {config['data']['dev']['gcc']['bmk_end']}"
            )
        else:
            comparison_str = (
                f"{config['data']['dev']['gcc']['bmk_start']} "
                + f"- {config['data']['dev']['gcc']['bmk_end']}"
                + f" Vs {config['data']['dev']['gchp']['bmk_start']} "
                + f"- {config['data']['dev']['gchp']['bmk_end']}"
            )

        # ==================================================================
        # GCHP vs GCC filepaths for StateMet collection data
        # ==================================================================
        refmet = get_filepath(gchp_vs_gcc_refdir, "StateMet", gcc_dev_date)
        devmet = get_filepath(
            gchp_vs_gcc_devdir,
            "StateMet",
            gchp_dev_date,
            is_gchp=True
        )

        # Get GCHP grid resolution from met collection file
        #ds_devmet = xr.open_dataset(devmet)
        #gchp_dev_res = str(get_input_res(ds_devmet)[0])

        # ==================================================================
        # GCHP vs GCC species concentration plots
        # ==================================================================
        if config["options"]["outputs"]["plot_conc"]:
            print("\n%%% Creating GCHP vs. GCC species concentration plots %%%")

            # Diagnostic collection files to read
            ref = get_filepath(gchp_vs_gcc_refdir, "SpeciesConc", gcc_dev_date)
            dev = get_filepath(
                gchp_vs_gcc_devdir,
                "SpeciesConc",
                gchp_dev_date,
                is_gchp=True
            )

            # Create plots
            make_benchmark_conc_plots(
                ref,
                gchp_vs_gcc_refstr,
                dev,
                gchp_vs_gcc_devstr,
                refmet=refmet,
                devmet=devmet,
                dst=gchp_vs_gcc_resultsdir,
                weightsdir=config["paths"]["weights_dir"],
                plot_by_spc_cat=config["options"]["outputs"]["plot_options"][
                    "by_spc_cat"],
                benchmark_type=config["options"]["bmk_type"],
                overwrite=True,
                sigdiff_files=gchp_vs_gcc_sigdiff,
                spcdb_dir=spcdb_dir,
                n_job=config["options"]["n_cores"]
            )

        # ==================================================================
        # GCHP vs. GCC Emissions plots
        # ==================================================================
        if config["options"]["outputs"]["plot_emis"]:
            print("\n%%% Creating GCHP vs. GCC emissions plots %%%")

            # Filepaths
            ref = get_filepath(gchp_vs_gcc_refdir, "Emissions", gcc_dev_date)
            dev = get_filepath(
                gchp_vs_gcc_devdir,
                "Emissions",
                gchp_dev_date,
                is_gchp=True
            )

            # Create emissions plots
            make_benchmark_emis_plots(
                ref,
                gchp_vs_gcc_refstr,
                dev,
                gchp_vs_gcc_devstr,
                dst=gchp_vs_gcc_resultsdir,
                weightsdir=config["paths"]["weights_dir"],
                plot_by_spc_cat=config["options"]["outputs"]["plot_options"][
                    "by_spc_cat"
                ],
                plot_by_hco_cat=config["options"]["outputs"]["plot_options"][
                    "by_hco_cat"
                ],
                benchmark_type=config["options"]["bmk_type"],
                overwrite=True,
                sigdiff_files=gchp_vs_gcc_sigdiff,
                spcdb_dir=spcdb_dir,
                n_job=config["options"]["n_cores"]
            )

        # ==================================================================
        # GCHP vs. GCC tables of emission and inventory totals
        # ==================================================================
        if config["options"]["outputs"]["emis_table"]:
            print("\n%%% Creating GCHP vs. GCC emissions/inventory tables %%%")

            # Filepaths
            ref = get_filepath(gchp_vs_gcc_refdir, "Emissions", gcc_dev_date)
            dev = get_filepath(
                gchp_vs_gcc_devdir,
                "Emissions",
                gchp_dev_date,
                is_gchp=True
            )

            # Create plots
            make_benchmark_emis_tables(
                ref,
                gchp_vs_gcc_refstr,
                dev,
                gchp_vs_gcc_devstr,
                dst=gchp_vs_gcc_resultsdir,
                ref_interval=[gcc_dev_sec_diff],
                dev_interval=[gchp_dev_sec_diff],
                benchmark_type=config["options"]["bmk_type"],
                overwrite=True,
                devmet=devmet,
                spcdb_dir=spcdb_dir,
            )

        # ==================================================================
        # GCHP vs. GCC J-values plots
        # ==================================================================
        if config["options"]["outputs"]["plot_jvalues"]:
            print("\n%%% Creating GCHP vs. GCC J-value plots %%%")

            # Filepaths
            ref = get_filepath(gchp_vs_gcc_refdir, "JValues", gcc_dev_date)
            dev = get_filepath(
                gchp_vs_gcc_devdir,
                "JValues",
                gchp_dev_date,
                is_gchp=True
            )

            # Create plots
            make_benchmark_jvalue_plots(
                ref,
                gchp_vs_gcc_refstr,
                dev,
                gchp_vs_gcc_devstr,
                dst=gchp_vs_gcc_resultsdir,
                weightsdir=config["paths"]["weights_dir"],
                overwrite=True,
                sigdiff_files=gchp_vs_gcc_sigdiff,
                spcdb_dir=spcdb_dir,
                n_job=config["options"]["n_cores"]
            )

        # ==================================================================
        # GCHP vs. GCC Budget plots
        # ==================================================================
        if config["options"]["outputs"]["plot_budget"]:
            print("\n%%% Creating GCHP vs. GCC Budget plots %%%")

            colname = "Budget"

            # Filepaths
            ref = get_filepath(gchp_vs_gcc_refdir, colname, gcc_dev_date)
            dev = get_filepath(
                gchp_vs_gcc_devdir,
                colname,
                gchp_dev_date,
                is_gchp=True
            )

            # Create plots
            make_benchmark_collection_2d_var_plots(
                ref,
                gchp_vs_gcc_refstr,
                dev,
                gchp_vs_gcc_devstr,
                colname=colname,
                var_prefix=colname,
                dst=gchp_vs_gcc_resultsdir,
                weightsdir=config["paths"]["weights_dir"],
                overwrite=True,
                spcdb_dir=spcdb_dir,
                n_job=config["options"]["n_cores"]
            )

        # ==================================================================
        # GCHP vs. GCC UVFlux plots
        # ==================================================================
        if config["options"]["outputs"]["plot_uvflux"]:
            print("\n%%% Creating GCHP vs. GCC UVFlux plots %%%")

            colname = "UVFlux"

            # Filepaths
            ref = get_filepath(gchp_vs_gcc_refdir, colname, gcc_dev_date)
            dev = get_filepath(
                gchp_vs_gcc_devdir,
                colname,
                gchp_dev_date,
                is_gchp=True
            )

            # Create plots
            make_benchmark_collection_3d_var_plots(
                ref,
                gchp_vs_gcc_refstr,
                dev,
                gchp_vs_gcc_devstr,
                colname=colname,
                dst=gchp_vs_gcc_resultsdir,
                weightsdir=config["paths"]["weights_dir"],
                overwrite=True,
                spcdb_dir=spcdb_dir,
                n_job=config["options"]["n_cores"]
            )

        # ==================================================================
        # GCHP vs. GCC 2D StateMet plots
        # ==================================================================
        if config["options"]["outputs"]["plot_2d_met"]:
            print("\n%%% Creating GCHP vs. GCC 2D StateMet plots %%%")

            # Create plots
            make_benchmark_collection_2d_var_plots(
                refmet,
                gchp_vs_gcc_refstr,
                devmet,
                gchp_vs_gcc_devstr,
                colname="StateMet",
                var_prefix='Met',
                dst=gchp_vs_gcc_resultsdir,
                weightsdir=config["paths"]["weights_dir"],
                overwrite=True,
                spcdb_dir=spcdb_dir,
                n_job=config["options"]["n_cores"]
            )

        # ==================================================================
        # GCHP vs. GCC 3D StateMet plots
        # ==================================================================
        if config["options"]["outputs"]["plot_3d_met"]:
            print("\n%%% Creating GCHP vs. GCC 3D StateMet plots %%%")

            # Create plots
            make_benchmark_collection_3d_var_plots(
                refmet,
                gchp_vs_gcc_refstr,
                devmet,
                gchp_vs_gcc_devstr,
                colname="StateMet",
                var_prefix='Met',
                dst=gchp_vs_gcc_resultsdir,
                weightsdir=config["paths"]["weights_dir"],
                overwrite=True,
                spcdb_dir=spcdb_dir,
                n_job=config["options"]["n_cores"]
            )

        # ==================================================================
        # GCHP vs GCC column AOD plots
        # ==================================================================
        if config["options"]["outputs"]["plot_aod"]:
            print("\n%%% Creating GCHP vs. GCC column AOD plots %%%")

            # Filepaths
            ref = get_filepath(gchp_vs_gcc_refdir, "Aerosols", gcc_dev_date)
            dev = get_filepath(
                gchp_vs_gcc_devdir,
                "Aerosols",
                gchp_dev_date,
                is_gchp=True
            )

            # Create plots
            make_benchmark_aod_plots(
                ref,
                gchp_vs_gcc_refstr,
                dev,
                gchp_vs_gcc_devstr,
                dst=gchp_vs_gcc_resultsdir,
                weightsdir=config["paths"]["weights_dir"],
                overwrite=True,
                sigdiff_files=gchp_vs_gcc_sigdiff,
                spcdb_dir=spcdb_dir,
                n_job=config["options"]["n_cores"]
            )

        # ==================================================================
        # GCHP vs GCC drydep plots
        # ==================================================================
        if config["options"]["outputs"]["plot_drydep"]:
            print("\n%%% Creating GCHP vs. GCC drydep plots %%%")

            # Filepaths
            ref = get_filepath(gchp_vs_gcc_refdir, "DryDep", gcc_ref_date)
            dev = get_filepath(
                gchp_vs_gcc_devdir,
                "DryDep",
                gchp_dev_date,
                is_gchp=True
            )

            # Create plots
            make_benchmark_drydep_plots(
                ref,
                gchp_vs_gcc_refstr,
                dev,
                gchp_vs_gcc_devstr,
                dst=gchp_vs_gcc_resultsdir,
                weightsdir=config["paths"]["weights_dir"],
                overwrite=True,
                sigdiff_files=gchp_vs_gcc_sigdiff,
                spcdb_dir=spcdb_dir,
                n_job=config["options"]["n_cores"],
                varlist=drydepvel_species(),
            )

        # ==================================================================
        # GCHP vs GCC global mass tables
        # ==================================================================
        if config["options"]["outputs"]["mass_table"]:
            print("\n%%% Creating GCHP vs. GCC global mass tables %%%")

            # Filepaths
            ref = get_filepath(
                gchp_vs_gcc_refrstdir,
                "Restart",
                gcc_end_dev_date
            )
            dev = get_filepath(
                gchp_vs_gcc_devrstdir,
                "Restart",
                gchp_end_dev_date,
                is_gchp=True,
                gchp_res=config["data"]["dev"]["gchp"]["resolution"],
                gchp_is_pre_14_0=config["data"]["dev"]["gchp"]["is_pre_14.0"]
            )

            # Create tables
            make_benchmark_mass_tables(
                ref,
                gchp_vs_gcc_refstr,
                dev,
                gchp_vs_gcc_devstr,
                dst=gchp_vs_gcc_tablesdir,
                overwrite=True,
                spcdb_dir=spcdb_dir,
            )

        # ==================================================================
        # GCHP vs GCC global mass accumulation tables
        # ==================================================================
        if config["options"]["outputs"]["mass_accum_table"]:
            print("\n%%% Creating GCHP vs. GCC mass accumulation tables %%%")

            # Filepaths for start and end restart files
            refs = get_filepath(
                gchp_vs_gcc_refrstdir,
                "Restart",
                gcc_dev_date
            )
            devs = get_filepath(
                gchp_vs_gcc_devrstdir,
                "Restart",
                gchp_dev_date,
                is_gchp=True,
                gchp_res=config["data"]["dev"]["gchp"]["resolution"],
                gchp_is_pre_14_0=config["data"]["dev"]["gchp"]["is_pre_14.0"]
            )
            refe = get_filepath(
                gchp_vs_gcc_refrstdir,
                "Restart",
                gcc_end_dev_date
            )
            deve = get_filepath(
                gchp_vs_gcc_devrstdir,
                "Restart",
                gchp_end_dev_date,
                is_gchp=True,
                gchp_res=config["data"]["dev"]["gchp"]["resolution"],
                gchp_is_pre_14_0=config["data"]["dev"]["gchp"]["is_pre_14.0"]
            )

            # Get period strings
            refs_str = np.datetime_as_string(gcc_dev_date, unit="s")
            devs_str = np.datetime_as_string(gchp_dev_date, unit="s")
            refe_str = np.datetime_as_string(gcc_end_dev_date, unit="s")
            deve_str = np.datetime_as_string(gchp_end_dev_date, unit="s")
            refperiod = refs_str + ' - ' + refe_str
            devperiod = devs_str + ' - ' + deve_str

            # Create tables
            make_benchmark_mass_accumulation_tables(
                refs,
                refe,
                config["data"]["dev"]["gcc"]["version"],
                refperiod,
                devs,
                deve,
                config["data"]["dev"]["gchp"]["version"],
                devperiod,
                overwrite=True,
                dst=gchp_vs_gcc_tablesdir,
                spcdb_dir=spcdb_dir,
            )

        # ==================================================================
        # GCHP vs GCC operations budgets tables
        # ==================================================================
        if config["options"]["outputs"]["ops_budget_table"]:
            print("\n%%% Creating GCHP vs GCC operations budget tables %%%")

            # Filepaths
            ref = get_filepath(gchp_vs_gcc_refdir, "Budget", gcc_dev_date)
            dev = get_filepath(
                gchp_vs_gcc_devdir,
                "Budget",
                gchp_dev_date,
                is_gchp=True
            )

            # Create table
            make_benchmark_operations_budget(
                config["data"]["dev"]["gcc"]["version"],
                ref,
                config["data"]["dev"]["gchp"]["version"],
                dev,
                gcc_dev_sec_diff,
                gchp_dev_sec_diff,
                benchmark_type=config["options"]["bmk_type"],
                label=comparison_str,
                operations=[
                    "Chemistry",
                    "Convection",
                    "EmisDryDep",
                    "Mixing",
                    "WetDep",
                ],
                compute_accum=False,
                dst=gchp_vs_gcc_tablesdir,
            )

        # ==================================================================
        # GCHP vs. GCC global mean OH, MCF Lifetime, CH4 Lifetime
        # ==================================================================
        if config["options"]["outputs"]["OH_metrics"]:
            print("\n%%% Creating GCHP vs GCC OH metrics table %%%")

            # Diagnostic collection files to read
            ref = get_filepath(gchp_vs_gcc_refdir, "Metrics", gcc_dev_date)
            dev = get_filepath(
                gchp_vs_gcc_devdir,
                "Metrics",
                gchp_dev_date,
                is_gchp=True
            )

            # Create table
            make_benchmark_oh_metrics(
                ref,
                gchp_vs_gcc_refstr,
                dev,
                gchp_vs_gcc_devstr,
                dst=gchp_vs_gcc_tablesdir,
                overwrite=True,
                spcdb_dir=spcdb_dir,
            )

        # ==================================================================
        # GCHP vs GCC Strat-Trop Exchange
        # ==================================================================
        if config["options"]["outputs"]["ste_table"]:
            title = "\n%%% Skipping GCHP vs. GCC Strat-Trop Exchange table %%%"
            print(title)

        # ==================================================================
        # GCHP vs. GCC summary table
        # ==================================================================
        if config["options"]["outputs"]["summary_table"]:
            print("\n%%% Creating GCHP vs. GCC summary table %%%")

            # Print summary of which collections are identical
            # between Ref & Dev, and which are not identical.
            create_benchmark_summary_table(
                gchp_vs_gcc_refdir,
                config["data"]["dev"]["gcc"]["version"],
                gcc_dev_date,
                gchp_vs_gcc_devdir,
                config["data"]["dev"]["gchp"]["version"],
                gchp_dev_date,
                collections=[
                    'AerosolMass',
                    'Aerosols',
                    'DryDep',
                    'Emissions',
                    'JValues',
                    'Metrics',
                    'SpeciesConc',
                    'StateMet',
                ],
                dst=gchp_vs_gcc_tablesdir,
                outfilename="Summary.txt",
                overwrite=True,
                verbose=False,
                dev_gchp=True
            )
            
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Create GCHP vs GCHP benchmark plots and tables
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if config["options"]["comparisons"]["gchp_vs_gchp"]["run"]:

        # ==================================================================
        # GCHP vs GCHP error checks
        # ==================================================================
        if not np.equal(gchp_ref_sec_diff, gchp_dev_sec_diff):
            print("Skipping GCHP vs GCHP emissions tables and operations")
            print("budget tables because months are different lengths")
            config["options"]["outputs"]["emis_table"] = False
            config["options"]["outputs"]["ops_budget_table"] = False

        # ==================================================================
        # GCHP vs GCHP string for month and year
        # (e.g.  "2019-07-01T00:00:00 - 2019-08-01T00:00:00")
        # ==================================================================
        if np.equal(gchp_ref_date, gchp_dev_date) and np.equal(
            gchp_end_ref_date, gchp_end_dev_date
        ):
            comparison_str = (
                f"{config['data']['dev']['gchp']['bmk_start']} "
                + f"- {config['data']['dev']['gchp']['bmk_end']}"
            )
        else:
            comparison_str = (
                f"{config['data']['ref']['gchp']['bmk_start']} "
                + f"- {config['data']['ref']['gchp']['bmk_end']}"
                + f" Vs {config['data']['dev']['gchp']['bmk_start']} "
                + f"- {config['data']['dev']['gchp']['bmk_end']}"
            )

        # ==================================================================
        # GCHP vs GCHP filepaths for StateMet collection data
        # ==================================================================
        refmet = get_filepath(
            gchp_vs_gchp_refdir,
            "StateMet",
            gchp_ref_date,
            is_gchp=True
        )
        devmet = get_filepath(
            gchp_vs_gchp_devdir,
            "StateMet",
            gchp_dev_date,
            is_gchp=True
        )

        # Get GCHP grid resolutions from met collection file
        #ds_refmet = xr.open_dataset(refmet)
        #ds_devmet = xr.open_dataset(devmet)
        #gchp_ref_res = str(get_input_res(ds_refmet)[0])
        #gchp_dev_res = str(get_input_res(ds_devmet)[0])

        # ==================================================================
        # GCHP vs GCHP species concentration plots
        # ==================================================================
        if config["options"]["outputs"]["plot_conc"]:
            print("\n%%% Creating GCHP vs. GCHP concentration plots %%%")

            # Filepaths
            ref = get_filepath(
                gchp_vs_gchp_refdir,
                "SpeciesConc",
                gchp_ref_date,
                is_gchp=True
            )
            dev = get_filepath(
                gchp_vs_gchp_devdir,
                "SpeciesConc",
                gchp_dev_date,
                is_gchp=True
            )

            # Create plots
            make_benchmark_conc_plots(
                ref,
                gchp_vs_gchp_refstr,
                dev,
                gchp_vs_gchp_devstr,
                refmet=refmet,
                devmet=devmet,
                dst=gchp_vs_gchp_resultsdir,
                weightsdir=config["paths"]["weights_dir"],
                plot_by_spc_cat=config["options"]["outputs"]["plot_options"][
                    "by_spc_cat"],
                benchmark_type=config["options"]["bmk_type"],
                overwrite=True,
                sigdiff_files=gchp_vs_gchp_sigdiff,
                spcdb_dir=spcdb_dir,
                n_job=config["options"]["n_cores"]
            )

        # ==================================================================
        # GCHP vs. GCHP Emissions plots
        # ==================================================================
        if config["options"]["outputs"]["plot_emis"]:
            print("\n%%% Creating GCHP vs. GCHP emissions plots %%%")

            # Filepaths
            ref = get_filepath(
                gchp_vs_gchp_refdir,
                "Emissions",
                gchp_ref_date,
                is_gchp=True
            )
            dev = get_filepath(
                gchp_vs_gchp_devdir,
                "Emissions",
                gchp_dev_date,
                is_gchp=True
            )

            # Create plots
            make_benchmark_emis_plots(
                ref,
                gchp_vs_gchp_refstr,
                dev,
                gchp_vs_gchp_devstr,
                dst=gchp_vs_gchp_resultsdir,
                weightsdir=config["paths"]["weights_dir"],
                plot_by_spc_cat=config["options"]["outputs"]["plot_options"][
                    "by_spc_cat"],
                plot_by_hco_cat=config["options"]["outputs"]["plot_options"][
                    "by_hco_cat"],
                benchmark_type=config["options"]["bmk_type"],
                overwrite=True,
                sigdiff_files=gchp_vs_gchp_sigdiff,
                spcdb_dir=spcdb_dir,
                n_job=config["options"]["n_cores"]
            )

        # ==================================================================
        # GCHP vs. GCHP tables of emission and inventory totals
        # ==================================================================
        if config["options"]["outputs"]["emis_table"]:
            print("\n%%% Creating GCHP vs. GCHP emissions/inventory tables %%%")

            # Filepaths
            ref = get_filepath(
                gchp_vs_gchp_refdir,
                "Emissions",
                gchp_ref_date,
                is_gchp=True
            )
            dev = get_filepath(
                gchp_vs_gchp_devdir,
                "Emissions",
                gchp_dev_date,
                is_gchp=True
            )

            # Create tables
            make_benchmark_emis_tables(
                ref,
                gchp_vs_gchp_refstr,
                dev,
                gchp_vs_gchp_devstr,
                dst=gchp_vs_gchp_resultsdir,
                ref_interval=[gchp_ref_sec_diff],
                dev_interval=[gchp_dev_sec_diff],
                benchmark_type=config["options"]["bmk_type"],
                overwrite=True,
                refmet=refmet,
                devmet=devmet,
                spcdb_dir=spcdb_dir,
            )

        # ==================================================================
        # GCHP vs. GCHP J-values plots
        # ==================================================================
        if config["options"]["outputs"]["plot_jvalues"]:
            print("\n%%% Creating GCHP vs. GCHP J-value plots %%%")

            # Filepaths
            ref = get_filepath(
                gchp_vs_gchp_refdir,
                "JValues",
                gchp_ref_date,
                is_gchp=True
            )
            dev = get_filepath(
                gchp_vs_gchp_devdir,
                "JValues",
                gchp_dev_date,
                is_gchp=True
            )

            # Create plots
            make_benchmark_jvalue_plots(
                ref,
                gchp_vs_gchp_refstr,
                dev,
                gchp_vs_gchp_devstr,
                dst=gchp_vs_gchp_resultsdir,
                weightsdir=config["paths"]["weights_dir"],
                overwrite=True,
                sigdiff_files=gchp_vs_gchp_sigdiff,
                spcdb_dir=spcdb_dir,
                n_job=config["options"]["n_cores"]
            )

        # ==================================================================
        # GCHP vs. GCHP Budget plots
        # ==================================================================
        if config["options"]["outputs"]["plot_budget"]:
            print("\n%%% Creating GCHP vs. GCHP Budget plots %%%")

            colname = "Budget"

            # Filepaths
            ref = get_filepath(
                gchp_vs_gchp_refdir,
                colname,
                gchp_ref_date,
                is_gchp=True
            )
            dev = get_filepath(
                gchp_vs_gchp_devdir,
                colname,
                gchp_dev_date,
                is_gchp=True
            )

            # Create plots
            make_benchmark_collection_2d_var_plots(
                ref,
                gchp_vs_gchp_refstr,
                dev,
                gchp_vs_gchp_devstr,
                colname=colname,
                var_prefix=colname,
                dst=gchp_vs_gchp_resultsdir,
                weightsdir=config["paths"]["weights_dir"],
                overwrite=True,
                spcdb_dir=spcdb_dir,
                n_job=config["options"]["n_cores"]
            )

        # ==================================================================
        # GCHP vs. GCHP UVFlux plots
        # ==================================================================
        if config["options"]["outputs"]["plot_uvflux"]:
            print("\n%%% Creating GCHP vs. GCHP UVFlux plots %%%")

            colname = "UVFlux"

            # Filepaths
            ref = get_filepath(
                gchp_vs_gchp_refdir,
                colname,
                gchp_ref_date,
                is_gchp=True
            )
            dev = get_filepath(
                gchp_vs_gchp_devdir,
                colname,
                gchp_dev_date,
                is_gchp=True
            )

            # Create plots
            make_benchmark_collection_3d_var_plots(
                ref,
                gchp_vs_gchp_refstr,
                dev,
                gchp_vs_gchp_devstr,
                colname=colname,
                dst=gchp_vs_gchp_resultsdir,
                weightsdir=config["paths"]["weights_dir"],
                overwrite=True,
                spcdb_dir=spcdb_dir,
                n_job=config["options"]["n_cores"]
            )

        # ==================================================================
        # GCHP vs. GCHP 2D StateMet plots
        # ==================================================================
        if config["options"]["outputs"]["plot_2d_met"]:
            print("\n%%% Creating GCHP vs. GCHP 2D StateMet plots %%%")

            # Create plots
            make_benchmark_collection_2d_var_plots(
                refmet,
                gchp_vs_gchp_refstr,
                devmet,
                gchp_vs_gchp_devstr,
                colname="StateMet",
                var_prefix='Met',
                dst=gchp_vs_gchp_resultsdir,
                weightsdir=config["paths"]["weights_dir"],
                overwrite=True,
                spcdb_dir=spcdb_dir,
                n_job=config["options"]["n_cores"]
            )

        # ==================================================================
        # GCHP vs. GCHP 3D StateMet plots
        # ==================================================================
        if config["options"]["outputs"]["plot_3d_met"]:
            print("\n%%% Creating GCHP vs. GCHP 3D StateMet plots %%%")

            # Create plots
            make_benchmark_collection_3d_var_plots(
                refmet,
                gchp_vs_gchp_refstr,
                devmet,
                gchp_vs_gchp_devstr,
                colname="StateMet",
                var_prefix='Met',
                dst=gchp_vs_gchp_resultsdir,
                weightsdir=config["paths"]["weights_dir"],
                overwrite=True,
                spcdb_dir=spcdb_dir,
                n_job=config["options"]["n_cores"]
            )

        # ==================================================================
        # GCHP vs GCHP column AOD plots
        # ==================================================================
        if config["options"]["outputs"]["plot_aod"]:
            print("\n%%% Creating GCHP vs. GCHP column AOD plots %%%")

            # Filepaths
            ref = get_filepath(
                gchp_vs_gchp_refdir,
                "Aerosols",
                gchp_ref_date,
                is_gchp=True
            )
            dev = get_filepath(
                gchp_vs_gchp_devdir,
                "Aerosols",
                gchp_dev_date,
                is_gchp=True
            )

            # Create plots
            make_benchmark_aod_plots(
                ref,
                gchp_vs_gchp_refstr,
                dev,
                gchp_vs_gchp_devstr,
                dst=gchp_vs_gchp_resultsdir,
                weightsdir=config["paths"]["weights_dir"],
                overwrite=True,
                sigdiff_files=gchp_vs_gchp_sigdiff,
                spcdb_dir=spcdb_dir,
                n_job=config["options"]["n_cores"]
            )

        # ==================================================================
        # GCHP vs GCHP drydep plots
        # ==================================================================
        if config["options"]["outputs"]["plot_drydep"]:
            print("\n%%% Creating GCHP vs. GCHP drydep plots %%%")

            # Filepaths
            ref = get_filepath(
                gchp_vs_gchp_refdir,
                "DryDep",
                gchp_ref_date,
                is_gchp=True
            )
            dev = get_filepath(
                gchp_vs_gchp_devdir,
                "DryDep",
                gchp_dev_date,
                is_gchp=True
            )

            # Create plots
            make_benchmark_drydep_plots(
                ref,
                gchp_vs_gchp_refstr,
                dev,
                gchp_vs_gchp_devstr,
                dst=gchp_vs_gchp_resultsdir,
                weightsdir=config["paths"]["weights_dir"],
                overwrite=True,
                sigdiff_files=gchp_vs_gchp_sigdiff,
                spcdb_dir=spcdb_dir,
                n_job=config["options"]["n_cores"],
                varlist=drydepvel_species()
            )

        # ==================================================================
        # GCHP vs GCHP global mass tables
        # ==================================================================
        if config["options"]["outputs"]["mass_table"]:
            print("\n%%% Creating GCHP vs. GCHP global mass tables %%%")

            # Filepaths
            ref = get_filepath(
                gchp_vs_gchp_refrstdir,
                "Restart",
                gchp_end_ref_date,
                is_gchp=True,
                gchp_res=config["data"]["ref"]["gchp"]["resolution"],
                gchp_is_pre_14_0=config["data"]["ref"]["gchp"]["is_pre_14.0"]
            )
            dev = get_filepath(
                gchp_vs_gchp_devrstdir,
                "Restart",
                gchp_end_dev_date,
                is_gchp=True,
                gchp_res=config["data"]["dev"]["gchp"]["resolution"],
                gchp_is_pre_14_0=config["data"]["dev"]["gchp"]["is_pre_14.0"]
            )

            # Create tables
            make_benchmark_mass_tables(
                ref,
                gchp_vs_gchp_refstr,
                dev,
                gchp_vs_gchp_devstr,
                dst=gchp_vs_gchp_tablesdir,
                overwrite=True,
                spcdb_dir=spcdb_dir,
            )

        # ==================================================================
        # GCHP vs GCHP global mass accumulation tables
        # ==================================================================
        if config["options"]["outputs"]["mass_accum_table"]:
            print("\n%%% Creating GCHP vs. GCHP mass accumulation tables %%%")

            # Filepaths for start and end restart files
            refs = get_filepath(
                gchp_vs_gchp_refrstdir,
                "Restart",
                gchp_ref_date,
                is_gchp=True,
                gchp_res=config["data"]["ref"]["gchp"]["resolution"],
                gchp_is_pre_14_0=config["data"]["ref"]["gchp"]["is_pre_14.0"]
            )
            devs = get_filepath(
                gchp_vs_gchp_devrstdir,
                "Restart",
                gchp_dev_date,
                is_gchp=True,
                gchp_res=config["data"]["dev"]["gchp"]["resolution"],
                gchp_is_pre_14_0=config["data"]["dev"]["gchp"]["is_pre_14.0"]
            )
            refe = get_filepath(
                gchp_vs_gchp_refrstdir,
                "Restart",
                gchp_end_ref_date,
                is_gchp=True,
                gchp_res=config["data"]["ref"]["gchp"]["resolution"],
                gchp_is_pre_14_0=config["data"]["ref"]["gchp"]["is_pre_14.0"]
            )
            deve = get_filepath(
                gchp_vs_gchp_devrstdir,
                "Restart",
                gchp_end_dev_date,
                is_gchp=True,
                gchp_res=config["data"]["dev"]["gchp"]["resolution"],
                gchp_is_pre_14_0=config["data"]["dev"]["gchp"]["is_pre_14.0"]
            )

            # Get period strings
            refs_str = np.datetime_as_string(gchp_ref_date, unit="s")
            devs_str = np.datetime_as_string(gchp_dev_date, unit="s")
            refe_str = np.datetime_as_string(gchp_end_ref_date, unit="s")
            deve_str = np.datetime_as_string(gchp_end_dev_date, unit="s")
            refperiod = refs_str + ' - ' + refe_str
            devperiod = devs_str + ' - ' + deve_str

            # Create tables
            make_benchmark_mass_accumulation_tables(
                refs,
                refe,
                config["data"]["ref"]["gchp"]["version"],
                refperiod,
                devs,
                deve,
                config["data"]["dev"]["gchp"]["version"],
                devperiod,
                overwrite=True,
                dst=gchp_vs_gchp_tablesdir,
                spcdb_dir=spcdb_dir,
            )

        # ==================================================================
        # GCHP vs GCHP operations budgets tables
        # ==================================================================
        if config["options"]["outputs"]["ops_budget_table"]:
            print("\n%%% Creating GCHP vs GCHP operations budget tables %%%")

            # Filepaths
            ref = get_filepath(
                gchp_vs_gchp_refdir,
                "Budget",
                gchp_ref_date,
                is_gchp=True
            )
            dev = get_filepath(
                gchp_vs_gchp_devdir,
                "Budget",
                gchp_dev_date,
                is_gchp=True
            )

            # Create tables
            make_benchmark_operations_budget(
                config["data"]["ref"]["gchp"]["version"],
                ref,
                config["data"]["dev"]["gchp"]["version"],
                dev,
                gchp_ref_sec_diff,
                gchp_dev_sec_diff,
                benchmark_type=config["options"]["bmk_type"],
                label=comparison_str,
                operations=[
                    "Chemistry",
                    "Convection",
                    "EmisDryDep",
                    "Mixing",
                    "WetDep",
                ],
                compute_accum=False,
                dst=gchp_vs_gchp_tablesdir,
            )

        # ==================================================================
        # GCHP vs. GCHP global mean OH, MCF Lifetime, CH4 Lifetime
        # ==================================================================
        if config["options"]["outputs"]["OH_metrics"]:
            print("\n%%% Creating GCHP vs GCHP OH metrics table %%%")

            # Diagnostic collection files to read
            ref = get_filepath(
                gchp_vs_gchp_refdir,
                "Metrics",
                gchp_ref_date,
                is_gchp=True
            )
            dev = get_filepath(
                gchp_vs_gchp_devdir,
                "Metrics",
                gchp_dev_date,
                is_gchp=True
            )

            # Create table
            make_benchmark_oh_metrics(
                ref,
                config["data"]["ref"]["gchp"]["version"],
                dev,
                config["data"]["dev"]["gchp"]["version"],
                dst=gchp_vs_gchp_tablesdir,
                overwrite=True,
                spcdb_dir=spcdb_dir,
            )

        # ==================================================================
        # GCHP vs GCHP Strat-Trop Exchange
        # ==================================================================
        if config["options"]["outputs"]["ste_table"]:
            print("\n%%% Skipping GCHP vs. GCHP Strat-Trop Exchange table %%%")

        # ==================================================================
        # GCHP vs. GCHP Benchmark Timing Table
        # ==================================================================
        if config["options"]["outputs"]["timing_table"]:
            print("\n%%% Creating GCHP vs. GCHP Benchmark Timing table %%%")

            # Filepaths
            ref = get_log_filepaths(
                gchp_vs_gchp_reflogdir,
                config["data"]["ref"]["gchp"]["logs_template"],
                gchp_ref_date,
            )
            dev = get_log_filepaths(
                gchp_vs_gchp_devlogdir,
                config["data"]["dev"]["gchp"]["logs_template"],
                gchp_dev_date,
            )

            # Create the table
            make_benchmark_gchp_timing_table(
                ref,
                config["data"]["ref"]["gchp"]["version"],
                dev,
                config["data"]["dev"]["gchp"]["version"],
                dst=gchp_vs_gchp_tablesdir,
                overwrite=True,
            )
            
        # ==================================================================
        # GCHP vs. GCHP summary table
        # ==================================================================
        if config["options"]["outputs"]["summary_table"]:
            print("\n%%% Creating GCHP vs. GCHP summary table %%%")

            # Print summary of which collections are identical
            # between Ref & Dev, and which are not identical.
            create_benchmark_summary_table(
                gchp_vs_gchp_refdir,
                config["data"]["ref"]["gchp"]["version"],
                gchp_ref_date,
                gchp_vs_gchp_devdir,
                config["data"]["dev"]["gchp"]["version"],
                gchp_dev_date,
                collections=[
                    'AerosolMass',
                    'Aerosols',
                    'DryDep',
                    'Emissions',
                    'JValues',
                    'Metrics',
                    'SpeciesConc',
                    'StateMet',
                ],
                dst=gchp_vs_gchp_tablesdir,
                outfilename="Summary.txt",
                overwrite=True,
                verbose=False,
                ref_gchp=True,
                dev_gchp=True,
            )

        # ==================================================================
        # GCHP vs. GCHP diagnostic sanity check table
        # ==================================================================
        if config["options"]["outputs"]["sanity_check_table"]:
            print("\n%%% Creating GCHP vs. GCHP sanity check table %%%")

            # Print summary of which collections are identical
            # between Ref & Dev, and which are not identical.
            create_benchmark_sanity_check_table(
                gchp_vs_gchp_devdir,
                config["data"]["dev"]["gchp"]["version"],
                gchp_dev_date,
                collections=[
                    'AerosolMass',
                    'Aerosols',
                    'DryDep',
                    'Emissions',
                    'JValues',
                    'Metrics',
                    'SpeciesConc',
                    'StateMet',
                ],
                dst=gchp_vs_gchp_tablesdir,
                is_gchp=True,
                outfilename="Diagnostic_Sanity_Check.txt",
                overwrite=True,
                verbose=False,
            )
            
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Create GCHP vs GCC difference of differences benchmark plots
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if config["options"]["comparisons"]["gchp_vs_gcc_diff_of_diffs"]["run"]:

        if config["options"]["outputs"]["plot_conc"]:
            print("\n%%% Creating GCHP vs. GCC diff-of-diffs conc plots %%%")

            # Filepaths
            gcc_ref = get_filepath(gcc_vs_gcc_refdir, "SpeciesConc", gcc_ref_date)
            gcc_dev = get_filepath(gcc_vs_gcc_devdir, "SpeciesConc", gcc_dev_date)
            gchp_ref = get_filepath(
                gchp_vs_gchp_refdir,
                "SpeciesConc",
                gchp_ref_date,
                is_gchp=True
            )
            gchp_dev = get_filepath(
                gchp_vs_gchp_devdir,
                "SpeciesConc",
                gchp_dev_date,
                is_gchp=True
            )

            # Create diff-of-diff plots for species concentrations
            # NOTE: for simplicity, do not convert aerosols to ug/m3.
            make_benchmark_conc_plots(
                gcc_ref,
                diff_of_diffs_refstr,
                gchp_ref,
                diff_of_diffs_devstr,
                dst=diff_of_diffs_resultsdir,
                weightsdir=config["paths"]["weights_dir"],
                benchmark_type=config["options"]["bmk_type"],
                overwrite=True,
                use_cmap_RdBu=True,
                second_ref=gcc_dev,
                second_dev=gchp_dev,
                cats_in_ugm3=None,
                spcdb_dir=spcdb_dir,
                n_job=config["options"]["n_cores"]
            )


    # ==================================================================
    # Print a message indicating that the benchmarks finished
    # ==================================================================
    print("\n%%%% All requested benchmark plots/tables created! %%%%")


def main(argv):
    """
    Driver program. Determines which benchmark script script to call
    for 1-hour, 1-day, 1-month, or 1-year benchmarks.
    """
    config_filename = argv[1] if len(argv) == 2 else "1mo_benchmark.yml"
    config = read_config_file(config_filename)
    config["configuration_file_name"] = config_filename
    choose_benchmark_type(config)


if __name__ == "__main__":
    main(sys.argv)
