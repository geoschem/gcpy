#!/usr/bin/env python
"""
run_1mo_benchmark.py: Driver script for creating benchmark plots and
                      testing gcpy 1-month benchmark capability

Run this script to generate benchmark comparisons between:

    (1) GCC (aka GEOS-Chem "Classic") vs. GCC
    (2) GCHP vs GCC
    (3) GCHP vs GCHP
    (4) GCHP vs GCC diff-of-diffs

You can customize this by editing the settings in the corresponding yaml
config file (eg. 1mo_benchmark.yml).

Calling sequence:

    ./run_1mo_benchmark.py <path-to-configuration-file>

To test gcpy, copy this script and the corresponding yaml config file
anywhere you want to run the test. Set gcpy_test to True at the top
of the script. Benchmark artifacts will be created locally in new folder
called Plots.

Remarks:

    By default, matplotlib will try to open an X window for plotting.
    If you are running this script in an environment where you do not have
    an active X display (such as in a computational queue), then you will
    need to use these commands to disable the X-window functionality.

        import os
        os.environ["QT_QPA_PLATFORM"]="offscreen"

    For more information, please see this issue posted at the ipython site:

        https://github.com/ipython/ipython/issues/10627

This script corresponds with GCPy 1.3.2. Edit this version ID if releasing
a new version of GCPy.
"""

# =====================================================================
# Imports and global settings (you should not need to edit these)
# ====================================q=================================

import os
import sys
from shutil import copyfile
import warnings
from datetime import datetime
import numpy as np
from gcpy.util import get_filepath, read_config_file
import gcpy.ste_flux as ste
import gcpy.oh_metrics as oh
import gcpy.benchmark as bmk
from gcpy.date_time import add_months, is_full_year
from modules.run_1yr_fullchem_benchmark import run_benchmark as run_1yr_benchmark
from modules.run_1yr_tt_benchmark import run_benchmark as run_1yr_tt_benchmark

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
    ):
        print(
            f"Error: invalid benchmark type {config['options']['bmk_type']}. "
            + "Please enter FullChemBenchmark or TransportTracersBenchmark."
        )
        sys.exit()

    start = np.datetime64(config["data"]["ref"]["gcc"]["bmk_start"])
    end = np.datetime64(config["data"]["ref"]["gcc"]["bmk_end"])
    # determine benchmark type and run relevant script
    if is_full_year(start, end):
        if config["options"]["bmk_type"] == "FullChemBenchmark":
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
    spcdb_dir = bmk.get_species_database_dir(config)

    # =====================================================================
    # Data directories
    # For gchp_vs_gcc_refdir use config["data"]["dev"]["gcc"]["version"],
    # not ref (mps, 6/27/19)
    # =====================================================================

    # Diagnostic file directory paths
    gcc_vs_gcc_refdir = os.path.join(
        config["paths"]["main_dir"],
        config["data"]["ref"]["gcc"]["dir"],
        config["data"]["ref"]["gcc"]["outputs_subdir"],
    )
    gcc_vs_gcc_devdir = os.path.join(
        config["paths"]["main_dir"],
        config["data"]["dev"]["gcc"]["dir"],
        config["data"]["dev"]["gcc"]["outputs_subdir"],
    )
    gchp_vs_gcc_refdir = os.path.join(
        config["paths"]["main_dir"],
        config["data"]["dev"]["gcc"]["dir"],
        config["data"]["dev"]["gcc"]["outputs_subdir"],
    )
    gchp_vs_gcc_devdir = os.path.join(
        config["paths"]["main_dir"],
        config["data"]["dev"]["gchp"]["dir"],
        config["data"]["dev"]["gchp"]["outputs_subdir"],
    )
    gchp_vs_gchp_refdir = os.path.join(
        config["paths"]["main_dir"],
        config["data"]["ref"]["gchp"]["dir"],
        config["data"]["ref"]["gchp"]["outputs_subdir"],
    )
    gchp_vs_gchp_devdir = os.path.join(
        config["paths"]["main_dir"],
        config["data"]["dev"]["gchp"]["dir"],
        config["data"]["dev"]["gchp"]["outputs_subdir"],
    )

    # Restart file directory paths
    gcc_vs_gcc_refrst = os.path.join(
        config["paths"]["main_dir"],
        config["data"]["ref"]["gcc"]["dir"],
        config["data"]["ref"]["gcc"]["restarts_subdir"]
    )
    gcc_vs_gcc_devrst = os.path.join(
        config["paths"]["main_dir"],
        config["data"]["dev"]["gcc"]["dir"],
        config["data"]["dev"]["gcc"]["restarts_subdir"]
    )
    gchp_vs_gcc_refrst = os.path.join(
        config["paths"]["main_dir"],
        config["data"]["dev"]["gcc"]["dir"],
        config["data"]["dev"]["gcc"]["restarts_subdir"]
    )
    gchp_vs_gcc_devrst = os.path.join(
        config["paths"]["main_dir"],
        config["data"]["dev"]["gchp"]["dir"],
        config["data"]["dev"]["gchp"]["restarts_subdir"]
    )
    gchp_vs_gchp_refrst = os.path.join(
        config["paths"]["main_dir"],
        config["data"]["ref"]["gchp"]["dir"],
        config["data"]["ref"]["gchp"]["restarts_subdir"]
    )
    gchp_vs_gchp_devrst = os.path.join(
        config["paths"]["main_dir"],
        config["data"]["dev"]["gchp"]["dir"],
        config["data"]["dev"]["gchp"]["restarts_subdir"]
    )

    # =====================================================================
    # Benchmark output directories
    # =====================================================================
    # Results directories
    if config["options"]["gcpy_test"]:
        mainresultsdir = os.path.join(".", config["paths"]["results_dir"])
        gcc_vs_gcc_resultsdir = os.path.join(
            mainresultsdir,
            config["options"]["comparisons"]["gcc_vs_gcc"]["dir"]
        )
        gchp_vs_gchp_resultsdir = os.path.join(
            mainresultsdir,
            config["options"]["comparisons"]["gchp_vs_gchp"]["dir"]
        )
        gchp_vs_gcc_resultsdir = os.path.join(
            mainresultsdir,
            "GCHP_GCC_comparison")
        diff_of_diffs_resultsdir = os.path.join(
            mainresultsdir,
            "GCHP_GCC_diff_of_diffs"
        )
        if not os.path.exists(mainresultsdir):
            os.mkdir(mainresultsdir)
        # Make copy of benchmark script in results directory
        curfile = os.path.realpath(__file__)
        dest = os.path.join(mainresultsdir, curfile.split("/")[-1])
        if not os.path.exists(dest):
            copyfile(curfile, dest)

    else:
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
        for resdir, plotting_type in zip(
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
                if resdir in [gcc_vs_gcc_resultsdir, base_gchp_resultsdir]:
                    # Make copy of benchmark script in results directory
                    curfile = os.path.realpath(__file__)
                    dest = os.path.join(resdir, curfile.split("/")[-1])
                    if os.path.exists(dest):
                        copyfile(curfile, dest)

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
    diff_of_diffs_refstr = (
        config["data"]["dev"]["gcc"]["version"]
        + " - "
        + config["data"]["ref"]["gcc"]["version"]
    )
    diff_of_diffs_devstr = (
        config["data"]["dev"]["gchp"]["version"]
        + " - "
        + config["data"]["ref"]["gchp"]["version"]
    )

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

    # TODO: remove is_pre_13.1 option with 14.0 release
    if config["data"]["ref"]["gchp"]["is_pre_13.1"]:
        if add_months(gchp_ref_date, 1) == gchp_end_ref_date:
            gchp_ref_date = np.datetime(
                config["data"]["ref"]["gchp"]["bmk_start"][0:8] + "16T12:00:00"
            )
        else:
            print("Error: `is_pre_13.1: True` option only supported for exactly 1 month and 1 year benchmarks")
            sys.exit()

    # Dev start used in diagnostic filename
    gcc_dev_date = np.datetime64(config["data"]["dev"]["gcc"]["bmk_start"])
    gchp_dev_date = np.datetime64(config["data"]["dev"]["gchp"]["bmk_start"])
    # Dev end used in restart filename
    gcc_end_dev_date = np.datetime64(config["data"]["dev"]["gcc"]["bmk_end"])
    gchp_end_dev_date = np.datetime64(config["data"]["dev"]["gchp"]["bmk_end"])

    # TODO: remove is_pre_13.1 option with 14.0 release
    if config["data"]["dev"]["gchp"]["is_pre_13.1"]:
        if add_months(gchp_dev_date, 1) == gchp_end_dev_date:
            gchp_dev_date = np.datetime(
                config["data"]["dev"]["gchp"]["bmk_start"][0:8] + "16T12:00:00"
            )
        else:
            print("Error: `is_pre_13.1: True` option only supported for exactly 1 month and 1 year benchmarks")
            sys.exit()


    # Seconds per month
    gcc_ref_sec_diff = (gcc_end_ref_date - gcc_ref_date).astype("float64")
    gchp_ref_sec_diff = (gchp_end_ref_date - gchp_ref_date).astype("float64")
    gcc_dev_sec_diff = (gcc_end_dev_date - gcc_dev_date).astype("float64")
    gchp_dev_sec_diff = (gchp_end_dev_date - gchp_dev_date).astype("float64")

    # Double gchp sec/month if mid-point timestamp in filename (legacy format)
    if config["data"]["ref"]["gchp"]["is_pre_13.1"]:
        gchp_ref_sec_diff = gchp_ref_sec_diff * 2
    if config["data"]["dev"]["gchp"]["is_pre_13.1"]:
        gchp_dev_sec_diff = gchp_dev_sec_diff * 2

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
    tmpstr = config["options"]["bmk_type"]
    print(
        f"The following plots and tables will be created for {tmpstr}"
    )
    if config["options"]["outputs"]["plot_conc"]:
        print(" - Concentration plots")
    if config["options"]["outputs"]["plot_emis"]:
        print(" - Emissions plots")
    if config["options"]["outputs"]["plot_jvalues"]:
        print(" - J-values (photolysis rates) plots")
    if config["options"]["outputs"]["plot_aod"]:
        print(" - Aerosol optical depth plots")
    if config["options"]["outputs"]["ops_budget_table"]:
        print(" - Operations budget tables")
    if config["options"]["outputs"]["emis_table"]:
        print(" - Table of emissions totals by spc and inventory")
    if config["options"]["outputs"]["mass_table"]:
        print(" - Table of species mass")
    if config["options"]["outputs"]["OH_metrics"]:
        print(" - Table of OH metrics")
    if config["options"]["outputs"]["ste_table"]:
        print(" - Table of strat-trop exchange")
    print("Comparisons will be made for the following combinations:")
    if config["options"]["comparisons"]["gcc_vs_gcc"]["run"]:
        print(" - GCC vs GCC")
    if config["options"]["comparisons"]["gchp_vs_gcc"]["run"]:
        print(" - GCHP vs GCC")
    if config["options"]["comparisons"]["gchp_vs_gchp"]["run"]:
        print(" - GCHP vs GCHP")
    if config["options"]["comparisons"]["gchp_vs_gcc_diff_of_diffs"]["run"]:
        print(" - GCHP vs GCC diff of diffs")

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
            bmk.make_benchmark_conc_plots(
                ref,
                gcc_vs_gcc_refstr,
                dev,
                gcc_vs_gcc_devstr,
                refmet=refmet,
                devmet=devmet,
                dst=gcc_vs_gcc_resultsdir,
                weightsdir=config["paths"]["weights_dir"],
                plot_by_spc_cat=config["options"]["outputs"]["plot_options"][
                    "by_spc_cat"
                ],
                overwrite=True,
                sigdiff_files=gcc_vs_gcc_sigdiff,
                spcdb_dir=spcdb_dir,
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
            bmk.make_benchmark_emis_plots(
                ref,
                gcc_vs_gcc_refstr,
                dev,
                gcc_vs_gcc_devstr,
                dst=gcc_vs_gcc_resultsdir,
                weightsdir=config["paths"]["weights_dir"],
                plot_by_spc_cat=config["options"]["outputs"]["plot_options"][
                    "by_spc_cat"
                ],
                plot_by_hco_cat=config["options"]["outputs"]["plot_options"][
                    "by_hco_cat"
                ],
                overwrite=True,
                sigdiff_files=gcc_vs_gcc_sigdiff,
                spcdb_dir=spcdb_dir,
            )

        # ==================================================================
        # GCC vs. GCC tables of emission and inventory totals
        # ==================================================================
        if config["options"]["outputs"]["emis_table"]:
            print("\n%%% Creating GCC vs. GCC emissions/inventory tables %%%")

            # Filepaths
            ref = get_filepath(gcc_vs_gcc_refdir, "Emissions", gcc_ref_date)
            dev = get_filepath(gcc_vs_gcc_devdir, "Emissions", gcc_dev_date)

            # Create tables
            bmk.make_benchmark_emis_tables(
                ref,
                gcc_vs_gcc_refstr,
                dev,
                gcc_vs_gcc_devstr,
                dst=gcc_vs_gcc_resultsdir,
                ref_interval=[gcc_ref_sec_diff],
                dev_interval=[gcc_dev_sec_diff],
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

            # Plot J-values
            bmk.make_benchmark_jvalue_plots(
                ref,
                gcc_vs_gcc_refstr,
                dev,
                gcc_vs_gcc_devstr,
                dst=gcc_vs_gcc_resultsdir,
                weightsdir=config["paths"]["weights_dir"],
                overwrite=True,
                sigdiff_files=gcc_vs_gcc_sigdiff,
                spcdb_dir=spcdb_dir,
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
            bmk.make_benchmark_aod_plots(
                ref,
                gcc_vs_gcc_refstr,
                dev,
                gcc_vs_gcc_devstr,
                dst=gcc_vs_gcc_resultsdir,
                weightsdir=config["paths"]["weights_dir"],
                overwrite=True,
                sigdiff_files=gcc_vs_gcc_sigdiff,
                spcdb_dir=spcdb_dir,
            )

        # ==================================================================
        # GCC vs GCC global mass tables
        # ==================================================================
        if config["options"]["outputs"]["mass_table"]:
            print("\n%%% Creating GCC vs. GCC global mass tables %%%")

            # Filepaths
            ref = get_filepath(gcc_vs_gcc_refrst, "Restart", gcc_end_ref_date)
            dev = get_filepath(gcc_vs_gcc_devrst, "Restart", gcc_end_dev_date)

            # Create tables
            bmk.make_benchmark_mass_tables(
                ref,
                config["data"]["ref"]["gcc"]["version"],
                dev,
                config["data"]["dev"]["gcc"]["version"],
                dst=gcc_vs_gcc_tablesdir,
                overwrite=True,
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
            bmk.make_benchmark_operations_budget(
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

            # Use this for benchmarks prior to GEOS-Chem 13.0.0
            #        # Diagnostic collection files to read
            #        col  = "ConcAfterChem"
            #        ref = get_filepath(gcc_vs_gcc_refdir, col, gcc_ref_date)
            #        dev = get_filepath(gcc_vs_gcc_devdir, col, gcc_dev_date)
            #
            #        # Meteorology data needed for calculations
            #        col = "StateMet"
            #        refmet = get_filepath(gcc_vs_gcc_refdir, col, gcc_ref_date)
            #        devmet = get_filepath(gcc_vs_gcc_devdir, col, gcc_dev_date)
            #
            #        # Print OH metrics
            #        bmk.make_benchmark_oh_metrics(
            #            ref,
            #            refmet,
            #            config["data"]["ref"]["gcc"]["version"],
            #            dev,
            #            devmet,
            #            config["data"]["dev"]["gcc"]["version"],
            #            dst=gcc_vs_gcc_tablesdir,
            #            overwrite=True
            #        )

            # Filepaths
            ref = get_filepath(gcc_vs_gcc_refdir, "Metrics", gcc_ref_date)
            dev = get_filepath(gcc_vs_gcc_devdir, "Metrics", gcc_dev_date)

            # Create table
            oh.make_benchmark_oh_metrics(
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
                ste.make_benchmark_ste_table(
                    config["data"]["dev"]["gcc"]["version"],
                    dev,
                    gcc_dev_date.astype(datetime).year,
                    bmk_type=config["options"]["bmk_type"],
                    dst=gcc_vs_gcc_tablesdir,
                    species=["O3"],
                    overwrite=True,
                    month=gcc_dev_date.astype(datetime).month,
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
            is_gchp=True,
            gchp_is_pre_13_1=config["data"]["dev"]["gchp"]["is_pre_13.1"],
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
                is_gchp=True,
                gchp_is_pre_13_1=config["data"]["dev"]["gchp"]["is_pre_13.1"],
            )

            # Create plots
            bmk.make_benchmark_conc_plots(
                ref,
                gchp_vs_gcc_refstr,
                dev,
                gchp_vs_gcc_devstr,
                refmet=refmet,
                devmet=devmet,
                dst=gchp_vs_gcc_resultsdir,
                weightsdir=config["paths"]["weights_dir"],
                plot_by_spc_cat=config["options"]["outputs"]["plot_options"][
                    "by_spc_cat"
                ],
                overwrite=True,
                sigdiff_files=gchp_vs_gcc_sigdiff,
                spcdb_dir=spcdb_dir,
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
                is_gchp=True,
                gchp_is_pre_13_1=config["data"]["dev"]["gchp"]["is_pre_13.1"],
            )

            # Create emissions plots
            bmk.make_benchmark_emis_plots(
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
                overwrite=True,
                sigdiff_files=gchp_vs_gcc_sigdiff,
                spcdb_dir=spcdb_dir,
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
                is_gchp=True,
                gchp_is_pre_13_1=config["data"]["dev"]["gchp"]["is_pre_13.1"],
            )

            # Create plots
            bmk.make_benchmark_emis_tables(
                ref,
                gchp_vs_gcc_refstr,
                dev,
                gchp_vs_gcc_devstr,
                dst=gchp_vs_gcc_resultsdir,
                ref_interval=[gcc_dev_sec_diff],
                dev_interval=[gchp_dev_sec_diff],
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
                is_gchp=True,
                gchp_is_pre_13_1=config["data"]["dev"]["gchp"]["is_pre_13.1"],
            )

            # Create plots
            bmk.make_benchmark_jvalue_plots(
                ref,
                gchp_vs_gcc_refstr,
                dev,
                gchp_vs_gcc_devstr,
                dst=gchp_vs_gcc_resultsdir,
                weightsdir=config["paths"]["weights_dir"],
                overwrite=True,
                sigdiff_files=gchp_vs_gcc_sigdiff,
                spcdb_dir=spcdb_dir,
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
                is_gchp=True,
                gchp_is_pre_13_1=config["data"]["dev"]["gchp"]["is_pre_13.1"],
            )

            # Create plots
            bmk.make_benchmark_aod_plots(
                ref,
                gchp_vs_gcc_refstr,
                dev,
                gchp_vs_gcc_devstr,
                dst=gchp_vs_gcc_resultsdir,
                weightsdir=config["paths"]["weights_dir"],
                overwrite=True,
                sigdiff_files=gchp_vs_gcc_sigdiff,
                spcdb_dir=spcdb_dir,
            )

        # ==================================================================
        # GCHP vs GCC global mass tables
        # ==================================================================
        if config["options"]["outputs"]["mass_table"]:
            print("\n%%% Creating GCHP vs. GCC global mass tables %%%")

            # Filepaths
            ref = get_filepath(
                gchp_vs_gcc_refrst,
                "Restart",
                gcc_end_dev_date
            )
            dev = get_filepath(
                gchp_vs_gcc_devrst,
                "Restart",
                gchp_end_dev_date,
                is_gchp=True,
                gchp_res=config["data"]["dev"]["gchp"]["resolution"],
                gchp_is_pre_14_0=config["data"]["dev"]["gchp"]["is_pre_14.0"]
            )

            # Create tables
            bmk.make_benchmark_mass_tables(
                ref,
                gchp_vs_gcc_refstr,
                dev,
                gchp_vs_gcc_devstr,
                dst=gchp_vs_gcc_tablesdir,
                overwrite=True,
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
                is_gchp=True,
                gchp_is_pre_13_1=config["data"]["dev"]["gchp"]["is_pre_13.1"],
            )

            # Create plots
            bmk.make_benchmark_operations_budget(
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
                is_gchp=True,
                gchp_is_pre_13_1=config["data"]["dev"]["gchp"]["is_pre_13.1"],
            )

            # Create table
            oh.make_benchmark_oh_metrics(
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
            is_gchp=True,
            gchp_is_pre_13_1=config["data"]["ref"]["gchp"]["is_pre_13.1"],
        )
        devmet = get_filepath(
            gchp_vs_gchp_devdir,
            "StateMet",
            gchp_dev_date,
            is_gchp=True,
            gchp_is_pre_13_1=config["data"]["dev"]["gchp"]["is_pre_13.1"],
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
                is_gchp=True,
                gchp_is_pre_13_1=config["data"]["ref"]["gchp"]["is_pre_13.1"],
            )
            dev = get_filepath(
                gchp_vs_gchp_devdir,
                "SpeciesConc",
                gchp_dev_date,
                is_gchp=True,
                gchp_is_pre_13_1=config["data"]["dev"]["gchp"]["is_pre_13.1"],
            )

            # Create plots
            bmk.make_benchmark_conc_plots(
                ref,
                gchp_vs_gchp_refstr,
                dev,
                gchp_vs_gchp_devstr,
                refmet=refmet,
                devmet=devmet,
                dst=gchp_vs_gchp_resultsdir,
                weightsdir=config["paths"]["weights_dir"],
                plot_by_spc_cat=config["options"]["outputs"]["plot_options"][
                    "by_spc_cat"
                ],
                overwrite=True,
                sigdiff_files=gchp_vs_gchp_sigdiff,
                spcdb_dir=spcdb_dir,
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
                is_gchp=True,
                gchp_is_pre_13_1=config["data"]["ref"]["gchp"]["is_pre_13.1"],
            )
            dev = get_filepath(
                gchp_vs_gchp_devdir,
                "Emissions",
                gchp_dev_date,
                is_gchp=True,
                gchp_is_pre_13_1=config["data"]["dev"]["gchp"]["is_pre_13.1"],
            )

            # Create plots
            bmk.make_benchmark_emis_plots(
                ref,
                gchp_vs_gchp_refstr,
                dev,
                gchp_vs_gchp_devstr,
                dst=gchp_vs_gchp_resultsdir,
                weightsdir=config["paths"]["weights_dir"],
                plot_by_spc_cat=config["options"]["outputs"]["plot_options"][
                    "by_spc_cat"
                ],
                plot_by_hco_cat=config["options"]["outputs"]["plot_options"][
                    "by_hco_cat"
                ],
                overwrite=True,
                sigdiff_files=gchp_vs_gchp_sigdiff,
                spcdb_dir=spcdb_dir,
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
                is_gchp=True,
                gchp_is_pre_13_1=config["data"]["ref"]["gchp"]["is_pre_13.1"],
            )
            dev = get_filepath(
                gchp_vs_gchp_devdir,
                "Emissions",
                gchp_dev_date,
                is_gchp=True,
                gchp_is_pre_13_1=config["data"]["ref"]["gchp"]["is_pre_13.1"],
            )

            # Create tables
            bmk.make_benchmark_emis_tables(
                ref,
                gchp_vs_gchp_refstr,
                dev,
                gchp_vs_gchp_devstr,
                dst=gchp_vs_gchp_resultsdir,
                ref_interval=[gchp_ref_sec_diff],
                dev_interval=[gchp_dev_sec_diff],
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
                is_gchp=True,
                gchp_is_pre_13_1=config["data"]["ref"]["gchp"]["is_pre_13.1"],
            )
            dev = get_filepath(
                gchp_vs_gchp_devdir,
                "JValues",
                gchp_dev_date,
                is_gchp=True,
                gchp_is_pre_13_1=config["data"]["dev"]["gchp"]["is_pre_13.1"],
            )

            # Create plots
            bmk.make_benchmark_jvalue_plots(
                ref,
                gchp_vs_gchp_refstr,
                dev,
                gchp_vs_gchp_devstr,
                dst=gchp_vs_gchp_resultsdir,
                weightsdir=config["paths"]["weights_dir"],
                overwrite=True,
                sigdiff_files=gchp_vs_gchp_sigdiff,
                spcdb_dir=spcdb_dir,
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
                is_gchp=True,
                gchp_is_pre_13_1=config["data"]["ref"]["gchp"]["is_pre_13.1"],
            )
            dev = get_filepath(
                gchp_vs_gchp_devdir,
                "Aerosols",
                gchp_dev_date,
                is_gchp=True,
                gchp_is_pre_13_1=config["data"]["dev"]["gchp"]["is_pre_13.1"],
            )

            # Create plots
            bmk.make_benchmark_aod_plots(
                ref,
                gchp_vs_gchp_refstr,
                dev,
                gchp_vs_gchp_devstr,
                dst=gchp_vs_gchp_resultsdir,
                weightsdir=config["paths"]["weights_dir"],
                overwrite=True,
                sigdiff_files=gchp_vs_gchp_sigdiff,
                spcdb_dir=spcdb_dir,
            )

        # ==================================================================
        # GCHP vs GCHP global mass tables
        # ==================================================================
        if config["options"]["outputs"]["mass_table"]:
            print("\n%%% Creating GCHP vs. GCHP global mass tables %%%")

            # Filepaths
            ref = get_filepath(
                gchp_vs_gchp_refrst,
                "Restart",
                gchp_end_ref_date,
                is_gchp=True,
                gchp_res=config["data"]["ref"]["gchp"]["resolution"],
                gchp_is_pre_14_0=config["data"]["ref"]["gchp"]["is_pre_14.0"]
            )
            dev = get_filepath(
                gchp_vs_gchp_devrst,
                "Restart",
                gchp_end_dev_date,
                is_gchp=True,
                gchp_res=config["data"]["dev"]["gchp"]["resolution"],
                gchp_is_pre_14_0=config["data"]["dev"]["gchp"]["is_pre_14.0"]
            )

            # Create tables
            bmk.make_benchmark_mass_tables(
                ref,
                gchp_vs_gchp_refstr,
                dev,
                gchp_vs_gchp_devstr,
                dst=gchp_vs_gchp_tablesdir,
                overwrite=True,
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
                is_gchp=True,
                gchp_is_pre_13_1=config["data"]["ref"]["gchp"]["is_pre_13.1"],
            )
            dev = get_filepath(
                gchp_vs_gchp_devdir,
                "Budget",
                gchp_dev_date,
                is_gchp=True,
                gchp_is_pre_13_1=config["data"]["dev"]["gchp"]["is_pre_13.1"],
            )

            # Create tables
            bmk.make_benchmark_operations_budget(
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
                is_gchp=True,
                gchp_is_pre_13_1=config["data"]["ref"]["gchp"]["is_pre_13.1"],
            )
            dev = get_filepath(
                gchp_vs_gchp_devdir,
                "Metrics",
                gchp_dev_date,
                is_gchp=True,
                gchp_is_pre_13_1=config["data"]["dev"]["gchp"]["is_pre_13.1"],
            )

            # Create table
            oh.make_benchmark_oh_metrics(
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
        # --------------------------------------------------------------
        if config["options"]["outputs"]["ste_table"]:
            print("\n%%% Skipping GCHP vs. GCHP Strat-Trop Exchange table %%%")

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
                is_gchp=True,
                gchp_is_pre_13_1=config["data"]["ref"]["gchp"]["is_pre_13.1"],
            )
            gchp_dev = get_filepath(
                gchp_vs_gchp_devdir,
                "SpeciesConc",
                gchp_dev_date,
                is_gchp=True,
                gchp_is_pre_13_1=config["data"]["dev"]["gchp"]["is_pre_13.1"],
            )

            # Create diff-of-diff plots for species concentrations
            # NOTE: for simplicity, do not convert aerosols to ug/m3.
            bmk.make_benchmark_conc_plots(
                gcc_ref,
                diff_of_diffs_refstr,
                gchp_ref,
                diff_of_diffs_devstr,
                dst=diff_of_diffs_resultsdir,
                weightsdir=config["paths"]["weights_dir"],
                overwrite=True,
                use_cmap_RdBu=True,
                second_ref=gcc_dev,
                second_dev=gchp_dev,
                cats_in_ugm3=None,
                spcdb_dir=spcdb_dir,
            )


    # ==================================================================
    # Print a message indicating that the benchmarks finished
    # ==================================================================
    print("\n %%%% All requested benchmark plots/tables created! %%%%")


def main():
    """
    Driver program. Determines which benchmark script script to call
    for 1-hour, 1-day, 1-month, or 1-year benchmarks.
    """
    config_filename = sys.argv[1] if len(sys.argv) == 2 else "1mo_benchmark.yml"
    config = read_config_file(config_filename)
    choose_benchmark_type(config)


if __name__ == "__main__":
    main()
