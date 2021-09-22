#!/usr/bin/env python
"""
run_1yr_fullchem_benchmark.py: Driver script for creating benchmark plots and
                               testing gcpy 1-year full-chemistry benchmark
                               capability.

Run this script to generate benchmark comparisons between:

    (1) GCC (aka GEOS-Chem "Classic") vs. GCC
    (2) GCHP vs GCC
    (3) GCHP vs GCHP

You can customize this by editing the settings in the corresponding yaml 
config file (eg. 1yr_fullchem_benchmark.yml).

Calling sequence:

    ./run_1yr_fullchem_benchmark.py <path-to-configuration-file>

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

This script corresponds with GCPy 1.1.0. Edit this version ID if releasing
a new version of GCPy.
"""

# =====================================================================
# Imports and global settings (you should not need to edit these)
# =====================================================================

import os
import sys
from os.path import join, exists
import warnings
from shutil import copyfile
from calendar import monthrange
import numpy as np
from joblib import Parallel, delayed
from gcpy.util import get_filepath, get_filepaths, read_config_file
import gcpy.ste_flux as ste
import gcpy.oh_metrics as oh
import gcpy.budget_ox as ox
from gcpy import benchmark as bmk


# Tell matplotlib not to look for an X-window
os.environ["QT_QPA_PLATFORM"] = "offscreen"
# Suppress annoying warning messages
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

def run_benchmark(config):
    """
    Runs 1 year benchmark with the given configuration settings.

    Args:
        config : dict
            Contains configuration for 1yr benchmark from yaml file.
    """

    # This script has a fixed benchmark type
    bmk_type = "FullChemBenchmark"
    bmk_year_ref = "2019"
    bmk_year_dev = "2019"
    bmk_mon_strs = ["Jan", "Apr", "Jul", "Oct"]
    bmk_mon_inds = [0, 3, 6, 9]
    bmk_n_months = len(bmk_mon_strs)

    ########################################################################
    ###           CONFIGURABLE SETTINGS: ***EDIT AS NEEDED ***           ###
    ########################################################################
    # Path to species_databse.yml
    spcdb_dir = join(
        config["paths"]["main_dir"], config["data"]["dev"]["gcc"]["version"]
    )

    # ======================================================================
    # Data directories
    # For gchp_vs_gcc_refdir use config["data"]["dev"]["gcc"]["version"], not ref (mps, 6/27/19)
    # ======================================================================

    # Diagnostics file directory paths
    gcc_vs_gcc_refdir = join(
        config["paths"]["main_dir"],
        config["data"]["ref"]["gcc"]["version"],
        config["data"]["ref"]["gcc"]["subdir"],
    )
    gcc_vs_gcc_devdir = join(
        config["paths"]["main_dir"],
        config["data"]["dev"]["gcc"]["version"],
        config["data"]["dev"]["gcc"]["subdir"],
    )
    gchp_vs_gcc_refdir = join(
        config["paths"]["main_dir"],
        config["data"]["dev"]["gcc"]["version"],
        config["data"]["dev"]["gcc"]["subdir"],
    )
    gchp_vs_gcc_devdir = join(
        config["paths"]["main_dir"],
        config["data"]["dev"]["gchp"]["version"],
        config["data"]["dev"]["gchp"]["subdir"],
    )
    gchp_vs_gchp_refdir = join(
        config["paths"]["main_dir"],
        config["data"]["ref"]["gchp"]["version"],
        config["data"]["ref"]["gchp"]["subdir"],
    )
    gchp_vs_gchp_devdir = join(
        config["paths"]["main_dir"],
        config["data"]["dev"]["gchp"]["version"],
        config["data"]["dev"]["gchp"]["subdir"],
    )

    # Restart file directory paths
    gcc_vs_gcc_refrstdir = join(
        config["paths"]["main_dir"], config["data"]["ref"]["gcc"]["version"], "restarts"
    )
    gcc_vs_gcc_devrstdir = join(
        config["paths"]["main_dir"], config["data"]["dev"]["gcc"]["version"], "restarts"
    )
    gchp_vs_gcc_refrstdir = join(
        config["paths"]["main_dir"], config["data"]["dev"]["gcc"]["version"], "restarts"
    )
    gchp_vs_gcc_devrstdir = join(
        config["paths"]["main_dir"], config["data"]["dev"]["gchp"]["version"]
    )
    gchp_vs_gchp_refrstdir = join(
        config["paths"]["main_dir"], config["data"]["ref"]["gchp"]["version"]
    )
    gchp_vs_gchp_devrstdir = join(
        config["paths"]["main_dir"], config["data"]["dev"]["gchp"]["version"]
    )

    # Log file directories -- GEOS-Chem "Classic" only
    gcc_vs_gcc_reflogdir = join(
        config["paths"]["main_dir"], config["data"]["ref"]["gcc"]["version"], "logs"
    )
    gcc_vs_gcc_devlogdir = join(
        config["paths"]["main_dir"], config["data"]["dev"]["gcc"]["version"], "logs"
    )

    # ======================================================================
    # Benchmark output directories
    # ======================================================================
    # Plot directories
    if config["options"]["gcpy_test"]:
        mainresultsdir = join(".", config["paths"]["results_dir"])
        gcc_vs_gcc_resultsdir = join(
            mainresultsdir, config["options"]["comparisons"]["gcc_vs_gcc"]["dir"]
        )
        gchp_vs_gchp_resultsdir = join(
            mainresultsdir, config["options"]["comparisons"]["gchp_vs_gchp"]["dir"]
        )
        gchp_vs_gcc_resultsdir = join(
            mainresultsdir, "GCHP_GCC_comparison"
        )
        if not exists(mainresultsdir):
            os.mkdir(mainresultsdir)
        # Make copy of benchmark script in results directory
        curfile = os.path.realpath(__file__)
        dest = join(mainresultsdir, curfile.split("/")[-1])
        if not exists(dest):
            copyfile(curfile, dest)

    else:
        gcc_vs_gcc_resultsdir = join(
            config["paths"]["main_dir"],
            config["data"]["dev"]["gcc"]["version"],
            config["paths"]["results_dir"],
        )
        gchp_vs_gchp_resultsdir = join(
            config["paths"]["main_dir"],
            config["data"]["dev"]["gchp"]["version"],
            config["paths"]["results_dir"],
            config["options"]["comparisons"]["gchp_vs_gchp"]["dir"],
        )
        gchp_vs_gcc_resultsdir = join(
            config["paths"]["main_dir"],
            config["data"]["dev"]["gchp"]["version"],
            config["paths"]["results_dir"],
            config["options"]["comparisons"]["gchp_vs_gcc"]["dir"],
        )
        base_gchp_resultsdir = join(
            config["paths"]["main_dir"],
            config["data"]["dev"]["gchp"]["version"],
            config["paths"]["results_dir"],
        )

        # make results directories that don't exist
        for resdir, plotting_type in zip(
            [
                gcc_vs_gcc_resultsdir,
                base_gchp_resultsdir,
                gchp_vs_gchp_resultsdir,
                gchp_vs_gcc_resultsdir,
            ],
            [
                config["options"]["comparisons"]["gcc_vs_gcc"]["run"],
                config["options"]["comparisons"]["gchp_vs_gcc"]["run"]
                or config["options"]["comparisons"]["gchp_vs_gchp"]["run"],
                config["options"]["comparisons"]["gchp_vs_gchp"]["run"],
                config["options"]["comparisons"]["gchp_vs_gcc"]["run"],
            ],
        ):
            if plotting_type and not exists(resdir):
                os.mkdir(resdir)
            if resdir in [gcc_vs_gcc_resultsdir, base_gchp_resultsdir]:
                # Make copy of benchmark script in results directory
                curfile = os.path.realpath(__file__)
                dest = join(resdir, curfile.split("/")[-1])
                if not exists(dest):
                    copyfile(curfile, dest)

    # Tables directories
    gcc_vs_gcc_tablesdir = join(
        gcc_vs_gcc_resultsdir,
        config["options"]["comparisons"]["gcc_vs_gcc"]["tables_subdir"],
    )
    gchp_vs_gcc_tablesdir = join(
        gchp_vs_gcc_resultsdir,
        config["options"]["comparisons"]["gchp_vs_gcc"]["tables_subdir"],
    )
    gchp_vs_gchp_tablesdir = join(
        gchp_vs_gchp_resultsdir,
        config["options"]["comparisons"]["gchp_vs_gchp"]["tables_subdir"],
    )

    # Budget directories
    gcc_vs_gcc_budgetdir = join(gcc_vs_gcc_resultsdir, "Budget")
    gchp_vs_gcc_budgetdir = join(gchp_vs_gcc_resultsdir, "Budget")
    gchp_vs_gchp_budgetdir = join(gchp_vs_gchp_resultsdir, "Budget")

    # ======================================================================
    # Plot title strings
    # For gchp_vs_gcc_refstr use config["data"]["dev"]["gcc"]["version"], not ref (mps, 6/27/19)
    # ======================================================================
    gcc_vs_gcc_refstr = config["data"]["ref"]["gcc"]["version"]
    gcc_vs_gcc_devstr = config["data"]["dev"]["gcc"]["version"]
    gchp_vs_gcc_refstr = config["data"]["dev"]["gcc"]["version"]
    gchp_vs_gcc_devstr = config["data"]["dev"]["gchp"]["version"]
    gchp_vs_gchp_refstr = config["data"]["ref"]["gchp"]["version"]
    gchp_vs_gchp_devstr = config["data"]["dev"]["gchp"]["version"]

    ########################################################################
    ###    THE REST OF THESE SETTINGS SHOULD NOT NEED TO BE CHANGED      ###
    ########################################################################

    def gchp_metname(prior_to_13):
        """
        Returns the proper name for the GCHP StateMet collection.
        """
        if prior_to_13:
            return "StateMet_avg"
        return "StateMet"

    # =====================================================================
    # Dates and times -- ref data
    # =====================================================================

    # Month/year strings for use in table subdirectories (e.g. Jan2016)
    bmk_mon_yr_strs_ref = [v + bmk_year_ref for v in bmk_mon_strs]

    # Get days per month and seconds per month for ref
    sec_per_month_ref = np.zeros(12)
    days_per_month_ref = np.zeros(12)
    for t in range(12):
        days_per_month_ref[t] = monthrange(int(bmk_year_ref), t + 1)[1]
        sec_per_month_ref[t] = days_per_month_ref[t] * 86400.0

    # Get all months array of start datetimes for benchmark year
    bmk_start_ref = np.datetime64(bmk_year_ref + "-01-01")
    bmk_end_ref = np.datetime64("{}-01-01".format(int(bmk_year_ref) + 1))
    all_months_ref = np.arange(
        bmk_start_ref, bmk_end_ref, step=np.timedelta64(1, "M"), dtype="datetime64[M]"
    )
    all_months_gchp_ref = all_months_ref

    # Reset all months datetime array if GCHP ref is legacy filename format.
    # Legacy format uses time-averaging period mid-point not start.
    if config["data"]["ref"]["gchp"]["is_legacy"]:
        all_months_gchp_ref = np.zeros(12, dtype="datetime64[h]")
        for t in range(12):
            middle_hr = int(days_per_month_ref[t] * 24 / 2)
            delta = np.timedelta64(middle_hr, "h")
            all_months_gchp_ref[t] = all_months_ref[t].astype("datetime64[h]") + delta

    # Get subset of month datetimes and seconds per month for only benchmark months
    bmk_mons_ref = all_months_ref[bmk_mon_inds]
    bmk_mons_gchp_ref = all_months_gchp_ref[bmk_mon_inds]
    bmk_sec_per_month_ref = sec_per_month_ref[bmk_mon_inds]

    # =====================================================================
    # Dates and times -- Dev data
    # =====================================================================

    # Month/year strings for use in table subdirectories (e.g. Jan2016)
    bmk_mon_yr_strs_dev = [v + bmk_year_dev for v in bmk_mon_strs]

    # Get days per month and seconds per month for dev
    sec_per_month_dev = np.zeros(12)
    days_per_month_dev = np.zeros(12)
    for t in range(12):
        days_per_month_dev[t] = monthrange(int(bmk_year_dev), t + 1)[1]
        sec_per_month_dev[t] = days_per_month_dev[t] * 86400.0

    # Get all months array of start datetimes for benchmark year
    bmk_start_dev = np.datetime64(bmk_year_dev + "-01-01")
    bmk_end_dev = np.datetime64("{}-01-01".format(int(bmk_year_dev) + 1))
    all_months_dev = np.arange(
        bmk_start_dev, bmk_end_dev, step=np.timedelta64(1, "M"), dtype="datetime64[M]"
    )
    all_months_gchp_dev = all_months_dev

    # Reset all months datetime array if GCHP dev is legacy filename format.
    # Legacy format uses time-averaging period mid-point not start.
    if config["data"]["dev"]["gchp"]["is_legacy"]:
        all_months_gchp_dev = np.zeros(12, dtype="datetime64[h]")
        for t in range(12):
            middle_hr = int(days_per_month_dev[t] * 24 / 2)
            delta = np.timedelta64(middle_hr, "h")
            all_months_gchp_dev[t] = all_months_dev[t].astype("datetime64[h]") + delta

    # Get subset of month datetimes and seconds per month for only benchmark months
    bmk_mons_dev = all_months_dev[bmk_mon_inds]
    bmk_mons_gchp_dev = all_months_gchp_dev[bmk_mon_inds]
    bmk_sec_per_month_dev = sec_per_month_dev[bmk_mon_inds]

    # ======================================================================
    # Print the list of plots & tables to the screen
    # ======================================================================

    print("The following plots and tables will be created for {}:".format(bmk_type))
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
    if config["options"]["outputs"]["aer_budget_table"]:
        print(" - Aerosol budget/burden tables")
    if config["options"]["outputs"]["emis_table"]:
        print(" - Table of emissions totals by species and inventory")
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

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Create GCC vs GCC benchmark plots and tables
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if config["options"]["comparisons"]["gcc_vs_gcc"]["run"]:

        # ==================================================================
        # GCC vs GCC filepaths for StateMet collection data
        # ==================================================================
        refmet = get_filepaths(gcc_vs_gcc_refdir, "StateMet", all_months_ref)[0]
        devmet = get_filepaths(gcc_vs_gcc_devdir, "StateMet", all_months_dev)[0]

        # ==================================================================
        # GCC vs GCC species concentration plots
        #
        # Includes lumped species and separates by category if plot_by_spc_cat
        # is true; otherwise excludes lumped species and writes to one file.
        # --------------------------------------------------------------
        if config["options"]["outputs"]["plot_conc"]:
            print("\n%%% Creating GCC vs. GCC concentration plots %%%")

            # --------------------------------------------------------------
            # GCC vs GCC species concentration plots: Annual mean
            # --------------------------------------------------------------

            # Filepaths
            ref = get_filepaths(gcc_vs_gcc_refdir, "SpeciesConc", all_months_ref)[0]
            dev = get_filepaths(gcc_vs_gcc_devdir, "SpeciesConc", all_months_dev)[0]

            # Create plots
            print("\nCreating plots for annual mean")
            bmk.make_benchmark_conc_plots(
                ref,
                gcc_vs_gcc_refstr,
                dev,
                gcc_vs_gcc_devstr,
                refmet=refmet,
                devmet=devmet,
                dst=gcc_vs_gcc_resultsdir,
                subdst="AnnualMean",
                time_mean=True,
                weightsdir=config["paths"]["weights_dir"],
                benchmark_type=bmk_type,
                overwrite=True,
                spcdb_dir=spcdb_dir,
            )

            # --------------------------------------------------------------
            # GCC vs GCC species concentration plots: Seasonal
            # --------------------------------------------------------------
            for t in range(bmk_n_months):
                print("\nCreating plots for {}".format(bmk_mon_strs[t]))

                # Create plots
                mon_ind = bmk_mon_inds[t]
                bmk.make_benchmark_conc_plots(
                    ref[mon_ind],
                    gcc_vs_gcc_refstr,
                    dev[mon_ind],
                    gcc_vs_gcc_devstr,
                    refmet=refmet[mon_ind],
                    devmet=devmet[mon_ind],
                    dst=gcc_vs_gcc_resultsdir,
                    subdst=bmk_mon_yr_strs_dev[t],
                    weightsdir=config["paths"]["weights_dir"],
                    benchmark_type=bmk_type,
                    plot_by_spc_cat=config["options"]["outputs"]["plot_options"][
                        "by_spc_cat"
                    ],
                    overwrite=True,
                    spcdb_dir=spcdb_dir,
                )

        # ==================================================================
        # GCC vs GCC emissions plots
        # ==================================================================
        if config["options"]["outputs"]["plot_emis"]:
            print("\n%%% Creating GCC vs. GCC emissions plots %%%")

            # --------------------------------------------------------------
            # GCC vs GCC emissions plots: Annual mean
            # --------------------------------------------------------------

            # Filepaths
            ref = get_filepaths(gcc_vs_gcc_refdir, "Emissions", all_months_ref)[0]
            dev = get_filepaths(gcc_vs_gcc_devdir, "Emissions", all_months_dev)[0]

            # Create plots
            print("\nCreating plots for annual mean")
            bmk.make_benchmark_emis_plots(
                ref,
                gcc_vs_gcc_refstr,
                dev,
                gcc_vs_gcc_devstr,
                dst=gcc_vs_gcc_resultsdir,
                subdst="AnnualMean",
                time_mean=True,
                weightsdir=config["paths"]["weights_dir"],
                plot_by_spc_cat=config["options"]["outputs"]["plot_options"][
                    "by_spc_cat"
                ],
                plot_by_hco_cat=config["options"]["outputs"]["plot_options"][
                    "by_hco_cat"
                ],
                overwrite=True,
                spcdb_dir=spcdb_dir,
            )

            # --------------------------------------------------------------
            # GCC vs GCC emissions plots: Seasonal
            # --------------------------------------------------------------
            for t in range(bmk_n_months):
                print("\nCreating plots for {}".format(bmk_mon_strs[t]))

                # Create plots
                mon_ind = bmk_mon_inds[t]
                bmk.make_benchmark_emis_plots(
                    ref[mon_ind],
                    gcc_vs_gcc_refstr,
                    dev[mon_ind],
                    gcc_vs_gcc_devstr,
                    dst=gcc_vs_gcc_resultsdir,
                    subdst=bmk_mon_yr_strs_dev[t],
                    weightsdir=config["paths"]["weights_dir"],
                    plot_by_spc_cat=config["options"]["outputs"]["plot_options"][
                        "by_spc_cat"
                    ],
                    plot_by_hco_cat=config["options"]["outputs"]["plot_options"][
                        "by_hco_cat"
                    ],
                    overwrite=True,
                    spcdb_dir=spcdb_dir,
                )

        # ==================================================================
        # GCC vs GCC tables of emission and inventory totals
        # ==================================================================
        if config["options"]["outputs"]["emis_table"]:
            print("\n%%% Creating GCC vs. GCC emissions & inventory totals %%%")

            # Filepaths
            ref = get_filepaths(gcc_vs_gcc_refdir, "Emissions", all_months_ref)[0]
            dev = get_filepaths(gcc_vs_gcc_devdir, "Emissions", all_months_dev)[0]

            # Create table
            bmk.make_benchmark_emis_tables(
                ref,
                gcc_vs_gcc_refstr,
                dev,
                gcc_vs_gcc_devstr,
                dst=gcc_vs_gcc_resultsdir,
                ref_interval=sec_per_month_ref,
                dev_interval=sec_per_month_dev,
                overwrite=True,
                spcdb_dir=spcdb_dir,
            )

        # ==================================================================
        # GCC vs GCC J-value plots
        # ==================================================================
        if config["options"]["outputs"]["plot_jvalues"]:
            print("\n%%% Creating GCC vs. GCC J-value plots %%%")

            # --------------------------------------------------------------
            # GCC vs GCC J-value plots: Annual mean
            # --------------------------------------------------------------

            # Filepaths
            ref = get_filepaths(gcc_vs_gcc_refdir, "JValues", all_months_ref)[0]
            dev = get_filepaths(gcc_vs_gcc_devdir, "JValues", all_months_dev)[0]

            # Create plots
            print("\nCreating plots for annual mean")
            bmk.make_benchmark_jvalue_plots(
                ref,
                gcc_vs_gcc_refstr,
                dev,
                gcc_vs_gcc_devstr,
                dst=gcc_vs_gcc_resultsdir,
                subdst="AnnualMean",
                time_mean=True,
                weightsdir=config["paths"]["weights_dir"],
                overwrite=True,
                spcdb_dir=spcdb_dir,
            )

            # --------------------------------------------------------------
            # GCC vs GCC J-value plots: Seasonal
            # --------------------------------------------------------------
            for t in range(bmk_n_months):
                print("\nCreating plots for {}".format(bmk_mon_strs[t]))

                # Create plots
                mon_ind = bmk_mon_inds[t]
                bmk.make_benchmark_jvalue_plots(
                    ref[mon_ind],
                    gcc_vs_gcc_refstr,
                    dev[mon_ind],
                    gcc_vs_gcc_devstr,
                    dst=gcc_vs_gcc_resultsdir,
                    subdst=bmk_mon_yr_strs_dev[t],
                    weightsdir=config["paths"]["weights_dir"],
                    overwrite=True,
                    spcdb_dir=spcdb_dir,
                )

        # ==================================================================
        # GCC vs. GCC column AOD plots
        # ==================================================================
        if config["options"]["outputs"]["plot_aod"]:
            print("\n%%% Creating GCC vs. GCC column AOD plots %%%")

            # --------------------------------------------------------------
            # GCC vs GCC column AOD plots: Annual mean
            # --------------------------------------------------------------

            # Filepaths
            ref = get_filepaths(gcc_vs_gcc_refdir, "Aerosols", all_months_ref)[0]
            dev = get_filepaths(gcc_vs_gcc_devdir, "Aerosols", all_months_dev)[0]

            # Create plots
            print("\nCreating plots for annual mean")
            bmk.make_benchmark_aod_plots(
                ref,
                gcc_vs_gcc_refstr,
                dev,
                gcc_vs_gcc_devstr,
                dst=gcc_vs_gcc_resultsdir,
                subdst="AnnualMean",
                time_mean=True,
                weightsdir=config["paths"]["weights_dir"],
                overwrite=True,
                spcdb_dir=spcdb_dir,
            )

            # --------------------------------------------------------------
            # GCC vs GCC column AOD plots: Seasonal
            # --------------------------------------------------------------
            for t in range(bmk_n_months):
                print("\nCreating plots for {}".format(bmk_mon_strs[t]))

                # Create plots
                mon_ind = bmk_mon_inds[t]
                bmk.make_benchmark_aod_plots(
                    ref[mon_ind],
                    gcc_vs_gcc_refstr,
                    dev[mon_ind],
                    gcc_vs_gcc_devstr,
                    dst=gcc_vs_gcc_resultsdir,
                    subdst=bmk_mon_yr_strs_dev[t],
                    weightsdir=config["paths"]["weights_dir"],
                    overwrite=True,
                    spcdb_dir=spcdb_dir,
                )

        # ==================================================================
        # GCC vs GCC mass tables
        # ==================================================================
        if config["options"]["outputs"]["mass_table"]:
            print("\n%%% Creating GCC vs. GCC mass tables %%%")

            def gcc_vs_gcc_mass_table(m):
                """
                Create mass table for each benchmark month m in parallel
                """

                # Filepaths
                refpath = get_filepath(gcc_vs_gcc_refrstdir, "Restart", bmk_mons_ref[m])
                devpath = get_filepath(gcc_vs_gcc_devrstdir, "Restart", bmk_mons_dev[m])

                # Create tables
                bmk.make_benchmark_mass_tables(
                    refpath,
                    gcc_vs_gcc_refstr,
                    devpath,
                    gcc_vs_gcc_devstr,
                    dst=gcc_vs_gcc_tablesdir,
                    subdst=bmk_mon_yr_strs_dev[m],
                    label="at 01{}".format(bmk_mon_yr_strs_dev[m]),
                    overwrite=True,
                    spcdb_dir=spcdb_dir,
                )

            # Run in parallel
            results = Parallel(n_jobs=-1)(
                delayed(gcc_vs_gcc_mass_table)(t) for t in range(bmk_n_months)
            )

        # ==================================================================
        # GCC vs GCC operations budgets tables
        # ==================================================================
        if config["options"]["outputs"]["ops_budget_table"]:
            print("\n%%% Creating GCC vs. GCC operations budget tables %%%")

            def gcc_vs_gcc_ops_budg(m):
                """
                Create budget table for each benchmark month m in parallel
                """

                # Filepaths
                refpath = get_filepath(gcc_vs_gcc_refdir, "Budget", bmk_mons_ref[m])
                devpath = get_filepath(gcc_vs_gcc_devdir, "Budget", bmk_mons_dev[m])

                # Create tables
                bmk.make_benchmark_operations_budget(
                    config["data"]["ref"]["gcc"]["version"],
                    refpath,
                    config["data"]["dev"]["gcc"]["version"],
                    devpath,
                    sec_per_month_ref[m],
                    sec_per_month_dev[m],
                    benchmark_type=bmk_type,
                    label="at 01{}".format(bmk_mon_yr_strs_dev[m]),
                    dst=gcc_vs_gcc_tablesdir,
                )

            # Run in parallel
            results = Parallel(n_jobs=-1)(
                delayed(gcc_vs_gcc_ops_budg)(t) for t in range(bmk_n_months)
            )

        # ==================================================================
        # GCC vs GCC aerosols budgets/burdens tables
        # ==================================================================
        if config["options"]["outputs"]["aer_budget_table"]:
            print("\n%%% Creating GCC vs. GCC aerosols budget tables %%%")

            # Filepaths
            devaero = get_filepaths(gcc_vs_gcc_devdir, "Aerosols", all_months_dev)[0]
            devspc = get_filepaths(gcc_vs_gcc_devdir, "SpeciesConc", all_months_dev)[0]

            # Compute tables
            bmk.make_benchmark_aerosol_tables(
                gcc_vs_gcc_devdir,
                devaero,
                devspc,
                devmet,
                config["data"]["dev"]["gcc"]["version"],
                bmk_year_dev,
                days_per_month_dev,
                dst=gcc_vs_gcc_tablesdir,
                overwrite=True,
                spcdb_dir=spcdb_dir,
            )

        # ==================================================================
        # GCC vs GCC Ox budget table
        # ==================================================================
        if config["options"]["outputs"]["Ox_budget_table"]:
            print("\n%%% Creating GCC vs. GCC Ox budget table %%%")
            ox.global_ox_budget(
                config["data"]["dev"]["gcc"]["version"],
                gcc_vs_gcc_devdir,
                gcc_vs_gcc_devrstdir,
                bmk_year_dev,
                dst=gcc_vs_gcc_tablesdir,
                overwrite=True,
                spcdb_dir=spcdb_dir,
            )

        # ==================================================================
        # GCC Strat-Trop Exchange
        # ==================================================================
        if config["options"]["outputs"]["ste_table"]:
            print("\n%%% Creating GCC vs. GCC Strat-Trop Exchange table %%%")

            # Filepaths
            dev = get_filepaths(gcc_vs_gcc_devdir, "AdvFluxVert", all_months_dev)[0]

            # Compute table
            ste.make_benchmark_ste_table(
                config["data"]["dev"]["gcc"]["version"],
                dev,
                bmk_year_dev,
                dst=gcc_vs_gcc_tablesdir,
                bmk_type=bmk_type,
                species=["O3"],
                overwrite=True,
            )

        # ==================================================================
        # GCC vs GCC Global mean OH, MCF Lifetime, CH4 Lifetime
        # ==================================================================
        if config["options"]["outputs"]["OH_metrics"]:
            print("\n%%% Creating GCC vs. GCC OH metrics %%%")

            # Filepaths
            ref = get_filepaths(gcc_vs_gcc_refdir, "Metrics", all_months_ref)[0]
            dev = get_filepaths(gcc_vs_gcc_devdir, "Metrics", all_months_dev)[0]

            # Create the OH Metrics table
            oh.make_benchmark_oh_metrics(
                ref,
                config["data"]["ref"]["gcc"]["version"],
                dev,
                config["data"]["dev"]["gcc"]["version"],
                dst=gcc_vs_gcc_tablesdir,
                overwrite=True,
                spcdb_dir=spcdb_dir,
            )

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Create GCHP vs GCC benchmark plots and tables
    #
    # (1) GCHP (Dev) and GCC (Ref) use the same benchmark year.
    # (2) The GCC version in "GCHP vs GCC" is the Dev of "GCC vs GCC".
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if config["options"]["comparisons"]["gchp_vs_gcc"]["run"]:

        # ==================================================================
        # GCHP vs GCC filepaths for StateMet collection data
        # ==================================================================
        refmet = get_filepaths(gchp_vs_gcc_refdir, "StateMet", all_months_dev)[0]
        devmet = get_filepaths(
            gchp_vs_gcc_devdir,
            gchp_metname(config["data"]["dev"]["gchp"]["prior_to_13"]),
            all_months_gchp_dev,
            is_gchp=True,
            gchp_format_is_legacy=config["data"]["dev"]["gchp"]["is_legacy"],
        )[0]

        # ==================================================================
        # GCHP vs GCC Concentration plots
        # ==================================================================
        if config["options"]["outputs"]["plot_conc"]:
            print("\n%%% Creating GCHP vs. GCC concentration plots %%%")

            # --------------------------------------------------------------
            # GCHP vs GCC species concentration plots: Annual Mean
            # --------------------------------------------------------------

            # Filepaths
            ref = get_filepaths(gchp_vs_gcc_refdir, "SpeciesConc", all_months_dev)[0]
            dev = get_filepaths(
                gchp_vs_gcc_devdir,
                "SpeciesConc",
                all_months_gchp_dev,
                is_gchp=True,
                gchp_format_is_legacy=config["data"]["dev"]["gchp"]["is_legacy"],
            )[0]

            # Create plots
            print("\nCreating plots for annual mean")
            bmk.make_benchmark_conc_plots(
                ref,
                gchp_vs_gcc_refstr,
                dev,
                gchp_vs_gcc_devstr,
                refmet=refmet,
                devmet=devmet,
                dst=gchp_vs_gcc_resultsdir,
                subdst="AnnualMean",
                time_mean=True,
                weightsdir=config["paths"]["weights_dir"],
                benchmark_type=bmk_type,
                overwrite=True,
                spcdb_dir=spcdb_dir,
            )

            # --------------------------------------------------------------
            # GCHP vs GCC species concentration plots: Seasonal
            # --------------------------------------------------------------
            for t in range(bmk_n_months):
                print("\nCreating plots for {}".format(bmk_mon_strs[t]))

                # Create plots
                mon_ind = bmk_mon_inds[t]
                bmk.make_benchmark_conc_plots(
                    ref[mon_ind],
                    gchp_vs_gcc_refstr,
                    dev[mon_ind],
                    gchp_vs_gcc_devstr,
                    refmet=refmet[mon_ind],
                    devmet=devmet[mon_ind],
                    dst=gchp_vs_gcc_resultsdir,
                    subdst=bmk_mon_yr_strs_dev[t],
                    weightsdir=config["paths"]["weights_dir"],
                    benchmark_type=bmk_type,
                    plot_by_spc_cat=config["options"]["outputs"]["plot_options"][
                        "by_spc_cat"
                    ],
                    overwrite=True,
                    spcdb_dir=spcdb_dir,
                )

        # ==============================================================
        # GCHP vs. GCC emissions plots
        # ==============================================================
        if config["options"]["outputs"]["plot_emis"]:
            print("\n%%% Creating GCHP vs. GCC emissions plots %%%")

            # --------------------------------------------------------------
            # GCHP vs GCC emissions plots: Annual Mean
            # --------------------------------------------------------------

            # Filepaths
            ref = get_filepaths(gchp_vs_gcc_refdir, "Emissions", all_months_dev)[0]
            dev = get_filepaths(
                gchp_vs_gcc_devdir,
                "Emissions",
                all_months_gchp_dev,
                is_gchp=True,
                gchp_format_is_legacy=config["data"]["dev"]["gchp"]["is_legacy"],
            )[0]

            # Create plots
            print("\nCreating plots for annual mean")
            bmk.make_benchmark_emis_plots(
                ref,
                gchp_vs_gcc_refstr,
                dev,
                gchp_vs_gcc_devstr,
                dst=gchp_vs_gcc_resultsdir,
                subdst="AnnualMean",
                time_mean=True,
                weightsdir=config["paths"]["weights_dir"],
                plot_by_spc_cat=config["options"]["outputs"]["plot_options"][
                    "by_spc_cat"
                ],
                plot_by_hco_cat=config["options"]["outputs"]["plot_options"][
                    "by_hco_cat"
                ],
                overwrite=True,
                spcdb_dir=spcdb_dir,
            )

            # --------------------------------------------------------------
            # GCHP vs GCC emissions plots: Seasonal
            # --------------------------------------------------------------
            for t in range(bmk_n_months):
                print("\nCreating plots for {}".format(bmk_mon_strs[t]))

                # Create plots
                mon_ind = bmk_mon_inds[t]
                bmk.make_benchmark_emis_plots(
                    ref[mon_ind],
                    gchp_vs_gcc_refstr,
                    dev[mon_ind],
                    gchp_vs_gcc_devstr,
                    dst=gchp_vs_gcc_resultsdir,
                    subdst=bmk_mon_yr_strs_dev[t],
                    weightsdir=config["paths"]["weights_dir"],
                    plot_by_spc_cat=config["options"]["outputs"]["plot_options"][
                        "by_spc_cat"
                    ],
                    plot_by_hco_cat=config["options"]["outputs"]["plot_options"][
                        "by_hco_cat"
                    ],
                    overwrite=True,
                    spcdb_dir=spcdb_dir,
                )

        # ==================================================================
        # GCHP vs. GCC tables of emission and inventory totals
        # ==================================================================
        if config["options"]["outputs"]["emis_table"]:
            print("\n%%% Creating GCHP vs. GCC emissions tables %%%")

            # Filepaths
            ref = get_filepaths(gchp_vs_gcc_refdir, "Emissions", all_months_dev)[0]
            dev = get_filepaths(
                gchp_vs_gcc_devdir,
                "Emissions",
                all_months_gchp_dev,
                is_gchp=True,
                gchp_format_is_legacy=config["data"]["dev"]["gchp"]["is_legacy"],
            )[0]

            # Create emissions table that spans entire year
            bmk.make_benchmark_emis_tables(
                ref,
                gchp_vs_gcc_refstr,
                dev,
                gchp_vs_gcc_devstr,
                devmet=devmet,
                dst=gchp_vs_gcc_resultsdir,
                ref_interval=sec_per_month_ref,
                dev_interval=sec_per_month_dev,
                overwrite=True,
                spcdb_dir=spcdb_dir,
            )

        # ==================================================================
        # GCHP vs. GCC J-values plots
        # ==================================================================
        if config["options"]["outputs"]["plot_jvalues"]:
            print("\n%%% Creating GCHP vs. GCC J-values plots %%%")

            # --------------------------------------------------------------
            # GCHP vs GCC J-values plots: Annual Mean
            # --------------------------------------------------------------

            # Filepaths
            ref = get_filepaths(gchp_vs_gcc_refdir, "JValues", all_months_dev)[0]
            dev = get_filepaths(
                gchp_vs_gcc_devdir,
                "JValues",
                all_months_gchp_dev,
                is_gchp=True,
                gchp_format_is_legacy=config["data"]["dev"]["gchp"]["is_legacy"],
            )[0]

            # Create plots
            print("\nCreating plots for annual mean")
            bmk.make_benchmark_jvalue_plots(
                ref,
                gchp_vs_gcc_refstr,
                dev,
                gchp_vs_gcc_devstr,
                dst=gchp_vs_gcc_resultsdir,
                subdst="AnnualMean",
                time_mean=True,
                weightsdir=config["paths"]["weights_dir"],
                overwrite=True,
                spcdb_dir=spcdb_dir,
            )

            # --------------------------------------------------------------
            # GCHP vs GCC J-values plots: Seasonal
            # --------------------------------------------------------------
            for t in range(bmk_n_months):
                print("\nCreating plots for {}".format(bmk_mon_strs[t]))

                # Create plots
                mon_ind = bmk_mon_inds[t]
                bmk.make_benchmark_jvalue_plots(
                    ref[mon_ind],
                    gchp_vs_gcc_refstr,
                    dev[mon_ind],
                    gchp_vs_gcc_devstr,
                    dst=gchp_vs_gcc_resultsdir,
                    subdst=bmk_mon_yr_strs_dev[t],
                    weightsdir=config["paths"]["weights_dir"],
                    overwrite=True,
                    spcdb_dir=spcdb_dir,
                )

        # ==================================================================
        # GCHP vs GCC column AOD plots
        # ==================================================================
        if config["options"]["outputs"]["plot_aod"]:
            print("\n%%% Creating GCHP vs. GCC AOD plots %%%")

            # --------------------------------------------------------------
            # GCHP vs GCC column AOD plots: Annual Mean
            # --------------------------------------------------------------

            # Filepaths
            ref = get_filepaths(gchp_vs_gcc_refdir, "Aerosols", all_months_dev)[0]
            dev = get_filepaths(
                gchp_vs_gcc_devdir,
                "Aerosols",
                all_months_gchp_dev,
                is_gchp=True,
                gchp_format_is_legacy=config["data"]["dev"]["gchp"]["is_legacy"],
            )[0]

            # Create plots
            print("\nCreating plots for annual mean")
            bmk.make_benchmark_aod_plots(
                ref,
                gchp_vs_gcc_refstr,
                dev,
                gchp_vs_gcc_devstr,
                dst=gchp_vs_gcc_resultsdir,
                subdst="AnnualMean",
                time_mean=True,
                weightsdir=config["paths"]["weights_dir"],
                overwrite=True,
                spcdb_dir=spcdb_dir,
            )

            # --------------------------------------------------------------
            # GCHP vs GCC column AOD plots: Seasonal
            # --------------------------------------------------------------
            for t in range(bmk_n_months):
                print("\nCreating plots for {}".format(bmk_mon_strs[t]))

                # Create plots
                mon_ind = bmk_mon_inds[t]
                bmk.make_benchmark_aod_plots(
                    ref[mon_ind],
                    gchp_vs_gcc_refstr,
                    dev[mon_ind],
                    gchp_vs_gcc_devstr,
                    dst=gchp_vs_gcc_resultsdir,
                    subdst=bmk_mon_yr_strs_dev[t],
                    weightsdir=config["paths"]["weights_dir"],
                    overwrite=True,
                    spcdb_dir=spcdb_dir,
                )

        # ==================================================================
        # GCHP vs GCC global mass tables
        # ==================================================================
        if config["options"]["outputs"]["mass_table"]:
            print("\n%%% Creating GCHP vs. GCC mass tables %%%")

            def gchp_vs_gcc_mass_table(m):
                """
                Create mass table for each benchmark month in parallel
                """

                # Filepaths
                refpath = get_filepath(
                    gchp_vs_gcc_refrstdir, "Restart", bmk_mons_dev[m]
                )
                devpath = get_filepath(
                    gchp_vs_gcc_devrstdir,
                    "Restart",
                    bmk_mons_dev[m],
                    is_gchp=True,
                    gchp_format_is_legacy=config["data"]["dev"]["gchp"]["is_legacy"],
                )

                # use initial restart if no checkpoint present (intended for
                # first month).  need to pass path of meteorology file with
                # area variable in this scenario
                dev_extra = ""
                if not os.path.isfile(devpath):
                    devpath = join(
                        gchp_vs_gcc_devrstdir,
                        "initial_GEOSChem_rst."
                        + config["data"]["dev"]["gchp"]["resolution"]
                        + "_benchmark.nc",
                    )
                    dev_extra = get_filepath(
                        gchp_vs_gcc_devrstdir,
                        "Restart",
                        bmk_mons_dev[m + 1],
                        is_gchp=True,
                        gchp_format_is_legacy=config["data"]["dev"]["gchp"][
                            "is_legacy"
                        ],
                    )

                # Create tables
                bmk.make_benchmark_mass_tables(
                    refpath,
                    gchp_vs_gcc_refstr,
                    devpath,
                    gchp_vs_gcc_devstr,
                    dst=gchp_vs_gcc_tablesdir,
                    subdst=bmk_mon_yr_strs_dev[m],
                    label="at 01{}".format(bmk_mon_yr_strs_dev[m]),
                    overwrite=True,
                    spcdb_dir=spcdb_dir,
                    dev_met_extra=dev_extra,
                )

            results = Parallel(n_jobs=-1)(
                delayed(gchp_vs_gcc_mass_table)(t) for t in range(bmk_n_months)
            )

        # ==================================================================
        # GCHP vs GCC operations budgets tables
        # ==================================================================
        if config["options"]["outputs"]["ops_budget_table"]:
            print("\n%%% Creating GCHP vs. GCC operations budget tables %%%")

            def gchp_vs_gcc_ops_budg(m):
                """
                Create operations budgets for each benchmark month m in parallel
                """

                # Filepaths
                refpath = get_filepath(gchp_vs_gcc_refdir, "Budget", bmk_mons_dev[m])
                devpath = get_filepath(
                    gchp_vs_gcc_devdir,
                    "Budget",
                    bmk_mons_gchp_dev[m],
                    is_gchp=True,
                    gchp_format_is_legacy=config["data"]["dev"]["gchp"]["is_legacy"],
                )

                # Create tables
                bmk.make_benchmark_operations_budget(
                    config["data"]["dev"]["gcc"]["version"],
                    refpath,
                    config["data"]["dev"]["gchp"]["version"],
                    devpath,
                    bmk_sec_per_month_dev[m],
                    bmk_sec_per_month_dev[m],
                    benchmark_type=bmk_type,
                    label="at 01{}".format(bmk_mon_yr_strs_dev[m]),
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

            results = Parallel(n_jobs=-1)(
                delayed(gchp_vs_gcc_ops_budg)(t) for t in range(bmk_n_months)
            )

        # ==================================================================
        # GCHP vs GCC aerosol budgets and burdens tables
        # ==================================================================
        if config["options"]["outputs"]["aer_budget_table"]:
            print("\n%%% Creating GCHP vs. GCC aerosol budget tables %%%")

            # Filepaths
            devaero = get_filepaths(
                gchp_vs_gcc_devdir,
                "Aerosols",
                all_months_gchp_dev,
                is_gchp=True,
                gchp_format_is_legacy=config["data"]["dev"]["gchp"]["is_legacy"],
            )[0]
            devspc = get_filepaths(
                gchp_vs_gcc_devdir,
                "SpeciesConc",
                all_months_gchp_dev,
                is_gchp=True,
                gchp_format_is_legacy=config["data"]["dev"]["gchp"]["is_legacy"],
            )[0]

            # Create tables
            bmk.make_benchmark_aerosol_tables(
                gchp_vs_gcc_devdir,
                devaero,
                devspc,
                devmet,
                config["data"]["dev"]["gchp"]["version"],
                bmk_year_dev,
                days_per_month_dev,
                dst=gchp_vs_gcc_tablesdir,
                overwrite=True,
                spcdb_dir=spcdb_dir,
                is_gchp=True,
            )

    # ==================================================================
    # GCHP vs GCC Ox budget tables
    # ==================================================================
    if config["options"]["outputs"]["Ox_budget_table"]:
        print("\n%%% Creating GCHP vs. GCC Ox budget tables %%%")

        # Compute Ox budget table for GCC
        ox.global_ox_budget(
            config["data"]["dev"]["gcc"]["version"],
            gcc_vs_gcc_devdir,
            gcc_vs_gcc_devrstdir,
            bmk_year_dev,
            dst=gcc_vs_gcc_tablesdir,
            overwrite=True,
            spcdb_dir=spcdb_dir
        )

        # Compute Ox budget table for GCHP
        ox.global_ox_budget(
            config["data"]["dev"]["gchp"]["version"],
            gchp_vs_gcc_devdir,
            gchp_vs_gcc_devrstdir,
            bmk_year_dev,
            dst=gchp_vs_gcc_tablesdir,
            overwrite=True,
            is_gchp=True,
            spcdb_dir=spcdb_dir
        )

        # ==================================================================
        # GCHP vs. GCC global mean OH, MCF Lifetime, CH4 Lifetime
        # ==================================================================
        if config["options"]["outputs"]["OH_metrics"]:
            print("\n%%% Creating GCHP vs. GCC OH metrics table %%%")

            # Filepaths
            ref = get_filepaths(gchp_vs_gcc_refdir, "Metrics", all_months_dev)[0]
            dev = get_filepaths(
                gchp_vs_gcc_devdir,
                "Metrics",
                all_months_gchp_dev,
                is_gchp=True,
                gchp_format_is_legacy=config["data"]["dev"]["gchp"]["is_legacy"],
            )[0]

            # Create table
            oh.make_benchmark_oh_metrics(
                ref,
                config["data"]["dev"]["gcc"]["version"],
                dev,
                config["data"]["dev"]["gchp"]["version"],
                dst=gchp_vs_gcc_tablesdir,
                overwrite=True,
                spcdb_dir=spcdb_dir,
            )

        # ==================================================================
        # GCHP Strat-Trop Exchange
        # -=================================================================
        if config["options"]["outputs"]["ste_table"]:
            print("\n%%% Skipping GCHP vs. GCC Strat-Trop Exchange table %%%")

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Create GCHP vs GCHP benchmark plots and tables
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if config["options"]["comparisons"]["gchp_vs_gchp"]["run"]:

        # ==================================================================
        # GCHP vs GCC filepaths for StateMet collection data
        # ==================================================================
        refmet = get_filepaths(
            gchp_vs_gchp_refdir,
            gchp_metname(config["data"]["ref"]["gchp"]["prior_to_13"]),
            all_months_gchp_ref,
            is_gchp=True,
            gchp_format_is_legacy=config["data"]["ref"]["gchp"]["is_legacy"],
        )[0]
        devmet = get_filepaths(
            gchp_vs_gcc_devdir,
            gchp_metname(config["data"]["dev"]["gchp"]["prior_to_13"]),
            all_months_gchp_dev,
            is_gchp=True,
            gchp_format_is_legacy=config["data"]["dev"]["gchp"]["is_legacy"],
        )[0]

        # ==================================================================
        # GCHP vs GCHP species concentration plots
        # ==================================================================
        if config["options"]["outputs"]["plot_conc"]:
            print("\n%%% Creating GCHP vs. GCHP concentration plots %%%")

            # --------------------------------------------------------------
            # GCHP vs GCHP species concentration plots: Annual Mean
            # --------------------------------------------------------------

            # Filepaths
            ref = get_filepaths(
                gchp_vs_gchp_refdir,
                "SpeciesConc",
                all_months_gchp_ref,
                is_gchp=True,
                gchp_format_is_legacy=config["data"]["ref"]["gchp"]["is_legacy"],
            )[0]
            dev = get_filepaths(
                gchp_vs_gchp_devdir,
                "SpeciesConc",
                all_months_gchp_dev,
                is_gchp=True,
                gchp_format_is_legacy=config["data"]["dev"]["gchp"]["is_legacy"],
            )[0]

            # Create plots
            print("\nCreating plots for annual mean")
            bmk.make_benchmark_conc_plots(
                ref,
                gchp_vs_gchp_refstr,
                dev,
                gchp_vs_gchp_devstr,
                refmet=refmet,
                devmet=devmet,
                dst=gchp_vs_gchp_resultsdir,
                subdst="AnnualMean",
                time_mean=True,
                weightsdir=config["paths"]["weights_dir"],
                benchmark_type=bmk_type,
                plot_by_spc_cat=config["options"]["outputs"]["plot_options"][
                    "by_spc_cat"
                ],
                overwrite=True,
                spcdb_dir=spcdb_dir,
            )

            # --------------------------------------------------------------
            # GCHP vs GCHP species concentration plots: Seasonal
            # --------------------------------------------------------------
            for t in range(bmk_n_months):
                print("\nCreating plots for {}".format(bmk_mon_strs[t]))

                # Create plots
                mon_ind = bmk_mon_inds[t]
                bmk.make_benchmark_conc_plots(
                    ref[mon_ind],
                    gchp_vs_gchp_refstr,
                    dev[mon_ind],
                    gchp_vs_gchp_devstr,
                    refmet=refmet[mon_ind],
                    devmet=devmet[mon_ind],
                    dst=gchp_vs_gchp_resultsdir,
                    subdst=bmk_mon_yr_strs_dev[t],
                    weightsdir=config["paths"]["weights_dir"],
                    benchmark_type=bmk_type,
                    plot_by_spc_cat=config["options"]["outputs"]["plot_options"][
                        "by_spc_cat"
                    ],
                    overwrite=True,
                    spcdb_dir=spcdb_dir,
                )

        # ==================================================================
        # GCHP vs. GCHP Emissions plots
        # ==================================================================
        if config["options"]["outputs"]["plot_emis"]:
            print("\n%%% Creating GCHP vs. GCHP emissions plots %%%")

            # --------------------------------------------------------------
            # GCHP vs GCHP species concentration plots: Annual Mean
            # --------------------------------------------------------------

            # Filepaths
            ref = get_filepaths(
                gchp_vs_gchp_refdir,
                "Emissions",
                all_months_gchp_ref,
                is_gchp=True,
                gchp_format_is_legacy=config["data"]["ref"]["gchp"]["is_legacy"],
            )[0]
            dev = get_filepaths(
                gchp_vs_gchp_devdir,
                "Emissions",
                all_months_gchp_dev,
                is_gchp=True,
                gchp_format_is_legacy=config["data"]["dev"]["gchp"]["is_legacy"],
            )[0]

            # Create plots
            print("\nCreating plots for annual mean")
            bmk.make_benchmark_emis_plots(
                ref,
                gchp_vs_gchp_refstr,
                dev,
                gchp_vs_gchp_devstr,
                dst=gchp_vs_gchp_resultsdir,
                subdst="AnnualMean",
                time_mean=True,
                weightsdir=config["paths"]["weights_dir"],
                plot_by_spc_cat=config["options"]["outputs"]["plot_options"][
                    "by_spc_cat"
                ],
                plot_by_hco_cat=config["options"]["outputs"]["plot_options"][
                    "by_hco_cat"
                ],
                overwrite=True,
                spcdb_dir=spcdb_dir,
            )

            # --------------------------------------------------------------
            # GCHP vs GCHP species concentration plots: Seasonal
            # --------------------------------------------------------------
            for t in range(bmk_n_months):
                print("\nCreating plots for {}".format(bmk_mon_strs[t]))

                # Create plots
                mon_ind = bmk_mon_inds[t]
                bmk.make_benchmark_emis_plots(
                    ref[mon_ind],
                    gchp_vs_gchp_refstr,
                    dev[mon_ind],
                    gchp_vs_gchp_devstr,
                    dst=gchp_vs_gchp_resultsdir,
                    subdst=bmk_mon_yr_strs_dev[t],
                    weightsdir=config["paths"]["weights_dir"],
                    plot_by_spc_cat=config["options"]["outputs"]["plot_options"][
                        "by_spc_cat"
                    ],
                    plot_by_hco_cat=config["options"]["outputs"]["plot_options"][
                        "by_hco_cat"
                    ],
                    overwrite=True,
                    spcdb_dir=spcdb_dir,
                )

        # ==================================================================
        # GCHP vs. GCHP tables of emission and inventory totals
        # ==================================================================
        if config["options"]["outputs"]["emis_table"]:
            print("\n%%% Creating GCHP vs. GCHP emissions tables %%%")

            # Filepaths
            ref = get_filepaths(
                gchp_vs_gchp_refdir,
                "Emissions",
                all_months_gchp_ref,
                is_gchp=True,
                gchp_format_is_legacy=config["data"]["ref"]["gchp"]["is_legacy"],
            )[0]
            dev = get_filepaths(
                gchp_vs_gchp_devdir,
                "Emissions",
                all_months_gchp_dev,
                is_gchp=True,
                gchp_format_is_legacy=config["data"]["dev"]["gchp"]["is_legacy"],
            )[0]

            # Create table
            bmk.make_benchmark_emis_tables(
                ref,
                gchp_vs_gchp_refstr,
                dev,
                gchp_vs_gchp_devstr,
                refmet=refmet,
                devmet=devmet,
                dst=gchp_vs_gchp_resultsdir,
                ref_interval=sec_per_month_ref,
                dev_interval=sec_per_month_dev,
                overwrite=True,
                spcdb_dir=spcdb_dir,
            )

        # ==================================================================
        # GCHP vs. GCHP J-values plots
        # ==================================================================
        if config["options"]["outputs"]["plot_jvalues"]:
            print("\n%%% Creating GCHP vs. GCHP J-values plots %%%")

            # --------------------------------------------------------------
            # GCHP vs GCHP J-values plots: Annual Mean
            # --------------------------------------------------------------

            # Filepaths
            ref = get_filepaths(
                gchp_vs_gchp_refdir,
                "JValues",
                all_months_gchp_ref,
                is_gchp=True,
                gchp_format_is_legacy=config["data"]["ref"]["gchp"]["is_legacy"],
            )[0]
            dev = get_filepaths(
                gchp_vs_gchp_devdir,
                "JValues",
                all_months_gchp_dev,
                is_gchp=True,
                gchp_format_is_legacy=config["data"]["dev"]["gchp"]["is_legacy"],
            )[0]

            # Create plots
            print("\nCreating plots for annual mean")
            bmk.make_benchmark_jvalue_plots(
                ref,
                gchp_vs_gchp_refstr,
                dev,
                gchp_vs_gchp_devstr,
                dst=gchp_vs_gchp_resultsdir,
                subdst='AnnualMean',
                time_mean=True,
                weightsdir=config["paths"]["weights_dir"],
                overwrite=True,
                spcdb_dir=spcdb_dir,
            )

            # --------------------------------------------------------------
            # GCHP vs GCHP J-values plots: Seasonal
            # --------------------------------------------------------------
            for t in range(bmk_n_months):
                print("\nCreating plots for {}".format(bmk_mon_strs[t]))

                # Create plots
                mon_ind = bmk_mon_inds[t]
                bmk.make_benchmark_jvalue_plots(
                    ref[mon_ind],
                    gchp_vs_gchp_refstr,
                    dev[mon_ind],
                    gchp_vs_gchp_devstr,
                    dst=gchp_vs_gchp_resultsdir,
                    subdst=bmk_mon_yr_strs_dev[t],
                    weightsdir=config["paths"]["weights_dir"],
                    overwrite=True,
                    spcdb_dir=spcdb_dir,
                )

        # ==================================================================
        # GCHP vs GCHP column AOD plots
        # ==================================================================
        if config["options"]["outputs"]["plot_aod"]:
            print("\n%%% Creating GCHP vs. GCHP AOD plots %%%")

            # --------------------------------------------------------------
            # GCHP vs GCHP column AOD plots: Annual Mean
            # --------------------------------------------------------------

            # Filepaths
            ref = get_filepaths(
                gchp_vs_gchp_refdir,
                "Aerosols",
                all_months_gchp_ref,
                is_gchp=True,
                gchp_format_is_legacy=config["data"]["ref"]["gchp"]["is_legacy"],
            )[0]
            dev = get_filepaths(
                gchp_vs_gchp_devdir,
                "Aerosols",
                all_months_gchp_dev,
                is_gchp=True,
                gchp_format_is_legacy=config["data"]["dev"]["gchp"]["is_legacy"],
            )[0]

            # Create plots
            print("\nCreating plots for annual mean")
            bmk.make_benchmark_aod_plots(
                ref,
                gchp_vs_gchp_refstr,
                dev,
                gchp_vs_gchp_devstr,
                dst=gchp_vs_gchp_resultsdir,
                subdst="AnnualMean",
                time_mean=True,
                weightsdir=config["paths"]["weights_dir"],
                overwrite=True,
                spcdb_dir=spcdb_dir,
            )

            # --------------------------------------------------------------
            # GCHP vs GCHP column AOD plots: Seasonal
            # --------------------------------------------------------------
            for t in range(bmk_n_months):
                print("\nCreating plots for {}".format(bmk_mon_strs[t]))

                # Create plots
                mon_ind = bmk_mon_inds[t]
                bmk.make_benchmark_aod_plots(
                    ref[mon_ind],
                    gchp_vs_gchp_refstr,
                    dev[mon_ind],
                    gchp_vs_gchp_devstr,
                    dst=gchp_vs_gchp_resultsdir,
                    subdst=bmk_mon_yr_strs_dev[t],
                    weightsdir=config["paths"]["weights_dir"],
                    overwrite=True,
                    spcdb_dir=spcdb_dir,
                )

        # ==================================================================
        # GCHP vs GCHP global mass tables
        # ==================================================================
        if config["options"]["outputs"]["mass_table"]:
            print("\n%%% Creating GCHP vs. GCHP mass tables %%%")

            def gchp_vs_gchp_mass_table(m):
                """
                Create mass table for each benchmark month m in parallel
                """

                # Ref filepaths
                refpath = get_filepath(
                    gchp_vs_gchp_refrstdir,
                    "Restart",
                    bmk_mons_ref[m],
                    is_gchp=True,
                    gchp_format_is_legacy=config["data"]["ref"]["gchp"]["is_legacy"],
                )

                # Use initial checkpoint if Ref restart is not present
                ref_extra = ""
                if not os.path.isfile(refpath):
                    refpath = join(
                        gchp_vs_gchp_refrstdir,
                        "initial_GEOSChem_rst."
                        + config["data"]["ref"]["gchp"]["resolution"]
                        + "_benchmark.nc",
                    )
                    ref_extra = get_filepath(
                        gchp_vs_gchp_refrstdir,
                        "Restart",
                        bmk_mons_ref[m + 1],
                        is_gchp=True,
                        gchp_format_is_legacy=config["data"]["ref"]["gchp"][
                            "is_legacy"
                        ],
                    )

                # Dev filepaths
                devpath = get_filepath(
                    gchp_vs_gchp_devrstdir,
                    "Restart",
                    bmk_mons_dev[m],
                    is_gchp=True,
                    gchp_format_is_legacy=config["data"]["dev"]["gchp"]["is_legacy"],
                )

                # Use initial checkpoint if Dev restart is not present
                dev_extra = ""
                if not os.path.isfile(devpath):
                    devpath = join(
                        gchp_vs_gchp_devrstdir,
                        "initial_GEOSChem_rst."
                        + config["data"]["dev"]["gchp"]["resolution"]
                        + "_benchmark.nc",
                    )
                    dev_extra = get_filepath(
                        gchp_vs_gchp_devrstdir,
                        "Restart",
                        bmk_mons_dev[m + 1],
                        is_gchp=True,
                        gchp_format_is_legacy=config["data"]["dev"]["gchp"][
                            "is_legacy"
                        ],
                    )

                # Create tables
                bmk.make_benchmark_mass_tables(
                    refpath,
                    gchp_vs_gchp_refstr,
                    devpath,
                    gchp_vs_gchp_devstr,
                    dst=gchp_vs_gchp_tablesdir,
                    subdst=bmk_mon_yr_strs_dev[m],
                    label="at 01{}".format(bmk_mon_yr_strs_dev[m]),
                    overwrite=True,
                    spcdb_dir=spcdb_dir,
                    ref_met_extra=ref_extra,
                    dev_met_extra=dev_extra,
                )

            # Run in parallel
            results = Parallel(n_jobs=-1)(
                delayed(gchp_vs_gchp_mass_table)(t) for t in range(bmk_n_months)
            )

        # ==================================================================
        # GCHP vs GCHP operations budgets tables
        # ==================================================================
        if config["options"]["outputs"]["ops_budget_table"]:
            print("\n%%% Creating GCHP vs. GCHP operations budget tables %%%")

            # Diagnostic collections to read
            def gchp_vs_gchp_ops_budg(m):
                """
                Creates operations budgets for each benchmark month m in parallel
                """

                # Filepaths
                refpath = get_filepath(
                    gchp_vs_gchp_refdir,
                    "Budget",
                    bmk_mons_gchp_ref[m],
                    is_gchp=True,
                    gchp_format_is_legacy=config["data"]["ref"]["gchp"]["is_legacy"],
                )
                devpath = get_filepath(
                    gchp_vs_gchp_devdir,
                    "Budget",
                    bmk_mons_gchp_dev[m],
                    is_gchp=True,
                    gchp_format_is_legacy=config["data"]["dev"]["gchp"]["is_legacy"],
                )

                # Compute tables
                bmk.make_benchmark_operations_budget(
                    config["data"]["ref"]["gchp"]["version"],
                    refpath,
                    config["data"]["dev"]["gchp"]["version"],
                    devpath,
                    bmk_sec_per_month_ref[m],
                    bmk_sec_per_month_dev[m],
                    benchmark_type=bmk_type,
                    label="at 01{}".format(bmk_mon_yr_strs_dev[m]),
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

            # Run in parallel
            results = Parallel(n_jobs=-1)(
                delayed(gchp_vs_gchp_ops_budg)(t) for t in range(bmk_n_months)
            )

        # ==================================================================
        # GCHP vs GCHP aerosol budgets and burdens tables
        # ==================================================================
        if config["options"]["outputs"]["aer_budget_table"]:
            print("\n%%% Creating GCHP vs. GCHP aerosol budget tables %%%")

            # Filepaths
            devaero = get_filepaths(
                gchp_vs_gchp_devdir,
                "Aerosols",
                all_months_gchp_dev,
                is_gchp=True,
                gchp_format_is_legacy=config["data"]["dev"]["gchp"]["is_legacy"]
            )[0]
            devspc = get_filepaths(
                gchp_vs_gchp_devdir,
                "SpeciesConc",
                all_months_gchp_dev,
                is_gchp=True,
                gchp_format_is_legacy=config["data"]["dev"]["gchp"]["is_legacy"]
            )[0]

            # Create tables
            bmk.make_benchmark_aerosol_tables(
                gchp_vs_gchp_devdir,
                devaero,
                devspc,
                devmet,
                config["data"]["dev"]["gchp"]["version"],
                bmk_year_dev,
                days_per_month_dev,
                dst=gchp_vs_gchp_tablesdir,
                overwrite=True,
                spcdb_dir=spcdb_dir,
                is_gchp=True,
            )

        # ==================================================================
        # GCHP vs GCHP Ox budget tables
        # ==================================================================
        if config["options"]["outputs"]["Ox_budget_table"]:
            print("\n%%% Creating GCHP Ox budget table %%%")

            # Compute Ox budget table for GCHP
            ox.global_ox_budget(
                config["data"]["dev"]["gchp"]["version"],
                gchp_vs_gchp_devdir,
                gchp_vs_gchp_devrstdir,
                bmk_year_dev,
                dst=gchp_vs_gchp_tablesdir,
                overwrite=True,
                is_gchp=True,
                spcdb_dir=spcdb_dir
            )

        # ==================================================================
        # GCHP vs. GCHP global mean OH, MCF Lifetime, CH4 Lifetime
        # ==================================================================
        if config["options"]["outputs"]["OH_metrics"]:
            print("\n%%% Creating GCHP vs. GCHP OH metrics table %%%")

            # Filepaths
            ref = get_filepaths(
                gchp_vs_gchp_refdir,
                "Metrics",
                all_months_gchp_ref,
                is_gchp=True,
                gchp_format_is_legacy=config["data"]["ref"]["gchp"]["is_legacy"],
            )[0]
            dev = get_filepaths(
                gchp_vs_gchp_devdir,
                "Metrics",
                all_months_gchp_dev,
                is_gchp=True,
                gchp_format_is_legacy=config["data"]["dev"]["gchp"]["is_legacy"],
            )[0]

            # Create the OH Metrics table
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
        # GCHP Strat-Trop Exchange
        # ==================================================================
        if config["options"]["outputs"]["ste_table"]:
            print("\n%%% Skipping GCHP vs. GCHP Strat-Trop Exchange table %%%")


def main():
    """
    Driver for extracting config information and running 1yr benchmark

    Args:
        accepts one optional argument pointing to the configuration file. Defaults to benchmarks.yml
    """
    config_filename = sys.argv[1] if len(sys.argv) == 2 else "1yr_fullchem_benchmark.yml"
    config = read_config_file(config_filename)
    run_benchmark(config)


if __name__ == "__main__":
    main()
