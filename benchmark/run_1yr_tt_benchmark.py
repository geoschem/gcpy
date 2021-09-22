#!/usr/bin/env python
"""
run_1yr_benchmark.py: Driver script for creating benchmark plots and testing
                      gcpy 1-year TransportTracers benchmark capability.

Run this script to generate benchmark comparisons between:

    (1) GCC (aka GEOS-Chem "Classic") vs. GCC
    (2) GCHP vs GCC (not yet tested)
    (3) GCHP vs GCHP (not yet tested)

You can customize this by editing the settings in the corresponding yaml 
config file (eg. 1yr_tt_benchmark.yml).

Calling sequence:

    ./run_1yr_tt_benchmark.py <path-to-configuration-file>

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

# ======================================================================
# Imports and global settings (you should not need to edit these)
# ======================================================================

import os
import sys
import warnings
from os.path import join, exists
from shutil import copyfile
from calendar import monthrange
import numpy as np
from gcpy.util import get_filepaths, read_config_file
from gcpy import benchmark as bmk
import gcpy.budget_tt as ttbdg
import gcpy.ste_flux as ste

# Tell matplotlib not to look for an X-window
os.environ["QT_QPA_PLATFORM"] = "offscreen"

# Suppress annoying warning messages
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)


def gchp_metname(prior_to_13):
    """
    Deterimines the correct collection name for GCHP StateMet data.
    """
    if prior_to_13:
        return "StateMet_avg"
    return "StateMet"


def run_benchmark(config):
    # This script has a fixed benchmark type, year, and months
    bmk_type = "TransportTracersBenchmark"
    bmk_year_ref = "2019"
    bmk_year_dev = "2019"
    bmk_mon_strs = ["Jan", "Apr", "Jul", "Oct"]
    bmk_mon_inds = [0, 3, 6, 9]
    bmk_n_months = len(bmk_mon_strs)

    ########################################################################
    ###           CONFIGURABLE SETTINGS: ***EDIT AS NEEDED***            ###
    ########################################################################

    # ======================================================================
    # Benchmark information
    # Note: When doing GCHP vs GCC comparisions gchp_dev will be compared
    # to gcc_dev (not gcc_ref!).
    # ======================================================================

    # Path to species_databse.yml
    spcdb_dir = join(
        config["paths"]["main_dir"], config["data"]["dev"]["gcc"]["version"]
    )

    # ======================================================================
    # Data directories
    # For gchp_vs_gcc_refdir use config["data"]["dev"]["gcc"]["version"], not ref (mps, 6/27/19)
    # ======================================================================

    # Diagnostic file directory paths
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

    # Plots directories
    if config["options"]["gcpy_test"]:
        mainresultsdir = join(".", config["paths"]["results_dir"])
        gcc_vs_gcc_resultsdir = join(
            mainresultsdir, config["options"]["comparisons"]["gcc_vs_gcc"]["dir"]
        )
        gchp_vs_gcc_resultsdir = join(
            mainresultsdir, config["options"]["comparisons"]["gchp_vs_gcc"]["dir"]
        )
        gchp_vs_gchp_resultsdir = join(
            mainresultsdir, config["options"]["comparisons"]["gcc_vs_gcc"]["dir"]
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
            config["options"]["comparisons"]["gcc_vs_gcc"]["dir"],
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

    # ======================================================================
    # Dates and times -- Ref data
    # ======================================================================

    # Month/year strings for use in tabl4e subdirectories (e.g. Jan2016)
    bmk_mon_yr_strs_ref = [v + bmk_year_ref for v in bmk_mon_strs]

    # Get all months array of start datetimes for benchmark year
    bmk_start_ref = np.datetime64(bmk_year_ref + "-01-01")
    bmk_end_ref = np.datetime64("{}-01-01".format(int(bmk_year_ref) + 1))
    all_months_ref = np.arange(
        bmk_start_ref, bmk_end_ref, step=np.timedelta64(1, "M"), dtype="datetime64[M]"
    )
    all_months_gchp_ref = all_months_ref

    # Overwrite all_months_gchp_ref if GCHP ref is legacy filename format.
    # Legacy format uses time-averaging period mid-point not start.
    if config["data"]["ref"]["gchp"]["is_legacy"]:
        sec_per_yr_ref = 0
        all_months_gchp_ref = np.zeros(12, dtype="datetime64[h]")
        for t in range(12):
            days_in_mon = monthrange(int(bmk_year_ref), t + 1)[1]
            sec_per_yr_ref += days_in_mon * 86400.0
            middle_hr = int(days_in_mon * 24 / 2)
            delta = np.timedelta64(middle_hr, "h")
            all_months_gchp_ref[t] = all_months_ref[t].astype("datetime64[h]") + delta

    # ======================================================================
    # Dates and times -- Dev data
    # ======================================================================

    # Month/year strings for use in table subdirectories (e.g. Jan2016)
    bmk_mon_yr_strs_dev = [v + bmk_year_dev for v in bmk_mon_strs]

    # Get all months array of start datetimes for benchmark year
    bmk_start_dev = np.datetime64(bmk_year_dev + "-01-01")
    bmk_end_dev = np.datetime64("{}-01-01".format(int(bmk_year_dev) + 1))
    all_months_dev = np.arange(
        bmk_start_dev, bmk_end_dev, step=np.timedelta64(1, "M"), dtype="datetime64[M]"
    )
    all_months_gchp_dev = all_months_dev

    # Overwrite all_months_gchp_ref if GCHP ref is legacy filename format.
    # Legacy format uses time-averaging period mid-point not start.
    if config["data"]["dev"]["gchp"]["is_legacy"]:
        sec_per_yr_dev = 0
        all_months_gchp_dev = np.zeros(12, dtype="datetime64[h]")
        for t in range(12):
            days_in_mon = monthrange(int(bmk_year_dev), t + 1)[1]
            sec_per_yr_dev += days_in_mon * 86400.0
            middle_hr = int(days_in_mon * 24 / 2)
            delta = np.timedelta64(middle_hr, "h")
            all_months_gchp_dev[t] = all_months_dev[t].astype("datetime64[h]") + delta

    # =======================================================================
    # Print the list of plots & tables to the screen
    # =======================================================================
    print("The following plots and tables will be created for {}:".format(bmk_type))
    if config["options"]["outputs"]["plot_conc"]:
        print(" - Concentration plots")
    if config["options"]["outputs"]["plot_wetdep"]:
        print(" - Convective and large-scale wet deposition plots")
    if config["options"]["outputs"]["rnpbbe_budget"]:
        print(" - Radionuclides budget table")
    if config["options"]["outputs"]["operations_budget"]:
        print(" - Operations budget table")
    if config["options"]["outputs"]["ste_table"]:
        print(" - Table of strat-trop exchange")
    if config["options"]["outputs"]["cons_table"]:
        print(" - Table of mass conservation")
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
        # ==================================================================
        if config["options"]["outputs"]["plot_conc"]:
            print("\n%%% Creating GCC vs. GCC concentration plots %%%")

            # Only plot concentration categories for TransportTracers
            restrict_cats = ["RnPbBeTracers", "PassiveTracers"]

            # --------------------------------------------------------------
            # GCC vs GCC species concentration plots: Annual mean
            # --------------------------------------------------------------

            # Filepaths
            ref = get_filepaths(gcc_vs_gcc_refdir, "SpeciesConc", all_months_ref)[0]
            dev = get_filepaths(gcc_vs_gcc_devdir, "SpeciesConc", all_months_dev)[0]

            # Create plots
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
                restrict_cats=restrict_cats,
                overwrite=True,
                spcdb_dir=spcdb_dir,
            )

            # --------------------------------------------------------------
            # GCC vs GCC species concentration plots: Seasonal
            # --------------------------------------------------------------
            for t in range(bmk_n_months):
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
                    restrict_cats=restrict_cats,
                    overwrite=True,
                    spcdb_dir=spcdb_dir,
                )

        # ==================================================================
        # GCC vs GCC wet deposition plots
        # ==================================================================
        if config["options"]["outputs"]["plot_wetdep"]:
            print("\n%%% Creating GCC vs. GCC wet deposition plots %%%")

            # Diagnostic collection files to read
            cols = ["WetLossConv", "WetLossLS"]

            # Loop over collections
            for col in cols:

                # ----------------------------------------------------------
                # GCC vs. GCC wet deposition plots: Annual mean
                # ----------------------------------------------------------

                # Filepaths
                ref = get_filepaths(gcc_vs_gcc_refdir, col, all_months_ref)[0]
                dev = get_filepaths(gcc_vs_gcc_devdir, col, all_months_dev)[0]

                # Create plots
                bmk.make_benchmark_wetdep_plots(
                    ref,
                    gcc_vs_gcc_refstr,
                    dev,
                    gcc_vs_gcc_devstr,
                    refmet=refmet,
                    devmet=devmet,
                    dst=gcc_vs_gcc_resultsdir,
                    datestr="AnnualMean",
                    time_mean=True,
                    weightsdir=config["paths"]["weights_dir"],
                    benchmark_type=bmk_type,
                    collection=col,
                    overwrite=True,
                    spcdb_dir=spcdb_dir,
                )

                # ----------------------------------------------------------
                # GCC vs GCC wet deposition plots: Seasonal
                # ----------------------------------------------------------
                for t in range(bmk_n_months):
                    mon_ind = bmk_mon_inds[t]
                    bmk.make_benchmark_wetdep_plots(
                        ref[mon_ind],
                        gcc_vs_gcc_refstr,
                        dev[mon_ind],
                        gcc_vs_gcc_devstr,
                        refmet=refmet[mon_ind],
                        devmet=devmet[mon_ind],
                        dst=gcc_vs_gcc_resultsdir,
                        datestr=bmk_mon_yr_strs_dev[t],
                        weightsdir=config["paths"]["weights_dir"],
                        benchmark_type=bmk_type,
                        collection=col,
                        overwrite=True,
                        spcdb_dir=spcdb_dir,
                    )

        # ==================================================================
        # GCC vs GCC radionuclides budget tables
        # ==================================================================
        if config["options"]["outputs"]["rnpbbe_budget"]:
            print("\n%%% Creating GCC vs. GCC radionuclides budget table %%%")
            ttbdg.transport_tracers_budgets(
                config["data"]["dev"]["gcc"]["version"],
                gcc_vs_gcc_devdir,
                gcc_vs_gcc_devrstdir,
                int(bmk_year_dev),
                dst=gcc_vs_gcc_tablesdir,
                overwrite=True,
                spcdb_dir=spcdb_dir,
            )

        # ==================================================================
        # GCC vs GCC operations budgets tables
        # ==================================================================
        if config["options"]["outputs"]["operations_budget"]:
            print("\n%%% Creating GCC vs. GCC operations budget tables %%%")

            # Filepaths
            refs = get_filepaths(gcc_vs_gcc_refdir, "Budget", all_months_ref)[0]
            devs = get_filepaths(gcc_vs_gcc_devdir, "Budget", all_months_dev)[0]

            # Create table
            bmk.make_benchmark_operations_budget(
                config["data"]["ref"]["gcc"]["version"],
                refs,
                config["data"]["dev"]["gcc"]["version"],
                devs,
                sec_per_yr_ref,
                sec_per_yr_dev,
                benchmark_type=bmk_type,
                label=bmk_year_dev,
                operations=[
                    "Chemistry",
                    "Convection",
                    "EmisDryDep",
                    "Mixing",
                    "WetDep",
                ],
                compute_accum=False,
                dst=gcc_vs_gcc_tablesdir,
            )

        # ==================================================================
        # GCC dev strat-trop exchange table
        # ==================================================================
        if config["options"]["outputs"]["ste_table"]:
            print("\n%%% Creating GCC vs. GCC Strat-Trop Exchange table %%%")

            # Diagnostic collection files to read (all 12 months)
            devs = get_filepaths(gcc_vs_gcc_devdir, "AdvFluxVert", all_months_dev)[0]

            # Make stat-trop exchange table for subset of species
            ste.make_benchmark_ste_table(
                config["data"]["dev"]["gcc"]["version"],
                devs,
                int(bmk_year_dev),
                dst=gcc_vs_gcc_tablesdir,
                bmk_type=bmk_type,
                species=["Pb210", "Be7", "Be10"],
                overwrite=True,
            )

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Create GCHP vs GCC benchmark plots and tables
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if config["options"]["comparisons"]["gchp_vs_gcc"]["run"]:

        # ==================================================================
        # GCHP vs GCC filepaths for StateMet collection data
        #
        # (1) GCHP (Dev) and GCC (Ref) use the same benchmark year.
        # (2) The GCC version in "GCHP vs GCC" is the Dev of "GCC vs GCC".
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
        # GCHP vs GCC species concentration plots
        # ==================================================================
        if config["options"]["outputs"]["plot_conc"]:
            print("\n%%% Creating GCHP vs. GCC concentration plots %%%")

            # Only plot concentration categories for TransportTracers
            restrict_cats = ["RnPbBeTracers", "PassiveTracers"]

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
                restrict_cats=restrict_cats,
                overwrite=True,
                spcdb_dir=spcdb_dir,
            )

            # --------------------------------------------------------------
            # GCHP vs GCC species concentration plots: Seasonal
            # --------------------------------------------------------------
            for t in range(bmk_n_months):
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
                    restrict_cats=restrict_cats,
                    overwrite=True,
                    spcdb_dir=spcdb_dir,
                )

        # ==================================================================
        # GCHP vs GCC wet deposition plots
        # ==================================================================
        if config["options"]["outputs"]["plot_wetdep"]:
            print("\n%%% Creating GCHP vs. GCC wet deposition plots %%%")

            # Create separate set of plots for each wetdep collection
            cols = ["WetLossConv", "WetLossLS"]

            # Create plots for each collection and benchmark month
            for col in cols:

                # ----------------------------------------------------------
                # GCHP vs GCC wet deposition plots: Annual mean
                # ----------------------------------------------------------

                # Filepaths
                ref = get_filepaths(gchp_vs_gcc_refdir, col, all_months_dev)[0]
                dev = get_filepaths(
                    gchp_vs_gcc_devdir,
                    col,
                    all_months_gchp_dev,
                    is_gchp=True,
                    gchp_format_is_legacy=config["data"]["dev"]["gchp"]["is_legacy"],
                )[0]

                # Create plots
                bmk.make_benchmark_wetdep_plots(
                    ref,
                    gchp_vs_gcc_refstr,
                    dev,
                    gchp_vs_gcc_devstr,
                    devmet=devmet,
                    collection=col,
                    dst=gchp_vs_gcc_resultsdir,
                    datestr="AnnualMean",
                    time_mean=True,
                    weightsdir=config["paths"]["weights_dir"],
                    overwrite=True,
                    benchmark_type=bmk_type,
                    normalize_by_area=True,
                    spcdb_dir=spcdb_dir,
                )

                # ----------------------------------------------------------
                # GCHP vs GCC wet deposition plots: Seasonal
                # ----------------------------------------------------------
                for t in range(bmk_n_months):
                    mon_ind = bmk_mon_inds[t]
                    bmk.make_benchmark_wetdep_plots(
                        ref[mon_ind],
                        gchp_vs_gcc_refstr,
                        dev[mon_ind],
                        gchp_vs_gcc_devstr,
                        refmet=refmet[mon_ind],
                        devmet=devmet[mon_ind],
                        collection=col,
                        dst=gchp_vs_gcc_resultsdir,
                        datestr=bmk_mon_yr_strs_dev[t],
                        weightsdir=config["paths"]["weights_dir"],
                        overwrite=True,
                        benchmark_type=bmk_type,
                        normalize_by_area=True,
                        spcdb_dir=spcdb_dir,
                    )

        # ==================================================================
        # GCHP vs GCC radionuclides budget tables
        # ==================================================================
        if config["options"]["outputs"]["rnpbbe_budget"]:
            print("\n%%% Creating GCHP vs. GCC radionuclides budget table %%%")
            ttbdg.transport_tracers_budgets(
                config["data"]["dev"]["gchp"]["version"],
                gchp_vs_gcc_devdir,
                gchp_vs_gcc_devrstdir,
                int(bmk_year_dev),
                dst=gchp_vs_gcc_tablesdir,
                is_gchp=True,
                overwrite=True,
                spcdb_dir=spcdb_dir,
            )

        # ==================================================================
        # GCHP vs GCC operations budgets tables
        # ==================================================================
        if config["options"]["outputs"]["operations_budget"]:
            print("\n%%% Creating GCHP vs. GCC operations budget tables %%%")

            # Filepaths
            col = "Budget"
            refs = get_filepaths(gchp_vs_gcc_refdir, "Budget", all_months_dev)[0]
            devs = get_filepaths(
                gchp_vs_gcc_devdir,
                col,
                all_months_gchp_dev,
                is_gchp=True,
                gchp_format_is_legacy=config["data"]["dev"]["gchp"]["is_legacy"],
            )[0]

            # Make operations budget table
            bmk.make_benchmark_operations_budget(
                config["data"]["dev"]["gcc"]["version"],
                refs,
                config["data"]["dev"]["gchp"]["version"],
                devs,
                sec_per_yr_ref,
                sec_per_yr_dev,
                benchmark_type=bmk_type,
                label=bmk_year_dev,
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

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Create GCHP vs GCHP benchmark plots and tables
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if config["options"]["comparisons"]["gchp_vs_gchp"]["run"]:

        # ==================================================================
        # GCHP vs GCHP filepaths for StateMet collection data
        # ==================================================================
        refmet = get_filepaths(
            gchp_vs_gchp_refdir,
            gchp_metname(config["data"]["ref"]["gchp"]["prior_to_13"]),
            all_months_gchp_ref,
            is_gchp=True,
            gchp_format_is_legacy=config["data"]["ref"]["gchp"]["is_legacy"],
        )[0]
        devmet = get_filepaths(
            gchp_vs_gchp_devdir,
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

            # Only plot concentration categories for TransportTracers
            restrict_cats = ["RnPbBeTracers", "PassiveTracers"]

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

            # Make concentration plots
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
                restrict_cats=restrict_cats,
                overwrite=True,
                spcdb_dir=spcdb_dir,
            )

            # --------------------------------------------------------------
            # GCHP vs GCHP species concentration plots: Seasonal
            # --------------------------------------------------------------
            for t in range(bmk_n_months):
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
                    restrict_cats=restrict_cats,
                    overwrite=True,
                    spcdb_dir=spcdb_dir,
                )

        # ==================================================================
        # GCHP vs GCHP wet deposition plots
        # ==================================================================
        if config["options"]["outputs"]["plot_wetdep"]:
            print("\n%%% Creating GCHP vs. GCHP wet deposition plots %%%")

            # Create separate set of plots for each wetdep collection
            cols = ["WetLossConv", "WetLossLS"]

            # Create plots for each collection and benchmark month
            for col in cols:

                # ----------------------------------------------------------
                # GCHP vs GCHP wet deposition plots: Annual Mean
                # ----------------------------------------------------------

                # Filepaths
                ref = get_filepaths(
                    gchp_vs_gchp_refdir,
                    col,
                    all_months_gchp_ref,
                    is_gchp=True,
                    gchp_format_is_legacy=config["data"]["ref"]["gchp"]["is_legacy"],
                )[0]
                dev = get_filepaths(
                    gchp_vs_gchp_devdir,
                    col,
                    all_months_gchp_dev,
                    is_gchp=True,
                    gchp_format_is_legacy=config["data"]["dev"]["gchp"]["is_legacy"],
                )[0]

                # Create plots
                bmk.make_benchmark_wetdep_plots(
                    ref,
                    gchp_vs_gchp_refstr,
                    dev,
                    gchp_vs_gchp_devstr,
                    refmet=refmet,
                    devmet=devmet,
                    collection=col,
                    dst=gchp_vs_gchp_resultsdir,
                    datestr="AnnualMean",
                    time_mean=True,
                    weightsdir=config["paths"]["weights_dir"],
                    overwrite=True,
                    benchmark_type=bmk_type,
                    normalize_by_area=True,
                    spcdb_dir=spcdb_dir,
                )

                # ----------------------------------------------------------
                # GCHP vs GCHP wet deposition plots: Seasonal
                # ----------------------------------------------------------
                for t in range(bmk_n_months):
                    mon_ind = bmk_mon_inds[t]
                    bmk.make_benchmark_wetdep_plots(
                        ref[mon_ind],
                        gchp_vs_gchp_refstr,
                        dev[mon_ind],
                        gchp_vs_gchp_devstr,
                        refmet=refmet[mon_ind],
                        devmet=devmet[mon_ind],
                        collection=col,
                        dst=gchp_vs_gchp_resultsdir,
                        datestr=bmk_mon_yr_strs_dev[t],
                        weightsdir=config["paths"]["weights_dir"],
                        overwrite=True,
                        benchmark_type=bmk_type,
                        normalize_by_area=True,
                        spcdb_dir=spcdb_dir,
                    )

        # ==================================================================
        # GCHP vs GCHP radionuclides budget table
        # ==================================================================
        if config["options"]["outputs"]["rnpbbe_budget"]:
            print("\n%%% Creating GCHP vs. GCHP radionuclides budget table %%%")
            ttbdg.transport_tracers_budgets(
                config["data"]["dev"]["gchp"]["version"],
                gchp_vs_gchp_devdir,
                gchp_vs_gchp_devrstdir,
                int(bmk_year_dev),
                dst=gchp_vs_gchp_tablesdir,
                is_gchp=True,
                overwrite=True,
                spcdb_dir=spcdb_dir,
            )

        # ==================================================================
        # GCHP vs GCHP operations budgets tables
        # ==================================================================
        if config["options"]["outputs"]["operations_budget"]:
            print("\n%%% Creating GCHP vs. GCHP operations budget tables %%%")

            # Filepaths
            refs = get_filepaths(
                gchp_vs_gchp_refdir,
                "Budget",
                all_months_gchp_ref,
                is_gchp=True,
                gchp_format_is_legacy=config["data"]["ref"]["gchp"]["is_legacy"],
            )[0]
            devs = get_filepaths(
                gchp_vs_gchp_devdir,
                "Budget",
                all_months_gchp_dev,
                is_gchp=True,
                gchp_format_is_legacy=config["data"]["dev"]["gchp"]["is_legacy"],
            )[0]

            # Create table
            bmk.make_benchmark_operations_budget(
                config["data"]["dev"]["gchp"]["version"],
                refs,
                config["data"]["dev"]["gchp"]["version"],
                devs,
                sec_per_yr_ref,
                sec_per_yr_dev,
                benchmark_type=bmk_type,
                label=bmk_year_dev,
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

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
    # Create mass conservations tables for GCC and GCHP
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if config["options"]["outputs"]["cons_table"]:

        # ======================================================================
        # Create mass conservation table for GCC_dev
        # ======================================================================
        if (
            config["options"]["comparisons"]["gcc_vs_gcc"]["run"]
            or config["options"]["comparisons"]["gchp_vs_gcc"]["run"]
        ):
            print("\n%%% Creating GCC dev mass conservation table %%%")

            # Filepaths
            datafiles = get_filepaths(gcc_vs_gcc_devrstdir, "Restart", all_months_dev)[
                0
            ]

            # Pick output folder
            if config["options"]["comparisons"]["gchp_vs_gcc"]["run"]:
                tablesdir = gchp_vs_gcc_tablesdir
            else:
                tablesdir = gcc_vs_gcc_tablesdir

            # Create table
            bmk.make_benchmark_mass_conservation_table(
                datafiles,
                config["data"]["dev"]["gcc"]["version"],
                dst=tablesdir,
                overwrite=True,
                spcdb_dir=spcdb_dir,
            )

        # ======================================================================
        # Create mass conservation table for GCHP_dev
        # ======================================================================
        if (
            config["options"]["comparisons"]["gchp_vs_gcc"]["run"]
            or config["options"]["comparisons"]["gchp_vs_gchp"]["run"]
        ):
            print("\n%%% Creating GCHP dev mass conservation table %%%")

            # Filepaths
            datafiles = get_filepaths(
                gchp_vs_gcc_devrstdir,
                "Restart",
                all_months_dev,
                is_gchp=True,
                gchp_format_is_legacy=config["data"]["dev"]["gchp"]["is_legacy"],
            )[0]

            # Pick output folder
            if config["options"]["comparisons"]["gchp_vs_gcc"]["run"]:
                tablesdir = gchp_vs_gcc_tablesdir
            else:
                tablesdir = gchp_vs_gchp_tablesdir

            # Create table
            bmk.make_benchmark_mass_conservation_table(
                datafiles,
                config["data"]["dev"]["gchp"]["version"],
                dst=tablesdir,
                overwrite=True,
                spcdb_dir=spcdb_dir,
            )


def main():
    """
    Driver for extracting config information and running 1yr tt benchmark

    Args:
        accepts one optional argument pointing to the configuration file. Defaults to benchmarks.yml
    """
    config_filename = sys.argv[1] if len(sys.argv) == 2 else "1yr_tt_benchmark.yml"
    config = read_config_file(config_filename)
    run_benchmark(config)


if __name__ == "__main__":
    main()
