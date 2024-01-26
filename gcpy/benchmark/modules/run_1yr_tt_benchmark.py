#!/usr/bin/env python
"""
run_1yr_tt_benchmark.py:
    Script containing code for creating benchmark plots and testing
    gcpy 1-year TransportTracers benchmark capability.

Run this script to generate benchmark comparisons between:

    (1) GCC (aka GEOS-Chem "Classic") vs. GCC
    (2) GCHP vs GCC
    (3) GCHP vs GCHP

You can customize this by editing the settings in the corresponding yaml
config file (eg. 1yr_tt_benchmark.yml).

To generate benchmark output:

    (1) Copy the gcpy/benchmark/run_benchmark.py script and the
        1yr_tt_benchmark.yml file anywhere you want to run the test.

    (2) Edit the 1yr_tt_benchmark.yml to point to the proper file paths
        on your disk space.

    (3) Make sure the /path/to/gcpy/benchmark is in your PYTHONPATH
        shell environment variable.

    (4) Type at the command line

        ./run_benchmark.py 1yr_tt_benchmark.yml

Remarks:

    By default, matplotlib will try to open an X window for plotting.
    If you are running this script in an environment where you do not have
    an active X display (such as in a computational queue), then you will
    need to use these commands to disable the X-window functionality.

        import os
        os.environ["QT_QPA_PLATFORM"]="offscreen"

    For more information, please see this issue posted at the ipython site:

        https://github.com/ipython/ipython/issues/10627

This script corresponds with GCPy 1.4.2. Edit this version ID if releasing
a new version of GCPy.
"""

# ======================================================================
# Imports and global settings (you should not need to edit these)
# ======================================================================

import os
import warnings
from shutil import copyfile
from calendar import monthrange
import numpy as np
from joblib import Parallel, delayed
from gcpy.util import get_filepath, get_filepaths
from gcpy import benchmark_funcs as bmk
import gcpy.budget_tt as ttbdg
import gcpy.ste_flux as ste

# Tell matplotlib not to look for an X-window
os.environ["QT_QPA_PLATFORM"] = "offscreen"

# Suppress annoying warning messages
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)


def run_benchmark(config, bmk_year_ref, bmk_year_dev):
    # This script has a fixed benchmark type, year, and months
    bmk_type = config["options"]["bmk_type"]
    bmk_mon_strs = ["Jan", "Apr", "Jul", "Oct"]
    bmk_mon_inds = [0, 3, 6, 9]
    bmk_n_months = len(bmk_mon_strs)

    # =====================================================================
    # Path to species_database.yml
    # =====================================================================
    spcdb_dir = bmk.get_species_database_dir(config)

    # ======================================================================
    # Data directories
    # For gchp_vs_gcc_refdir use config["data"]["dev"]["gcc"]["version"], not ref (mps, 6/27/19)
    # ======================================================================

    # Diagnostic file directory paths
    gcc_vs_gcc_refdir = os.path.join(
        config["paths"]["main_dir"],
        config["data"]["ref"]["gcc"]["dir"],
        config["data"]["ref"]["gcc"]["outputs_subdir"]
    )
    gcc_vs_gcc_devdir = os.path.join(
        config["paths"]["main_dir"],
        config["data"]["dev"]["gcc"]["dir"],
        config["data"]["dev"]["gcc"]["outputs_subdir"]
    )
    gchp_vs_gcc_refdir = os.path.join(
        config["paths"]["main_dir"],
        config["data"]["dev"]["gcc"]["dir"],
        config["data"]["dev"]["gcc"]["outputs_subdir"]
    )
    gchp_vs_gcc_devdir = os.path.join(
        config["paths"]["main_dir"],
        config["data"]["dev"]["gchp"]["dir"],
        config["data"]["dev"]["gchp"]["outputs_subdir"]
    )
    gchp_vs_gchp_refdir = os.path.join(
        config["paths"]["main_dir"],
        config["data"]["ref"]["gchp"]["dir"],
        config["data"]["ref"]["gchp"]["outputs_subdir"]
    )
    gchp_vs_gchp_devdir = os.path.join(
        config["paths"]["main_dir"],
        config["data"]["dev"]["gchp"]["dir"],
        config["data"]["dev"]["gchp"]["outputs_subdir"]
    )

    # Diagnostic file directory paths
    gcc_vs_gcc_refrstdir = os.path.join(
        config["paths"]["main_dir"],
        config["data"]["ref"]["gcc"]["dir"],
        config["data"]["ref"]["gcc"]["restarts_subdir"]
    )
    gchp_vs_gcc_refrstdir = os.path.join(
        config["paths"]["main_dir"],
        config["data"]["ref"]["gchp"]["dir"],
        config["data"]["ref"]["gchp"]["restarts_subdir"]
    )
    gchp_vs_gchp_refrstdir = os.path.join(
        config["paths"]["main_dir"],
        config["data"]["ref"]["gchp"]["dir"],
        config["data"]["ref"]["gchp"]["restarts_subdir"]
    )
    gcc_vs_gcc_devrstdir = os.path.join(
        config["paths"]["main_dir"],
        config["data"]["dev"]["gcc"]["dir"],
        config["data"]["dev"]["gcc"]["restarts_subdir"]
    )
    gchp_vs_gcc_devrstdir = os.path.join(
        config["paths"]["main_dir"],
        config["data"]["dev"]["gchp"]["dir"],
        config["data"]["dev"]["gchp"]["restarts_subdir"]
    )
    gchp_vs_gchp_devrstdir = os.path.join(
        config["paths"]["main_dir"],
        config["data"]["dev"]["gchp"]["dir"],
        config["data"]["dev"]["gchp"]["restarts_subdir"]
    )

    # Directories where plots & tables will be created
    mainresultsdir = os.path.join(
        config["paths"]["results_dir"]
    )
    gcc_vs_gcc_resultsdir = os.path.join(
        mainresultsdir,
        config["options"]["comparisons"]["gcc_vs_gcc"]["dir"]
    )
    gchp_vs_gcc_resultsdir = os.path.join(
        mainresultsdir,
        config["options"]["comparisons"]["gchp_vs_gcc"]["dir"]
    )
    gchp_vs_gchp_resultsdir = os.path.join(
        mainresultsdir,
        config["options"]["comparisons"]["gchp_vs_gchp"]["dir"]
    )

    # Create the main results directory
    if not os.path.exists(mainresultsdir):
        os.mkdir(mainresultsdir)

    # Make copy of benchmark script in results directory
    curfile = os.path.realpath(__file__)
    dest = os.path.join(mainresultsdir, curfile.split("/")[-1])
    if not os.path.exists(dest):
        copyfile(curfile, dest)

    # Create results directories that don't exist,
    # and place a copy of this file in each results directory
    resdir_list = [
        gcc_vs_gcc_resultsdir,
        gchp_vs_gcc_resultsdir,
        gchp_vs_gchp_resultsdir
    ]
    comparisons_list =  [
        config["options"]["comparisons"]["gcc_vs_gcc"]["run"],
        config["options"]["comparisons"]["gchp_vs_gcc"]["run"],
        config["options"]["comparisons"]["gchp_vs_gchp"]["run"]
    ]
    for resdir, plotting_type in zip(resdir_list, comparisons_list):
        if plotting_type and not os.path.exists(resdir):
            os.mkdir(resdir)
            if resdir in resdir_list:
                curfile = os.path.realpath(__file__)
                dest = os.path.join(resdir, curfile.split("/")[-1])
                if not os.path.exists(dest):
                    copyfile(curfile, dest)

    # Tables directories
    gcc_vs_gcc_tablesdir = os.path.join(
        gcc_vs_gcc_resultsdir,
        config["options"]["comparisons"]["gcc_vs_gcc"]["tables_subdir"],
    )
    gchp_vs_gcc_tablesdir = os.path.join(
        gchp_vs_gcc_resultsdir,
        config["options"]["comparisons"]["gchp_vs_gcc"]["tables_subdir"],
    )
    gchp_vs_gchp_tablesdir = os.path.join(
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

    # Get days per month and seconds per month for ref
    sec_per_month_ref = np.zeros(12)
    days_per_month_ref = np.zeros(12)
    for mon in range(12):
        days_per_month_ref[mon] = monthrange(int(bmk_year_ref), mon + 1)[1]
        sec_per_month_ref[mon] = days_per_month_ref[mon] * 86400.0

    # Get all months array of start datetimes for benchmark year
    bmk_start_ref = np.datetime64(bmk_year_ref + "-01-01")
    bmk_end_ref = np.datetime64(f"{int(bmk_year_ref) + 1}-01-01")
    all_months_ref = np.arange(
        bmk_start_ref,
        bmk_end_ref,
        step=np.timedelta64(1, "M"),
        dtype="datetime64[M]"
    )
    all_months_gchp_ref = all_months_ref

    # Get subset of month datetimes and seconds per month
    # for only benchmark months
    bmk_mons_ref = all_months_ref[bmk_mon_inds]
    bmk_mons_gchp_ref = all_months_gchp_ref[bmk_mon_inds]
    bmk_sec_per_month_ref = sec_per_month_ref[bmk_mon_inds]

    # Compute seconds in the Ref year
    sec_per_yr_ref = 0
    for t in range(12):
        days_in_mon = monthrange(int(bmk_year_ref), t + 1)[1]
        sec_per_yr_ref += days_in_mon * 86400.0

    # ======================================================================
    # Dates and times -- Dev data
    # ======================================================================

    # Month/year strings for use in table subdirectories (e.g. Jan2016)
    bmk_mon_yr_strs_dev = [v + bmk_year_dev for v in bmk_mon_strs]

    # Get days per month and seconds per month for dev
    sec_per_month_dev = np.zeros(12)
    days_per_month_dev = np.zeros(12)
    for mon in range(12):
        days_per_month_dev[mon] = monthrange(int(bmk_year_dev), mon + 1)[1]
        sec_per_month_dev[mon] = days_per_month_dev[mon] * 86400.0

    # Get all months array of start datetimes for benchmark year
    bmk_start_dev = np.datetime64(bmk_year_dev + "-01-01")
    bmk_end_dev = np.datetime64(f"{int(bmk_year_dev) + 1}-01-01")
    all_months_dev = np.arange(
        bmk_start_dev,
        bmk_end_dev,
        step=np.timedelta64(1, "M"),
        dtype="datetime64[M]"
    )
    all_months_gchp_dev = all_months_dev

    bmk_mons_dev = all_months_dev[bmk_mon_inds]
    bmk_mons_gchp_dev = all_months_gchp_dev[bmk_mon_inds]
    bmk_sec_per_month_dev = sec_per_month_dev[bmk_mon_inds]

    # Compute seconds in the Dev year
    sec_per_yr_dev = 0
    for t in range(12):
        days_in_mon = monthrange(int(bmk_year_dev), t + 1)[1]
        sec_per_yr_dev += days_in_mon * 86400.0

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
    if config["options"]["outputs"]["mass_table"]:
        print(" - Table of species mass")
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
        refmet = get_filepaths(
            gcc_vs_gcc_refdir,
            "StateMet",
            all_months_ref
        )[0]
        devmet = get_filepaths(
            gcc_vs_gcc_devdir,
            "StateMet",
            all_months_dev
        )[0]

        # ==================================================================
        # GCC vs GCC species concentration plots
        # ==================================================================
        if config["options"]["outputs"]["plot_conc"]:
            print("\n%%% Creating GCC vs. GCC concentration plots %%%")

            # Only plot concentration categories for TransportTracers
            restrict_cats = ["RnPbBeTracers", "TransportTracers"]

            # --------------------------------------------------------------
            # GCC vs GCC species concentration plots: Annual mean
            # --------------------------------------------------------------

            # Filepaths
            collection = "SpeciesConc"
            ref = get_filepaths(
                gcc_vs_gcc_refdir,
                collection,
                all_months_ref
            )[0]
            dev = get_filepaths(
                gcc_vs_gcc_devdir,
                collection,
                all_months_dev
            )[0]

            # Create plots
            bmk.make_benchmark_conc_plots(
                ref,
                gcc_vs_gcc_refstr,
                dev,
                gcc_vs_gcc_devstr,
                collection=collection,
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
                n_job=config["options"]["n_cores"]
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
                    collection=collection,
                    refmet=refmet[mon_ind],
                    devmet=devmet[mon_ind],
                    dst=gcc_vs_gcc_resultsdir,
                    subdst=bmk_mon_yr_strs_dev[t],
                    weightsdir=config["paths"]["weights_dir"],
                    benchmark_type=bmk_type,
                    restrict_cats=restrict_cats,
                    overwrite=True,
                    spcdb_dir=spcdb_dir,
                    n_job=config["options"]["n_cores"]
                )

        # ==================================================================
        # GCC vs GCC wet deposition plots
        # ==================================================================
        if config["options"]["outputs"]["plot_wetdep"]:
            print("\n%%% Creating GCC vs. GCC wet deposition plots %%%")

            # Diagnostic collection files to read
            wetdep_collections = ["WetLossConv", "WetLossLS"]

            # Loop over collections
            for collection in wetdep_collections:

                # ----------------------------------------------------------
                # GCC vs. GCC wet deposition plots: Annual mean
                # ----------------------------------------------------------

                # Filepaths
                ref = get_filepaths(
                    gcc_vs_gcc_refdir,
                    collection,
                    all_months_ref
                )[0]
                dev = get_filepaths(
                    gcc_vs_gcc_devdir,
                    collection,
                    all_months_dev
                )[0]

                # Create plots
                bmk.make_benchmark_wetdep_plots(
                    ref,
                    gcc_vs_gcc_refstr,
                    dev,
                    gcc_vs_gcc_devstr,
                    collection=collection,
                    refmet=refmet,
                    devmet=devmet,
                    dst=gcc_vs_gcc_resultsdir,
                    datestr="AnnualMean",
                    time_mean=True,
                    weightsdir=config["paths"]["weights_dir"],
                    benchmark_type=bmk_type,
                    overwrite=True,
                    spcdb_dir=spcdb_dir,
                    n_job=config["options"]["n_cores"]
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
                        collection=collection,
                        refmet=refmet[mon_ind],
                        devmet=devmet[mon_ind],
                        dst=gcc_vs_gcc_resultsdir,
                        datestr=bmk_mon_yr_strs_dev[t],
                        weightsdir=config["paths"]["weights_dir"],
                        benchmark_type=bmk_type,
                        overwrite=True,
                        spcdb_dir=spcdb_dir,
                        n_job=config["options"]["n_cores"]
                    )


        # ==================================================================
        # GCC vs GCC radionuclides budget tables
        # ==================================================================
        if config["options"]["outputs"]["rnpbbe_budget"]:
            print("\n%%% Creating GCC vs. GCC radionuclides budget table %%%")
            ttbdg.transport_tracers_budgets(
                config["data"]["dev"]["gcc"]["dir"],
                gcc_vs_gcc_devdir,
                gcc_vs_gcc_devrstdir,
                int(bmk_year_dev),
                dst=gcc_vs_gcc_tablesdir,
                overwrite=True,
                spcdb_dir=spcdb_dir,
            )

        # ==================================================================
        # GCC vs GCC mass tables
        # ==================================================================
        if config["options"]["outputs"]["mass_table"]:
            print("\n%%% Creating GCC vs. GCC mass tables %%%")

            def gcc_vs_gcc_mass_table(mon):
                """
                Create mass table for each benchmark month mon in parallel
                """

                # Filepaths
                refpath = get_filepath(
                    gcc_vs_gcc_refrstdir,
                    "Restart",
                    bmk_mons_ref[mon]
                )
                devpath = get_filepath(
                    gcc_vs_gcc_devrstdir,
                    "Restart",
                    bmk_mons_dev[mon]
                )

                # Create tables
                bmk.make_benchmark_mass_tables(
                    refpath,
                    gcc_vs_gcc_refstr,
                    devpath,
                    gcc_vs_gcc_devstr,
                    dst=gcc_vs_gcc_tablesdir,
                    subdst=bmk_mon_yr_strs_dev[mon],
                    label=f"at 01{bmk_mon_yr_strs_dev[mon]}",
                    overwrite=True,
                    spcdb_dir=spcdb_dir,
                )

            # Create tables in parallel
            # Turn off parallelization if n_jobs==1
            if config["options"]["n_cores"] != 1:
                results = Parallel(n_jobs=config["options"]["n_cores"])(
                    delayed(gcc_vs_gcc_mass_table)(mon) \
                    for mon in range(bmk_n_months)
                )
            else:
                results = []
                for mon in range(bmk_n_months):
                    results.append(gcc_vs_gcc_mass_table(mon))

        # ==================================================================
        # GCC vs GCC operations budgets tables
        # ==================================================================
        if config["options"]["outputs"]["operations_budget"]:
            print("\n%%% Creating GCC vs. GCC operations budget tables %%%")

            # Filepaths
            col = "Budget"
            refs = get_filepaths(
                gcc_vs_gcc_refdir,
                col,
                all_months_ref
            )[0]
            devs = get_filepaths(
                gcc_vs_gcc_devdir,
                col,
                all_months_dev
            )[0]

            # Create table
            bmk.make_benchmark_operations_budget(
                config["data"]["ref"]["gcc"]["dir"],
                refs,
                config["data"]["dev"]["gcc"]["dir"],
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
            col = "AdvFluxVert"
            devs = get_filepaths(
                gcc_vs_gcc_devdir,
                col,
                all_months_dev
            )[0]

            # Make stat-trop exchange table for subset of species
            ste.make_benchmark_ste_table(
                config["data"]["dev"]["gcc"]["dir"],
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
        refmet = get_filepaths(
            gchp_vs_gcc_refdir,
            "StateMet",
            all_months_dev
        )[0]
        devmet = get_filepaths(
            gchp_vs_gcc_devdir,
            "StateMet",
            all_months_gchp_dev,
            is_gchp=True
        )[0]

        # ==================================================================
        # GCHP vs GCC species concentration plots
        # ==================================================================
        if config["options"]["outputs"]["plot_conc"]:
            print("\n%%% Creating GCHP vs. GCC concentration plots %%%")

            # Only plot concentration categories for TransportTracers
            restrict_cats = ["RnPbBeTracers", "TransportTracers"]

            # --------------------------------------------------------------
            # GCHP vs GCC species concentration plots: Annual Mean
            # --------------------------------------------------------------

            # Filepaths
            collection = "SpeciesConc"
            ref = get_filepaths(
                gchp_vs_gcc_refdir,
                collection,
                all_months_dev
            )[0]
            dev = get_filepaths(
                gchp_vs_gcc_devdir,
                collection,
                all_months_gchp_dev,
                is_gchp=True
            )[0]

            # Create plots
            bmk.make_benchmark_conc_plots(
                ref,
                gchp_vs_gcc_refstr,
                dev,
                gchp_vs_gcc_devstr,
                collection=collection,
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
                n_job=config["options"]["n_cores"]
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
                    collection=collection,
                    refmet=refmet[mon_ind],
                    devmet=devmet[mon_ind],
                    dst=gchp_vs_gcc_resultsdir,
                    subdst=bmk_mon_yr_strs_dev[t],
                    weightsdir=config["paths"]["weights_dir"],
                    benchmark_type=bmk_type,
                    restrict_cats=restrict_cats,
                    overwrite=True,
                    spcdb_dir=spcdb_dir,
                    n_job=config["options"]["n_cores"]
                )

        # ==================================================================
        # GCHP vs GCC wet deposition plots
        # ==================================================================
        if config["options"]["outputs"]["plot_wetdep"]:
            print("\n%%% Creating GCHP vs. GCC wet deposition plots %%%")

            # Create separate set of plots for each wetdep collection
            wetdep_collections = ["WetLossConv", "WetLossLS"]

            # Create plots for each collection and benchmark month
            for collection in wetdep_collections:

                # ----------------------------------------------------------
                # GCHP vs GCC wet deposition plots: Annual mean
                # ----------------------------------------------------------

                # Filepaths
                ref = get_filepaths(
                    gchp_vs_gcc_refdir,
                    collection,
                    all_months_dev
                )[0]
                dev = get_filepaths(
                    gchp_vs_gcc_devdir,
                    collection,
                    all_months_gchp_dev,
                    is_gchp=True
                )[0]

                # Create plots
                bmk.make_benchmark_wetdep_plots(
                    ref,
                    gchp_vs_gcc_refstr,
                    dev,
                    gchp_vs_gcc_devstr,
                    collection=collection,
                    devmet=devmet,
                    dst=gchp_vs_gcc_resultsdir,
                    datestr="AnnualMean",
                    time_mean=True,
                    weightsdir=config["paths"]["weights_dir"],
                    overwrite=True,
                    benchmark_type=bmk_type,
                    normalize_by_area=True,
                    spcdb_dir=spcdb_dir,
                    n_job=config["options"]["n_cores"]
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
                        collection=collection,
                        refmet=refmet[mon_ind],
                        devmet=devmet[mon_ind],
                        dst=gchp_vs_gcc_resultsdir,
                        datestr=bmk_mon_yr_strs_dev[t],
                        weightsdir=config["paths"]["weights_dir"],
                        overwrite=True,
                        benchmark_type=bmk_type,
                        normalize_by_area=True,
                        spcdb_dir=spcdb_dir,
                        n_job=config["options"]["n_cores"]
                    )

        # ==================================================================
        # GCHP vs GCC radionuclides budget tables
        # ==================================================================
        if config["options"]["outputs"]["rnpbbe_budget"]:
            print("\n%%% Creating GCHP vs. GCC radionuclides budget table %%%")
            ttbdg.transport_tracers_budgets(
                config["data"]["dev"]["gchp"]["dir"],
                gchp_vs_gcc_devdir,
                gchp_vs_gcc_devrstdir,
                int(bmk_year_dev),
                dst=gchp_vs_gcc_tablesdir,
                is_gchp=True,
                gchp_res=config["data"]["dev"]["gchp"]["resolution"],
                gchp_is_pre_14_0=config["data"]["dev"]["gchp"]["is_pre_14.0"],
                overwrite=True,
                spcdb_dir=spcdb_dir,
            )

        # ==================================================================
        # GCHP vs GCC global mass tables
        # ==================================================================
        if config["options"]["outputs"]["mass_table"]:
            print("\n%%% Creating GCHP vs. GCC mass tables %%%")

            def gchp_vs_gcc_mass_table(mon):
                """
                Create mass table for each benchmark month in parallel
                """

                # Filepaths
                refpath = get_filepath(
                    gchp_vs_gcc_refrstdir,
                    "Restart",
                    bmk_mons_dev[mon]
                )
                devpath = get_filepath(
                    gchp_vs_gcc_devrstdir,
                    "Restart",
                    bmk_mons_dev[mon],
                    is_gchp=True,
                    gchp_res=config["data"]["dev"]["gchp"]["resolution"],
                    gchp_is_pre_14_0=config["data"]["dev"]["gchp"][
                        "is_pre_14.0"]
                )

                # KLUDGE: ewl, bmy, 13 Oct 2022
                # Use last GCHP restart file, which has correct area values
                devareapath = get_filepath(
                    gchp_vs_gcc_devrstdir,
                    "Restart",
                    bmk_end_dev,
                    is_gchp=True,
                    gchp_res=config["data"]["dev"]["gchp"]["resolution"],
                    gchp_is_pre_14_0=config["data"]["dev"]["gchp"][
                        "is_pre_14.0"]
                )

            # Create tables in parallel
            # Turn off parallelization if n_jobs==1
            if config["options"]["n_cores"] != 1:
                results = Parallel(n_jobs=config["options"]["n_cores"])(
                    delayed(gchp_vs_gcc_mass_table)(mon) \
                    for mon in range(bmk_n_months)
                )
            else:
                results = []
                for mon in range(bmk_n_months):
                    results.append(gchp_vs_gcc_mass_table(mon))

        # ==================================================================
        # GCHP vs GCC operations budgets tables
        # ==================================================================
        if config["options"]["outputs"]["operations_budget"]:
            print("\n%%% Creating GCHP vs. GCC operations budget tables %%%")

            # Filepaths
            collection = "Budget"
            refs = get_filepaths(
                gchp_vs_gcc_refdir,
                collection,
                all_months_dev
            )[0]
            devs = get_filepaths(
                gchp_vs_gcc_devdir,
                collection,
                all_months_gchp_dev,
                is_gchp=True
            )[0]

            # Make operations budget table
            bmk.make_benchmark_operations_budget(
                config["data"]["dev"]["gcc"]["dir"],
                refs,
                config["data"]["dev"]["gchp"]["dir"],
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

        # Option to specify grid resolution of comparison plots
        # This is a temporary hack until cs->cs regridding in GCPy is fixed
        cmpres="1x1.25"

        # ==================================================================
        # GCHP vs GCHP filepaths for StateMet collection data
        # ==================================================================
        refmet = get_filepaths(
            gchp_vs_gchp_refdir,
            "StateMet",
            all_months_gchp_ref,
            is_gchp=True
        )[0]
        devmet = get_filepaths(
            gchp_vs_gchp_devdir,
            "StateMet",
            all_months_gchp_dev,
            is_gchp=True
        )[0]

        # ==================================================================
        # GCHP vs GCHP species concentration plots
        # ==================================================================
        if config["options"]["outputs"]["plot_conc"]:
            print("\n%%% Creating GCHP vs. GCHP concentration plots %%%")

            # Only plot concentration categories for TransportTracers
            restrict_cats = ["RnPbBeTracers", "TransportTracers"]

            # --------------------------------------------------------------
            # GCHP vs GCHP species concentration plots: Annual Mean
            # --------------------------------------------------------------

            # Filepaths
            collection = "SpeciesConc"
            ref = get_filepaths(
                gchp_vs_gchp_refdir,
                collection,
                all_months_gchp_ref,
                is_gchp=True
            )[0]
            dev = get_filepaths(
                gchp_vs_gchp_devdir,
                collection,
                all_months_gchp_dev,
                is_gchp=True
            )[0]

            # Make concentration plots
            bmk.make_benchmark_conc_plots(
                ref,
                gchp_vs_gchp_refstr,
                dev,
                gchp_vs_gchp_devstr,
                collection=collection,
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
                cmpres=cmpres,
                n_job=config["options"]["n_cores"]
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
                    collection=collection,
                    refmet=refmet[mon_ind],
                    devmet=devmet[mon_ind],
                    dst=gchp_vs_gchp_resultsdir,
                    subdst=bmk_mon_yr_strs_dev[t],
                    weightsdir=config["paths"]["weights_dir"],
                    benchmark_type=bmk_type,
                    restrict_cats=restrict_cats,
                    overwrite=True,
                    spcdb_dir=spcdb_dir,
                    cmpres=cmpres,
                    n_job=config["options"]["n_cores"]
                )

        # ==================================================================
        # GCHP vs GCHP wet deposition plots
        # ==================================================================
        if config["options"]["outputs"]["plot_wetdep"]:
            print("\n%%% Creating GCHP vs. GCHP wet deposition plots %%%")

            # Create separate set of plots for each wetdep collection
            wetdep_collections = ["WetLossConv", "WetLossLS"]

            # Create plots for each collection and benchmark month
            for collection in wetdep_collections:

                # ----------------------------------------------------------
                # GCHP vs GCHP wet deposition plots: Annual Mean
                # ----------------------------------------------------------

                # Filepaths
                ref = get_filepaths(
                    gchp_vs_gchp_refdir,
                    collection,
                    all_months_gchp_ref,
                    is_gchp=True
                )[0]
                dev = get_filepaths(
                    gchp_vs_gchp_devdir,
                    collection,
                    all_months_gchp_dev,
                    is_gchp=True
                )[0]

                # Create plots
                bmk.make_benchmark_wetdep_plots(
                    ref,
                    gchp_vs_gchp_refstr,
                    dev,
                    gchp_vs_gchp_devstr,
                    collection=collection,
                    refmet=refmet,
                    devmet=devmet,
                    dst=gchp_vs_gchp_resultsdir,
                    datestr="AnnualMean",
                    time_mean=True,
                    weightsdir=config["paths"]["weights_dir"],
                    overwrite=True,
                    benchmark_type=bmk_type,
                    normalize_by_area=True,
                    spcdb_dir=spcdb_dir,
                    cmpres=cmpres,
                    n_job=config["options"]["n_cores"]
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
                        collection=collection,
                        refmet=refmet[mon_ind],
                        devmet=devmet[mon_ind],
                        dst=gchp_vs_gchp_resultsdir,
                        datestr=bmk_mon_yr_strs_dev[t],
                        weightsdir=config["paths"]["weights_dir"],
                        overwrite=True,
                        benchmark_type=bmk_type,
                        normalize_by_area=True,
                        spcdb_dir=spcdb_dir,
                        cmpres=cmpres,
                        n_job=config["options"]["n_cores"]
                    )

        # ==================================================================
        # GCHP vs GCHP radionuclides budget table
        # ==================================================================
        if config["options"]["outputs"]["rnpbbe_budget"]:
            print("\n%%% Creating GCHP vs. GCHP radionuclides budget table %%%")
            ttbdg.transport_tracers_budgets(
                config["data"]["dev"]["gchp"]["dir"],
                gchp_vs_gchp_devdir,
                gchp_vs_gchp_devrstdir,
                int(bmk_year_dev),
                dst=gchp_vs_gchp_tablesdir,
                is_gchp=True,
                gchp_res=config["data"]["dev"]["gchp"]["resolution"],
                gchp_is_pre_14_0=config["data"]["dev"]["gchp"]["is_pre_14.0"],
                overwrite=True,
                spcdb_dir=spcdb_dir
            )

        # ==================================================================
        # GCHP vs GCHP global mass tables
        # ==================================================================
        if config["options"]["outputs"]["mass_table"]:
            print("\n%%% Creating GCHP vs. GCHP mass tables %%%")

            def gchp_vs_gchp_mass_table(mon):
                """
                Create mass table for each benchmark month m in parallel
                """

                # Ref filepaths
                refpath = get_filepath(
                    gchp_vs_gchp_refrstdir,
                    "Restart",
                    bmk_mons_ref[mon],
                    is_gchp=True,
                    gchp_res=config["data"]["ref"]["gchp"]["resolution"],
                    gchp_is_pre_14_0=config["data"]["ref"]["gchp"][
                        "is_pre_14.0"]
                )

                # Dev filepaths
                devpath = get_filepath(
                    gchp_vs_gchp_devrstdir,
                    "Restarts",
                    bmk_mons_dev[mon],
                    is_gchp=True,
                    gchp_res=config["data"]["dev"]["gchp"]["resolution"],
                    gchp_is_pre_14_0=config["data"]["dev"]["gchp"][
                        "is_pre_14.0"]
                )

                # KLUDGE: ewl, bmy, 13 Oct 2022
                # Use last GCHP restart file, which has correct area values
                refareapath = get_filepath(
                    gchp_vs_gchp_refrstdir,
                    "Restart",
                    bmk_end_ref,
                    is_gchp=True,
                    gchp_res=config["data"]["dev"]["gchp"]["resolution"],
                    gchp_is_pre_14_0=config["data"]["dev"]["gchp"][
                        "is_pre_14.0"]
                )
                devareapath = get_filepath(
                    gchp_vs_gchp_devrstdir,
                    "Restart",
                    bmk_end_dev,
                    is_gchp=True,
                    gchp_res=config["data"]["dev"]["gchp"]["resolution"],
                    gchp_is_pre_14_0=config["data"]["dev"]["gchp"][
                        "is_pre_14.0"]
                )

                # Create tables
                bmk.make_benchmark_mass_tables(
                    refpath,
                    gchp_vs_gchp_refstr,
                    devpath,
                    gchp_vs_gchp_devstr,
                    dst=gchp_vs_gchp_tablesdir,
                    subdst=bmk_mon_yr_strs_dev[mon],
                    label=f"at 01{bmk_mon_yr_strs_dev[mon]}",
                    overwrite=True,
                    spcdb_dir=spcdb_dir,
                    ref_met_extra=refareapath,
                    dev_met_extra=devareapath
                )

            # Create tables in parallel
            # Turn off parallelization if n_jobs==1
            if config["options"]["n_cores"] != 1:
                results = Parallel(n_jobs=config["options"]["n_cores"])(
                    delayed(gchp_vs_gchp_mass_table)(mon) \
                    for mon in range(bmk_n_months)
                )
            else:
                results = []
                for mon in range(bmk_n_months):
                    results.append(gchp_vs_gchp_mass_table(mon))

        # ==================================================================
        # GCHP vs GCHP operations budgets tables
        # ==================================================================
        if config["options"]["outputs"]["operations_budget"]:
            print("\n%%% Creating GCHP vs. GCHP operations budget tables %%%")

            # Filepaths
            col = "Budget"
            refs = get_filepaths(
                gchp_vs_gchp_refdir,
                col,
                all_months_gchp_ref,
                is_gchp=True
            )[0]
            devs = get_filepaths(
                gchp_vs_gchp_devdir,
                col,
                all_months_gchp_dev,
                is_gchp=True
            )[0]

            # Create table
            bmk.make_benchmark_operations_budget(
                config["data"]["dev"]["gchp"]["dir"],
                refs,
                config["data"]["dev"]["gchp"]["dir"],
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

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Create mass conservations tables for GCC and GCHP
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if config["options"]["outputs"]["cons_table"]:

        # ===================================================================
        # Create mass conservation table for GCC_dev
        # ===================================================================
        if (
            config["options"]["comparisons"]["gcc_vs_gcc"]["run"]
            or config["options"]["comparisons"]["gchp_vs_gcc"]["run"]
        ):
            print("\n%%% Creating GCC dev mass conservation table %%%")

            # Filepaths
            datafiles = get_filepaths(
                gcc_vs_gcc_devrstdir,
                "Restart",
                all_months_dev
            )[0]

            # Pick output folder
            if config["options"]["comparisons"]["gchp_vs_gcc"]["run"]:
                tablesdir = gchp_vs_gcc_tablesdir
            else:
                tablesdir = gcc_vs_gcc_tablesdir

            # Create table
            bmk.make_benchmark_mass_conservation_table(
                datafiles,
                config["data"]["dev"]["gcc"]["dir"],
                dst=tablesdir,
                overwrite=True,
                spcdb_dir=spcdb_dir,
            )

        # =====================================================================
        # Create mass conservation table for GCHP_dev
        # =====================================================================
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
                gchp_res=config["data"]["dev"]["gchp"]["resolution"],
                gchp_is_pre_14_0=config["data"]["dev"]["gchp"]["is_pre_14.0"],
            )[0]

            # KLUDGE: ewl, bmy, 14 Oct 2022
            # Read the AREA from the restart file at the end of the
            # simulation.  Due to a GCHP bug, intermediate restarts
            # have AREA with all zeroes.
            areapath = get_filepath(
                gchp_vs_gcc_devrstdir,
                "Restart",
                bmk_end_dev,
                is_gchp=True,
                gchp_res=config["data"]["dev"]["gchp"]["resolution"],
                gchp_is_pre_14_0=config["data"]["dev"]["gchp"]["is_pre_14.0"],
            )

            # Pick output folder
            if config["options"]["comparisons"]["gchp_vs_gcc"]["run"]:
                tablesdir = gchp_vs_gcc_tablesdir
            else:
                tablesdir = gchp_vs_gchp_tablesdir

            # KLUDGE: ewl, bmy, 14 Oct 2022
            # Read the AREA from the restart file at the end of the
            # simulation.  Due to a GCHP bug, intermediate restarts
            # have AREA with all zeroes.
            bmk.make_benchmark_mass_conservation_table(
                datafiles,
                config["data"]["dev"]["gchp"]["dir"],
                dst=tablesdir,
                overwrite=True,
                spcdb_dir=spcdb_dir,
                areapath=areapath
            )

    # ==================================================================
    # Print a message indicating that the benchmarks finished
    # ==================================================================
    print("\n %%%% All requested benchmark plots/tables created! %%%%")
