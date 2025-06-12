#!/usr/bin/env python
"""
run_1yr_fullchem_benchmark.py:
    Driver script for creating benchmark plots and testing gcpy
    1-year full-chemistry benchmark capability.

Run this script to generate benchmark comparisons between:

    (1) GCC (aka GEOS-Chem "Classic") vs. GCC
    (2) GCHP vs GCC
    (3) GCHP vs GCHP

You can customize this by editing the settings in the corresponding yaml
config file (eg. 1yr_fullchem_benchmark.yml).

To generate benchmark output:

    (1) Copy the file gcpy/benchmark/config/1yr_fullchem_benchmark.yml
        to a folder of your choice.

    (2) Edit the 1yr_fullchem_benchmark.yml to select the desired options
        and to point to the proper file paths on your system.

    (3) Run the command:

        $ python -m gcpy.benchmark.run_benchmark.py 1yr_fullchem_benchmark.yml

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

# =====================================================================
# Imports and global settings (you should not need to edit these)
# =====================================================================

import os
import warnings
from calendar import monthrange
import numpy as np
from joblib import Parallel, delayed
from gcpy.util import copy_file_to_dir, get_filepath, get_filepaths
from gcpy.benchmark.modules.ste_flux import make_benchmark_ste_table
from gcpy.benchmark.modules.oh_metrics import make_benchmark_oh_metrics
from gcpy.benchmark.modules.budget_ox import global_ox_budget
#TODO: Peel out routines from benchmark_funcs.py into smaller
# routines in the gcpy/benchmark/modules folder, such as these:
from gcpy.benchmark.modules.benchmark_funcs import \
    diff_of_diffs_toprow_title, get_species_database_dir, \
    make_benchmark_conc_plots, make_benchmark_emis_plots, \
    make_benchmark_emis_tables, make_benchmark_jvalue_plots, \
    make_benchmark_aod_plots, make_benchmark_mass_tables, \
    make_benchmark_operations_budget, make_benchmark_aerosol_tables
from gcpy.benchmark.modules.benchmark_utils import \
    gcc_vs_gcc_dirs, gchp_vs_gcc_dirs, gchp_vs_gchp_dirs, \
    get_log_filepaths, print_benchmark_info
from gcpy.benchmark.modules.benchmark_models_vs_obs \
    import make_benchmark_models_vs_obs_plots
from gcpy.benchmark.modules.benchmark_models_vs_sondes \
    import make_benchmark_models_vs_sondes_plots
from gcpy.benchmark.modules.benchmark_drydep \
    import drydepvel_species, make_benchmark_drydep_plots
from gcpy.benchmark.modules.benchmark_scrape_gcclassic_timers import \
    make_benchmark_gcclassic_timing_table
from gcpy.benchmark.modules.benchmark_scrape_gchp_timers import \
    make_benchmark_gchp_timing_table

# Tell matplotlib not to look for an X-window
os.environ["QT_QPA_PLATFORM"] = "offscreen"

# Suppress annoying warning messages
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

def run_benchmark(config, bmk_year_ref, bmk_year_dev):
    """
    Runs 1 year benchmark with the given configuration settings.

    Args:
        config : dict
            Contains configuration for 1yr benchmark from yaml file.
    """

    # This script has a fixed benchmark type
    bmk_type = config["options"]["bmk_type"]
    bmk_mon_strs = ["Jan", "Apr", "Jul", "Oct"]
    bmk_mon_inds = [0, 3, 6, 9]
    bmk_n_months = len(bmk_mon_strs)

    # Folder in which the species_database.yml file is found
    spcdb_dir = get_species_database_dir(config)

    # ======================================================================
    # Data directories
    # ======================================================================

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

    # Directories where plots & tables will be created
    mainresultsdir = os.path.join(config["paths"]["results_dir"])
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
    diff_of_diffs_resultsdir = os.path.join(
        mainresultsdir,
        "GCHP_GCC_diff_of_diffs"
    )

    # Create the main results directory first
    if not os.path.exists(mainresultsdir):
        os.mkdir(mainresultsdir)

    # Create results directories that don't exist.  Also place a copy of
    # this file plus the YAML configuration file in each results directory.
    results_list = [
        gcc_vs_gcc_resultsdir,
        gchp_vs_gchp_resultsdir,
        gchp_vs_gcc_resultsdir,
        diff_of_diffs_resultsdir
    ]
    comparisons_list = [
        config["options"]["comparisons"]["gcc_vs_gcc"]["run"],
        config["options"]["comparisons"]["gchp_vs_gchp"]["run"],
        config["options"]["comparisons"]["gchp_vs_gcc"]["run"],
        config["options"]["comparisons"]["gchp_vs_gcc_diff_of_diffs"]["run"]
    ]
    for (resdir, plotting_type) in zip(results_list, comparisons_list):
        if plotting_type and not os.path.exists(resdir):
            os.mkdir(resdir)
            if resdir in results_list:
                copy_file_to_dir(__file__, resdir)
                copy_file_to_dir(config["configuration_file_name"], resdir)

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

    ## Budget directories
    #gcc_vs_gcc_budgetdir = os.path.join(gcc_vs_gcc_resultsdir, "Budget")
    #gchp_vs_gcc_budgetdir = os.path.join(gchp_vs_gcc_resultsdir, "Budget")
    #gchp_vs_gchp_budgetdir = os.path.join(gchp_vs_gchp_resultsdir, "Budget")

    # Models vs. observations directories
    s = "ModelVsObs"
    gcc_vs_gcc_models_vs_obs_dir = os.path.join(gcc_vs_gcc_resultsdir, s)
    gchp_vs_gcc_models_vs_obs_dir = os.path.join(gchp_vs_gcc_resultsdir, s)
    gchp_vs_gchp_models_vs_obs_dir = os.path.join(gchp_vs_gchp_resultsdir, s)

    # ======================================================================
    # Plot title strings
    # ======================================================================
    gcc_vs_gcc_refstr = config["data"]["ref"]["gcc"]["version"]
    gcc_vs_gcc_devstr = config["data"]["dev"]["gcc"]["version"]
    gchp_vs_gcc_refstr = config["data"]["dev"]["gcc"]["version"]
    gchp_vs_gcc_devstr = config["data"]["dev"]["gchp"]["version"]
    gchp_vs_gchp_refstr = config["data"]["ref"]["gchp"]["version"]
    gchp_vs_gchp_devstr = config["data"]["dev"]["gchp"]["version"]
    diff_of_diffs_refstr = diff_of_diffs_toprow_title(config, "gcc")
    diff_of_diffs_devstr = diff_of_diffs_toprow_title(config, "gchp")

    # ======================================================================
    # Observational data files
    # ======================================================================
    sondes_data_file = os.path.join(
        config["paths"]["obs_data"]["sondes"]["data_dir"],
        config["paths"]["obs_data"]["sondes"]["data_file"],
    )
    sondes_site_file = os.path.join(
        config["paths"]["obs_data"]["sondes"]["data_dir"],
        config["paths"]["obs_data"]["sondes"]["site_file"],
    )

    ########################################################################
    ###    THE REST OF THESE SETTINGS SHOULD NOT NEED TO BE CHANGED      ###
    ########################################################################

    # =====================================================================
    # Dates and times -- ref data
    # =====================================================================

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
        bmk_start_ref, bmk_end_ref, step=np.timedelta64(1, "M"), dtype="datetime64[M]"
    )
    all_months_gchp_ref = all_months_ref

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
    for mon in range(12):
        days_per_month_dev[mon] = monthrange(int(bmk_year_dev), mon + 1)[1]
        sec_per_month_dev[mon] = days_per_month_dev[mon] * 86400.0

    # Get all months array of start datetimes for benchmark year
    bmk_start_dev = np.datetime64(bmk_year_dev + "-01-01")
    bmk_end_dev = np.datetime64(f"{int(bmk_year_dev) + 1}-01-01")
    all_months_dev = np.arange(
        bmk_start_dev, bmk_end_dev, step=np.timedelta64(1, "M"), dtype="datetime64[M]"
    )
    all_months_gchp_dev = all_months_dev

    # Get subset of month datetimes and seconds per month for only benchmark months
    bmk_mons_dev = all_months_dev[bmk_mon_inds]
    bmk_mons_gchp_dev = all_months_gchp_dev[bmk_mon_inds]
    bmk_sec_per_month_dev = sec_per_month_dev[bmk_mon_inds]

    # ======================================================================
    # Print the list of plots & tables being generated
    # ======================================================================
    print_benchmark_info(config)

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
            ref = get_filepaths(
                gcc_vs_gcc_refdir,
                "SpeciesConc",
                all_months_ref
            )[0]
            dev = get_filepaths(
                gcc_vs_gcc_devdir,
                "SpeciesConc",
                all_months_dev
            )[0]

            # Create plots
            print("\nCreating plots for annual mean")
            make_benchmark_conc_plots(
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
                n_job=config["options"]["n_cores"]
            )

            # --------------------------------------------------------------
            # GCC vs GCC species concentration plots: Seasonal
            # --------------------------------------------------------------
            for mon in range(bmk_n_months):
                print(f"\nCreating plots for {bmk_mon_strs[mon]}")

                # Create plots
                mon_ind = bmk_mon_inds[mon]
                make_benchmark_conc_plots(
                    ref[mon_ind],
                    gcc_vs_gcc_refstr,
                    dev[mon_ind],
                    gcc_vs_gcc_devstr,
                    refmet=refmet[mon_ind],
                    devmet=devmet[mon_ind],
                    dst=gcc_vs_gcc_resultsdir,
                    subdst=bmk_mon_yr_strs_dev[mon],
                    weightsdir=config["paths"]["weights_dir"],
                    benchmark_type=bmk_type,
                    plot_by_spc_cat=config["options"]["outputs"][
                        "plot_options"]["by_spc_cat"],
                    overwrite=True,
                    spcdb_dir=spcdb_dir,
                    n_job=config["options"]["n_cores"]
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
            ref = get_filepaths(
                gcc_vs_gcc_refdir,
                "Emissions",
                all_months_ref
            )[0]
            dev = get_filepaths(
                gcc_vs_gcc_devdir,
                "Emissions",
                all_months_dev
            )[0]

            # Create plots
            print("\nCreating plots for annual mean")
            make_benchmark_emis_plots(
                ref,
                gcc_vs_gcc_refstr,
                dev,
                gcc_vs_gcc_devstr,
                dst=gcc_vs_gcc_resultsdir,
                subdst="AnnualMean",
                time_mean=True,
                weightsdir=config["paths"]["weights_dir"],
                plot_by_spc_cat=config["options"]["outputs"][
                    "plot_options"]["by_spc_cat"],
                plot_by_hco_cat=config["options"]["outputs"][
                    "plot_options"]["by_hco_cat"],
                benchmark_type=bmk_type,
                overwrite=True,
                spcdb_dir=spcdb_dir,
                n_job=config["options"]["n_cores"]
            )

            # --------------------------------------------------------------
            # GCC vs GCC emissions plots: Seasonal
            # --------------------------------------------------------------
            for mon in range(bmk_n_months):
                print(f"\nCreating plots for {bmk_mon_strs[mon]}")

                # Create plots
                mon_ind = bmk_mon_inds[mon]
                make_benchmark_emis_plots(
                    ref[mon_ind],
                    gcc_vs_gcc_refstr,
                    dev[mon_ind],
                    gcc_vs_gcc_devstr,
                    dst=gcc_vs_gcc_resultsdir,
                    subdst=bmk_mon_yr_strs_dev[mon],
                    weightsdir=config["paths"]["weights_dir"],
                    plot_by_spc_cat=config["options"]["outputs"][
                        "plot_options"]["by_spc_cat"],
                    plot_by_hco_cat=config["options"]["outputs"][
                        "plot_options"]["by_hco_cat"],
                    benchmark_type=bmk_type,
                    overwrite=True,
                    spcdb_dir=spcdb_dir,
                    n_job=config["options"]["n_cores"]
                )

        # ==================================================================
        # GCC vs GCC tables of emission and inventory totals
        # ==================================================================
        if config["options"]["outputs"]["emis_table"]:
            print("\n%%% Creating GCC vs. GCC emissions & inventory totals %%%")

            # Filepaths
            ref = get_filepaths(
                gcc_vs_gcc_refdir,
                "Emissions",
                all_months_ref
            )[0]
            dev = get_filepaths(
                gcc_vs_gcc_devdir,
                "Emissions",
                all_months_dev
            )[0]

            # Create table
            make_benchmark_emis_tables(
                ref,
                gcc_vs_gcc_refstr,
                dev,
                gcc_vs_gcc_devstr,
                dst=gcc_vs_gcc_resultsdir,
                benchmark_type=bmk_type,
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
            ref = get_filepaths(
                gcc_vs_gcc_refdir,
                "JValues",
                all_months_ref
            )[0]
            dev = get_filepaths(
                gcc_vs_gcc_devdir,
                "JValues",
                all_months_dev
            )[0]

            # Create plots
            print("\nCreating plots for annual mean")
            make_benchmark_jvalue_plots(
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
                n_job=config["options"]["n_cores"]
            )

            # --------------------------------------------------------------
            # GCC vs GCC J-value plots: Seasonal
            # --------------------------------------------------------------
            for mon in range(bmk_n_months):
                print(f"\nCreating plots for {bmk_mon_strs[mon]}")

                # Create plots
                mon_ind = bmk_mon_inds[mon]
                make_benchmark_jvalue_plots(
                    ref[mon_ind],
                    gcc_vs_gcc_refstr,
                    dev[mon_ind],
                    gcc_vs_gcc_devstr,
                    dst=gcc_vs_gcc_resultsdir,
                    subdst=bmk_mon_yr_strs_dev[mon],
                    weightsdir=config["paths"]["weights_dir"],
                    overwrite=True,
                    spcdb_dir=spcdb_dir,
                    n_job=config["options"]["n_cores"]
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
            ref = get_filepaths(
                gcc_vs_gcc_refdir,
                "Aerosols",
                all_months_ref
            )[0]
            dev = get_filepaths(
                gcc_vs_gcc_devdir,
                "Aerosols",
                all_months_dev
            )[0]

            # Create plots
            print("\nCreating plots for annual mean")
            make_benchmark_aod_plots(
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
                n_job=config["options"]["n_cores"]
            )

            # --------------------------------------------------------------
            # GCC vs GCC column AOD plots: Seasonal
            # --------------------------------------------------------------
            for mon in range(bmk_n_months):
                print(f"\nCreating plots for {bmk_mon_strs[mon]}")

                # Create plots
                mon_ind = bmk_mon_inds[mon]
                make_benchmark_aod_plots(
                    ref[mon_ind],
                    gcc_vs_gcc_refstr,
                    dev[mon_ind],
                    gcc_vs_gcc_devstr,
                    dst=gcc_vs_gcc_resultsdir,
                    subdst=bmk_mon_yr_strs_dev[mon],
                    weightsdir=config["paths"]["weights_dir"],
                    overwrite=True,
                    spcdb_dir=spcdb_dir,
                    n_job=config["options"]["n_cores"]
                )

        # ==================================================================
        # GCC vs. GCC drydep plots
        # ==================================================================
        if config["options"]["outputs"]["plot_drydep"]:
            print("\n%%% Creating GCC vs. GCC drydep plots %%%")

            # --------------------------------------------------------------
            # GCC vs GCC drydep plots: Annual mean
            # --------------------------------------------------------------

            # Filepaths
            ref = get_filepaths(
                gcc_vs_gcc_refdir,
                "DryDep",
                all_months_ref
            )[0]
            dev = get_filepaths(
                gcc_vs_gcc_devdir,
                "DryDep",
                all_months_dev
            )[0]

            # Create plots
            print("\nCreating plots for annual mean")
            make_benchmark_drydep_plots(
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
                n_job=config["options"]["n_cores"],
                varlist=drydepvel_species()
            )

            # --------------------------------------------------------------
            # GCC vs GCC drydep plots: Seasonal
            # --------------------------------------------------------------
            for mon in range(bmk_n_months):
                print(f"\nCreating plots for {bmk_mon_strs[mon]}")

                # Create plots
                mon_ind = bmk_mon_inds[mon]
                make_benchmark_drydep_plots(
                    ref[mon_ind],
                    gcc_vs_gcc_refstr,
                    dev[mon_ind],
                    gcc_vs_gcc_devstr,
                    dst=gcc_vs_gcc_resultsdir,
                    subdst=bmk_mon_yr_strs_dev[mon],
                    weightsdir=config["paths"]["weights_dir"],
                    overwrite=True,
                    spcdb_dir=spcdb_dir,
                    n_job=config["options"]["n_cores"],
                    varlist=drydepvel_species()
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
                make_benchmark_mass_tables(
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
                    delayed(gcc_vs_gcc_mass_table)(mon)
                    for mon in range(bmk_n_months)
                )
            else:
                results = []
                for mon in range(bmk_n_months):
                    results.append(gcc_vs_gcc_mass_table(mon))

        # ==================================================================
        # GCC vs GCC operations budgets tables
        # ==================================================================
        if config["options"]["outputs"]["ops_budget_table"]:
            print("\n%%% Creating GCC vs. GCC operations budget tables %%%")

            def gcc_vs_gcc_ops_budg(mon):
                """
                Create budget table for each benchmark month m in parallel
                """

                # Filepaths
                refpath = get_filepath(
                    gcc_vs_gcc_refdir,
                    "Budget",
                    bmk_mons_ref[mon]
                )
                devpath = get_filepath(
                    gcc_vs_gcc_devdir,
                    "Budget",
                    bmk_mons_dev[mon]
                )

                # Create tables
                make_benchmark_operations_budget(
                    config["data"]["ref"]["gcc"]["version"],
                    refpath,
                    config["data"]["dev"]["gcc"]["version"],
                    devpath,
                    sec_per_month_ref[mon],
                    sec_per_month_dev[mon],
                    benchmark_type=bmk_type,
                    label=f"at 01{bmk_mon_yr_strs_dev[mon]}",
                    dst=gcc_vs_gcc_tablesdir,
                )

            # Create tables in parallel
            # Turn off parallelization if n_jobs==1
            if config["options"]["n_cores"] != 1:
                results = Parallel(n_jobs=config["options"]["n_cores"])(
                    delayed(gcc_vs_gcc_ops_budg)(mon) \
                    for mon in range(bmk_n_months)
                )
            else:
                results = []
                for mon in range(bmk_n_months):
                    results.append(gcc_vs_gcc_ops_budg(mon))

        # ==================================================================
        # GCC vs GCC aerosols budgets/burdens tables
        # ==================================================================
        if config["options"]["outputs"]["aer_budget_table"]:
            print("\n%%% Creating GCC vs. GCC aerosols budget tables %%%")

            # Filepaths
            devaero = get_filepaths(
                gcc_vs_gcc_devdir,
                "Aerosols",
                all_months_dev
            )[0]
            devspc = get_filepaths(
                gcc_vs_gcc_devdir,
                "SpeciesConc",
                all_months_dev
            )[0]

            # Compute tables
            make_benchmark_aerosol_tables(
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
            global_ox_budget(
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
            dev = get_filepaths(
                gcc_vs_gcc_devdir,
                "AdvFluxVert",
                all_months_dev
            )[0]

            # Compute table
            make_benchmark_ste_table(
                config["data"]["dev"]["gcc"]["version"],
                dev,
                bmk_year_dev,
                dst=gcc_vs_gcc_tablesdir,
                bmk_type=bmk_type,
                species=["O3"],
                overwrite=True,
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
                all_months_ref
            )
            dev = get_log_filepaths(
                gcc_vs_gcc_devlogdir,
                config["data"]["dev"]["gcc"]["logs_template"],
                all_months_dev
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
        # GCC vs GCC Global mean OH, MCF Lifetime, CH4 Lifetime
        # ==================================================================
        if config["options"]["outputs"]["OH_metrics"]:
            print("\n%%% Creating GCC vs. GCC OH metrics %%%")

            # Filepaths
            ref = get_filepaths(
                gcc_vs_gcc_refdir,
                "Metrics",
                all_months_ref
            )[0]
            dev = get_filepaths(
                gcc_vs_gcc_devdir,
                "Metrics",
                all_months_dev
            )[0]

            # Create the OH Metrics table
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
        # GCC vs GCC Model vs. Observations plots
        # ==================================================================
        if config["options"]["outputs"]["plot_models_vs_obs"]:
            print("\n%%% Creating GCC vs. GCC models vs. obs. plots %%%")

            # Filepaths
            ref = get_filepaths(
                gcc_vs_gcc_refdir,
                "SpeciesConc",
                all_months_ref
            )[0]
            dev = get_filepaths(
                gcc_vs_gcc_devdir,
                "SpeciesConc",
                all_months_dev
            )[0]

            # Plot models vs. observations
            make_benchmark_models_vs_obs_plots(
                config["paths"]["obs_data"]["ebas_o3"]["data_dir"],
                config["paths"]["obs_data"]["ebas_o3"]["data_label"],
                ref,
                config["data"]["ref"]["gcc"]["version"],
                dev,
                config["data"]["dev"]["gcc"]["version"],
                dst=gcc_vs_gcc_models_vs_obs_dir,
                overwrite=True,
                verbose=False
            )

            # Plot models vs. sondes
            make_benchmark_models_vs_sondes_plots(
                sondes_data_file,
                sondes_site_file,
                ref,
                config["data"]["ref"]["gcc"]["version"],
                dev,
                config["data"]["dev"]["gcc"]["version"],
                dst=gcc_vs_gcc_models_vs_obs_dir,
                overwrite=True,
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
        refmet = get_filepaths(
            gchp_vs_gcc_refdir,
            "StateMet",
            all_months_dev
        )[0]
        devmet = get_filepaths(
            gchp_vs_gcc_devdir,
            "StateMet",
            all_months_gchp_dev,
            is_gchp=True,
        )[0]

        # Get GCHP grid resolution from met collection file
        #ds_devmet = xr.open_dataset(devmet[0])
        #gchp_dev_res = str(get_input_res(ds_devmet)[0])

        # ==================================================================
        # GCHP vs GCC Concentration plots
        # ==================================================================
        if config["options"]["outputs"]["plot_conc"]:
            print("\n%%% Creating GCHP vs. GCC concentration plots %%%")

            # --------------------------------------------------------------
            # GCHP vs GCC species concentration plots: Annual Mean
            # --------------------------------------------------------------

            # Filepaths
            ref = get_filepaths(
                gchp_vs_gcc_refdir,
                "SpeciesConc",
                all_months_dev
            )[0]
            dev = get_filepaths(
                gchp_vs_gcc_devdir,
                "SpeciesConc",
                all_months_gchp_dev,
                is_gchp=True
            )[0]

            # Create plots
            print("\nCreating plots for annual mean")
            make_benchmark_conc_plots(
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
                n_job=config["options"]["n_cores"]
            )

            # --------------------------------------------------------------
            # GCHP vs GCC species concentration plots: Seasonal
            # --------------------------------------------------------------
            for mon in range(bmk_n_months):
                print(f"\nCreating plots for {bmk_mon_strs[mon]}")

                # Create plots
                mon_ind = bmk_mon_inds[mon]
                make_benchmark_conc_plots(
                    ref[mon_ind],
                    gchp_vs_gcc_refstr,
                    dev[mon_ind],
                    gchp_vs_gcc_devstr,
                    refmet=refmet[mon_ind],
                    devmet=devmet[mon_ind],
                    dst=gchp_vs_gcc_resultsdir,
                    subdst=bmk_mon_yr_strs_dev[mon],
                    weightsdir=config["paths"]["weights_dir"],
                    benchmark_type=bmk_type,
                    plot_by_spc_cat=config["options"]["outputs"][
                        "plot_options"]["by_spc_cat"],
                    overwrite=True,
                    spcdb_dir=spcdb_dir,
                    n_job=config["options"]["n_cores"]
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
            ref = get_filepaths(
                gchp_vs_gcc_refdir,
                "Emissions",
                all_months_dev
            )[0]
            dev = get_filepaths(
                gchp_vs_gcc_devdir,
                "Emissions",
                all_months_gchp_dev,
                is_gchp=True
            )[0]

            # Create plots
            print("\nCreating plots for annual mean")
            make_benchmark_emis_plots(
                ref,
                gchp_vs_gcc_refstr,
                dev,
                gchp_vs_gcc_devstr,
                dst=gchp_vs_gcc_resultsdir,
                subdst="AnnualMean",
                time_mean=True,
                weightsdir=config["paths"]["weights_dir"],
                plot_by_spc_cat=config["options"]["outputs"][
                    "plot_options"]["by_spc_cat"],
                plot_by_hco_cat=config["options"]["outputs"][
                    "plot_options"]["by_hco_cat"],
                benchmark_type=bmk_type,
                overwrite=True,
                spcdb_dir=spcdb_dir,
                n_job=config["options"]["n_cores"]
            )

            # --------------------------------------------------------------
            # GCHP vs GCC emissions plots: Seasonal
            # --------------------------------------------------------------
            for mon in range(bmk_n_months):
                print(f"\nCreating plots for {bmk_mon_strs[mon]}")

                # Create plots
                mon_ind = bmk_mon_inds[mon]
                make_benchmark_emis_plots(
                    ref[mon_ind],
                    gchp_vs_gcc_refstr,
                    dev[mon_ind],
                    gchp_vs_gcc_devstr,
                    dst=gchp_vs_gcc_resultsdir,
                    subdst=bmk_mon_yr_strs_dev[mon],
                    weightsdir=config["paths"]["weights_dir"],
                    plot_by_spc_cat=config["options"]["outputs"][
                        "plot_options"]["by_spc_cat"],
                    plot_by_hco_cat=config["options"]["outputs"][
                        "plot_options"]["by_hco_cat"],
                    benchmark_type=bmk_type,
                    overwrite=True,
                    spcdb_dir=spcdb_dir,
                    n_job=config["options"]["n_cores"]
                )

        # ==================================================================
        # GCHP vs. GCC tables of emission and inventory totals
        # ==================================================================
        if config["options"]["outputs"]["emis_table"]:
            print("\n%%% Creating GCHP vs. GCC emissions tables %%%")

            # Filepaths
            ref = get_filepaths(
                gchp_vs_gcc_refdir,
                "Emissions",
                all_months_dev
            )[0]
            dev = get_filepaths(
                gchp_vs_gcc_devdir,
                "Emissions",
                all_months_gchp_dev,
                is_gchp=True
            )[0]

            # Create emissions table that spans entire year
            make_benchmark_emis_tables(
                ref,
                gchp_vs_gcc_refstr,
                dev,
                gchp_vs_gcc_devstr,
                devmet=devmet,
                dst=gchp_vs_gcc_resultsdir,
                ref_interval=sec_per_month_ref,
                dev_interval=sec_per_month_dev,
                benchmark_type=bmk_type,
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
            ref = get_filepaths(
                gchp_vs_gcc_refdir,
                "JValues",
                all_months_dev
            )[0]
            dev = get_filepaths(
                gchp_vs_gcc_devdir,
                "JValues",
                all_months_gchp_dev,
                is_gchp=True
            )[0]

            # Create plots
            print("\nCreating plots for annual mean")
            make_benchmark_jvalue_plots(
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
                n_job=config["options"]["n_cores"]
            )

            # --------------------------------------------------------------
            # GCHP vs GCC J-values plots: Seasonal
            # --------------------------------------------------------------
            for mon in range(bmk_n_months):
                print(f"\nCreating plots for {bmk_mon_strs[mon]}")

                # Create plots
                mon_ind = bmk_mon_inds[mon]
                make_benchmark_jvalue_plots(
                    ref[mon_ind],
                    gchp_vs_gcc_refstr,
                    dev[mon_ind],
                    gchp_vs_gcc_devstr,
                    dst=gchp_vs_gcc_resultsdir,
                    subdst=bmk_mon_yr_strs_dev[mon],
                    weightsdir=config["paths"]["weights_dir"],
                    overwrite=True,
                    spcdb_dir=spcdb_dir,
                    n_job=config["options"]["n_cores"]
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
            ref = get_filepaths(
                gchp_vs_gcc_refdir,
                "Aerosols",
                all_months_dev
            )[0]
            dev = get_filepaths(
                gchp_vs_gcc_devdir,
                "Aerosols",
                all_months_gchp_dev,
                is_gchp=True
            )[0]

            # Create plots
            print("\nCreating plots for annual mean")
            make_benchmark_aod_plots(
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
                n_job=config["options"]["n_cores"]
            )

            # --------------------------------------------------------------
            # GCHP vs GCC column AOD plots: Seasonal
            # --------------------------------------------------------------
            for mon in range(bmk_n_months):
                print(f"\nCreating plots for {bmk_mon_strs[mon]}")

                # Create plots
                mon_ind = bmk_mon_inds[mon]
                make_benchmark_aod_plots(
                    ref[mon_ind],
                    gchp_vs_gcc_refstr,
                    dev[mon_ind],
                    gchp_vs_gcc_devstr,
                    dst=gchp_vs_gcc_resultsdir,
                    subdst=bmk_mon_yr_strs_dev[mon],
                    weightsdir=config["paths"]["weights_dir"],
                    overwrite=True,
                    spcdb_dir=spcdb_dir,
                    n_job=config["options"]["n_cores"]
                )

        # ==================================================================
        # GCHP vs. GCC drydep plots
        # ==================================================================
        if config["options"]["outputs"]["plot_drydep"]:
            print("\n%%% Creating GCHP vs. GCC drydep plots %%%")

            # --------------------------------------------------------------
            # GCHP vs GCC drydep plots: Annual mean
            # --------------------------------------------------------------

            # Filepaths
            ref = get_filepaths(
                gchp_vs_gcc_refdir,
                "DryDep",
                all_months_dev
            )[0]
            dev = get_filepaths(
                gchp_vs_gcc_devdir,
                "DryDep",
                all_months_gchp_dev,
                is_gchp=True
            )[0]

            # Create plots
            print("\nCreating plots for annual mean")
            make_benchmark_drydep_plots(
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
                n_job=config["options"]["n_cores"],
                varlist=drydepvel_species()
            )

            # --------------------------------------------------------------
            # GCHP vs GCC drydep plots: Seasonal
            # --------------------------------------------------------------
            for mon in range(bmk_n_months):
                print(f"\nCreating plots for {bmk_mon_strs[mon]}")

                # Create plots
                mon_ind = bmk_mon_inds[mon]
                make_benchmark_drydep_plots(
                    ref[mon_ind],
                    gchp_vs_gcc_refstr,
                    dev[mon_ind],
                    gchp_vs_gcc_devstr,
                    dst=gchp_vs_gcc_resultsdir,
                    subdst=bmk_mon_yr_strs_dev[mon],
                    weightsdir=config["paths"]["weights_dir"],
                    overwrite=True,
                    spcdb_dir=spcdb_dir,
                    n_job=config["options"]["n_cores"],
                    varlist=drydepvel_species()
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

                # Create tables
                make_benchmark_mass_tables(
                    refpath,
                    gchp_vs_gcc_refstr,
                    devpath,
                    gchp_vs_gcc_devstr,
                    dst=gchp_vs_gcc_tablesdir,
                    subdst=bmk_mon_yr_strs_dev[mon],
                    label=f"at 01{bmk_mon_yr_strs_dev[mon]}",
                    overwrite=True,
                    spcdb_dir=spcdb_dir,
                    dev_met_extra=devareapath
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
        if config["options"]["outputs"]["ops_budget_table"]:
            print("\n%%% Creating GCHP vs. GCC operations budget tables %%%")

            def gchp_vs_gcc_ops_budg(mon):
                """
                Create operations budgets for each benchmark month m in parallel
                """

                # Filepaths
                refpath = get_filepath(
                    gchp_vs_gcc_refdir,
                    "Budget",
                    bmk_mons_dev[mon]
                )
                devpath = get_filepath(
                    gchp_vs_gcc_devdir,
                    "Budget",
                    bmk_mons_gchp_dev[mon],
                    is_gchp=True
                )

                # Create tables
                make_benchmark_operations_budget(
                    config["data"]["dev"]["gcc"]["version"],
                    refpath,
                    config["data"]["dev"]["gchp"]["version"],
                    devpath,
                    bmk_sec_per_month_dev[mon],
                    bmk_sec_per_month_dev[mon],
                    benchmark_type=bmk_type,
                    label=f"at 01{bmk_mon_yr_strs_dev[mon]}",
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

            # Create tables in parallel
            # Turn off parallelization if n_jobs==1
            if config["options"]["n_cores"] != 1:
                results = Parallel(n_jobs=config["options"]["n_cores"])(
                    delayed(gchp_vs_gcc_ops_budg)(mon) \
                    for mon in range(bmk_n_months)
                )
            else:
                results = []
                for mon in range(bmk_n_months):
                    results.append(gchp_vs_gcc_ops_budg(mon))

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
                is_gchp=True
            )[0]
            devspc = get_filepaths(
                gchp_vs_gcc_devdir,
                "SpeciesConc",
                all_months_gchp_dev,
                is_gchp=True
            )[0]

            # Create tables
            make_benchmark_aerosol_tables(
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
            global_ox_budget(
                config["data"]["dev"]["gcc"]["version"],
                gcc_vs_gcc_devdir,
                gcc_vs_gcc_devrstdir,
                bmk_year_dev,
                dst=gcc_vs_gcc_tablesdir,
                overwrite=True,
                spcdb_dir=spcdb_dir
            )

            # Compute Ox budget table for GCHP
            global_ox_budget(
                config["data"]["dev"]["gchp"]["version"],
                gchp_vs_gcc_devdir,
                gchp_vs_gcc_devrstdir,
                bmk_year_dev,
                dst=gchp_vs_gcc_tablesdir,
                overwrite=True,
                spcdb_dir=spcdb_dir,
                is_gchp=True,
                gchp_res=config["data"]["dev"]["gchp"]["resolution"],
                gchp_is_pre_14_0=config["data"]["dev"]["gchp"]["is_pre_14.0"]
            )

        # ==================================================================
        # GCHP vs. GCC global mean OH, MCF Lifetime, CH4 Lifetime
        # ==================================================================
        if config["options"]["outputs"]["OH_metrics"]:
            print("\n%%% Creating GCHP vs. GCC OH metrics table %%%")

            # Filepaths
            ref = get_filepaths(
                gchp_vs_gcc_refdir,
                "Metrics",
                all_months_dev
            )[0]
            dev = get_filepaths(
                gchp_vs_gcc_devdir,
                "Metrics",
                all_months_gchp_dev,
                is_gchp=True
            )[0]

            # Create table
            make_benchmark_oh_metrics(
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
        # ==================================================================
        if config["options"]["outputs"]["ste_table"]:
            print("\n%%% Skipping GCHP vs. GCC Strat-Trop Exchange table %%%")

        # ==================================================================
        # GCHP vs GCC Model vs. Observations plots
        # ==================================================================
        if config["options"]["outputs"]["plot_models_vs_obs"]:
            print("\n%%% Creating GCHP vs. GCC models vs. obs. plots %%%")

            # Filepaths
            ref = get_filepaths(
                gchp_vs_gcc_refdir,
                "SpeciesConc",
                all_months_dev
            )[0]
            dev = get_filepaths(
                gchp_vs_gcc_devdir,
                "SpeciesConc",
                all_months_gchp_dev,
                is_gchp=True
            )[0]

            # Plot models vs. observations
            make_benchmark_models_vs_obs_plots(
                config["paths"]["obs_data"]["ebas_o3"]["data_dir"],
                config["paths"]["obs_data"]["ebas_o3"]["data_label"],
                ref,
                config["data"]["dev"]["gcc"]["version"],
                dev,
                config["data"]["dev"]["gchp"]["version"],
                dst=gchp_vs_gcc_models_vs_obs_dir,
                overwrite=True,
                verbose=False
            )

            # Plot models vs. sondes
            make_benchmark_models_vs_sondes_plots(
                sondes_data_file,
                sondes_site_file,
                ref,
                config["data"]["dev"]["gcc"]["version"],
                dev,
                config["data"]["dev"]["gchp"]["version"],
                dst=gchp_vs_gcc_models_vs_obs_dir,
                overwrite=True,
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
            "StateMet",
            all_months_gchp_ref,
            is_gchp=True
        )[0]
        devmet = get_filepaths(
            gchp_vs_gcc_devdir,
            "StateMet",
            all_months_gchp_dev,
            is_gchp=True
        )[0]

        # Get GCHP grid resolutions from met collection file
        #ds_refmet = xr.open_dataset(refmet[0])
        #ds_devmet = xr.open_dataset(devmet[0])
        #gchp_ref_res = str(get_input_res(ds_refmet)[0])
        #gchp_dev_res = str(get_input_res(ds_devmet)[0])

        # Option to specify grid resolution of comparison plots
        # This is a temporary hack until cs->cs regridding in GCPy is fixed
        cmpres="1x1.25"

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
                is_gchp=True
            )[0]
            dev = get_filepaths(
                gchp_vs_gchp_devdir,
                "SpeciesConc",
                all_months_gchp_dev,
                is_gchp=True
            )[0]

            # Create plots
            print("\nCreating plots for annual mean")
            make_benchmark_conc_plots(
                ref,
                gchp_vs_gchp_refstr,
                dev,
                gchp_vs_gchp_devstr,
                refmet=refmet,
                devmet=devmet,
                cmpres=cmpres,
                dst=gchp_vs_gchp_resultsdir,
                subdst="AnnualMean",
                time_mean=True,
                weightsdir=config["paths"]["weights_dir"],
                benchmark_type=bmk_type,
                plot_by_spc_cat=config["options"]["outputs"][
                    "plot_options"]["by_spc_cat"],
                overwrite=True,
                spcdb_dir=spcdb_dir,
                n_job=config["options"]["n_cores"]
            )

            # --------------------------------------------------------------
            # GCHP vs GCHP species concentration plots: Seasonal
            # --------------------------------------------------------------
            for mon in range(bmk_n_months):
                print(f"\nCreating plots for {bmk_mon_strs[mon]}")

                # Create plots
                mon_ind = bmk_mon_inds[mon]
                make_benchmark_conc_plots(
                    ref[mon_ind],
                    gchp_vs_gchp_refstr,
                    dev[mon_ind],
                    gchp_vs_gchp_devstr,
                    cmpres=cmpres,
                    refmet=refmet[mon_ind],
                    devmet=devmet[mon_ind],
                    dst=gchp_vs_gchp_resultsdir,
                    subdst=bmk_mon_yr_strs_dev[mon],
                    weightsdir=config["paths"]["weights_dir"],
                    benchmark_type=bmk_type,
                    plot_by_spc_cat=config["options"]["outputs"][
                        "plot_options"]["by_spc_cat"],
                    overwrite=True,
                    spcdb_dir=spcdb_dir,
                    n_job=config["options"]["n_cores"]
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
                is_gchp=True
            )[0]
            dev = get_filepaths(
                gchp_vs_gchp_devdir,
                "Emissions",
                all_months_gchp_dev,
                is_gchp=True
            )[0]

            # Create plots
            print("\nCreating plots for annual mean")
            make_benchmark_emis_plots(
                ref,
                gchp_vs_gchp_refstr,
                dev,
                gchp_vs_gchp_devstr,
                dst=gchp_vs_gchp_resultsdir,
                subdst="AnnualMean",
                cmpres=cmpres,
                time_mean=True,
                weightsdir=config["paths"]["weights_dir"],
                plot_by_spc_cat=config["options"]["outputs"][
                    "plot_options"]["by_spc_cat"],
                plot_by_hco_cat=config["options"]["outputs"][
                    "plot_options"]["by_hco_cat"],
                benchmark_type=bmk_type,
                overwrite=True,
                spcdb_dir=spcdb_dir,
                n_job=config["options"]["n_cores"]
            )

            # --------------------------------------------------------------
            # GCHP vs GCHP species concentration plots: Seasonal
            # --------------------------------------------------------------
            for mon in range(bmk_n_months):
                print(f"\nCreating plots for {bmk_mon_strs[mon]}")

                # Create plots
                mon_ind = bmk_mon_inds[mon]
                make_benchmark_emis_plots(
                    ref[mon_ind],
                    gchp_vs_gchp_refstr,
                    dev[mon_ind],
                    gchp_vs_gchp_devstr,
                    dst=gchp_vs_gchp_resultsdir,
                    cmpres=cmpres,
                    subdst=bmk_mon_yr_strs_dev[mon],
                    weightsdir=config["paths"]["weights_dir"],
                    plot_by_spc_cat=config["options"]["outputs"][
                        "plot_options"]["by_spc_cat"],
                    plot_by_hco_cat=config["options"]["outputs"][
                        "plot_options"]["by_hco_cat"],
                    benchmark_type=bmk_type,
                    overwrite=True,
                    spcdb_dir=spcdb_dir,
                    n_job=config["options"]["n_cores"]
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
                is_gchp=True
            )[0]
            dev = get_filepaths(
                gchp_vs_gchp_devdir,
                "Emissions",
                all_months_gchp_dev,
                is_gchp=True
            )[0]

            # Create table
            make_benchmark_emis_tables(
                ref,
                gchp_vs_gchp_refstr,
                dev,
                gchp_vs_gchp_devstr,
                refmet=refmet,
                devmet=devmet,
                dst=gchp_vs_gchp_resultsdir,
                ref_interval=sec_per_month_ref,
                dev_interval=sec_per_month_dev,
                benchmark_type=bmk_type,
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
                is_gchp=True
            )[0]
            dev = get_filepaths(
                gchp_vs_gchp_devdir,
                "JValues",
                all_months_gchp_dev,
                is_gchp=True
            )[0]

            # Create plots
            print("\nCreating plots for annual mean")
            make_benchmark_jvalue_plots(
                ref,
                gchp_vs_gchp_refstr,
                dev,
                gchp_vs_gchp_devstr,
                dst=gchp_vs_gchp_resultsdir,
                subdst='AnnualMean',
                cmpres=cmpres,
                time_mean=True,
                weightsdir=config["paths"]["weights_dir"],
                overwrite=True,
                spcdb_dir=spcdb_dir,
                n_job=config["options"]["n_cores"]
            )

            # --------------------------------------------------------------
            # GCHP vs GCHP J-values plots: Seasonal
            # --------------------------------------------------------------
            for mon in range(bmk_n_months):
                print(f"\nCreating plots for {bmk_mon_strs[mon]}")

                # Create plots
                mon_ind = bmk_mon_inds[mon]
                make_benchmark_jvalue_plots(
                    ref[mon_ind],
                    gchp_vs_gchp_refstr,
                    dev[mon_ind],
                    gchp_vs_gchp_devstr,
                    dst=gchp_vs_gchp_resultsdir,
                    subdst=bmk_mon_yr_strs_dev[mon],
                    cmpres=cmpres,
                    weightsdir=config["paths"]["weights_dir"],
                    overwrite=True,
                    spcdb_dir=spcdb_dir,
                    n_job=config["options"]["n_cores"]
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
                is_gchp=True
            )[0]
            dev = get_filepaths(
                gchp_vs_gchp_devdir,
                "Aerosols",
                all_months_gchp_dev,
                is_gchp=True
            )[0]

            # Create plots
            print("\nCreating plots for annual mean")
            make_benchmark_aod_plots(
                ref,
                gchp_vs_gchp_refstr,
                dev,
                gchp_vs_gchp_devstr,
                dst=gchp_vs_gchp_resultsdir,
                subdst="AnnualMean",
                cmpres=cmpres,
                time_mean=True,
                weightsdir=config["paths"]["weights_dir"],
                overwrite=True,
                spcdb_dir=spcdb_dir,
                n_job=config["options"]["n_cores"]
            )

            # --------------------------------------------------------------
            # GCHP vs GCHP column AOD plots: Seasonal
            # --------------------------------------------------------------
            for mon in range(bmk_n_months):
                print(f"\nCreating plots for {bmk_mon_strs[mon]}")

                # Create plots
                mon_ind = bmk_mon_inds[mon]
                make_benchmark_aod_plots(
                    ref[mon_ind],
                    gchp_vs_gchp_refstr,
                    dev[mon_ind],
                    gchp_vs_gchp_devstr,
                    dst=gchp_vs_gchp_resultsdir,
                    subdst=bmk_mon_yr_strs_dev[mon],
                    cmpres=cmpres,
                    weightsdir=config["paths"]["weights_dir"],
                    overwrite=True,
                    spcdb_dir=spcdb_dir,
                    n_job=config["options"]["n_cores"]
                )

        # ==================================================================
        # GCHP vs GCHP drydep plots
        # ==================================================================
        if config["options"]["outputs"]["plot_drydep"]:
            print("\n%%% Creating GCHP vs. GCHP drydep plots %%%")

            # --------------------------------------------------------------
            # GCHP vs GCHP drydep: Annual Mean
            # --------------------------------------------------------------

            # Filepaths
            ref = get_filepaths(
                gchp_vs_gchp_refdir,
                "DryDep",
                all_months_gchp_ref,
                is_gchp=True
            )[0]
            dev = get_filepaths(
                gchp_vs_gchp_devdir,
                "DryDep",
                all_months_gchp_dev,
                is_gchp=True
            )[0]

            # Create plots
            print("\nCreating plots for annual mean")
            make_benchmark_drydep_plots(
                ref,
                gchp_vs_gchp_refstr,
                dev,
                gchp_vs_gchp_devstr,
                dst=gchp_vs_gchp_resultsdir,
                subdst="AnnualMean",
                cmpres=cmpres,
                time_mean=True,
                weightsdir=config["paths"]["weights_dir"],
                overwrite=True,
                spcdb_dir=spcdb_dir,
                n_job=config["options"]["n_cores"],
                varlist=drydepvel_species()
            )

            # --------------------------------------------------------------
            # GCHP vs GCHP drydep plots: Seasonal
            # --------------------------------------------------------------
            for mon in range(bmk_n_months):
                print(f"\nCreating plots for {bmk_mon_strs[mon]}")

                # Create plots
                mon_ind = bmk_mon_inds[mon]
                make_benchmark_drydep_plots(
                    ref[mon_ind],
                    gchp_vs_gchp_refstr,
                    dev[mon_ind],
                    gchp_vs_gchp_devstr,
                    dst=gchp_vs_gchp_resultsdir,
                    subdst=bmk_mon_yr_strs_dev[mon],
                    cmpres=cmpres,
                    weightsdir=config["paths"]["weights_dir"],
                    overwrite=True,
                    spcdb_dir=spcdb_dir,
                    n_job=config["options"]["n_cores"],
                    varlist=drydepvel_species()
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
                    gchp_res=config["data"]["ref"]["gchp"]["resolution"],
                    gchp_is_pre_14_0=config["data"]["ref"]["gchp"][
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
                make_benchmark_mass_tables(
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
        if config["options"]["outputs"]["ops_budget_table"]:
            print("\n%%% Creating GCHP vs. GCHP operations budget tables %%%")

            # Diagnostic collections to read
            def gchp_vs_gchp_ops_budg(mon):
                """
                Creates operations budgets for each benchmark month m in parallel
                """

                # Filepaths
                refpath = get_filepath(
                    gchp_vs_gchp_refdir,
                    "Budget",
                    bmk_mons_gchp_ref[mon],
                    is_gchp=True
                )
                devpath = get_filepath(
                    gchp_vs_gchp_devdir,
                    "Budget",
                    bmk_mons_gchp_dev[mon],
                    is_gchp=True
                )

                # Compute tables
                make_benchmark_operations_budget(
                    config["data"]["ref"]["gchp"]["version"],
                    refpath,
                    config["data"]["dev"]["gchp"]["version"],
                    devpath,
                    bmk_sec_per_month_ref[mon],
                    bmk_sec_per_month_dev[mon],
                    benchmark_type=bmk_type,
                    label=f"at 01{bmk_mon_yr_strs_dev[mon]}",
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

            # Create tables in parallel
            # Turn off parallelization if n_jobs==1
            if config["options"]["n_cores"] != 1:
                results = Parallel(n_jobs=config["options"]["n_cores"])(
                    delayed(gchp_vs_gchp_ops_budg)(mon) \
                    for mon in range(bmk_n_months)
                )
            else:
                results = []
                for mon in range(bmk_n_months):
                    results.append(gchp_vs_gchp_ops_budg(mon))

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
                is_gchp=True
            )[0]
            devspc = get_filepaths(
                gchp_vs_gchp_devdir,
                "SpeciesConc",
                all_months_gchp_dev,
                is_gchp=True
            )[0]

            # Create tables
            make_benchmark_aerosol_tables(
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
            global_ox_budget(
                config["data"]["dev"]["gchp"]["version"],
                gchp_vs_gchp_devdir,
                gchp_vs_gchp_devrstdir,
                bmk_year_dev,
                dst=gchp_vs_gchp_tablesdir,
                overwrite=True,
                spcdb_dir=spcdb_dir,
                is_gchp=True,
                gchp_res=config["data"]["dev"]["gchp"]["resolution"],
                gchp_is_pre_14_0=config["data"]["dev"]["gchp"]["is_pre_14.0"]
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
                is_gchp=True
            )[0]
            dev = get_filepaths(
                gchp_vs_gchp_devdir,
                "Metrics",
                all_months_gchp_dev,
                is_gchp=True
            )[0]

            # Create the OH Metrics table
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
        # GCHP Strat-Trop Exchange
        # ==================================================================
        if config["options"]["outputs"]["ste_table"]:
            print("\n%%% Skipping GCHP vs. GCHP Strat-Trop Exchange table %%%")

        # ==================================================================
        # GCHP vs. GCHP Benchmark Timing Table
        # ==================================================================
        if config["options"]["outputs"]["timing_table"]:
            print("\n%%% Creating GCHP vs. GCHP Benchmark Timing table %%%")

            # Filepaths
            # NOTE: Usually the GCHP 1-yr benchmark is run as
            # one job, so we only need to take the 1st log file.
            ref = get_log_filepaths(
                gchp_vs_gchp_reflogdir,
                config["data"]["ref"]["gchp"]["logs_template"],
                all_months_gchp_ref,
            )[0]
            dev = get_log_filepaths(
                gchp_vs_gchp_devlogdir,
                config["data"]["dev"]["gchp"]["logs_template"],
                all_months_gchp_dev,
            )[0]

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
        # GCHP vs GCHP Model vs. Observations plots
        # ==================================================================
        if config["options"]["outputs"]["plot_models_vs_obs"]:
            print("\n%%% Creating GCHP vs. GCHP models vs. obs. plots %%%")

            # Filepaths
            # NOTE: If the GCHP benchmark is done in one-shot
            # then you need the [0] after the call to get_filepaths.
            ref = get_filepaths(
                gchp_vs_gchp_refdir,
                "SpeciesConc",
                all_months_gchp_ref,
                is_gchp=True
            )[0]
            dev = get_filepaths(
                gchp_vs_gchp_devdir,
                "SpeciesConc",
                all_months_gchp_dev,
                is_gchp=True
            )[0]

            # Plot models vs. observations
            make_benchmark_models_vs_obs_plots(
                config["paths"]["obs_data"]["ebas_o3"]["data_dir"],
                config["paths"]["obs_data"]["ebas_o3"]["data_label"],
                ref,
                config["data"]["ref"]["gchp"]["version"],
                dev,
                config["data"]["dev"]["gchp"]["version"],
                dst=gchp_vs_gchp_models_vs_obs_dir,
                overwrite=True,
                verbose=False
            )

            # Plot models vs. sondes
            make_benchmark_models_vs_sondes_plots(
                sondes_data_file,
                sondes_site_file,
                ref,
                config["data"]["ref"]["gchp"]["version"],
                dev,
                config["data"]["dev"]["gchp"]["version"],
                dst=gchp_vs_gchp_models_vs_obs_dir,
                overwrite=True,
            )

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Create GCHP vs GCC difference of differences benchmark plots
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if config["options"]["comparisons"]["gchp_vs_gcc_diff_of_diffs"]["run"]:

        if config["options"]["outputs"]["plot_conc"]:
            print("\n%%% Creating GCHP vs. GCC diff-of-diffs conc plots %%%")


            # --------------------------------------------------------------
            # GCHP vs GCC diff-of-diff species concentration plots:
            # Annual Mean
            # --------------------------------------------------------------

            # Filepaths
            gcc_ref = get_filepaths(
                gcc_vs_gcc_refdir,
                "SpeciesConc",
                all_months_ref
            )[0]
            gcc_dev = get_filepaths(
                gcc_vs_gcc_devdir,
                "SpeciesConc",
                all_months_dev
            )[0]
            gchp_ref = get_filepaths(
                gchp_vs_gchp_refdir,
                "SpeciesConc",
                all_months_gchp_ref,
                is_gchp=True
            )[0]
            gchp_dev = get_filepaths(
                gchp_vs_gchp_devdir,
                "SpeciesConc",
                all_months_gchp_dev,
                is_gchp=True
            )[0]

            # Create plots
            print("\nCreating plots for annual mean")
            make_benchmark_conc_plots(
                gcc_ref,
                diff_of_diffs_refstr,
                gchp_ref,
                diff_of_diffs_devstr,
                dst=diff_of_diffs_resultsdir,
                subdst="AnnualMean",
                time_mean=True,
                weightsdir=config["paths"]["weights_dir"],
                benchmark_type=bmk_type,
                plot_by_spc_cat=config["options"]["outputs"][
                    "plot_options"]["by_spc_cat"],
                overwrite=True,
                use_cmap_RdBu=True,
                second_ref=gcc_dev,
                second_dev=gchp_dev,
                cats_in_ugm3=None,
                spcdb_dir=spcdb_dir,
                n_job=config["options"]["n_cores"]
            )

            # --------------------------------------------------------------
            # GCHP vs GCC diff-of-diff species concentration plots: Seasonal
            # --------------------------------------------------------------
            for mon in range(bmk_n_months):
                print(f"\nCreating plots for {bmk_mon_strs[mon]}")

                # Create plots
                mon_ind = bmk_mon_inds[mon]
                make_benchmark_conc_plots(
                    gcc_ref[mon_ind],
                    diff_of_diffs_refstr,
                    gchp_ref[mon_ind],
                    diff_of_diffs_devstr,
                    dst=diff_of_diffs_resultsdir,
                    subdst=bmk_mon_yr_strs_dev[mon],
                    weightsdir=config["paths"]["weights_dir"],
                    benchmark_type=bmk_type,
                    plot_by_spc_cat=config["options"]["outputs"][
                        "plot_options"]["by_spc_cat"],
                    overwrite=True,
                    use_cmap_RdBu=True,
                    second_ref=gcc_dev[mon_ind],
                    second_dev=gchp_dev[mon_ind],
                    cats_in_ugm3=None,
                    spcdb_dir=spcdb_dir,
                    n_job=config["options"]["n_cores"]
                )

    # ==================================================================
    # Print a message indicating that the benchmarks finished
    # ==================================================================
    print("\n%%%% All requested benchmark plots/tables created! %%%%")
