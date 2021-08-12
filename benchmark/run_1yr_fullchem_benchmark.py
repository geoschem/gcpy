#!/usr/bin/env python
"""
run_1yr_fullchem_benchmark.py: Driver script for creating benchmark plots and
                               testing gcpy 1-year full-chemistry benchmark
                               capability.

Run this script to generate benchmark comparisons between:

    (1) GCC (aka GEOS-Chem "Classic") vs. GCC
    (2) GCHP vs GCC
    (3) GCHP vs GCHP

You can customize this by editing the following settings in the
"Configurables" section below:

    (1) Edit the path variables so that they point to folders w/ model data
    (2) Edit the version strings for each benchmark simulation
    (3) Edit the switches that turn on/off creating of plots and tables
    (4) If necessary, edit labels for the dev and ref versions

Calling sequence:

    ./run_1yr_fullchem_benchmark.py

To test gcpy, copy this script anywhere you want to run the test and
set gcpy_test to True at the top of the script. Benchmark artifacts will
be created locally in new folder called Plots.

Remarks:

    By default, matplotlib will try to open an X window for plotting.
    If you are running this script in an environment where you do not have
    an active X display (such as in a computational queue), then you will
    need to use these commands to disable the X-window functionality.

        import os
        os.environ["QT_QPA_PLATFORM"]="offscreen"

    For more information, please see this issue posted at the ipython site:

        https://github.com/ipython/ipython/issues/10627

This script corresponds with GCPy 1.0.4. Edit this version ID if releasing
a new version of GCPy.
"""

# =====================================================================
# Imports and global settings (you should not need to edit these)
# =====================================================================

import os
from os.path import join, exists
import warnings
from shutil import copyfile
from calendar import monthrange
import numpy as np
from joblib import Parallel, delayed
from gcpy.util import get_filepath, get_filepaths
import gcpy.ste_flux as ste
import gcpy.oh_metrics as oh
import gcpy.budget_ox as ox
#import gcpy.mean_oh_from_logs as moh  # NOTE: to be removed after 13.0.0
from gcpy import benchmark as bmk

# Tell matplotlib not to look for an X-window
os.environ["QT_QPA_PLATFORM"] = "offscreen"

# Suppress annoying warning messages
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

# This script has a fixed benchmark type
bmk_type     = "FullChemBenchmark"
bmk_year_ref = '2019'
bmk_year_dev = '2019'
bmk_mon_strs = ["Jan", "Apr", "Jul", "Oct"]
bmk_mon_inds = [0, 3, 6, 9]
bmk_n_months = len(bmk_mon_strs)

########################################################################
###           CONFIGURABLE SETTINGS: ***EDIT AS NEEDED ***           ###
########################################################################

# =====================================================================
# Benchmark information
# Note: When doing GCHP vs GCC comparisions gchp_dev will be compared
# to gcc_dev (not gcc_ref!).
# =====================================================================

# High-level directory containing subdirectories with data
maindir  = "/n/holyscratch01/external_repos/GEOS-CHEM/gcgrid/geos-chem/validation/gcpy_test_data/1yr_fullchem"

# Version strings
# NOTE: these will be used in some filenames and so should not have spaces
# or other characters not appropriate for a filename.
gcc_ref_version = "GCC_ref"
gcc_dev_version = "GCC_dev"
gchp_ref_version = "GCHP_ref"
gchp_dev_version = "GCHP_dev"

# Name to be used for directory of output from this script
results_dir = "BenchmarkResults"

# Path to regridding weights
weightsdir = "/n/holylfs/EXTERNAL_REPOS/GEOS-CHEM/gcgrid/gcdata/ExtData/GCHP/RegriddingWeights"

# Path to species_databse.yml
spcdb_dir = join(maindir, gcc_dev_version)

# GCHP initial restart resolution (for mass tables)
gchp_ref_res = 'c48'
gchp_dev_res = 'c48'

# Kludge: Set switches that will pick the proper StateMet collection
# for GCHP.  Versions prior to 13.0.0 used StateMet_avg.
gchp_ref_prior_to_13 = False
gchp_dev_prior_to_13 = False

# Whether GCHP files are legacy (pre-13.1) format
gchp_ref_is_legacy=True
gchp_dev_is_legacy=False

# ======================================================================
# Specify if this is a gcpy test validation run
# ======================================================================
gcpy_test = True

# ======================================================================
# Comparisons to run
# ======================================================================
gcc_vs_gcc   = True
gchp_vs_gcc  = True
gchp_vs_gchp = True
# GCHP vs GCC diff of diffs not included in 1-yr full chemistry benchmark

# ======================================================================
# Output to generate (plots/tables will be created in this order):
# ======================================================================
plot_conc        = True
plot_emis        = True
emis_table       = True
plot_jvalues     = True
plot_aod         = True
mass_table       = True
ops_budget_table = False
aer_budget_table = True
Ox_budget_table  = True
ste_table        = True # GCC only
OH_metrics       = True

# Plot concentrations and emissions by category?
plot_by_spc_cat = True
plot_by_hco_cat = True

# ======================================================================
# Data directories
# For gchp_vs_gcc_refdir use gcc_dev_version, not ref (mps, 6/27/19)
# ======================================================================

# Diagnostics file directory paths
gcc_vs_gcc_refdir = join(maindir, gcc_ref_version, "OutputDir")
gcc_vs_gcc_devdir = join(maindir, gcc_dev_version, "OutputDir")
gchp_vs_gcc_refdir = join(maindir, gcc_dev_version, "OutputDir")
gchp_vs_gcc_devdir = join(maindir, gchp_dev_version, "OutputDir")
gchp_vs_gchp_refdir = join(maindir, gchp_ref_version, "OutputDir")
gchp_vs_gchp_devdir = join(maindir, gchp_dev_version, "OutputDir")

# Restart file directory paths
gcc_vs_gcc_refrstdir = join(maindir, gcc_ref_version, "restarts")
gcc_vs_gcc_devrstdir = join(maindir, gcc_dev_version, "restarts")
gchp_vs_gcc_refrstdir = join(maindir, gcc_dev_version, "restarts")
gchp_vs_gcc_devrstdir = join(maindir, gchp_dev_version)
gchp_vs_gchp_refrstdir = join(maindir, gchp_ref_version)
gchp_vs_gchp_devrstdir = join(maindir, gchp_dev_version)

# Log file directories -- GEOS-Chem "Classic" only
gcc_vs_gcc_reflogdir = join(maindir, gcc_ref_version, "logs")
gcc_vs_gcc_devlogdir = join(maindir, gcc_dev_version, "logs")

# ======================================================================
# Benchmark output directories
# ======================================================================
# Plot directories
if gcpy_test:
    mainresultsdir = join(
        '.', 
        results_dir
    )
    gcc_vs_gcc_resultsdir = join(
        mainresultsdir,
        'GCC_version_comparison'
    )
    gchp_vs_gchp_resultsdir = join(
        mainresultsdir,
        'GCHP_version_comparison'
    )
    gchp_vs_gcc_resultsdir = join(
        mainresultsdir,
        'GCHP_GCC_comparison'
    )
    if not exists(mainresultsdir):
        os.mkdir(mainresultsdir)
    # Make copy of benchmark script in results directory
    curfile = os.path.realpath(__file__)
    dest = join(mainresultsdir,curfile.split('/')[-1])
    if not exists(dest):
        copyfile(curfile, dest)

else:
    gcc_vs_gcc_resultsdir = join(
        maindir,
        gcc_dev_version,
        results_dir
    )
    gchp_vs_gchp_resultsdir = join(
        maindir, 
        gchp_dev_version,
        results_dir, 
        "GCHP_version_comparison"
    )
    gchp_vs_gcc_resultsdir = join(
        maindir,
        gchp_dev_version,
        results_dir, 
        "GCHP_GCC_comparison"
    )
    base_gchp_resultsdir = join(
        maindir, 
        gchp_dev_version, 
        results_dir
    )

    #make results directories that don't exist
    for resdir, plotting_type in zip(
            [
                gcc_vs_gcc_resultsdir,
                base_gchp_resultsdir,
                gchp_vs_gchp_resultsdir,
                gchp_vs_gcc_resultsdir
            ],
            [
                gcc_vs_gcc,
                gchp_vs_gcc or gchp_vs_gchp,
                gchp_vs_gchp,
                gchp_vs_gcc
            ]
    ):
        if plotting_type and not exists(resdir):
            os.mkdir(resdir)
        if resdir in [gcc_vs_gcc_resultsdir, base_gchp_resultsdir]:
            # Make copy of benchmark script in results directory
            curfile = os.path.realpath(__file__)
            dest = join(resdir,curfile.split('/')[-1])
            if not exists(dest):
                copyfile(curfile, dest)

# Tables directories
gcc_vs_gcc_tablesdir = join(gcc_vs_gcc_resultsdir, "Tables")
gchp_vs_gcc_tablesdir = join(gchp_vs_gcc_resultsdir, "Tables")
gchp_vs_gchp_tablesdir = join(gchp_vs_gchp_resultsdir, "Tables")

# Budget directories
gcc_vs_gcc_budgetdir = join(gcc_vs_gcc_resultsdir, "Budget")
gchp_vs_gcc_budgetdir = join(gchp_vs_gcc_resultsdir, "Budget")
gchp_vs_gchp_budgetdir = join(gchp_vs_gchp_resultsdir, "Budget")

# ======================================================================
# Plot title strings
# For gchp_vs_gcc_refstr use gcc_dev_version, not ref (mps, 6/27/19)
# ======================================================================
gcc_vs_gcc_refstr    = gcc_ref_version
gcc_vs_gcc_devstr    = gcc_dev_version
gchp_vs_gcc_refstr   = gcc_dev_version
gchp_vs_gcc_devstr   = gchp_dev_version
gchp_vs_gchp_refstr  = gchp_ref_version
gchp_vs_gchp_devstr  = gchp_dev_version

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
bmk_end_ref = np.datetime64("{}-01-01".format(int(bmk_year_ref)+1))
all_months_ref = np.arange(bmk_start_ref,
                           bmk_end_ref,
                           step=np.timedelta64(1, "M"),
                           dtype="datetime64[M]")
all_months_gchp_ref = all_months_ref

# Reset all months datetime array if GCHP ref is legacy filename format.
# Legacy format uses time-averaging period mid-point not start.
if gchp_ref_is_legacy:
    all_months_gchp_ref = np.zeros(12, dtype="datetime64[h]")
    for t in range(12):
        middle_hr = int(days_per_month_ref[t] * 24 / 2)
        delta = np.timedelta64(middle_hr, 'h')
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
bmk_end_dev = np.datetime64("{}-01-01".format(int(bmk_year_dev)+1))
all_months_dev = np.arange(bmk_start_dev,
                           bmk_end_dev,
                           step=np.timedelta64(1, "M"),
                           dtype="datetime64[M]")
all_months_gchp_dev = all_months_dev

# Reset all months datetime array if GCHP dev is legacy filename format.
# Legacy format uses time-averaging period mid-point not start.
if gchp_dev_is_legacy:
    all_months_gchp_dev = np.zeros(12, dtype="datetime64[h]")
    for t in range(12):
        middle_hr = int(days_per_month_dev[t] * 24 / 2)
        delta = np.timedelta64(middle_hr, 'h')
        all_months_gchp_dev[t] = all_months_dev[t].astype("datetime64[h]") + delta

# Get subset of month datetimes and seconds per month for only benchmark months
bmk_mons_dev = all_months_dev[bmk_mon_inds]
bmk_mons_gchp_dev = all_months_gchp_dev[bmk_mon_inds]
bmk_sec_per_month_dev = sec_per_month_dev[bmk_mon_inds]

# ======================================================================
# Print the list of plots & tables to the screen
# ======================================================================

print("The following plots and tables will be created for {}:".format(
    bmk_type
))
if plot_conc:
    print(" - Concentration plots")
if plot_emis:
    print(" - Emissions plots")
if plot_jvalues:
    print(" - J-values (photolysis rates) plots")
if plot_aod:
    print(" - Aerosol optical depth plots")
if ops_budget_table:
    print(" - Operations budget tables")
if aer_budget_table:
    print(" - Aerosol budget/burden tables")
if emis_table:
    print(" - Table of emissions totals by species and inventory")
if mass_table:
    print(" - Table of species mass")
if OH_metrics:
    print(" - Table of OH metrics")
if ste_table:
    print(" - Table of strat-trop exchange")
print("Comparisons will be made for the following combinations:")
if gcc_vs_gcc:
    print(" - GCC vs GCC")
if gchp_vs_gcc:
    print(" - GCHP vs GCC")
if gchp_vs_gchp:
    print(" - GCHP vs GCHP")

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Create GCC vs GCC benchmark plots and tables
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if gcc_vs_gcc:

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
    if plot_conc:
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
        bmk.make_benchmark_conc_plots(
            ref,
            gcc_vs_gcc_refstr,
            dev,
            gcc_vs_gcc_devstr,
            refmet=refmet,
            devmet=devmet,
            dst=gcc_vs_gcc_resultsdir,
            subdst='AnnualMean',
            time_mean=True,
            weightsdir=weightsdir,
            benchmark_type=bmk_type,
            overwrite=True,
            spcdb_dir=spcdb_dir
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
                weightsdir=weightsdir,
                benchmark_type=bmk_type,
                plot_by_spc_cat=plot_by_spc_cat,
                overwrite=True,
                spcdb_dir=spcdb_dir
            )

    # ==================================================================
    # GCC vs GCC emissions plots
    # ==================================================================
    if plot_emis:
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
        bmk.make_benchmark_emis_plots(
            ref,
            gcc_vs_gcc_refstr,
            dev,
            gcc_vs_gcc_devstr,
            dst=gcc_vs_gcc_resultsdir,
            subdst='AnnualMean',
            time_mean=True,
            weightsdir=weightsdir,
            plot_by_spc_cat=plot_by_spc_cat,
            plot_by_hco_cat=plot_by_hco_cat,
            overwrite=True,
            spcdb_dir=spcdb_dir
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
                weightsdir=weightsdir,
                plot_by_spc_cat=plot_by_spc_cat,
                plot_by_hco_cat=plot_by_hco_cat,
                overwrite=True,
                spcdb_dir=spcdb_dir
            )

    # ==================================================================
    # GCC vs GCC tables of emission and inventory totals
    # ==================================================================
    if emis_table:
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
        bmk.make_benchmark_emis_tables(
            ref,
            gcc_vs_gcc_refstr,
            dev,
            gcc_vs_gcc_devstr,
            dst=gcc_vs_gcc_resultsdir,
            ref_interval=sec_per_month_ref,
            dev_interval=sec_per_month_dev,
            overwrite=True,
            spcdb_dir=spcdb_dir
        )

    # ==================================================================
    # GCC vs GCC J-value plots
    # ==================================================================
    if plot_jvalues:
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
        bmk.make_benchmark_jvalue_plots(
            ref,
            gcc_vs_gcc_refstr,
            dev,
            gcc_vs_gcc_devstr,
            dst=gcc_vs_gcc_resultsdir,
            subdst='AnnualMean',
            time_mean=True,
            weightsdir=weightsdir,
            overwrite=True,
            spcdb_dir=spcdb_dir
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
                weightsdir=weightsdir,
                overwrite=True,
                spcdb_dir=spcdb_dir
            )

    # ==================================================================
    # GCC vs. GCC column AOD plots
    # ==================================================================
    if plot_aod:
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
        bmk.make_benchmark_aod_plots(
            ref,
            gcc_vs_gcc_refstr,
            dev,
            gcc_vs_gcc_devstr,
            dst=gcc_vs_gcc_resultsdir,
            subdst='AnnualMean',
            time_mean=True,
            weightsdir=weightsdir,
            overwrite=True,
            spcdb_dir=spcdb_dir
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
                weightsdir=weightsdir,
                overwrite=True,
                spcdb_dir=spcdb_dir
            )

    # ==================================================================
    # GCC vs GCC mass tables
    # ==================================================================
    if mass_table:
        print("\n%%% Creating GCC vs. GCC mass tables %%%")

        def gcc_vs_gcc_mass_table(m):
            """
            Create mass table for each benchmark month m in parallel
            """

            # Filepaths
            refpath = get_filepath(
                gcc_vs_gcc_refrstdir,
                "Restart",
                bmk_mons_ref[m]
            )
            devpath = get_filepath(
                gcc_vs_gcc_devrstdir,
                "Restart",
                bmk_mons_dev[m]
            )

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
                spcdb_dir=spcdb_dir
            )

        # Run in parallel
        results = Parallel(n_jobs=-1)\
            (delayed(gcc_vs_gcc_mass_table)(t) for t in range(bmk_n_months))

    # ==================================================================
    # GCC vs GCC operations budgets tables
    # ==================================================================
    if ops_budget_table:
        print("\n%%% Creating GCC vs. GCC operations budget tables %%%")

        def gcc_vs_gcc_ops_budg(m):
            """
            Create budget table for each benchmark month m in parallel
            """

            # Filepaths
            refpath = get_filepath(
                gcc_vs_gcc_refdir,
                "Budget",
                bmk_mons_ref[m]
            )
            devpath = get_filepath(
                gcc_vs_gcc_devdir,
                "Budget",
                bmk_mons_dev[m]
            )

            # Create tables
            bmk.make_benchmark_operations_budget(
                gcc_ref_version,
                refpath,
                gcc_dev_version,
                devpath,
                sec_per_month_ref[m],
                sec_per_month_dev[m],
                benchmark_type=bmk_type,
                label="at 01{}".format(bmk_mon_yr_strs_dev[m]),
                dst=gcc_vs_gcc_tablesdir
            )

        # Run in parallel
        results = Parallel(n_jobs=-1)\
            (delayed(gcc_vs_gcc_ops_budg)(t) for t in range(bmk_n_months))

    # ==================================================================
    # GCC vs GCC aerosols budgets/burdens tables
    # ==================================================================
    if aer_budget_table:
        print("\n%%% Creating GCC vs. GCC aerosols budget tables %%%")

        # Filepaths
        devaero = get_filepaths(
            gcc_vs_gcc_devdir,
            "Aerosols",
            all_months_dev
        )
        devspc = get_filepaths(
            gcc_vs_gcc_devdir,
            "SpeciesConc",
            all_months_dev
        )

        # Compute tables
        bmk.make_benchmark_aerosol_tables(
            gcc_vs_gcc_devdir,
            devaero,
            devspc,
            devmet,
            gcc_dev_version,
            bmk_year_dev,
            days_per_month_dev,
            dst=gcc_vs_gcc_tablesdir,
            overwrite=True,
            spcdb_dir=spcdb_dir
        )

    # ==================================================================
    # GCC vs GCC Ox budget table
    # ==================================================================
    if Ox_budget_table:
        print("\n%%% Creating GCC vs. GCC Ox budget table %%%")
        ox.global_ox_budget(
            gcc_dev_version,
            gcc_vs_gcc_devdir,
            gcc_vs_gcc_devrstdir,
            bmk_year_dev,
            dst=gcc_vs_gcc_tablesdir,
            overwrite=True,
            spcdb_dir=spcdb_dir
        )

    # ==================================================================
    # GCC Strat-Trop Exchange
    # ==================================================================
    if ste_table:
        print("\n%%% Creating GCC vs. GCC Strat-Trop Exchange table %%%")

        # Filepaths
        dev = get_filepaths(
            gcc_vs_gcc_devdir,
            "AdvFluxVert",
            all_months_dev
        )[0]

        # Compute table
        ste.make_benchmark_ste_table(
            gcc_dev_version,
            dev,
            bmk_year_dev,
            dst=gcc_vs_gcc_tablesdir,
            bmk_type=bmk_type,
            species=['O3'],
            overwrite=True
        )

    # ==================================================================
    # GCC vs GCC Global mean OH, MCF Lifetime, CH4 Lifetime
    # ==================================================================
    if OH_metrics:
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
        oh.make_benchmark_oh_metrics(
            ref,
            gcc_ref_version,
            dev,
            gcc_dev_version,
            dst=gcc_vs_gcc_tablesdir,
            overwrite=True,
            spcdb_dir=spcdb_dir
        )

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Create GCHP vs GCC benchmark plots and tables
#
# (1) GCHP (Dev) and GCC (Ref) use the same benchmark year.
# (2) The GCC version in "GCHP vs GCC" is the Dev of "GCC vs GCC".
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if gchp_vs_gcc:

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
        gchp_metname(gchp_dev_prior_to_13),
        all_months_gchp_dev,
        is_gchp=True,
        gchp_format_is_legacy=gchp_dev_is_legacy
    )[0]

    # ==================================================================
    # GCHP vs GCC Concentration plots
    # ==================================================================
    if plot_conc:
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
            is_gchp=True,
            gchp_format_is_legacy=gchp_dev_is_legacy
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
            subdst='AnnualMean',
            time_mean=True,
            weightsdir=weightsdir,
            benchmark_type=bmk_type,
            overwrite=True,
            spcdb_dir=spcdb_dir
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
                subdst=datestr,
                weightsdir=weightsdir,
                benchmark_type=bmk_type,
                plot_by_spc_cat=plot_by_spc_cat,
                overwrite=True,
                spcdb_dir=spcdb_dir
            )

    # ==============================================================
    # GCHP vs. GCC emissions plots
    # ==============================================================
    if plot_emis:
        print("\n%%% Creating GCHP vs. GCC emissions plots %%%")

        # Diagnostic collections to read
        col = "Emissions"

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
            is_gchp=True,
            gchp_format_is_legacy=gchp_dev_is_legacy
        )[0]

        # Create plots
        print("\nCreating plots for annual mean")
        bmk.make_benchmark_emis_plots(
            ref,
            gchp_vs_gcc_refstr,
            dev,
            gchp_vs_gcc_devstr,
            dst=gchp_vs_gcc_resultsdir,
            subdst='AnnualMean',
            time_mean=True,
            weightsdir=weightsdir,
            plot_by_spc_cat=plot_by_spc_cat,
            plot_by_hco_cat=plot_by_hco_cat,
            overwrite=True,
            spcdb_dir=spcdb_dir
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
                weightsdir=weightsdir,
                plot_by_spc_cat=plot_by_spc_cat,
                plot_by_hco_cat=plot_by_hco_cat,
                overwrite=True,
                spcdb_dir=spcdb_dir
            )

    # ==================================================================
    # GCHP vs. GCC tables of emission and inventory totals
    # ==================================================================
    if emis_table:
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
            is_gchp=True,
            gchp_format_is_legacy=gchp_dev_is_legacy
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
            spcdb_dir=spcdb_dir
        )

    # ==================================================================
    # GCHP vs. GCC J-values plots
    # ==================================================================
    if plot_jvalues:
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
            is_gchp=True,
            gchp_format_is_legacy=gchp_dev_is_legacy
        )[0]

        # Create plots
        print("\nCreating plots for annual mean")
        bmk.make_benchmark_jvalue_plots(
            ref,
            gchp_vs_gcc_refstr,
            dev,
            gchp_vs_gcc_devstr,
            dst=gchp_vs_gcc_resultsdir,
            subdst='AnnualMean',
            time_mean=True,
            weightsdir=weightsdir,
            overwrite=True,
            spcdb_dir=spcdb_dir
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
                weightsdir=weightsdir,
                overwrite=True,
                spcdb_dir=spcdb_dir
            )

    # ==================================================================
    # GCHP vs GCC column AOD plots
    # ==================================================================
    if plot_aod:
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
            is_gchp=True,
            gchp_format_is_legacy=gchp_dev_is_legacy

        )[0]

        # Create plots
        print("\nCreating plots for annual mean")
        bmk.make_benchmark_aod_plots(
            ref,
            gchp_vs_gcc_refstr,
            dev,
            gchp_vs_gcc_devstr,
            dst=gchp_vs_gcc_resultsdir,
            subdst='AnnualMean',
            time_mean=True,
            weightsdir=weightsdir,
            overwrite=True,
            spcdb_dir=spcdb_dir
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
                weightsdir=weightsdir,
                overwrite=True,
                spcdb_dir=spcdb_dir
            )

    # ==================================================================
    # GCHP vs GCC global mass tables
    # ==================================================================
    if mass_table:
        print("\n%%% Creating GCHP vs. GCC mass tables %%%")

        def gchp_vs_gcc_mass_table(m):
            """
            Create mass table for each benchmark month in parallel
            """

            # Filepaths
            refpath = get_filepath(
                gchp_vs_gcc_refrstdir,
                "Restart",
                bmk_mons_dev[m]
            )
            devpath = get_filepath(
                gchp_vs_gcc_devrstdir,
                "Restart",
                bmk_mons_dev[m],
                is_gchp=True,
                gchp_format_is_legacy=gchp_dev_is_legacy
            )

            # use initial restart if no checkpoint present (intended for
            # first month).  need to pass path of meteorology file with
            # area variable in this scenario
            dev_extra = ''
            if not os.path.isfile(devpath):
                devpath = join(
                    gchp_vs_gcc_devrstdir,
                    'initial_GEOSChem_rst.' + gchp_dev_res + '_benchmark.nc'
                )
                dev_extra = get_filepath(
                    gchp_vs_gcc_devrstdir,
                    "Restart",
                    bmk_mons_dev[m+1],
                    is_gchp=True,
                    gchp_format_is_legacy=gchp_dev_is_legacy
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
                dev_met_extra=dev_extra
            )

        results = Parallel(n_jobs=-1)\
            (delayed(gchp_vs_gcc_mass_table)(t) for t in range(bmk_n_months))

    # ==================================================================
    # GCHP vs GCC operations budgets tables
    # ==================================================================
    if ops_budget_table:
        print("\n%%% Creating GCHP vs. GCC operations budget tables %%%")

        def gchp_vs_gcc_ops_budg(m):
            """
            Create operations budgets for each benchmark month m in parallel
            """

            # Filepaths
            refpath = get_filepath(
                gchp_vs_gcc_refdir,
                "Budget",
                bmk_mons_dev[m]
            )
            devpath = get_filepath(
                gchp_vs_gcc_devdir,
                "Budget",
                bmk_mons_gchp_dev[m],
                is_gchp=True,
                gchp_format_is_legacy=gchp_dev_is_legacy
            )

            # Create tables
            bmk.make_benchmark_operations_budget(
                gcc_dev_version,
                refpath,
                gchp_dev_version,
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
                    "WetDep"
                ],
                compute_accum=False,
                dst=gchp_vs_gcc_tablesdir
            )

        results = Parallel(n_jobs=-1)\
            (delayed(gchp_vs_gcc_ops_budg)(t) for t in range(bmk_n_months))

    # ==================================================================
    # GCHP vs GCC aerosol budgets and burdens tables
    # ==================================================================
    if aer_budget_table:
        print("\n%%% Creating GCHP vs. GCC aerosol budget tables %%%")

        # Filepaths
        devaero = get_filepaths(
            gchp_vs_gcc_devdir,
            "Aerosols",
            all_months_gchp_dev,
            is_gchp=True,
            gchp_format_is_legacy=gchp_dev_is_legacy
        )[0]
        devspc = get_filepaths(
            gchp_vs_gcc_devdir,
            "SpeciesConc",
            all_months_gchp_dev,
            is_gchp=True,
            gchp_format_is_legacy=gchp_dev_is_legacy
        )[0]

        # Create tables
        bmk.make_benchmark_aerosol_tables(
            gchp_vs_gcc_devdir,
            devaero,
            devspc,
            devmet,
            gchp_dev_version,
            bmk_year_dev,
            days_per_month_dev,
            dst=gchp_vs_gcc_tablesdir,
            overwrite=True,
            spcdb_dir=spcdb_dir,
            is_gchp=True,
        )

    # Comment out the budget tables until we are sure that GCHP
    # benchmarks archive wetdep fields for HNO3 (bmy, 4/1/21)
    ## ==================================================================
    ## GCHP vs GCC Ox budget tables
    ## ==================================================================
    #if Ox_budget_table:
    #    print("\n%%% Creating GCHP vs. GCC Ox budget tables %%%")
    #
    #    # Compute Ox budget table for GCC
    #    ox.global_ox_budget(
    #        gcc_dev_version,
    #        gcc_vs_gcc_devdir,
    #        gcc_vs_gcc_devrstdir,
    #        bmk_year_dev,
    #        dst=gcc_vs_gcc_tablesdir,
    #        overwrite=True,
    #        spcdb_dir=spcdb_dir
    #    )
    #
    #    # Compute Ox budget table for GCHP
    #    ox.global_ox_budget(
    #        gchp_dev_version,
    #        gchp_vs_gcc_devdir,
    #        gchp_vs_gcc_devrstdir,
    #        bmk_year_dev,
    #        dst=gchp_vs_gcc_tablesdir,
    #        overwrite=True,
    #        is_gchp=True,
    #        spcdb_dir=spcdb_dir
    #    )

    # ==================================================================
    # GCHP vs. GCC global mean OH, MCF Lifetime, CH4 Lifetime
    # ==================================================================
    if OH_metrics:
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
            is_gchp=True,
            gchp_format_is_legacy=gchp_dev_is_legacy
        )[0]

        # Create table
        oh.make_benchmark_oh_metrics(
            ref,
            gcc_dev_version,
            dev,
            gchp_dev_version,
            dst=gchp_vs_gcc_tablesdir,
            overwrite=True,
            spcdb_dir=spcdb_dir
        )

    # ==================================================================
    # GCHP Strat-Trop Exchange
    # -=================================================================
    if ste_table:
        print("\n%%% Skipping GCHP vs. GCC Strat-Trop Exchange table %%%")

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Create GCHP vs GCHP benchmark plots and tables
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if gchp_vs_gchp:

    # ==================================================================
    # GCHP vs GCC filepaths for StateMet collection data
    # ==================================================================
    refmet = get_filepaths(
        gchp_vs_gchp_refdir,
        gchp_metname(gchp_ref_prior_to_13),
        all_months_gchp_ref,
        is_gchp=True,
        gchp_format_is_legacy=gchp_ref_is_legacy
    )[0]
    devmet = get_filepaths(
        gchp_vs_gcc_devdir,
        gchp_metname(gchp_dev_prior_to_13),
        all_months_gchp_dev,
        is_gchp=True,
        gchp_format_is_legacy=gchp_dev_is_legacy
    )[0]

    # ==================================================================
    # GCHP vs GCHP species concentration plots
    # ==================================================================
    if plot_conc:
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
            gchp_format_is_legacy=gchp_ref_is_legacy
        )
        dev = get_filepaths(
            gchp_vs_gchp_devdir,
            col,
            all_months_gchp_dev,
            is_gchp=True,
            gchp_format_is_legacy=gchp_dev_is_legacy
        )

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
            subdst='AnnualMean',
            time_mean=True,
            weightsdir=weightsdir,
            benchmark_type=bmk_type,
            plot_by_spc_cat=plot_by_spc_cat,
            overwrite=True,
            spcdb_dir=spcdb_dir
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
                weightsdir=weightsdir,
                benchmark_type=bmk_type,
                plot_by_spc_cat=plot_by_spc_cat,
                overwrite=True,
                spcdb_dir=spcdb_dir
            )

    # ==================================================================
    # GCHP vs. GCHP Emissions plots
    # ==================================================================
    if plot_emis:
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
            gchp_format_is_legacy=gchp_ref_is_legacy
        )[0]
        dev = get_filepaths(
            gchp_vs_gchp_devdir,
            "Emissions",
            all_months_gchp_dev,
            is_gchp=True,
            gchp_format_is_legacy=gchp_dev_is_legacy
        )[0]

        # Create plots
        print("\nCreating plots for annual mean")
        bmk.make_benchmark_emis_plots(
            ref,
            gchp_vs_gchp_refstr,
            dev,
            gchp_vs_gchp_devstr,
            dst=gchp_vs_gchp_resultsdir,
            subdst='AnnualMean',
            time_mean=True,
            weightsdir=weightsdir,
            plot_by_spc_cat=plot_by_spc_cat,
            plot_by_hco_cat=plot_by_hco_cat,
            overwrite=True,
            spcdb_dir=spcdb_dir
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
                weightsdir=weightsdir,
                plot_by_spc_cat=plot_by_spc_cat,
                plot_by_hco_cat=plot_by_hco_cat,
                overwrite=True,
                spcdb_dir=spcdb_dir
            )

    # ==================================================================
    # GCHP vs. GCHP tables of emission and inventory totals
    # ==================================================================
    if emis_table:
        print("\n%%% Creating GCHP vs. GCHP emissions tables %%%")

        # Filepaths
        ref = get_filepaths(
            gchp_vs_gchp_refdir,
            "Emissions",
            all_months_gchp_ref,
            is_gchp=True,
            gchp_format_is_legacy=gchp_ref_is_legacy
        )
        dev = get_filepaths(
            gchp_vs_gchp_devdir,
            "Emissions",
            all_months_gchp_dev,
            is_gchp=True,
            gchp_format_is_legacy=gchp_dev_is_legacy
        )

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
            spcdb_dir=spcdb_dir
        )

    # ==================================================================
    # GCHP vs. GCHP J-values plots
    # ==================================================================
    if plot_jvalues:
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
            gchp_format_is_legacy=gchp_ref_is_legacy
        )
        dev = get_filepaths(
            gchp_vs_gchp_devdir,
            "JValues",
            all_months_gchp_dev,
            is_gchp=True,
            gchp_format_is_legacy=gchp_dev_is_legacy
        )

        # Create plots
        print("\nCreating plots for annual mean")
        bmk.make_benchmark_jvalue_plots(
            ref,
            gchp_vs_gchp_refstr,
            dev,
            gchp_vs_gchp_devstr,
            dst=gchp_vs_gchp_resultsdir,
            subdst=bmk_mon_yr_strs_dev[t],
            weightsdir=weightsdir,
            overwrite=True,
            spcdb_dir=spcdb_dir
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
                weightsdir=weightsdir,
                overwrite=True,
                spcdb_dir=spcdb_dir
            )

    # ==================================================================
    # GCHP vs GCHP column AOD plots
    # ==================================================================
    if plot_aod:
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
            gchp_format_is_legacy=gchp_ref_is_legacy
        )
        dev = get_filepaths(
            gchp_vs_gchp_devdir,
            "Aerosols",
            all_months_gchp_dev,
            is_gchp=True,
            gchp_format_is_legacy=gchp_dev_is_legacy
        )

        # Create plots
        print("\nCreating plots for annual mean")
        bmk.make_benchmark_aod_plots(
            ref,
            gchp_vs_gchp_refstr,
            dev,
            gchp_vs_gchp_devstr,
            dst=gchp_vs_gchp_resultsdir,
            subdst='AnnualMean',
            time_mean=True,
            weightsdir=weightsdir,
            overwrite=True,
            spcdb_dir=spcdb_dir
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
                weightsdir=weightsdir,
                overwrite=True,
                spcdb_dir=spcdb_dir
            )


    # ==================================================================
    # GCHP vs GCHP global mass tables
    # ==================================================================
    if mass_table:
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
                gchp_format_is_legacy=gchp_ref_is_legacy
            )

            # Use initial checkpoint if Ref restart is not present
            ref_extra = ''
            if not os.path.isfile(refpath):
                refpath = join(
                    gchp_vs_gchp_refrstdir,
                    'initial_GEOSChem_rst.' + gchp_ref_res + '_benchmark.nc'
                )
                ref_extra = get_filepath(
                    gchp_vs_gchp_refrstdir,
                    "Restart",
                    bmk_mons_ref[m+1],
                    is_gchp=True,
                    gchp_format_is_legacy=gchp_ref_is_legacy
                )

            # Dev filepaths
            devpath = get_filepath(
                gchp_vs_gchp_devrstdir,
                "Restart",
                bmk_mons_dev[m],
                is_gchp=True,
                gchp_format_is_legacy=gchp_dev_is_legacy
            )

            # Use initial checkpoint if Dev restart is not present
            dev_extra = ''
            if not os.path.isfile(devpath):
                devpath = join(
                    gchp_vs_gchp_devrstdir,
                    'initial_GEOSChem_rst.' + gchp_dev_res + '_benchmark.nc'
                )
                dev_extra = get_filepath(
                    gchp_vs_gchp_devrstdir,
                    "Restart",
                    bmk_mons_dev[m+1],
                    is_gchp=True,
                    gchp_format_is_legacy=gchp_dev_is_legacy
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
                dev_met_extra=dev_extra
            )

        # Run in parallel
        results = Parallel(n_jobs=-1)\
            (delayed(gchp_vs_gchp_mass_table)(t) for t in range(bmk_n_months))

    # ==================================================================
    # GCHP vs GCHP operations budgets tables
    # ==================================================================
    if ops_budget_table:
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
                gchp_format_is_legacy=gchp_ref_is_legacy
            )
            devpath = get_filepath(
                gchp_vs_gchp_devdir,
                "Budget",
                bmk_mons_gchp_dev[m],
                is_gchp=True,
                gchp_format_is_legacy=gchp_dev_is_legacy
            )

            # Compute tables
            bmk.make_benchmark_operations_budget(
                gchp_ref_version,
                refpath,
                gchp_dev_version,
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
                    "WetDep"
                ],
                compute_accum=False,
                dst=gchp_vs_gchp_tablesdir
            )

        # Run in parallel
        results = Parallel(n_jobs=-1)\
            (delayed(gchp_vs_gchp_ops_budg)(t) for t in range(bmk_n_months))

    # ==================================================================
    # GCHP vs GCHP aerosol budgets and burdens tables
    # ==================================================================
    if aer_budget_table:
        print("\n%%% Skipping GCHP vs. GCHP aerosol budget tables: Redundant%%%")

    # Comment out the budget tables until we are sure that GCHP
    # benchmarks archive wetdep fields for HNO3 (bmy, 4/1/21)
    ## ==================================================================
    ## GCHP vs GCHP Ox budget tables
    ## ==================================================================
    #if Ox_budget_table:
    #    print("\n%%% Creating GCHP Ox budget table %%%")
    #
    #    # Compute Ox budget table for GCHP
    #    ox.global_ox_budget(
    #        gchp_dev_version,
    #        gchp_vs_gchp_devdir,
    #        gchp_vs_gchp_devrstdir,
    #        bmk_year_dev,
    #        dst=gchp_vs_gchp_tablesdir,
    #        overwrite=True,
    #        is_gchp=True,
    #        spcdb_dir=spcdb_dir
    #    )

    # ==================================================================
    # GCHP vs. GCHP global mean OH, MCF Lifetime, CH4 Lifetime
    # ==================================================================
    if OH_metrics:
        print("\n%%% Creating GCHP vs. GCHP OH metrics table %%%")

        # Filepaths
        ref = get_filepaths(
            gchp_vs_gchp_refdir,
            "Metrics",
            all_months_gchp_ref,
            is_gchp=True,
            gchp_format_is_legacy=gchp_ref_is_legacy
        )[0]
        dev = get_filepaths(
            gchp_vs_gchp_devdir,
            "Metrics",
            all_months_gchp_dev,
            is_gchp=True,
            gchp_format_is_legacy=gchp_dev_is_legacy
        )[0]

        # Create the OH Metrics table
        oh.make_benchmark_oh_metrics(
            ref,
            gchp_ref_version,
            dev,
            gchp_dev_version,
            dst=gchp_vs_gchp_tablesdir,
            overwrite=True,
            spcdb_dir=spcdb_dir
        )

    # ==================================================================
    # GCHP Strat-Trop Exchange
    # ==================================================================
    if ste_table:
        print("\n%%% Skipping GCHP vs. GCHP Strat-Trop Exchange table %%%")
