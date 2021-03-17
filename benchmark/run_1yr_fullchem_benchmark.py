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

    This issue might be fixed in matplotlib 3.0.

This script corresponds with GCPy 1.0.2. Edit this version ID if releasing
a new version of GCPy.
"""

# =====================================================================
# Imports and global settings (you should not need to edit these)
# =====================================================================

import os
from os.path import join
from shutil import copyfile
import warnings

from calendar import monthrange
import numpy as np
import xarray as xr

from gcpy import benchmark as bmk
from gcpy.util import get_filepath, get_filepaths, get_area_from_dataset
import gcpy.ste_flux as ste
import gcpy.oh_metrics as oh
import gcpy.mean_oh_from_logs as moh  # NOTE: to be removed after 13.0.0
from joblib import Parallel, delayed, cpu_count, parallel_backend

# Tell matplotlib not to look for an X-window
os.environ["QT_QPA_PLATFORM"]="offscreen"

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
spcdb_dir   = join(maindir, gcc_dev_version)

# GCHP initial restart resolution (for mass tables)
gchp_ref_res = 'c48'
gchp_dev_res = 'c48'

# Kludge: Set switches that will pick the proper StateMet collection
# for GCHP.  Versions prior to 13.0.0 used StateMet_avg.
gchp_ref_prior_to_13 = False
gchp_dev_prior_to_13 = False

# =====================================================================
# Specify if this is a gcpy test validation run
# =====================================================================
gcpy_test = True

# =====================================================================
# Comparisons to run
# =====================================================================
gcc_vs_gcc   = True
gchp_vs_gcc  = True
gchp_vs_gchp = True
# GCHP vs GCC diff of diffs not included in 1-yr full chemistry benchmark

# =====================================================================
# Output to generate (plots/tables will be created in this order):
# =====================================================================
plot_conc        = True
plot_emis        = True
emis_table       = True
plot_jvalues     = True
plot_aod         = True
mass_table       = True
ops_budget_table = True
aer_budget_table = True
ste_table        = True # GCC only
OH_metrics       = True

# Plot concentrations and emissions by category?
plot_by_spc_cat = True
plot_by_hco_cat = True

# =====================================================================
# Data directories
# For gchp_vs_gcc_refdir use gcc_dev_version, not ref (mps, 6/27/19)
# =====================================================================

# Diagnostics file directory paths
gcc_vs_gcc_refdir      = join(maindir, gcc_ref_version,  "OutputDir")
gcc_vs_gcc_devdir      = join(maindir, gcc_dev_version,  "OutputDir")
gchp_vs_gcc_refdir     = join(maindir, gcc_dev_version,  "OutputDir")
gchp_vs_gcc_devdir     = join(maindir, gchp_dev_version, "OutputDir")
gchp_vs_gchp_refdir    = join(maindir, gchp_ref_version, "OutputDir")
gchp_vs_gchp_devdir    = join(maindir, gchp_dev_version, "OutputDir")

# Restart file directory paths
gcc_vs_gcc_refrstdir   = join(maindir, gcc_ref_version,  "restarts")
gcc_vs_gcc_devrstdir   = join(maindir, gcc_dev_version,  "restarts")
gchp_vs_gcc_refrstdir  = join(maindir, gcc_dev_version,  "restarts")
gchp_vs_gcc_devrstdir  = join(maindir, gchp_dev_version)
gchp_vs_gchp_refrstdir = join(maindir, gchp_ref_version)
gchp_vs_gchp_devrstdir = join(maindir, gchp_dev_version)

# Log file directories -- GEOS-Chem "Classic" only
gcc_vs_gcc_reflogdir   = join(maindir, gcc_ref_version,  "logs")
gcc_vs_gcc_devlogdir   = join(maindir, gcc_dev_version,  "logs")

# =====================================================================
# Benchmark output directories
# =====================================================================
# Plot directories
if gcpy_test:
    mainresultsdir           = join('.', results_dir)
    gcc_vs_gcc_resultsdir    = join(mainresultsdir,'GCC_version_comparison')
    gchp_vs_gchp_resultsdir  = join(mainresultsdir,'GCHP_version_comparison')
    gchp_vs_gcc_resultsdir   = join(mainresultsdir,'GCHP_GCC_comparison')
    if not os.path.exists(mainresultsdir): os.mkdir(mainresultsdir)
    # Make copy of benchmark script in results directory
    curfile = os.path.realpath(__file__)
    copyfile(curfile, join(mainresultsdir,curfile.split('/')[-1]))

else:
    gcc_vs_gcc_resultsdir    = join(maindir, gcc_dev_version, results_dir)
    gchp_vs_gchp_resultsdir  = join(maindir, gchp_dev_version,
                                    results_dir, "GCHP_version_comparison")
    gchp_vs_gcc_resultsdir   = join(maindir, gchp_dev_version,
                                    results_dir, "GCHP_GCC_comparison")
    base_gchp_resultsdir     = join(maindir, gchp_dev_version, results_dir)
    #make results directories that don't exist
    for resdir, plotting_type in zip([gcc_vs_gcc_resultsdir, base_gchp_resultsdir,
                                      gchp_vs_gchp_resultsdir, gchp_vs_gcc_resultsdir],
                                     [gcc_vs_gcc, gchp_vs_gcc or gchp_vs_gchp,
                                      gchp_vs_gchp, gchp_vs_gcc]):
        if plotting_type and not os.path.exists(resdir): os.mkdir(resdir)
        if resdir == gcc_vs_gcc_resultsdir or resdir == base_gchp_resultsdir:
            # Make copy of benchmark script in results directory
            curfile = os.path.realpath(__file__)
            copyfile(curfile, join(resdir,curfile.split('/')[-1]))

# Tables directories
gcc_vs_gcc_tablesdir   = join(gcc_vs_gcc_resultsdir,"Tables")
gchp_vs_gcc_tablesdir  = join(gchp_vs_gcc_resultsdir,"Tables")
gchp_vs_gchp_tablesdir = join(gchp_vs_gchp_resultsdir,"Tables")

# Budget directories
gcc_vs_gcc_budgetdir   = join(gcc_vs_gcc_resultsdir,"Budget")
gchp_vs_gcc_budgetdir  = join(gchp_vs_gcc_resultsdir,"Budget")
gchp_vs_gchp_budgetdir = join(gchp_vs_gchp_resultsdir,"Budget")

# =====================================================================
# Plot title strings
# For gchp_vs_gcc_refstr use gcc_dev_version, not ref (mps, 6/27/19)
# =====================================================================
gcc_vs_gcc_refstr    = gcc_ref_version
gcc_vs_gcc_devstr    = gcc_dev_version
gchp_vs_gcc_refstr   = gcc_dev_version
gchp_vs_gcc_devstr   = gchp_dev_version
gchp_vs_gchp_refstr  = gchp_ref_version
gchp_vs_gchp_devstr  = gchp_dev_version

########################################################################
###    THE REST OF THESE SETTINGS SHOULD NOT NEED TO BE CHANGED      ###
########################################################################

def gchp_ref_metname():
    if gchp_ref_prior_to_13:
        return "StateMet_avg"
    else:
        return "StateMet"

def gchp_dev_metname():
    if gchp_dev_prior_to_13:
        return "StateMet_avg"
    else:
        return "StateMet"

# =====================================================================
# Dates and times -- ref data
# =====================================================================

# Month/year strings for use in table subdirectories (e.g. Jan2016)
bmk_mon_yr_strs_ref = [v + bmk_year_ref for v in bmk_mon_strs]

# Get all months array of start datetimes for benchmark year
bmk_start_ref = np.datetime64(bmk_year_ref + "-01-01")
bmk_end_ref = np.datetime64("{}-01-01".format(int(bmk_year_ref)+1))
all_months_ref = np.arange(bmk_start_ref,
                           bmk_end_ref,
                           step=np.timedelta64(1, "M"),
                           dtype="datetime64[M]")

# Get all months array of mid-point datetime per month for benchmark year,
# and # sec and # days per month
# NOTE: GCHP time-averaged files have time in the middle of the month
sec_per_month_ref = np.zeros(12)
days_per_month_ref = np.zeros(12)
all_months_mid_ref = np.zeros(12, dtype="datetime64[h]")
for t in range(12):
    days_per_month_ref[t] = monthrange(int(bmk_year_ref), t + 1)[1]
    sec_per_month_ref[t] = days_per_month_ref[t] * 86400.0
    middle_hr = int(days_per_month_ref[t] * 24 / 2)
    delta = np.timedelta64(middle_hr, 'h')
    all_months_mid_ref[t] = all_months_ref[t].astype("datetime64[h]") + delta

# Get subset of month datetimes for only benchmark months
bmk_mons_ref = all_months_ref[bmk_mon_inds]
bmk_mons_mid_ref = all_months_mid_ref[bmk_mon_inds]
bmk_sec_per_month_ref = sec_per_month_ref[bmk_mon_inds]

# =====================================================================
# Dates and times -- Dev data
# =====================================================================

# Month/year strings for use in table subdirectories (e.g. Jan2016)
bmk_mon_yr_strs_dev = [v + bmk_year_dev for v in bmk_mon_strs]

# Get all months array of start datetimes for benchmark year
bmk_start_dev = np.datetime64(bmk_year_dev + "-01-01")
bmk_end_dev = np.datetime64("{}-01-01".format(int(bmk_year_dev)+1))
all_months_dev = np.arange(bmk_start_dev,
                           bmk_end_dev,
                           step=np.timedelta64(1, "M"),
                           dtype="datetime64[M]")

# Get all months array of mid-point datetime per month for benchmark year,
# and # sec and # days per month
# NOTE: GCHP time-averaged files have time in the middle of the month
sec_per_month_dev = np.zeros(12)
days_per_month_dev = np.zeros(12)
all_months_mid_dev = np.zeros(12, dtype="datetime64[h]")
for t in range(12):
    days_per_month_dev[t] = monthrange(int(bmk_year_dev), t + 1)[1]
    sec_per_month_dev[t] = days_per_month_dev[t] * 86400.0
    middle_hr = int(days_per_month_dev[t] * 24 / 2)
    delta = np.timedelta64(middle_hr, 'h')
    all_months_mid_dev[t] = all_months_dev[t].astype("datetime64[h]") + delta

# Get subset of month datetimes for only benchmark months
bmk_mons_dev = all_months_dev[bmk_mon_inds]
bmk_mons_mid_dev = all_months_mid_dev[bmk_mon_inds]
bmk_sec_per_month_dev = sec_per_month_dev[bmk_mon_inds]

# ======================================================================
# Print the list of plots & tables to the screen
# ======================================================================

print("The following plots and tables will be created for {}:".format(bmk_type))
if plot_conc:        print(" - Concentration plots")
if plot_emis:        print(" - Emissions plots")
if plot_jvalues:     print(" - J-values (photolysis rates) plots")
if plot_aod:         print(" - Aerosol optical depth plots")
if ops_budget_table: print(" - Operations budget tables")
if aer_budget_table: print(" - Aerosol budget/burden tables")
if emis_table:       print(" - Table of emissions totals by species and inventory")
if mass_table:       print(" - Table of species mass")
if OH_metrics:       print(" - Table of OH metrics")
if ste_table:        print(" - Table of strat-trop exchange")
print("Comparisons will be made for the following combinations:")
if gcc_vs_gcc:   print(" - GCC vs GCC")
if gchp_vs_gcc:  print(" - GCHP vs GCC")
if gchp_vs_gchp: print(" - GCHP vs GCHP")

# ======================================================================
# Create GCC vs GCC benchmark plots and tables
# ======================================================================
if gcc_vs_gcc:

    # --------------------------------------------------------------
    # GCC vs GCC Concentration plots
    #
    # Includes lumped species and separates by category if plot_by_spc_cat
    # is true; otherwise excludes lumped species and writes to one file.
    # --------------------------------------------------------------
    if plot_conc:
        print("\n%%% Creating GCC vs. GCC concentration plots %%%")

        # Diagnostic collections to read
        col = "SpeciesConc"
        colmet = "StateMet"

        # Create concentration plots for each benchmark month
        for t in range(bmk_n_months):

            # Time & date quantities
            reftime = bmk_mons_ref[t]
            devtime = bmk_mons_dev[t]
            datestr = bmk_mon_yr_strs_dev[t]

            # Read seasonal data
            ref = get_filepath(gcc_vs_gcc_refdir, col, reftime)
            dev = get_filepath(gcc_vs_gcc_devdir, col, devtime)
            refmet = get_filepath(gcc_vs_gcc_refdir, colmet, reftime)
            devmet = get_filepath(gcc_vs_gcc_devdir, colmet, devtime)

            # Create plots
            bmk.make_benchmark_conc_plots(
                ref,
                gcc_vs_gcc_refstr,
                dev,
                gcc_vs_gcc_devstr,
                refmet=refmet,
                devmet=devmet,
                dst=gcc_vs_gcc_resultsdir,
                subdst=datestr,
                weightsdir=weightsdir,
                benchmark_type=bmk_type,
                plot_by_spc_cat=plot_by_spc_cat,
                overwrite=True,
                spcdb_dir=spcdb_dir
            )

    # --------------------------------------------------------------
    # GCC vs GCC emissions plots
    # --------------------------------------------------------------
    if plot_emis:
        print("\n%%% Creating GCC vs. GCC emissions plots %%%")

        # Diagnostic collections to read
        col = "Emissions"

        # Create concentration plots for each benchmark month
        for t in range(bmk_n_months):

            # Time & date quantities
            reftime = bmk_mons_ref[t]
            devtime = bmk_mons_dev[t]
            datestr = bmk_mon_yr_strs_dev[t]

            # Read seasonal data
            ref = get_filepath(gcc_vs_gcc_refdir, col, reftime)
            dev = get_filepath(gcc_vs_gcc_devdir, col, devtime)

            # Create plots
            bmk.make_benchmark_emis_plots(
                ref,
                gcc_vs_gcc_refstr,
                dev,
                gcc_vs_gcc_devstr,
                dst=gcc_vs_gcc_resultsdir,
                subdst=datestr,
                weightsdir=weightsdir,
                plot_by_spc_cat=plot_by_spc_cat,
                plot_by_hco_cat=plot_by_hco_cat,
                overwrite=True,
                spcdb_dir=spcdb_dir
            )

    # --------------------------------------------------------------
    # GCC vs GCC tables of emission and inventory totals
    # --------------------------------------------------------------
    if emis_table:
        print("\n%%% Creating GCC vs. GCC emissions & inventory totals %%%")

        # Read all months of data
        col = "Emissions"
        ref = get_filepaths(gcc_vs_gcc_refdir, col, all_months_ref)
        dev = get_filepaths(gcc_vs_gcc_devdir, col, all_months_dev)

        # Create emissions table that spans entire year
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

    # --------------------------------------------------------------
    # GCC vs GCC J-value plots
    # --------------------------------------------------------------
    if plot_jvalues:
        print("\n%%% Creating GCC vs. GCC J-value plots %%%")

        # Diagnostic collections to read
        col = "JValues"

        # Create J-value plots for each benchmark month
        for t in range(bmk_n_months):

            # Time & date quantities
            reftime = bmk_mons_ref[t]
            devtime = bmk_mons_dev[t]
            datestr = bmk_mon_yr_strs_dev[t]

            # Read seasonal data
            ref = get_filepath(gcc_vs_gcc_refdir, col, reftime)
            dev = get_filepath(gcc_vs_gcc_devdir, col, devtime)

            # Make plots
            bmk.make_benchmark_jvalue_plots(
                ref,
                gcc_vs_gcc_refstr,
                dev,
                gcc_vs_gcc_devstr,
                dst=gcc_vs_gcc_resultsdir,
                subdst=datestr,
                weightsdir=weightsdir,
                overwrite=True,
                spcdb_dir=spcdb_dir
            )

    # --------------------------------------------------------------
    # GCC vs. GCC column AOD plots
    # --------------------------------------------------------------
    if plot_aod:
        print("\n%%% Creating GCC vs. GCC column AOD plots %%%")

        # Diagnostic collections to read
        col = "Aerosols"

        # Create AOD plots for each benchmark month
        for t in range(bmk_n_months):

            # Time & date quantities
            reftime = bmk_mons_ref[t]
            devtime = bmk_mons_dev[t]
            datestr = bmk_mon_yr_strs_dev[t]

            # Read seasonal data
            ref = get_filepath(gcc_vs_gcc_refdir, col, reftime)
            dev = get_filepath(gcc_vs_gcc_devdir, col, devtime)

            # Create plots
            bmk.make_benchmark_aod_plots(
                ref,
                gcc_vs_gcc_refstr,
                dev,
                gcc_vs_gcc_devstr,
                dst=gcc_vs_gcc_resultsdir,
                subdst=datestr,
                weightsdir=weightsdir,
                overwrite=True,
                spcdb_dir=spcdb_dir
            )

    # --------------------------------------------------------------
    # GCC vs GCC mass tables
    # --------------------------------------------------------------
    if mass_table:
        print("\n%%% Creating GCC vs. GCC mass tables %%%")

        # Diagnostic collections to read
        col = "Restart"
        # Create mass table for each benchmark month
        def parallel_mass_table(t):

            # Time & date quantities
            reftime = bmk_mons_ref[t]
            devtime = bmk_mons_dev[t]
            datestr = bmk_mon_yr_strs_dev[t]
            label ="at 01{}".format(datestr)

            # Read seasonal data
            ref = get_filepath(gcc_vs_gcc_refrstdir, col, reftime)
            dev = get_filepath(gcc_vs_gcc_devrstdir, col, devtime)

            # Create tables
            bmk.make_benchmark_mass_tables(
                ref,
                gcc_vs_gcc_refstr,
                dev,
                gcc_vs_gcc_devstr,
                dst=gcc_vs_gcc_tablesdir,
                subdst=datestr,
                label=label,
                overwrite=True,
                spcdb_dir=spcdb_dir
            )

        # Run in parallel
        results = Parallel(n_jobs=-1)\
            (delayed(parallel_mass_table) (t) for t in range(bmk_n_months))

    # --------------------------------------------------------------
    # GCC vs GCC operations budgets tables
    # --------------------------------------------------------------
    if ops_budget_table:
        print("\n%%% Creating GCC vs. GCC operations budget tables %%%")

        # Diagnostic collections to read
        col = "Budget"

        # Create budget table for each benchmark month (ewl??)
        def parallel_ops_budg(t):

            # Time & date quantities
            reftime = bmk_mons_ref[t]
            devtime = bmk_mons_dev[t]
            datestr = bmk_mon_yr_strs_dev[t]
            label ="at 01{}".format(datestr)
            sec_ref = sec_per_month_ref[t]
            sec_dev = sec_per_month_dev[t]

            # Read seasonal data
            ref = get_filepath(gcc_vs_gcc_refdir, col, reftime)
            dev = get_filepath(gcc_vs_gcc_devdir, col, devtime)

            # Create tables
            bmk.make_benchmark_operations_budget(
                gcc_ref_version,
                ref,
                gcc_dev_version,
                dev,
                sec_ref,
                sec_dev,
                benchmark_type=bmk_type,
                label=label,
                dst=gcc_vs_gcc_tablesdir
            )

        # Run in parallel
        results = Parallel(n_jobs=-1)\
            (delayed(parallel_ops_budg) (t) for t in range(bmk_n_months))

    # --------------------------------------------------------------
    # GCC vs GCC aerosols budgets/burdens tables
    # --------------------------------------------------------------
    if aer_budget_table:
        print("\n%%% Creating GCC vs. GCC aerosols budget tables %%%")

        # Diagnostic collections to read
        col_aero = "Aerosols"
        col_spc = "SpeciesConc"
        col_met = "StateMet"
        dev_aero = get_filepaths(gcc_vs_gcc_devdir, col_aero, all_months_dev)
        dev_spc = get_filepaths(gcc_vs_gcc_devdir, col_spc, all_months_dev)
        dev_met = get_filepaths(gcc_vs_gcc_devdir, col_met, all_months_dev)

        # Compute global aerosol budgets and burdens
        bmk.make_benchmark_aerosol_tables(
            gcc_vs_gcc_devdir,
            dev_aero,
            dev_spc,
            dev_met,
            gcc_dev_version,
            bmk_year_dev,
            days_per_month_dev,
            dst=gcc_vs_gcc_tablesdir,
            overwrite=True,
            spcdb_dir=spcdb_dir
        )

    # --------------------------------------------------------------
    # GCC Strat-Trop Exchange
    # --------------------------------------------------------------
    if ste_table:
        print("\n%%% Creating GCC vs. GCC Strat-Trop Exchange table %%%")

        # Diagnostic collections to read (all 12 months)
        col = "AdvFluxVert"
        dev = get_filepaths(gcc_vs_gcc_devdir, col, all_months_dev)[0]

        # Compute monthly and annual average strat-trop exchange of O3
        ste.make_benchmark_ste_table(
            gcc_dev_version,
            dev,
            bmk_year_dev,
            dst=gcc_vs_gcc_tablesdir,
            bmk_type=bmk_type,
            species=['O3'],
            overwrite=True
        )

    # --------------------------------------------------------------
    # GCC vs GCC Global mean OH, MCF Lifetime, CH4 Lifetime
    # --------------------------------------------------------------
    if OH_metrics:
        print("\n%%% Creating GCC vs. GCC OH metrics %%%")

# NOTE: Use this if the benchmark is prior to 13.0.0!
#        # Compute mean OH from the log files
#        # NOTE: Only works for GEOS-Chem "Classic" benchmarks!
#        moh.make_benchmark_oh_from_logs(
#            gcc_vs_gcc_reflogdir,
#            gcc_vs_gcc_refstr,
#           gcc_vs_gcc_devlogdir,
#            gcc_vs_gcc_devstr,
#            bmk_year_ref,
#            dst=gcc_vs_gcc_tablesdir,
#            overwrite=True
#        )

        # Diagnostic collection files to read
        col = "Metrics"
        ref = get_filepaths(gcc_vs_gcc_refdir, col, all_months_ref)[0]
        dev = get_filepaths(gcc_vs_gcc_devdir, col, all_months_dev)[0]

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

# ======================================================================
# Create GCHP vs GCC benchmark plots and tables
# ======================================================================
if gchp_vs_gcc:

    #---------------------------------------------------------------
    # GCHP vs GCC Concentration plots
    #---------------------------------------------------------------
    if plot_conc:
        print("\n%%% Creating GCHP vs. GCC concentration plots %%%")

        # Diagnostic collections to read
        col = "SpeciesConc"
        col_met_gcc = "StateMet"
        col_met_gchp = gchp_dev_metname()

        # Create concentration plots for each benchmark month
        for t in range(bmk_n_months):

            # Time & date quantities
            reftime = bmk_mons_dev[t]
            devtime = bmk_mons_mid_dev[t]
            datestr = bmk_mon_yr_strs_dev[t]

            # SpeciesConc
            ref = get_filepath(gchp_vs_gcc_refdir, col, reftime)
            dev = get_filepath(gchp_vs_gcc_devdir, col, devtime, is_gchp=True)

            # StateMet
            refmet = get_filepath(gchp_vs_gcc_refdir, col_met_gcc, reftime)
            devmet = get_filepath(gchp_vs_gcc_devdir, col_met_gchp, devtime,
                                  is_gchp=True)

            bmk.make_benchmark_conc_plots(
                ref,
                gchp_vs_gcc_refstr,
                dev,
                gchp_vs_gcc_devstr,
                refmet=refmet,
                devmet=devmet,
                dst=gchp_vs_gcc_resultsdir,
                subdst=datestr,
                weightsdir=weightsdir,
                benchmark_type=bmk_type,
                plot_by_spc_cat=plot_by_spc_cat,
                overwrite=True,
                spcdb_dir=spcdb_dir
            )

    #---------------------------------------------------------------
    # GCHP vs. GCC Emissions plots
    #---------------------------------------------------------------
    if plot_emis:
        print("\n%%% Creating GCHP vs. GCC emissions plots %%%")

        # Diagnostic collections to read
        col = "Emissions"

        # Create concentration plots for each benchmark month
        for t in range(bmk_n_months):

            # Time & date quantities
            reftime = bmk_mons_dev[t]
            devtime = bmk_mons_mid_dev[t]
            datestr = bmk_mon_yr_strs_dev[t]

            # Read data
            ref = get_filepath(gchp_vs_gcc_refdir, col, reftime)
            dev = get_filepath(gchp_vs_gcc_devdir, col, devtime, is_gchp=True)

            # Create plots
            bmk.make_benchmark_emis_plots(
                ref,
                gchp_vs_gcc_refstr,
                dev,
                gchp_vs_gcc_devstr,
                dst=gchp_vs_gcc_resultsdir,
                subdst=datestr,
                weightsdir=weightsdir,
                plot_by_spc_cat=plot_by_spc_cat,
                plot_by_hco_cat=plot_by_hco_cat,
                overwrite=True,
                spcdb_dir=spcdb_dir
            )

    #---------------------------------------------------------------
    # GCHP vs. GCC tables of emission and inventory totals
    #---------------------------------------------------------------
    if emis_table:
        print("\n%%% Creating GCHP vs. GCC emissions tables %%%")

        # Diagnostic collections to read
        col = "Emissions"
        col_met_gchp = gchp_dev_metname()

        # Read data
        ref = get_filepaths(gchp_vs_gcc_refdir, col, all_months_dev)
        dev = get_filepaths(gchp_vs_gcc_devdir, col, all_months_mid_dev,
                            is_gchp=True)
        devmet = get_filepaths(gchp_vs_gcc_devdir, col_met_gchp,
                               all_months_mid_dev, is_gchp=True)

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

    #---------------------------------------------------------------
    # GCHP vs. GCC J-values plots
    #---------------------------------------------------------------
    if plot_jvalues:
        print("\n%%% Creating GCHP vs. GCC J-values plots %%%")

        # Diagnostic collections to read
        col = "JValues"

        # Create J-value plots for each benchmark month
        for t in range(bmk_n_months):

            # Time & date quantities
            reftime = bmk_mons_dev[t]
            devtime = bmk_mons_mid_dev[t]
            datestr = bmk_mon_yr_strs_dev[t]

            # Read data
            ref = get_filepath(gchp_vs_gcc_refdir, col, reftime)
            dev = get_filepath(gchp_vs_gcc_devdir, col, devtime, is_gchp=True)

            # Create plots
            bmk.make_benchmark_jvalue_plots(
                ref,
                gchp_vs_gcc_refstr,
                dev,
                gchp_vs_gcc_devstr,
                dst=gchp_vs_gcc_resultsdir,
                subdst=datestr,
                weightsdir=weightsdir,
                overwrite=True,
                spcdb_dir=spcdb_dir
            )

    #---------------------------------------------------------------
    # GCHP vs GCC column AOD plots
    #---------------------------------------------------------------
    if plot_aod:
        print("\n%%% Creating GCHP vs. GCC AOD plots %%%")

        # Diagnostic collections to read
        col = "Aerosols"

        # Create AOD plots for each benchmark month
        for t in range(bmk_n_months):

            # Time & date quantities
            reftime = bmk_mons_dev[t]
            devtime = bmk_mons_mid_dev[t]
            datestr = bmk_mon_yr_strs_dev[t]

            # Read data
            ref = get_filepath(gchp_vs_gcc_refdir, col, reftime)
            dev = get_filepath(gchp_vs_gcc_devdir, col, devtime, is_gchp=True)

            # Create plots
            bmk.make_benchmark_aod_plots(
                ref,
                gchp_vs_gcc_refstr,
                dev,
                gchp_vs_gcc_devstr,
                dst=gchp_vs_gcc_resultsdir,
                subdst=datestr,
                weightsdir=weightsdir,
                overwrite=True,
                spcdb_dir=spcdb_dir
            )

    #---------------------------------------------------------------
    # GCHP vs GCC global mass tables
    #---------------------------------------------------------------
    if mass_table:
        print("\n%%% Creating GCHP vs. GCC mass tables %%%")

        # Diagnostic collections to read
        col = "Restart"

        # Create mass table for each benchmark month in parallel
        def parallel_mass_table(t):

            # Time & date quantities
            reftime = bmk_mons_ref[t]
            devtime = bmk_mons_dev[t]
            datestr = bmk_mon_yr_strs_dev[t]
            label = "at 01{}".format(datestr)

            # Read data
            ref = get_filepath(gchp_vs_gcc_refrstdir, col, reftime)
            dev = get_filepath(gchp_vs_gcc_devrstdir, col, devtime,
                               is_gchp=True)

            # use initial restart if no checkpoint present (intended for
            # first month).  need to pass path of meteorology file with
            # area variable in this scenario
            dev_extra=''
            if not os.path.isfile(dev):
                dev = join(gchp_vs_gcc_devrstdir,
                           'initial_GEOSChem_rst.' + gchp_dev_res +
                           '_benchmark.nc')
                extratime = bmk_mons_dev[t+1]
                dev_extra = get_filepath(gchp_vs_gcc_devrstdir, col,
                                         extratime, is_gchp=True)
            # Create tables
            bmk.make_benchmark_mass_tables(
                ref,
                gchp_vs_gcc_refstr,
                dev,
                gchp_vs_gcc_devstr,
                dst=gchp_vs_gcc_tablesdir,
                subdst=datestr,
                label=label,
                overwrite=True,
                spcdb_dir=spcdb_dir,
                dev_met_extra=dev_extra
            )

        results = Parallel(n_jobs=-1)\
            (delayed(parallel_mass_table) (t) for t in range(bmk_n_months))

    #---------------------------------------------------------------
    # GCHP vs GCC operations budgets tables
    #---------------------------------------------------------------
    if ops_budget_table:
        print("\n%%% Creating GCHP vs. GCC operations budget tables %%%")

        # Diagnostic collections to read
        col = "Budget"

        # Create operations table for each benchmark month in parallel
        def parallel_ops_budg(t):

            # Time & date quantities
            reftime = bmk_mons_ref[t]
            devtime = bmk_mons_mid_dev[t]
            datestr = bmk_mon_yr_strs_dev[t]
            label = "at 01{}".format(datestr)
            sec_ref = bmk_sec_per_month_ref[t]
            sec_dev = bmk_sec_per_month_dev[t]

            # Read data
            ref = get_filepath(gchp_vs_gcc_refdir, col, reftime)
            dev = get_filepath(gchp_vs_gcc_devdir, col, devtime, is_gchp=True)

            # Create tables
            bmk.make_benchmark_operations_budget(
                gcc_dev_version,
                ref,
                gchp_dev_version,
                dev,
                sec_ref,
                sec_dev,
                benchmark_type=bmk_type,
                label=datestr,
                operations=["Chemistry", "Convection", "EmisDryDep",
                            "Mixing", "WetDep"],
                compute_accum=False,
                dst=gchp_vs_gcc_tablesdir
            )

        results = Parallel(n_jobs=-1)\
            (delayed(parallel_ops_budg) (t) for t in range(bmk_n_months))

    #---------------------------------------------------------------
    # GCHP vs GCC aerosol budgets and burdens tables
    #---------------------------------------------------------------
    if aer_budget_table:
        print("\n%%% Creating GCHP vs. GCC aerosol budget tables %%%")

        # Compute annual mean AOD budgets and aerosol burdens
        # Diagnostic collections to read
        col_aero = "Aerosols"
        col_spc = "SpeciesConc"
        col_met = gchp_dev_metname()

        # Read data from all months
        dev_aero = get_filepaths(gchp_vs_gcc_devdir, col_aero,
                                 all_months_mid_dev, is_gchp=True)
        dev_spc = get_filepaths(gchp_vs_gcc_devdir, col_spc,
                                all_months_mid_dev, is_gchp=True)
        dev_met = get_filepaths(gchp_vs_gcc_devdir, col_met,
                                all_months_mid_dev, is_gchp=True)

        # Compute global aerosol budgets and burdens
        bmk.make_benchmark_aerosol_tables(
            gchp_vs_gcc_devdir,
            dev_aero,
            dev_spc,
            dev_met,
            gchp_dev_version,
            bmk_year_dev,
            days_per_month_dev,
            dst=gchp_vs_gcc_tablesdir,
            overwrite=True,
            spcdb_dir=spcdb_dir,
            is_gchp=True
        )

    #---------------------------------------------------------------
    # GCHP vs. GCC global mean OH, MCF Lifetime, CH4 Lifetime
    #---------------------------------------------------------------
    if OH_metrics:
        print("\n%%% Creating GCHP vs. GCC OH metrics table %%%")

        # Diagnostic collection files to read
        col = "Metrics"
        ref = get_filepaths(gchp_vs_gcc_refdir, col, all_months_dev)[0]
        dev = get_filepaths(gchp_vs_gcc_devdir, col, all_months_mid_dev,
                            is_gchp=True)[0]

        # Create the OH Metrics table
        oh.make_benchmark_oh_metrics(
            ref,
            gcc_dev_version,
            dev,
            gchp_dev_version,
            dst=gchp_vs_gcc_tablesdir,
            overwrite=True,
            spcdb_dir=spcdb_dir
        )

    # --------------------------------------------------------------
    # GCHP Strat-Trop Exchange
    # --------------------------------------------------------------
    if ste_table:
        print("\n%%% Skipping GCHP vs. GCC Strat-Trop Exchange table %%%")

# =====================================================================
# Create GCHP vs GCHP benchmark plots and tables
# =====================================================================
if gchp_vs_gchp:

    if plot_conc:
        print("\n%%% Creating GCHP vs. GCHP concentration plots %%%")

        # Diagnostic collections to read
        col = "SpeciesConc"
        col_met_ref = gchp_ref_metname()
        col_met_dev = gchp_dev_metname()

        # Create concentration plots for each benchmark month
        for t in range(bmk_n_months):

            # SpeciesConc
            ref = get_filepath(gchp_vs_gchp_refdir, col, bmk_mons_mid_ref[t],
                               is_gchp=True)
            dev = get_filepath(gchp_vs_gchp_devdir, col, bmk_mons_mid_dev[t],
                               is_gchp=True)

            # StateMet
            refmet = get_filepath(gchp_vs_gchp_refdir, col_met_ref,
                                  bmk_mons_mid_ref[t], is_gchp=True)
            devmet = get_filepath(gchp_vs_gchp_devdir, col_met_dev,
                                  bmk_mons_mid_dev[t], is_gchp=True)

            bmk.make_benchmark_conc_plots(
                ref,
                gchp_vs_gchp_refstr,
                dev,
                gchp_vs_gchp_devstr,
                refmet=refmet,
                devmet=devmet,
                dst=gchp_vs_gchp_resultsdir,
                subdst=bmk_mon_yr_strs_dev[t],
                weightsdir=weightsdir,
                benchmark_type=bmk_type,
                plot_by_spc_cat=plot_by_spc_cat,
                overwrite=True,
                spcdb_dir=spcdb_dir
            )

    #---------------------------------------------------------------
    # GCHP vs. GCHP Emissions plots
    #---------------------------------------------------------------
    if plot_emis:
        print("\n%%% Creating GCHP vs. GCHP emissions plots %%%")

        # Diagnostic collections to read
        col = "Emissions"

        # Create concentration plots for each benchmark month
        for t in range(bmk_n_months):

            ref = get_filepath(gchp_vs_gchp_refdir, col,
                               bmk_mons_mid_ref[t],  is_gchp=True)
            dev = get_filepath(gchp_vs_gchp_devdir, col,
                               bmk_mons_mid_dev[t],  is_gchp=True)

            # Create the plots
            bmk.make_benchmark_emis_plots(
                ref,
                gchp_vs_gchp_refstr,
                dev,
                gchp_vs_gchp_devstr,
                dst=gchp_vs_gchp_resultsdir,
                subdst=bmk_mon_yr_strs_dev[t],
                weightsdir=weightsdir,
                plot_by_spc_cat=plot_by_spc_cat,
                plot_by_hco_cat=plot_by_hco_cat,
                overwrite=True,
                spcdb_dir=spcdb_dir
            )

    #---------------------------------------------------------------
    # GCHP vs. GCHP tables of emission and inventory totals
    #---------------------------------------------------------------
    if emis_table:
        print("\n%%% Creating GCHP vs. GCHP emissions tables %%%")

        # Diagnostic collections to read
        col = "Emissions"
        col_met_ref = gchp_ref_metname()
        col_met_dev = gchp_dev_metname()

        # Read data
        ref = get_filepaths(gchp_vs_gchp_refdir, col,
                            all_months_mid_ref, is_gchp=True)
        dev = get_filepaths(gchp_vs_gchp_devdir, col,
                            all_months_mid_dev, is_gchp=True)
        refmet = get_filepaths(gchp_vs_gchp_refdir, col_met_ref,
                               all_months_mid_ref, is_gchp=True)
        devmet = get_filepaths(gchp_vs_gchp_devdir, col_met_dev,
                               all_months_mid_dev, is_gchp=True)

        # Create emissions table that spans entire year
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

    #---------------------------------------------------------------
    # GCHP vs. GCHP J-values plots
    #---------------------------------------------------------------
    if plot_jvalues:
        print("\n%%% Creating GCHP vs. GCHP J-values plots %%%")

        # Diagnostic collections to read
        col = "JValues"

        # Create J-value plots for each benchmark month
        for t in range(bmk_n_months):
            ref = get_filepath(gchp_vs_gchp_refdir, col,
                               bmk_mons_mid_ref[t], is_gchp=True)
            dev = get_filepath(gchp_vs_gchp_devdir, col,
                               bmk_mons_mid_dev[t], is_gchp=True)
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

    #---------------------------------------------------------------
    # GCHP vs GCHP column AOD plots
    #---------------------------------------------------------------
    if plot_aod:
        print("\n%%% Creating GCHP vs. GCHP AOD plots %%%")

        # Diagnostic collections to read
        col = "Aerosols"

        # Create AOD plots for each benchmark month
        for t in range(bmk_n_months):
            ref = get_filepath(gchp_vs_gchp_refdir, col,
                               bmk_mons_mid_ref[t], is_gchp=True)
            dev = get_filepath(gchp_vs_gchp_devdir, col,
                               bmk_mons_mid_dev[t], is_gchp=True)
            bmk.make_benchmark_aod_plots(
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

    #---------------------------------------------------------------
    # GCHP vs GCHP global mass tables
    #---------------------------------------------------------------
    if mass_table:
        print("\n%%% Creating GCHP vs. GCHP mass tables %%%")

        # Diagnostic collections to read
        col = "Restart"
        # Create mass table for each benchmark month in parallel
        def parallel_mass_table(t):

            # Time & date quantities
            reftime = bmk_mons_ref[t]
            devtime = bmk_mons_dev[t]
            datestr = bmk_mon_yr_strs_dev[t]
            label = "at 01{}".format(datestr)

            # Read Ref data
            ref = get_filepath(gchp_vs_gchp_refrstdir, col, reftime,
                               is_gchp=True)
            ref_extra=''
            if not os.path.isfile(ref):
                ref = join(gchp_vs_gchp_refrstdir,
                           'initial_GEOSChem_rst.' + gchp_ref_res
                           + '_benchmark.nc')
                extraref = bmk_mons_ref[t+1]
                ref_extra = get_filepath(gchp_vs_gchp_refrstdir, col,
                                         extraref, is_gchp=True)

            # Read Dev data
            dev = get_filepath(gchp_vs_gchp_devrstdir, col, devtime,
                               is_gchp=True)
            dev_extra=''
            if not os.path.isfile(dev):
                dev = join(gchp_vs_gchp_devrstdir,
                           'initial_GEOSChem_rst.' + gchp_dev_res
                           + '_benchmark.nc')
                extradev = bmk_mons_dev[t+1]
                dev_extra = get_filepath(gchp_vs_gchp_devrstdir, col,
                                         extradev, is_gchp=True)

            # Create tables
            bmk.make_benchmark_mass_tables(
                ref,
                gchp_vs_gchp_refstr,
                dev,
                gchp_vs_gchp_devstr,
                dst=gchp_vs_gchp_tablesdir,
                subdst=datestr,
                label=label,
                overwrite=True,
                spcdb_dir=spcdb_dir,
                ref_met_extra=ref_extra,
                dev_met_extra=dev_extra
            )

        # Run in parallel
        results = Parallel(n_jobs=-1)\
            (delayed(parallel_mass_table) (t) for t in range(bmk_n_months))

    #---------------------------------------------------------------
    # GCHP vs GCHP operations budgets tables
    #---------------------------------------------------------------
    if ops_budget_table:
        print("\n%%% Creating GCHP vs. GCHP operations budget tables %%%")

        # Diagnostic collections to read
        col = "Budget"
        def parallel_ops_budg(t):

            # Time & date quantities
            reftime = bmk_mons_mid_ref[t]
            devtime = bmk_mons_mid_dev[t]
            datestr = bmk_mon_yr_strs_dev[t]
            label = "at 01{}".format(datestr)
            sec_ref = bmk_sec_per_month_ref[t]
            sec_dev = bmk_sec_per_month_dev[t]
            # Create budget table for each benchmark month (ewl??)
            ref = get_filepath(gchp_vs_gchp_refdir, col, bmk_mons_mid_ref[t],
                               is_gchp=True)
            dev = get_filepath(gchp_vs_gchp_devdir, col, bmk_mons_mid_dev[t],
                               is_gchp=True)
            plot_dir = join(gchp_vs_gchp_budgetdir, bmk_mon_yr_strs_dev[t])
            bmk.make_benchmark_operations_budget(
                gchp_ref_version,
                ref,
                gchp_dev_version,
                dev,
                sec_ref,
                sec_dev,
                benchmark_type=bmk_type,
                label=datestr,
                operations=["Chemistry", "Convection", "EmisDryDep",
                            "Mixing", "WetDep"],
                compute_accum=False,
                dst=gchp_vs_gchp_tablesdir
            )

        # Run in parallel
        results = Parallel(n_jobs=-1)\
            (delayed(parallel_ops_budg) (t) for t in range(bmk_n_months))

    #---------------------------------------------------------------
    # GCHP vs GCHP aerosol budgets and burdens tables
    #---------------------------------------------------------------
    if aer_budget_table:
        print("\n%%% Skipping GCHP vs. GCHP aerosol budget tables: Redundant%%%")
        '''
        # Compute annual mean AOD budgets and aerosol burdens
        aerbdg.aerosol_budgets_and_burdens(
            gchp_dev_version,
            gchp_vs_gchp_devdir,
            bmk_year,
            dst=gchp_vs_gchp_tablesdir,
            overwrite=True,
            spcdb_dir=spcdb_dir
        )
        '''
    #---------------------------------------------------------------
    # GCHP vs. GCHP global mean OH, MCF Lifetime, CH4 Lifetime
    #---------------------------------------------------------------
    if OH_metrics:
        print("\n%%% Creating GCHP vs. GCHP OH metrics table %%%")

        # Diagnostic collection files to read
        col = "Metrics"
        ref = get_filepaths(gchp_vs_gchp_refdir, col, all_months_mid_ref, is_gchp=True)[0]
        dev = get_filepaths(gchp_vs_gchp_devdir, col, all_months_mid_dev, is_gchp=True)[0]

        # Create the OH Metrics table
        oh.make_benchmark_oh_metrics(
            ref,
            gchp_ref_version,
            dev,
            gchp_dev_version,
            dst=gchp_vs_gcc_tablesdir,
            overwrite=True,
            spcdb_dir=spcdb_dir
        )

    # --------------------------------------------------------------
    # GCHP Strat-Trop Exchange
    # --------------------------------------------------------------
    if ste_table:
        print("\n%%% Skipping GCHP vs. GCHP Strat-Trop Exchange table %%%")
