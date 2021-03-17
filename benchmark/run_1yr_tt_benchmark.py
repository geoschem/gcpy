#!/usr/bin/env python
"""
run_1yr_benchmark.py: Driver script for creating benchmark plots and testing
                      gcpy 1-year TransportTracers benchmark capability.

Run this script to generate benchmark comparisons between:

    (1) GCC (aka GEOS-Chem "Classic") vs. GCC
    (2) GCHP vs GCC (not yet tested)
    (3) GCHP vs GCHP (not yet tested)

You can customize this script by editing the following settings in the
"Configurables" section below:

    (1) Edit the path variables so that they point to folders w/ model data
    (2) Edit the version strings for each benchmark simulation
    (3) Edit the switches that turn on/off creating of plots and tables
    (4) If necessary, edit labels for the dev and ref versions

Calling sequence:

    ./run_1yr_tt_benchmark.py

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
from gcpy.util import get_filepath, get_filepaths
import gcpy.budget_tt as ttbdg
import gcpy.ste_flux as ste

# Tell matplotlib not to look for an X-window
os.environ["QT_QPA_PLATFORM"]="offscreen"

# Suppress annoying warning messages
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

# This script has a fixed benchmark type, year, and months
bmk_type     = "TransportTracersBenchmark"
bmk_year_ref = '2019'
bmk_year_dev = '2019'
bmk_mon_strs = ["Jan", "Apr", "Jul", "Oct"]
bmk_mon_inds = [0, 3, 6, 9]
bmk_n_months = len(bmk_mon_strs)

########################################################################
###           CONFIGURABLE SETTINGS: ***EDIT AS NEEDED***            ###
########################################################################

# =====================================================================
# Benchmark information
# Note: When doing GCHP vs GCC comparisions gchp_dev will be compared
# to gcc_dev (not gcc_ref!).
# =====================================================================

# High-level directory containing subdirectories with data
maindir  = "/n/holyscratch01/external_repos/GEOS-CHEM/gcgrid/geos-chem/validation/gcpy_test_data/1yr_transporttracer"

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
weightsdir = "/n/holyscratch01/external_repos/GEOS-CHEM/gcgrid/gcdata/ExtData/GCHP/RegriddingWeights"

# Path to species_databse.yml
spcdb_dir   = join(maindir, gcc_dev_version)

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
# GCHP vs GCC diff of diffs not included in transport tracer benchmark

# =====================================================================
# Output to generate (plots/tables will be created in this order):
# =====================================================================
plot_conc         = True
plot_wetdep       = True
rnpbbe_budget     = True
operations_budget = True
ste_table         = True  # GCC only
cons_table        = True

# =====================================================================
# Data directories
# For gchp_vs_gcc_refdir use gcc_dev_version, not ref (mps, 6/27/19)
# =====================================================================

# Diagnostic file directory paths
gcc_vs_gcc_refdir   = join(maindir, gcc_ref_version,  "OutputDir")
gcc_vs_gcc_devdir   = join(maindir, gcc_dev_version,  "OutputDir")
gchp_vs_gcc_refdir  = join(maindir, gcc_dev_version,  "OutputDir")
gchp_vs_gcc_devdir  = join(maindir, gchp_dev_version, "OutputDir")
gchp_vs_gchp_refdir = join(maindir, gchp_ref_version, "OutputDir")
gchp_vs_gchp_devdir = join(maindir, gchp_dev_version, "OutputDir")

# Restart file directory paths
gcc_vs_gcc_refrstdir   = join(maindir, gcc_ref_version, "restarts")
gcc_vs_gcc_devrstdir   = join(maindir, gcc_dev_version, "restarts")
gchp_vs_gcc_refrstdir  = join(maindir, gcc_dev_version, "restarts")
gchp_vs_gcc_devrstdir  = join(maindir, gchp_dev_version)
gchp_vs_gchp_refrstdir = join(maindir, gchp_ref_version)
gchp_vs_gchp_devrstdir = join(maindir, gchp_dev_version)

# Plots directories
if gcpy_test:
    mainresultsdir           = join('.', results_dir)
    gcc_vs_gcc_resultsdir    = join(mainresultsdir,'GCC_version_comparison')
    gchp_vs_gcc_resultsdir   = join(mainresultsdir,'GCHP_GCC_comparison')
    gchp_vs_gchp_resultsdir  = join(mainresultsdir,'GCHP_version_comparison')
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

# =====================================================================
# Dates and times -- Ref data
# =====================================================================

# Month/year strings for use in tabl4e subdirectories (e.g. Jan2016)
bmk_mon_yr_strs_ref = [v + bmk_year_ref for v in bmk_mon_strs]

# Get all months array of start datetimes for benchmark year
bmk_start_ref = np.datetime64(bmk_year_ref + "-01-01")
bmk_end_ref = np.datetime64("{}-01-01".format(int(bmk_year_ref)+1))
all_months_ref = np.arange(bmk_start_ref, bmk_end_ref,
                           step=np.timedelta64(1, "M"),
                           dtype="datetime64[M]")

# Get all months array of mid-point datetime per month for benchmark year
# and # sec in year
# NOTE: GCHP time-averaged files have time in the middle of the month
sec_per_yr_ref = 0
all_months_mid_ref = np.zeros(12, dtype="datetime64[h]")
for t in range(12):
    days_in_mon = monthrange(int(bmk_year_ref), t + 1)[1]
    sec_per_yr_ref += days_in_mon * 86400.0
    middle_hr = int(days_in_mon * 24 / 2)
    delta = np.timedelta64(middle_hr, 'h')
    all_months_mid_ref[t] = all_months_ref[t].astype("datetime64[h]") + delta

# Get subset of month datetimes for only benchmark months
bmk_mons_ref = all_months_ref[bmk_mon_inds]
bmk_mons_mid_ref = all_months_mid_ref[bmk_mon_inds]

# =====================================================================
# Dates and times -- Dev data
# =====================================================================

# Month/year strings for use in table subdirectories (e.g. Jan2016)
bmk_mon_yr_strs_dev = [v + bmk_year_dev for v in bmk_mon_strs]

# Get all months array of start datetimes for benchmark year
bmk_start_dev = np.datetime64(bmk_year_dev + "-01-01")
bmk_end_dev = np.datetime64("{}-01-01".format(int(bmk_year_dev)+1))
all_months_dev = np.arange(bmk_start_dev, bmk_end_dev,
                           step=np.timedelta64(1, "M"),
                           dtype="datetime64[M]")

# Get all months array of mid-point datetime per month for benchmark year
# and # sec in year
# NOTE: GCHP time-averaged files have time in the middle of the month
sec_per_yr_dev = 0
all_months_mid_dev = np.zeros(12, dtype="datetime64[h]")
for t in range(12):
    days_in_mon = monthrange(int(bmk_year_dev), t + 1)[1]
    sec_per_yr_dev += days_in_mon * 86400.0
    middle_hr = int(days_in_mon* 24 / 2)
    delta = np.timedelta64(middle_hr, 'h')
    all_months_mid_dev[t] = all_months_dev[t].astype("datetime64[h]") + delta

# Get subset of month datetimes for only benchmark months
bmk_mons_dev = all_months_dev[bmk_mon_inds]
bmk_mons_mid_dev = all_months_mid_dev[bmk_mon_inds]

# ======================================================================
# Print the list of plots & tables to the screen
# ======================================================================

print("The following plots and tables will be created for {}:".format(bmk_type))
if plot_conc:    print(" - Concentration plots")
if plot_wetdep:  print(" - Convective and large-scale wet deposition plots")
if rnpbbe_budget: print(" - Radionuclides budget table")
if operations_budget: print(" - Operations budget table")
if ste_table:    print(" - Table of strat-trop exchange")
if cons_table:   print(" - Table of mass conservation")
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
    # Separates RnPbBe tracers and passive tracers into separate files.
    # --------------------------------------------------------------
    if plot_conc:
        print("\n%%% Creating GCC vs. GCC concentration plots %%%")

        # Only plot concentration categories for TransportTracers
        restrict_cats = ["RnPbBeTracers", "PassiveTracers"]

        # Diagnostic collections to read
        col = "SpeciesConc"
        colmet = "StateMet"
        colmet_gchp="StateMet_avg" # Use this for benchmarks prior to 13.0

        # Create concentration plots for each benchmark month
        for t in range(bmk_n_months):

            # Time & date quantities
            reftime = bmk_mons_ref[t]
            devtime = bmk_mons_dev[t]
            datestr = bmk_mon_yr_strs_dev[t]

            # Seasonal diagnostic collection files to read
            ref = get_filepath(gcc_vs_gcc_refdir, col, reftime)
            dev = get_filepath(gcc_vs_gcc_devdir, col, devtime)
            refmet = get_filepath(gcc_vs_gcc_refdir, colmet, reftime)
            devmet = get_filepath(gcc_vs_gcc_devdir, colmet, devtime)

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
                restrict_cats=restrict_cats,
                overwrite=True,
                spcdb_dir=spcdb_dir
            )

    # --------------------------------------------------------------
    # GCC vs GCC wet deposition plots
    # --------------------------------------------------------------
    if plot_wetdep:
        print("\n%%% Creating GCC vs. GCC wet deposition plots %%%")

        # Diagnostic collection files to read
        cols = ["WetLossConv", "WetLossLS"]
        colmet = "StateMet"
        colmet_gchp="StateMet_avg" # Use this for benchmarks prior to 13.0

        # Loop over wet deposition collections and benchmark months
        for col in cols:
            for t in range(bmk_n_months):

                # Time & date quantities
                reftime = bmk_mons_ref[t]
                devtime = bmk_mons_dev[t]
                datestr = bmk_mon_yr_strs_dev[t]

                 # Seasonal diagnostic collection files to read
                ref = get_filepath(gcc_vs_gcc_refdir, col, reftime)
                dev = get_filepath(gcc_vs_gcc_devdir, col, devtime)
                refmet = get_filepath(gcc_vs_gcc_refdir, colmet, reftime)
                devmet = get_filepath(gcc_vs_gcc_devdir, colmet, devtime)

                # Make wet deposition plots
                bmk.make_benchmark_wetdep_plots(
                    ref,
                    gcc_vs_gcc_refstr,
                    dev,
                    gcc_vs_gcc_devstr,
                    refmet=refmet,
                    devmet=devmet,
                    dst=gcc_vs_gcc_resultsdir,
                    datestr=datestr,
                    weightsdir=weightsdir,
                    benchmark_type=bmk_type,
                    collection=col,
                    overwrite=True,
                    spcdb_dir=spcdb_dir
                )

    # --------------------------------------------------------------
    # GCC vs GCC radionuclides budget tables
    # --------------------------------------------------------------
    if rnpbbe_budget:
        print("\n%%% Creating GCC vs. GCC radionuclides budget table %%%")

        # Make radionuclides budget table
        ttbdg.transport_tracers_budgets(
            gcc_dev_version,
            gcc_vs_gcc_devdir,
            gcc_vs_gcc_devrstdir,
            int(bmk_year_dev),
            dst=gcc_vs_gcc_tablesdir,
            overwrite=True,
            spcdb_dir=spcdb_dir
        )

    # --------------------------------------------------------------
    # GCC vs GCC operations budgets tables
    # --------------------------------------------------------------
    if operations_budget:
        print("\n%%% Creating GCC vs. GCC operations budget tables %%%")

        # Diagnostic collection files to read (all 12 months)
        col = "Budget"
        refs = get_filepaths(gcc_vs_gcc_refdir, col, all_months_ref)
        devs = get_filepaths(gcc_vs_gcc_devdir, col, all_months_dev)

        # Make operations budget table
        bmk.make_benchmark_operations_budget(
            gcc_ref_version,
            refs,
            gcc_dev_version,
            devs,
            sec_per_yr_ref,
            sec_per_yr_dev,
            benchmark_type=bmk_type,
            label=bmk_year_dev,
            compute_accum=False,
            dst=gcc_vs_gcc_tablesdir
        )

    # --------------------------------------------------------------
    # GCC dev strat-trop exchange table
    # --------------------------------------------------------------
    if ste_table:
        print("\n%%% Creating GCC vs. GCC Strat-Trop Exchange table %%%")

        # Diagnostic collection files to read (all 12 months)
        col = "AdvFluxVert"
        devs = get_filepaths(gcc_vs_gcc_devdir, col, all_months_dev)[0]

        # Make stat-trop exchange table for subset of species
        species = ["Pb210","Be7","Be10"]
        ste.make_benchmark_ste_table(
            gcc_dev_version,
            devs,
            int(bmk_year_dev),
            dst=gcc_vs_gcc_tablesdir,
            bmk_type=bmk_type,
            species=species,
            overwrite=True
        )

# ======================================================================
# Create GCHP vs GCC benchmark plots and tables
# ======================================================================

if gchp_vs_gcc:

    # --------------------------------------------------------------
    # GCHP vs GCC Concentration plots
    # --------------------------------------------------------------
    if plot_conc:
        print("\n%%% Creating GCHP vs. GCC concentration plots %%%")

        # Only plot concentration categories for TransportTracers
        restrict_cats = ["RnPbBeTracers", "PassiveTracers"]

        # Diagnostic collections to read
        col = "SpeciesConc"
        colmet = "StateMet"
        colmet_gchp = "StateMet_avg" # Use this for benchmarks prior to 13.0

        # Create concentration plots for each benchmark month
        for t in range(bmk_n_months):

            # Time & date quantities
            reftime = bmk_mons_dev[t]
            devtime = bmk_mons_mid_dev[t]
            datestr = bmk_mon_yr_strs_dev[t]

            # Seasonal diagnostic collection files to read
            ref = get_filepath(gchp_vs_gcc_refdir, col, reftime)
            dev = get_filepath(gchp_vs_gcc_devdir, col, devtime,
                               is_gchp=True)
            devmet = get_filepath(gchp_vs_gcc_devdir, colmet, devtime,
                                  is_gchp=True)

            # Create plots
            bmk.make_benchmark_conc_plots(
                ref,
                gchp_vs_gcc_refstr,
                dev,
                gchp_vs_gcc_devstr,
                devmet=devmet,
                dst=gchp_vs_gcc_resultsdir,
                subdst=datestr,
                weightsdir=weightsdir,
                benchmark_type=bmk_type,
                restrict_cats=restrict_cats,
                overwrite=True,
                spcdb_dir=spcdb_dir
            )

    # --------------------------------------------------------------
    # GCHP vs GCC wet deposition plots
    # --------------------------------------------------------------
    if plot_wetdep:
        print("\n%%% Creating GCHP vs. GCC wet deposition plots %%%")

        # Create separate set of plots for each wetdep collection
        cols = ["WetLossConv", "WetLossLS"]
        colmet = "StateMet"
        colmet_gchp = "StateMet_avg" # Use this for benchmarks prior to 13.0

        # Create plots for each collection and benchmark month
        for col in cols:
            for t in range(bmk_n_months):

                # Time & date quantities
                reftime = bmk_mons_dev[t]
                devtime = bmk_mons_mid_dev[t]
                datestr = bmk_mon_yr_strs_dev[t]

                # Seasonal diagnostic quantities to read
                ref = get_filepath(gchp_vs_gcc_refdir, col, reftime)
                dev = get_filepath(gchp_vs_gcc_devdir, col, devtime,
                                   is_gchp=True)
                devmet = get_filepath(gchp_vs_gcc_devdir, colmet, devtime,
                                      is_gchp=True)

                # Create plots
                bmk.make_benchmark_wetdep_plots(
                    ref,
                    gchp_vs_gcc_refstr,
                    dev,
                    gchp_vs_gcc_devstr,
                    devmet=devmet,
                    collection=col,
                    dst=gchp_vs_gcc_resultsdir,
                    datestr=datestr,
                    weightsdir=weightsdir,
                    overwrite=True,
                    benchmark_type=bmk_type,
                    normalize_by_area=True,
                    spcdb_dir=spcdb_dir
                )

    # --------------------------------------------------------------
    # GCHP vs GCC radionuclides budget tables
    # --------------------------------------------------------------
    if rnpbbe_budget:
        print("\n%%% Creating GCHP vs. GCC radionuclides budget table %%%")
        # Make radionuclides budget table
        ttbdg.transport_tracers_budgets(
            gchp_dev_version,
            gchp_vs_gcc_devdir,
            gchp_vs_gcc_devrstdir,
            int(bmk_year_dev),
            dst=gchp_vs_gcc_tablesdir,
            is_gchp=True,
            overwrite=True,
            spcdb_dir=spcdb_dir
        )

    # --------------------------------------------------------------
    # GCHP vs GCC operations budgets tables
    # --------------------------------------------------------------
    if operations_budget:
        print("\n%%% Creating GCHP vs. GCC operations budget tables %%%")

        # Diagnostic collection files to read (all 12 months)
        col = "Budget"
        refs = get_filepaths(gchp_vs_gcc_refdir, col, all_months_dev)
        devs = get_filepaths(gchp_vs_gcc_devdir, col, all_months_mid_dev,
                             is_gchp=True)

        # Make operations budget table
        bmk.make_benchmark_operations_budget(
            gcc_dev_version,
            refs,
            gchp_dev_version,
            devs,
            sec_per_yr_ref,
            sec_per_yr_dev,
            benchmark_type=bmk_type,
            label=bmk_year_dev,
            operations=["Chemistry", "Convection", "EmisDryDep",
                        "Mixing", "WetDep"],
            compute_accum=False,
            dst=gchp_vs_gcc_tablesdir
        )

# =====================================================================
# Create GCHP vs GCHP benchmark plots and tables
# =====================================================================

if gchp_vs_gchp:

    # --------------------------------------------------------------
    # GCHP vs GCHP Concentration plots
    # --------------------------------------------------------------
    if plot_conc:
        print("\n%%% Creating GCHP vs. GCHP concentration plots %%%")

        # Only plot concentration categories for TransportTracers
        restrict_cats = ["RnPbBeTracers", "PassiveTracers"]

        # Diagnostic collections to read
        col = "SpeciesConc"
        colmet = "StateMet"
        colmet_gchp = "StateMet_avg" # Use this for benchmarks prior to 13.0

        # Create concentration plots for each benchmark month
        for t in range(bmk_n_months):

            # Time & date quantities
            reftime = bmk_mons_mid_ref[t]
            devtime = bmk_mons_mid_dev[t]
            datestr = bmk_mon_yr_strs_dev[t]

            # Seasonal diagnostic collection files to read
            ref = get_filepath(gchp_vs_gchp_refdir, col, reftime,
                               is_gchp=True)
            dev = get_filepath(gchp_vs_gchp_devdir, col, devtime,
                               is_gchp=True)
            refmet = get_filepath(gchp_vs_gchp_refdir, colmet, devtime,
                                  is_gchp=True)
            devmet = get_filepath(gchp_vs_gchp_devdir, colmet, devtime,
                                  is_gchp=True)

            # Make concentration plots
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
                restrict_cats=restrict_cats,
                overwrite=True,
                spcdb_dir=spcdb_dir
            )


    # --------------------------------------------------------------
    # GCHP vs GCHP wet deposition plots
    # --------------------------------------------------------------
    if plot_wetdep:
        print("\n%%% Creating GCHP vs. GCHP wet deposition plots %%%")

        # Create separate set of plots for each wetdep collection
        cols = ["WetLossConv", "WetLossLS"]
        colmet = "StateMet"
        colmet_gchp = "StateMet_avg" # Use this for benchmarks prior to 13.0

        # Create plots for each collection and benchmark month
        for col in cols:
            for t in range(bmk_n_months):

                # Time & date quantities
                reftime = bmk_mons_mid_ref[t]
                devtime = bmk_mons_mid_dev[t]
                datestr = bmk_mon_yr_strs_dev[t]

                # Seasonal diagnostic quantity files to read
                ref = get_filepath(gchp_vs_gchp_refdir, col, reftime,
                                   is_gchp=True)
                dev = get_filepath(gchp_vs_gchp_devdir, col, devtime,
                                   is_gchp=True)
                refmet = get_filepath(gchp_vs_gchp_refdir, colmet, reftime,
                                      is_gchp=True)
                devmet = get_filepath(gchp_vs_gchp_devdir, colmet, devtime,
                                      is_gchp=True)

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
                    datestr=datestr,
                    weightsdir=weightsdir,
                    overwrite=True,
                    benchmark_type=bmk_type,
                    normalize_by_area=True,
                    spcdb_dir=spcdb_dir
                )

    # --------------------------------------------------------------
    # GCHP vs GCHP radionuclides budget table
    # --------------------------------------------------------------
    if rnpbbe_budget:
        print("\n%%% Creating GCHP vs. GCHP radionuclides budget table %%%")

        # Make radionuclides budget table
        ttbdg.transport_tracers_budgets(
            gchp_dev_version,
            gchp_vs_gchp_devdir,
            gchp_vs_gchp_devrstdir,
            int(bmk_year_dev),
            dst=gchp_vs_gchp_tablesdir,
            is_gchp=True,
            overwrite=True,
            spcdb_dir=spcdb_dir
        )

    # --------------------------------------------------------------
    # GCHP vs GCHP operations budgets tables
    # --------------------------------------------------------------
    if operations_budget:
        print("\n%%% Creating GCHP vs. GCHP operations budget tables %%%")

        # Diagnostic collection files to read (all 12 months)
        col = "Budget"
        refs = get_filepaths(gchp_vs_gchp_refdir, col,
                             all_months_mid_ref, is_gchp=True)
        devs = get_filepaths(gchp_vs_gchp_devdir, col,
                             all_months_mid_dev, is_gchp=True)

        # Make operations budget table
        bmk.make_benchmark_operations_budget(
            gchp_dev_version,
            refs,
            gchp_dev_version,
            devs,
            sec_per_yr_ref,
            sec_per_yr_dev,
            benchmark_type=bmk_type,
            label=bmk_year_dev,
            operations=["Chemistry", "Convection", "EmisDryDep",
                        "Mixing", "WetDep"],
            compute_accum=False,
            dst=gchp_vs_gchp_tablesdir
        )

# ======================================================================
# Create mass conservations tables for GCC and GCHP
# ======================================================================

if cons_table:

    col = "Restart"

    # Create mass conservation table for gcc_ref
    if gcc_vs_gcc:        
        print("\n%%% Creating GCC ref mass conservation table %%%")
        
        # Get monthly restart files in the gcc refrst directory
        datafiles = get_filepaths(gcc_vs_gcc_refrstdir, col, all_months_ref,
                                  is_gchp=False)[0]
        
        # Make mass conservation table
        bmk.make_benchmark_mass_conservation_table(
            datafiles,
            gcc_ref_version,
            dst=gcc_vs_gcc_tablesdir,
            overwrite=True,
            spcdb_dir=spcdb_dir
        )

    if gcc_vs_gcc or gchp_vs_gcc:        
        print("\n%%% Creating GCC dev mass conservation table %%%")
        
        # Get monthly restart files in the gcc devrst directory
        datafiles = get_filepaths(gcc_vs_gcc_devrstdir, col, all_months_dev,
                                  is_gchp=False)[0]
        
        if gchp_vs_gcc:
            tablesdir=gchp_vs_gcc_tablesdir
        else:
            tablesdir=gcc_vs_gcc_tablesdir

        # Make mass conservation table
        bmk.make_benchmark_mass_conservation_table(
            datafiles,
            gcc_dev_version,
            dst=tablesdir,
            overwrite=True,
            spcdb_dir=spcdb_dir
        )
    
    if gchp_vs_gcc or gchp_vs_gchp:
        print("\n%%% Creating GCHP dev mass conservation table %%%")
        
        # Get monthly restart files in the gcc devrst directory
        datafiles = get_filepaths(gchp_vs_gcc_devrstdir, col, all_months_dev,
                                  is_gchp=True)[0]
        
        if gchp_vs_gcc:
            tablesdir=gchp_vs_gcc_tablesdir
        else:
            tablesdir=gchp_vs_gchp_tablesdir

        # Make mass conservation table
        bmk.make_benchmark_mass_conservation_table(
            datafiles,
            gchp_dev_version,
            dst=tablesdir,
            overwrite=True,
            spcdb_dir=spcdb_dir
        )

    if gchp_vs_gchp:
        print("\n%%% Creating GCHP ref mass conservation table %%%")
        
        # Get monthly restart files in the gcc devrst directory
        datafiles = get_filepaths(gchp_vs_gchp_refrstdir, col, all_months_ref,
                                  is_gchp=True)[0]
        
        # Make mass conservation table
        bmk.make_benchmark_mass_conservation_table(
            datafiles,
            gchp_ref_version,
            dst=gchp_vs_gchp_tablesdir,
            overwrite=True,
            spcdb_dir=spcdb_dir
        )
