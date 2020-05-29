#!/usr/bin/env python
"""
run_1yr_fullchem_benchmark.py: Driver script for creating benchmark plots and
                               testing gcpy 1-year full-chemistry benchmark
                               capability.

Run this script to generate benchmark comparisons between:

    (1) GCC (aka GEOS-Chem "Classic") vs. GCC
    (2) GCHP vs GCC (not yet tested)
    (3) GCHP vs GCHP (not yet tested)

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
"""

# =====================================================================
# Imports and global settings (you should not need to edit these)
# =====================================================================

from calendar import monthrange
import os
from os.path import join
import xarray as xr
from gcpy import benchmark as bmk
from gcpy.util import get_filepaths
import gcpy.ste_flux as ste
import gcpy.budget_aer as aerbdg
import gcpy.budget_ops as opbdg
import gcpy.mean_oh_from_logs as moh
import numpy as np
import warnings

# Tell matplotlib not to look for an X-window
os.environ["QT_QPA_PLATFORM"]="offscreen"

# Suppress harmless run-time warnings (mostly about underflow in division)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# This script has a fixed benchmark type
bmk_type     = "FullChemBenchmark"

########################################################################
###           CONFIGURABLE SETTINGS: EDIT THESE ACCORDINGLY          ###
########################################################################

# =====================================================================
# Benchmark information (**EDIT AS NEEDED**)
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

# Path to regridding weights (*MUST EDIT*)
weightsdir = "/n/holyscratch01/external_repos/GEOS-CHEM/gcgrid/gcdata/ExtData/GCHP/RegriddingWeights"

# =====================================================================
# Specify if this is a gcpy test validation run
# =====================================================================
gcpy_test = True

# =====================================================================
# Comparisons to run (**EDIT AS NEEDED**)
# =====================================================================
gcc_vs_gcc   = True
gchp_vs_gcc  = False # not yet tested
gchp_vs_gchp = False # not yet tested
# GCHP vs GCC diff of diffs not included in 1-yr full chemistry benchmark

# =====================================================================
# Output to generate (**EDIT AS NEEDED**)
# Plots/tables will be created in this order:
# =====================================================================
plot_conc    = True
plot_emis    = True
emis_table   = True
plot_jvalues = True
plot_aod     = True
mass_table   = True
budget_table = True
ste_table    = True
OH_metrics   = True

# Plot concentrations and emissions by category?
plot_by_spc_cat = True
plot_by_hco_cat = True

# =====================================================================
# Data directories (**EDIT AS NEEDED**)
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

# Plots directories
if gcpy_test:
    mainplotsdir          = './Plots'
    gcc_vs_gcc_plotsdir    = join(mainplotsdir,'GCC_version_comparison')
    gchp_vs_gchp_plotsdir  = join(mainplotsdir,'GCHP_version_comparison')
    gchp_vs_gcc_plotsdir   = join(mainplotsdir,'GCHP_GCC_comparison')
    if not os.path.exists(mainplotsdir): os.mkdir(mainplotsdir)
else:
    gcc_vs_gcc_plotsdir    = join(maindir, gcc_dev_version, "Plots")
    gchp_vs_gchp_plotsdir  = join(maindir, gchp_dev_version,
                              "Plots", "GCHP_version_comparison")
    gchp_vs_gcc_plotsdir   = join(maindir, gchp_dev_version,
                              "Plots", "GCHP_GCC_comparison")

# Tables directories
gcc_vs_gcc_tablesdir   = join(gcc_vs_gcc_plotsdir,"Tables") 
gchp_vs_gcc_tablesdir  = join(gchp_vs_gcc_plotsdir,"Tables") 
gchp_vs_gchp_tablesdir = join(gchp_vs_gchp_plotsdir,"Tables")

# =====================================================================
# Plot title strings (edit as needed)
# For gchp_vs_gcc_refstr use gcc_dev_version, not ref (mps, 6/27/19)
# =====================================================================
gcc_vs_gcc_refstr    = "{}".format(gcc_ref_version)
gcc_vs_gcc_devstr    = "{}".format(gcc_dev_version)
gchp_vs_gcc_refstr   = "{}".format(gcc_dev_version)
gchp_vs_gcc_devstr   = "{}".format(gchp_dev_version)
gchp_vs_gchp_refstr  = "{}".format(gchp_ref_version)
gchp_vs_gchp_devstr  = "{}".format(gchp_dev_version)

########################################################################
###    THE REST OF THESE SETTINGS SHOULD NOT NEED TO BE CHANGED      ###
########################################################################

# =====================================================================
# Dates and times
# =====================================================================

# Start and end of the benchmark year (as numpy.datetime64 dates)
bmk_year = 2016
bmk_start_str = "{}-01-01".format(str(bmk_year))
bmk_end_str = "{}-01-01".format(str(bmk_year+1))
bmk_start = np.datetime64(bmk_start_str)
bmk_end = np.datetime64(bmk_end_str)

# Monthly array of dates
bmk_delta_1m = np.timedelta64(1, "M")
bmk_months = np.arange(bmk_start, bmk_end,
                       step=bmk_delta_1m, dtype="datetime64[M]")

# Get the benchmark year from the datetime64 object
bmk_year = bmk_months[0].astype("datetime64[Y]").astype(int) + 1970

# Seasonal array of dates
bmk_delta_3m = np.timedelta64(3, "M")
bmk_seasons = np.arange(bmk_start, bmk_end,
                        step=bmk_delta_3m, dtype="datetime64[M]")
bmk_nseasons = len(bmk_seasons)

# Names for each season (e.g. Jan2016, Apr2016, Jul2016, Oct2016)
bmk_seasons_names = ["Jan", "Apr", "Jul", "Oct"]
bmk_seasons_names = [v + str(bmk_year) for v in bmk_seasons_names]

# Seconds in each month of the benchmark year
sec_per_month = np.zeros(12)
for m in range(1, 13):
    month_info = monthrange(bmk_year, m)
    sec_per_month[m-1] = month_info[1] * 86400.0

## Timestamps for GCHP (these are in the middle of the month)
#gchp_months = np.zeros(12, dtype="datetime64[h]")
#for m in range(12):
#    if days_per_month[m] == 31:
#        delta = np.timedelta64(((15 * 24) + 12), 'h')
#    elif days_per_month[m] == 30:
#        delta = np.timedelta64((15 * 24), 'h')
#    elif days_per_month[m] == 29:
#        delta = np.timedelta64(((14 * 24) + 12), 'h')
#    else:
#        delta = np.timedelta64((14 * 24), 'h')
#    gchp_months[m] = bmk_months[m].astype("datetime64[h]") + delta
#gchp_seasons = gchp_months[[0, 3, 6, 9]]

# ======================================================================
# Echo the list of plots & tables that will be made to the screen
# ======================================================================

print("The following plots and tables will be created for {}:".format(bmk_type))
if plot_conc:    print(" - Concentration plots")
if plot_emis:    print(" - Emissions plots")
if plot_jvalues: print(" - J-values (photolysis rates) plots")
if plot_aod:     print(" - Aerosol optical depth plots")
if budget_table: print(" - Budget tables")
if emis_table:   print(" - Table of emissions totals by species and inventory")
if mass_table:   print(" - Table of species mass")
if OH_metrics:   print(" - Table of OH metrics")
if ste_table:    print(" - Table of strat-trop exchange")
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
    # is true; otherwise excludes lumped species and writes to one file
    # --------------------------------------------------------------
    if plot_conc:
        print("\n%%% Creating GCC vs. GCC concentration plots %%%")

        # File lists for emissions data (seasonal)
        collection = "SpeciesConc"
        gcc_vs_gcc_refspc = get_filepaths(gcc_vs_gcc_refdir, collection,
                                          bmk_seasons, is_gcc=True)
        gcc_vs_gcc_devspc = get_filepaths(gcc_vs_gcc_devdir, collection,
                                          bmk_seasons, is_gcc=True)

        # Create seasonal concentration plots
        for s in range(bmk_nseasons):
            mon_yr_str = bmk_seasons_names[s]
            bmk.make_benchmark_plots(gcc_vs_gcc_refspc[s],
                                     gcc_vs_gcc_refstr,
                                     gcc_vs_gcc_devspc[s],
                                     gcc_vs_gcc_devstr,
                                     dst=gcc_vs_gcc_plotsdir,
                                     subdst=mon_yr_str,
                                     weightsdir=weightsdir,
                                     benchmark_type=bmk_type,
                                     collection=collection,
                                     plot_by_spc_cat=plot_by_spc_cat,
                                     overwrite=True)

    # --------------------------------------------------------------
    # GCC vs GCC emissions plots
    # --------------------------------------------------------------
    if plot_emis:
        print("\n%%% Creating GCC vs. GCC emissions plots %%%")

        # File lists for emissions data (seasonal)
        gcc_vs_gcc_refhco = get_filepaths(gcc_vs_gcc_refdir, "Emissions",
                                          bmk_seasons, is_gcc=True)
        gcc_vs_gcc_devhco = get_filepaths(gcc_vs_gcc_devdir, "Emissions",
                                          bmk_seasons, is_gcc=True)

        # Create seasonal emissions plots
        for s in range(bmk_nseasons):
            mon_yr_str = bmk_seasons_names[s]
            bmk.make_benchmark_emis_plots(gcc_vs_gcc_refhco[s],
                                          gcc_vs_gcc_refstr,
                                          gcc_vs_gcc_devhco[s],
                                          gcc_vs_gcc_devstr,
                                          dst=gcc_vs_gcc_plotsdir,
                                          subdst=mon_yr_str,
                                          weightsdir=weightsdir,
                                          plot_by_spc_cat=plot_by_spc_cat,
                                          plot_by_hco_cat=plot_by_hco_cat,
                                          overwrite=True)

    # --------------------------------------------------------------
    # GCC vs GCC tables of emission and inventory totals
    # --------------------------------------------------------------
    if emis_table:
        print("\n%%% Creating GCC vs. GCC emissions & inventory totals %%%")

        # File lists for J-values data (monthly)
        gcc_vs_gcc_refhco = get_filepaths(gcc_vs_gcc_refdir, "Emissions",
                                          bmk_months, is_gcc=True)
        gcc_vs_gcc_devhco = get_filepaths(gcc_vs_gcc_devdir, "Emissions",
                                          bmk_months, is_gcc=True)

        # Create emission tables
        bmk.make_benchmark_emis_tables(gcc_vs_gcc_refhco,
                                       gcc_vs_gcc_refstr,
                                       gcc_vs_gcc_devhco,
                                       gcc_vs_gcc_devstr,
                                       dst=gcc_vs_gcc_plotsdir,
                                       interval=sec_per_month,
                                       overwrite=True)

    # --------------------------------------------------------------
    # GCC vs GCC J-value plots
    # --------------------------------------------------------------
    if plot_jvalues:
        print("\n%%% Creating GCC vs. GCC J-value plots %%%")

        # Paths to J-value data (seasonal)
        gcc_vs_gcc_refjv = get_filepaths(gcc_vs_gcc_refdir, "JValues",
                                          bmk_seasons, is_gcc=True)
        gcc_vs_gcc_devjv = get_filepaths(gcc_vs_gcc_devdir, "JValues",
                                          bmk_seasons, is_gcc=True)

        # Create seasonal J-values plots
        for s in range(bmk_nseasons):
            mon_yr_str = bmk_seasons_names[s]
            bmk.make_benchmark_jvalue_plots(gcc_vs_gcc_refjv[s],
                                            gcc_vs_gcc_refstr,
                                            gcc_vs_gcc_devjv[s],
                                            gcc_vs_gcc_devstr,
                                            dst=gcc_vs_gcc_plotsdir,
                                            subdst=mon_yr_str,
                                            weightsdir=weightsdir,
                                            overwrite=True)

    # --------------------------------------------------------------
    # GCC vs. GCC column AOD plots
    # --------------------------------------------------------------
    if plot_aod:
        print("\n%%% Creating GCC vs. GCC column AOD plots %%%")

        ## Paths to aerosol optical depth data
        gcc_vs_gcc_refaod = get_filepaths(gcc_vs_gcc_refdir, "Aerosols",
                                          bmk_seasons, is_gcc=True)
        gcc_vs_gcc_devaod = get_filepaths(gcc_vs_gcc_devdir, "Aerosols",
                                          bmk_seasons, is_gcc=True)

        # Create seasonal column AOD plots
        for s in range(bmk_nseasons):
            mon_yr_str = bmk_seasons_names[s]
            bmk.make_benchmark_aod_plots(gcc_vs_gcc_refaod[s],
                                         gcc_vs_gcc_refstr,
                                         gcc_vs_gcc_devaod[s],
                                         gcc_vs_gcc_devstr,
                                         dst=gcc_vs_gcc_plotsdir,
                                         subdst=mon_yr_str,
                                         weightsdir=weightsdir,
                                         overwrite=True)

    # --------------------------------------------------------------
    # GCC vs GCC mass tables
    # --------------------------------------------------------------
    if mass_table:
        print("\n%%% Creating GCC vs. GCC mass tables %%%")

        ## Paths to restart files
        gcc_vs_gcc_refrst = get_filepaths(gcc_vs_gcc_refrstdir, "Restart",
                                          bmk_seasons, is_gcc=True)
        gcc_vs_gcc_devrst = get_filepaths(gcc_vs_gcc_devrstdir, "Restart",
                                          bmk_seasons, is_gcc=True)

        # Create seasonal budget tables (mass at end of each season month)
        for s in range(bmk_nseasons):
            mon_yr_str = bmk_seasons_names[s]
            label = "at 01{}".format(mon_yr_str)
            plot_dir = join(gcc_vs_gcc_plotsdir, "Tables", mon_yr_str)
            bmk.make_benchmark_mass_tables(gcc_vs_gcc_refrst[s],
                                           gcc_vs_gcc_refstr,
                                           gcc_vs_gcc_devrst[s],
                                           gcc_vs_gcc_devstr,
                                           dst=plot_dir,
                                           label=label,
                                           overwrite=True,
                                           subdst=mon_yr_str)

    # --------------------------------------------------------------
    # GCC vs GCC operations budgets tables
    # --------------------------------------------------------------
    if budget_table:
        print("\n%%% Creating GCC vs. GCC budget tables %%%")

        # Change in mass of species after each operation
        print("-- Change in mass after each operation")
        collection = "Budget"
        gcc_vs_gcc_reflist = get_filepaths(gcc_vs_gcc_refdir, collection,
                                           bmk_months, is_gcc=True)
        gcc_vs_gcc_devlist = get_filepaths(gcc_vs_gcc_devdir, collection,
                                           bmk_months, is_gcc=True)

        # Compute change in mass after each operation
        for s in range(bmk_nseasons):
            mon_yr_str = bmk_seasons_names[s]
            plot_dir = join(gcc_vs_gcc_plotsdir, 'Budget', mon_yr_str)
            opbdg.make_operations_budget_table(gcc_ref_version,
                                               gcc_vs_gcc_reflist[s],
                                               gcc_dev_version,
                                               gcc_vs_gcc_devlist[s],
                                               bmk_type,
                                               dst=plot_dir,
                                               label=mon_yr_str,
                                               interval=sec_per_month[s*3],
                                               overwrite=True)

        # Compute annual mean AOD budgets and aerosol burdens
        print("-- AOD budgets and aerosol burdens")
        plot_dir = join(gcc_vs_gcc_plotsdir, "Tables")
        aerbdg.aerosol_budgets_and_burdens(gcc_dev_version,
                                           gcc_vs_gcc_devdir,
                                           bmk_year,
                                           dst=plot_dir,
                                           overwrite=True)

    # --------------------------------------------------------------
    # GCC Strat-Trop Exchange
    # --------------------------------------------------------------
    if ste_table:
        print("\n%%% Creating GCC vs. GCC Strat-Trop Exchange table %%%")

        # Compute monthly and annual average strat-trop exchange of O3
        gcc_vs_gcc_devflx = get_filepaths(gcc_vs_gcc_devdir, "AdvFluxVert",
                                          bmk_months, is_gcc=True)
        plot_dir = join(gcc_vs_gcc_plotsdir, 'Tables')
        ste.make_benchmark_ste_table(gcc_dev_version,
                                     gcc_vs_gcc_devflx,
                                     bmk_year,
                                     dst=plot_dir,
                                     bmk_type=bmk_type,
                                     species=['O3'],
                                     overwrite=True)

    # --------------------------------------------------------------
    # GCC vs GCC Global mean OH, MCF Lifetime, CH4 Lifetime
    # --------------------------------------------------------------
    if OH_metrics:
        print("\n%%% Creating GCC vs. GCC OH metrics %%%")

        ####################################################################
        # NOTE: Need to better validate this routine
        # for now, use the mean OH from the log files (bmy, 3/12/20)
        ## Paths to data files
        #collections = ["ConcAfterChem", "StateMet"]
        #gcc_vs_gcc_reflist = get_filepaths(gcc_vs_gcc_refdir, collections,
        #                                   bmk_months, is_gcc=True)
        #gcc_vs_gcc_devlist = get_filepaths(gcc_vs_gcc_devdir, collections,
        #                                   bmk_months, is_gcc=True)
        #
        ## Create OH metrics table
        #plot_dir = join(gcc_vs_gcc_plotsdir, "Tables")
        #bmk.make_benchmark_oh_metrics(gcc_vs_gcc_reflist,
        #                              gcc_vs_gcc_refstr,
        #                              gcc_vs_gcc_devlist,
        #                              gcc_vs_gcc_devstr,
        #                              dst=plot_dir,
        #                              overwrite=True)
        #####################################################################

        # Compute mean OH from the log files
        # NOTE: Only works for GEOS-Chem "Classic" benchmarks!
        plot_dir = join(gcc_vs_gcc_plotsdir, "Tables")
        moh.make_benchmark_oh_from_logs(gcc_vs_gcc_reflogdir,
                                        gcc_vs_gcc_refstr,
                                        gcc_vs_gcc_devlogdir,
                                        gcc_vs_gcc_devstr,
                                        bmk_year,
                                        dst=plot_dir,
                                        overwrite=True)

# ======================================================================
# Create GCHP vs GCC benchmark plots and tables
# ======================================================================

if gchp_vs_gcc:

    # --------------------------------------------------------------
    # GCHP vs GCC Concentration plots
    # --------------------------------------------------------------
    if plot_conc:
        print("\n%%% Skipping GCHP vs. GCC concentration plots %%%")
#        # Concentration plots
#        # (includes lumped species and separates by category)
#        print("\n%%% Creating GCHP vs. GCC concentration plots %%%")
#        bmk.make_benchmark_conc_plots(gchp_vs_gcc_refspc,
#                                      gchp_vs_gcc_refstr,
#                                      gchp_vs_gcc_devspc,
#                                      gchp_vs_gcc_devstr,
#                                      dst=gchp_vs_gcc_plotsdir,
#                                      plot_by_spc_cat=plot_by_spc_cat,
#                                      overwrite=True)

    # --------------------------------------------------------------
    # GCHP vs GCC emissions plots
    # --------------------------------------------------------------
    if plot_emis:
        print("\n%%% Skipping GCHP vs. GCC emissions plots %%%")
#        # Emissions plots
#        print("\n%%% Creating GCHP vs. GCC emissions plots %%%")
#        bmk.make_benchmark_emis_plots(gchp_vs_gcc_refhco,
#                                      gchp_vs_gcc_refstr,
#                                      gchp_vs_gcc_devhco,
#                                      gchp_vs_gcc_devstr,
#                                      dst=gchp_vs_gcc_plotsdir,
#                                      plot_by_spc_cat=plot_by_spc_cat,
#                                      plot_by_hco_cat=plot_by_hco_cat,
#                                      overwrite=True,
#                                      flip_dev=True)

    if emis_table:
        print("\n%%% Skipping GCHP vs. GCC emissions tables %%%")
#        # Tables of emissions and inventory totals
#        print("\n%%% Creating GCHP vs. GCC emissions and inventory tables %%%")
#        gchp_vs_gcc_reflist = [gchp_vs_gcc_refhco]
#        gchp_vs_gcc_devlist = [gchp_vs_gcc_devhco, gchp_vs_gcc_devmet]
#        bmk.make_benchmark_emis_tables(gchp_vs_gcc_reflist,
#                                       gchp_vs_gcc_refstr,
#                                       gchp_vs_gcc_devlist,
#                                       gchp_vs_gcc_devstr,
#                                       dst=gchp_vs_gcc_plotsdir,
#                                       overwrite=True)

    if plot_jvalues:
        print("\n%%% Skipping GCHP vs. GCC J-values plots %%%")
#        # Local noon J-values plots
#        print("\n%%% Creating GCHP vs. GCC J-value plots %%%")
#        bmk.make_benchmark_jvalue_plots(gchp_vs_gcc_refjv,
#                                        gchp_vs_gcc_refstr,
#                                        gchp_vs_gcc_devjv,
#                                        gchp_vs_gcc_devstr,
#                                        dst=gchp_vs_gcc_plotsdir,
#                                        overwrite=True)

    if plot_aod:
        print("\n%%% Skipping GCHP vs. GCC AOD plots %%%")
#        # Column AOD plots
#        print("\n%%% Creating GCHP vs. GCC column AOD plots %%%")
#        bmk.make_benchmark_aod_plots(gchp_vs_gcc_refaod,
#                                     gchp_vs_gcc_refstr,
#                                     gchp_vs_gcc_devaod,
#                                     gchp_vs_gcc_devstr,
#                                     dst=gchp_vs_gcc_plotsdir,
#                                     overwrite=True)

    if mass_table:
        print("\n%%% Skipping GCHP vs. GCC mass tables %%%")
#        # Global mass tables
#        print("\n%%% Creating GCHP vs. GCC global mass tables %%%")
#        gchp_vs_gcc_reflist = [gchp_vs_gcc_refrst]
#        gchp_vs_gcc_devlist = [gchp_vs_gcc_devrst, gchp_vs_gcc_devmetinst]
#        bmk.make_benchmark_mass_tables(gchp_vs_gcc_reflist,
#                                       gchp_vs_gcc_refstr,
#                                       gchp_vs_gcc_devlist,
#                                       gchp_vs_gcc_devstr,
#                                       dst=gchp_vs_gcc_plotsdir,
#                                       overwrite=True)

    if OH_metrics:
        print("\n%%% Skipping GCHP vs. GCC OH metrics %%%")
#        # Global mean OH, MCF Lifetime, CH4 Lifetime
#        print("\n%%% Creating GCHP vs. GCC OH metrics %%%")
#        gchp_vs_gcc_reflist = [gchp_vs_gcc_refcac, gchp_vs_gcc_refmet]
#        gchp_vs_gcc_devlist = [gchp_vs_gcc_devcac, gchp_vs_gcc_devmet]
#        bmk.make_benchmark_oh_metrics(gchp_vs_gcc_reflist,
#                                      gchp_vs_gcc_refstr,
#                                      gchp_vs_gcc_devlist,
#                                      gchp_vs_gcc_devstr,
#                                      dst=gchp_vs_gcc_plotsdir,
#                                      overwrite=True)

# =====================================================================
# Create GCHP vs GCHP benchmark plots and tables
# =====================================================================

if gchp_vs_gchp:

    # --------------------------------------------------------------
    # GCHP vs GCHP Concentration plots
    # --------------------------------------------------------------
    if plot_conc:
        print("\n%%% Skipping GCHP vs. GCHP concentration plots %%%")
#        # Concentration plots
#        # (includes lumped species and separates by category)
#        print("\n%%% Creating GCHP vs. GCHP concentration plots %%%")
#        bmk.make_benchmark_conc_plots(gchp_vs_gchp_refspc,
#                                      gchp_vs_gchp_refstr,
#                                      gchp_vs_gchp_devspc,
#                                      gchp_vs_gchp_devstr,
#                                      dst=gchp_vs_gchp_plotsdir,
#                                      plot_by_spc_cat=plot_by_spc_cat,
#                                      overwrite=True)


    # --------------------------------------------------------------
    # GCHP vs GCHP emissions plots
    # --------------------------------------------------------------
    if plot_emis:
        print("\n%%% Skipping GCHP vs. GCHP emissions plots %%%")
#        # Emissions plots
#        print("\n%%% Creating GCHP vs. GCHP emissions plots %%%")
#        bmk.make_benchmark_emis_plots(gchp_vs_gchp_refhco,
#                                      gchp_vs_gchp_refstr,
#                                      gchp_vs_gchp_devhco,
#                                      gchp_vs_gchp_devstr,
#                                      dst=gchp_vs_gchp_plotsdir,
#                                      plot_by_spc_cat=plot_by_spc_cat,
#                                      plot_by_hco_cat=plot_by_hco_cat,
#                                      overwrite=True,
#                                      flip_ref=True,
#                                      flip_dev=True)

    if emis_table:
        print("\n%%% Skipping GCHP vs. GCHP emissions tables %%%")
#        # Tables of emissions and inventory totals
#        print("\n%%% Creating GCHP vs. GCHP emissions and inventory tables %%%")
#        gchp_vs_gchp_reflist = [gchp_vs_gchp_refhco, gchp_vs_gchp_refmet]
#        gchp_vs_gchp_devlist = [gchp_vs_gchp_devhco, gchp_vs_gchp_devmet]
#        bmk.make_benchmark_emis_tables(gchp_vs_gchp_reflist,
#                                       gchp_vs_gchp_refstr,
#                                       gchp_vs_gchp_devlist,
#                                       gchp_vs_gchp_devstr,
#                                       dst=gchp_vs_gchp_plotsdir,
#                                       overwrite=True)

    if plot_jvalues:
        print("\n%%% Skipping GCHP vs. GCHP J-values plots %%%")
#        # Local noon J-values plots
#        print("\n%%% Creating GCHP vs. GCHP J-value plots %%%")
#        bmk.make_benchmark_jvalue_plots(gchp_vs_gchp_refjv,
#                                        gchp_vs_gchp_refstr,
#                                        gchp_vs_gchp_devjv,
#                                        gchp_vs_gchp_devstr,
#                                        dst=gchp_vs_gchp_plotsdir,
#                                        overwrite=True)

    if plot_aod:
        print("\n%%% Skipping GCHP vs. GCHP AOD plots %%%")
#        # Column AOD plots
#        print("\n%%% Creating GCHP vs. GCHP column AOD plots %%%")
#        bmk.make_benchmark_aod_plots(gchp_vs_gchp_refaod,
#                                     gchp_vs_gchp_refstr,
#                                     gchp_vs_gchp_devaod,
#                                     gchp_vs_gchp_devstr,
#                                     dst=gchp_vs_gchp_plotsdir,
#                                     overwrite=True)

    if mass_table:
        print("\n%%% Skipping GCHP vs. GCHP mass tables %%%")
#        # Global mass tables
#        print("\n%%% Creating GCHP vs. GCHP global mass tables %%%")
#        gchp_vs_gchp_reflist = [gchp_vs_gchp_refrst, gchp_vs_gchp_refmetinst]
#        gchp_vs_gchp_devlist = [gchp_vs_gchp_devrst, gchp_vs_gchp_devmetinst]
#        bmk.make_benchmark_mass_tables(gchp_vs_gchp_reflist,
#                                       gchp_vs_gchp_refstr,
#                                       gchp_vs_gchp_devlist,
#                                       gchp_vs_gchp_devstr,
#                                       dst=gchp_vs_gchp_plotsdir,
#                                       overwrite=True)

    if OH_metrics:
        print("\n%%% Skipping GCHP vs. GCHP OH metrics %%%")
#        # Global mean OH, MCF Lifetime, CH4 Lifetime
#        print("\n%%% Creating GCHP vs. GCHP OH metrics %%%")
#        gchp_vs_gcc_reflist = [gchp_vs_gchp_refcac, gchp_vs_gchp_refmet]
#        gchp_vs_gcc_devlist = [gchp_vs_gchp_devcac, gchp_vs_gchp_devmet]
#        bmk.make_benchmark_oh_metrics(gchp_vs_gchp_reflist,
#                                      gchp_vs_gchp_refstr,
#                                      gchp_vs_gchp_devlist,
#                                      gchp_vs_gchp_devstr,
#                                      dst=gchp_vs_gchp_plotsdir,
#                                      overwrite=True)

