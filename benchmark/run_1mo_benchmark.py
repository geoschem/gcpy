#!/usr/bin/env python
"""
run_1mo_benchmark.py: Driver script for creating benchmark plots and
                      testing gcpy 1-month benchmark capability

Run this script to generate benchmark comparisons between:

    (1) GCC (aka GEOS-Chem "Classic") vs. GCC
    (2) GCHP vs GCC
    (3) GCHP vs GCHP
    (4) GCHP vs GCC diff-of-diffs

You can customize this by editing the following settings in the
"Configurables" section below:

    (1) Edit the path variables so that they point to folders w/ model data
    (2) Edit the version strings for each benchmark simulation
    (3) Edit the switches that turn on/off creating of plots and tables
    (4) If necessary, edit labels for the dev and ref versions

Calling sequence:

    ./run_1mo_benchmark.py

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
# ====================================q=================================

import os
from os.path import join

import warnings
import calendar
import numpy as np
import xarray as xr
from gcpy.util import get_filepath
from gcpy import benchmark as bmk
import gcpy.ste_flux as ste

# Tell matplotlib not to look for an X-window
os.environ["QT_QPA_PLATFORM"]="offscreen"

# Suppress harmless run-time warnings (mostly about underflow in division)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# This script has a fixed benchmark type
bmk_type     = "FullChemBenchmark"
bmk_year     = '2016'
bmk_mon      = '7'

########################################################################
###           CONFIGURABLE SETTINGS: ***EDIT AS NEEDED***            ###
########################################################################

# =====================================================================
# Benchmark information (**EDIT AS NEEDED**)
# Note: When doing GCHP vs GCC comparisions gchp_dev will be compared
# to gcc_dev (not gcc_ref!). This ensures consistency in version names
# when doing GCHP vs GCC diff-of-diffs (mps, 6/27/19)
# =====================================================================

# High-level directory containing subdirectories with data
maindir  = "/n/holyscratch01/external_repos/GEOS-CHEM/gcgrid/geos-chem/validation/gcpy_test_data/1mon"

# Version strings
# NOTE: these will be used in some filenames and so should not have spaces
# or other characters not appropriate for a filename.
gcc_ref_version = "GCC_ref"
gcc_dev_version = "GCC_dev"
gchp_ref_version = "GCHP_ref"
gchp_dev_version = "GCHP_dev"

# Name to be used for directory of output from this script
results_dir = "Results"

# Path to regridding weights
weightsdir = "/n/holyscratch01/external_repos/GEOS-CHEM/gcgrid/gcdata/ExtData/GCHP/RegriddingWeights"

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
gchp_vs_gcc_diff_of_diffs = True

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
OH_metrics       = True # GCC only
ste_table        = True # GCC only

# Plot concentrations and emissions by category?
plot_by_spc_cat = True
plot_by_hco_cat = True

# =====================================================================
# Data directories
# For gchp_vs_gcc_refdir use gcc_dev_version, not ref (mps, 6/27/19)
# =====================================================================

# Directory names (edit if not same as version strings)
gcc_ref_dir = gcc_ref_version
gcc_dev_dir = gcc_dev_version
gchp_ref_dir = gchp_ref_version
gchp_dev_dir = gchp_dev_version

# Diagnostic file directory paths
gcc_vs_gcc_refdir   = join(maindir, gcc_ref_dir,  "OutputDir")
gcc_vs_gcc_devdir   = join(maindir, gcc_dev_dir,  "OutputDir")
gchp_vs_gcc_refdir  = join(maindir, gcc_dev_dir,  "OutputDir")
gchp_vs_gcc_devdir  = join(maindir, gchp_dev_dir, "OutputDir")
gchp_vs_gchp_refdir = join(maindir, gchp_ref_dir, "OutputDir")
gchp_vs_gchp_devdir = join(maindir, gchp_dev_dir, "OutputDir")

# Restart file directory paths
gcc_vs_gcc_refrst   = join(maindir, gcc_ref_dir )
gcc_vs_gcc_devrst   = join(maindir, gcc_dev_dir )
gchp_vs_gcc_refrst  = join(maindir, gcc_dev_dir )
gchp_vs_gcc_devrst  = join(maindir, gchp_dev_dir)
gchp_vs_gchp_refrst = join(maindir, gchp_ref_dir)
gchp_vs_gchp_devrst = join(maindir, gchp_dev_dir)


# =====================================================================
# Path to species_databse.yml
# =====================================================================
spcdb_dir   = join(maindir, gcc_dev_dir)

# =====================================================================
# Benchmark output directories
# =====================================================================
# Results directories
if gcpy_test:
    mainresultsdir           = join(".", results_dir)
    gcc_vs_gcc_resultsdir    = join(mainresultsdir,'GCC_version_comparison')
    gchp_vs_gchp_resultsdir  = join(mainresultsdir,'GCHP_version_comparison')
    gchp_vs_gcc_resultsdir   = join(mainresultsdir,'GCHP_GCC_comparison')
    diff_of_diffs_resultsdir = join(mainresultsdir,'GCHP_GCC_diff_of_diffs')
    if not os.path.exists(mainresultsdir): os.mkdir(mainresultsdir)
else:
    gcc_vs_gcc_resultsdir    = join(maindir, gcc_dev_dir, results_dir)
    gchp_vs_gchp_resultsdir  = join(maindir, gchp_dev_dir,
                              results_dir, "GCHP_version_comparison")
    gchp_vs_gcc_resultsdir   = join(maindir, gchp_dev_dir,
                              results_dir, "GCHP_GCC_comparison")
    diff_of_diffs_resultsdir = join(maindir, gchp_dev_dir,
                              results_dir, "GCHP_GCC_diff_of_diffs")
gcc_vs_gcc_tablesdir    = join(gcc_vs_gcc_resultsdir, "Tables")   
gchp_vs_gchp_tablesdir  = join(gchp_vs_gchp_resultsdir, "Tables")
gchp_vs_gcc_tablesdir   = join(gchp_vs_gcc_resultsdir, "Tables")

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
diff_of_diffs_refstr = [gcc_ref_version, gcc_dev_version]
diff_of_diffs_devstr = [gchp_ref_version, gchp_dev_version]

########################################################################
###    THE REST OF THESE SETTINGS SHOULD NOT NEED TO BE CHANGED      ###
########################################################################

# =====================================================================
# Dates and times
# =====================================================================

# Start and end months of the benchmark
b_start   = (int(bmk_year), int(bmk_mon))
b_stop    = (int(bmk_year), int(bmk_mon) + 1)

# Convert to strings
s_start   = (str(b_start[0]), str(b_start[1]).zfill(2))
s_stop    = (str(b_stop[0]),  str(b_stop[1]).zfill(2))

# Timestamps for files
gcc_date  = np.datetime64( "{}-{}-01T00:00:00".format(s_start[0], s_start[1]))
gchp_date = np.datetime64("{}-{}-16T12:00:00".format(s_start[0], s_start[1]))
end_date  = np.datetime64("{}-{}-01T00:00:00".format(s_stop[0], s_stop[1]))

# Seconds per month
sec_in_bmk_month = (end_date - gcc_date).astype("float64")

# String for month and year (e.g. "Jul2016")
mon_yr_str = calendar.month_abbr[b_start[1]] + s_start[0]

# ======================================================================
# Significant difference filenames
# ======================================================================

vstr = "{}_vs_{}".format(gcc_ref_version, gcc_dev_version)
gcc_vs_gcc_sigdiff = [
    join(gcc_vs_gcc_resultsdir, "{}_sig_diffs_sfc.txt".format(vstr)),
    join(gcc_vs_gcc_resultsdir, "{}_sig_diffs_500hpa.txt".format(vstr)),
    join(gcc_vs_gcc_resultsdir, "{}_sig_diffs_zonalmean.txt".format(vstr)),
    join(gcc_vs_gcc_resultsdir, "{}_sig_diffs_emissions.txt".format(vstr))]

vstr = "{}_vs_{}".format(gcc_dev_version, gchp_dev_version)
gchp_vs_gcc_sigdiff = [
    join(gchp_vs_gcc_resultsdir, "{}_sig_diffs_sfc.txt".format(vstr)),
    join(gchp_vs_gcc_resultsdir, "{}_sig_diffs_500hpa.txt".format(vstr)),
    join(gchp_vs_gcc_resultsdir, "{}_sig_diffs_zonalmean.txt".format(vstr)),
    join(gchp_vs_gcc_resultsdir, "{}_sig_diffs_emissions.txt".format(vstr))]

vstr = "{}_vs_{}".format(gchp_ref_version, gchp_dev_version)
gchp_vs_gchp_sigdiff = [
    join(gchp_vs_gchp_resultsdir, "{}_sig_diffs_sfc.txt".format(vstr)),
    join(gchp_vs_gchp_resultsdir, "{}_sig_diffs_500hpa.txt".format(vstr)),
    join(gchp_vs_gchp_resultsdir, "{}_sig_diffs_zonalmean.txt".format(vstr)),
    join(gchp_vs_gchp_resultsdir, "{}_sig_diffs_emissions.txt").format(vstr)]

# ======================================================================
# Print the list of plots & tables to the screen
# ======================================================================
print("The following plots and tables will be created for {}:".format(bmk_type))
if plot_conc:        print(" - Concentration plots")
if plot_emis:        print(" - Emissions plots")
if plot_jvalues:     print(" - J-values (photolysis rates) plots")
if plot_aod:         print(" - Aerosol optical depth plots")
if ops_budget_table: print(" - Operations budget tables")
if emis_table:       print(" - Table of emissions totals by spc and inventory")
if mass_table:       print(" - Table of species mass")
if OH_metrics:       print(" - Table of OH metrics")
if ste_table:        print(" - Table of strat-trop exchange")
print("Comparisons will be made for the following combinations:")
if gcc_vs_gcc:   print(" - GCC vs GCC")
if gchp_vs_gcc:  print(" - GCHP vs GCC")
if gchp_vs_gchp: print(" - GCHP vs GCHP")
if gchp_vs_gcc_diff_of_diffs: print(" - GCHP vs GCC diff of diffs")

# ======================================================================
# Create GCC vs GCC benchmark plots and tables
# ======================================================================
if gcc_vs_gcc:

    #---------------------------------------------------------------
    # GCC vs GCC Concentration plots
    #
    # Includes lumped species and separates by category if plot_by_spc_cat
    # is true; otherwise excludes lumped species and writes to one file
    #---------------------------------------------------------------
    if plot_conc:
        title = "\n%%% Creating GCC vs. GCC concentration plots %%%"

        # Diagnostic collection files to read
        col = "SpeciesConc"
        ref = get_filepath(gcc_vs_gcc_refdir, col, gcc_date)
        dev = get_filepath(gcc_vs_gcc_devdir, col, gcc_date)

        # Meteorology data needed for calculations
        colmet = "StateMet"
        refmet = get_filepath(gcc_vs_gcc_refdir, colmet, gcc_date)
        devmet = get_filepath(gcc_vs_gcc_devdir, colmet, gcc_date)

        # Make concentration plots
        bmk.make_benchmark_conc_plots(
            ref,
            gcc_vs_gcc_refstr,
            dev,
            gcc_vs_gcc_devstr,
            refmet=refmet,
            devmet=devmet,
            dst=gcc_vs_gcc_resultsdir,
            weightsdir=weightsdir,
            plot_by_spc_cat=plot_by_spc_cat,
            overwrite=True,
            sigdiff_files=gcc_vs_gcc_sigdiff,
            spcdb_dir=spcdb_dir
        )

    #---------------------------------------------------------------
    # GCC vs. GCC emissions plots
    #---------------------------------------------------------------
    if plot_emis:
        print("\n%%% Creating GCC vs. GCC emissions plots %%%")

        # Diagnostic collection files to read
        col = "Emissions"
        ref = get_filepath(gcc_vs_gcc_refdir, col, gcc_date)
        dev = get_filepath(gcc_vs_gcc_devdir, col, gcc_date)

        # Create emissions plots
        bmk.make_benchmark_emis_plots(
            ref,
            gcc_vs_gcc_refstr,
            dev,
            gcc_vs_gcc_devstr,
            dst=gcc_vs_gcc_resultsdir,
            weightsdir=weightsdir,
            plot_by_spc_cat=plot_by_spc_cat,
            plot_by_hco_cat=plot_by_hco_cat,
            overwrite=True,
            sigdiff_files=gcc_vs_gcc_sigdiff,
            spcdb_dir=spcdb_dir
        )

    #---------------------------------------------------------------
    # GCC vs. GCC tables of emission and inventory totals
    #---------------------------------------------------------------
    if emis_table:
        print("\n%%% Creating GCC vs. GCC emissions/inventory tables %%%")

        # Diagnostic collection files to read
        col = "Emissions"
        ref = get_filepath(gcc_vs_gcc_refdir, col, gcc_date)
        dev = get_filepath(gcc_vs_gcc_devdir, col, gcc_date)

        # Print emisisons and inventory tables
        bmk.make_benchmark_emis_tables(
            ref,
            gcc_vs_gcc_refstr,
            dev,
            gcc_vs_gcc_devstr,
            dst=gcc_vs_gcc_resultsdir,
            interval=[sec_in_bmk_month],
            overwrite=True,
            spcdb_dir=spcdb_dir
        )

    # --------------------------------------------------------------
    # GCC vs GCC J-value plots
    # --------------------------------------------------------------
    if plot_jvalues:
        print("\n%%% Creating GCC vs. GCC J-value plots %%%")

        # Diagnostic collection files to read
        col = "JValues"
        ref = get_filepath(gcc_vs_gcc_refdir, col, gcc_date)
        dev = get_filepath(gcc_vs_gcc_devdir, col, gcc_date)

        # Plot J-values
        bmk.make_benchmark_jvalue_plots(
            ref,
            gcc_vs_gcc_refstr,
            dev,
            gcc_vs_gcc_devstr,
            dst=gcc_vs_gcc_resultsdir,
            weightsdir=weightsdir,
            overwrite=True,
            sigdiff_files=gcc_vs_gcc_sigdiff,
            spcdb_dir=spcdb_dir
        )

    #---------------------------------------------------------------
    # GCC vs GCC column AOD plots
    #---------------------------------------------------------------
    if plot_aod:
        print("\n%%% Creating GCC vs. GCC column AOD plots %%%")

        # Diagnostic collection files to read
        col = "Aerosols"
        ref = get_filepath(gcc_vs_gcc_refdir, col, gcc_date)
        dev = get_filepath(gcc_vs_gcc_devdir, col, gcc_date)

        # Plot AODs
        bmk.make_benchmark_aod_plots(
            ref,
            gcc_vs_gcc_refstr,
            dev,
            gcc_vs_gcc_devstr,
            dst=gcc_vs_gcc_resultsdir,
            weightsdir=weightsdir,
            overwrite=True,
            sigdiff_files=gcc_vs_gcc_sigdiff,
            spcdb_dir=spcdb_dir
        )

    #---------------------------------------------------------------
    # GCC vs GCC global mass tables
    #---------------------------------------------------------------
    if mass_table:
        print("\n%%% Creating GCC vs. GCC global mass tables %%%")

        # Diagnostic collection files to read
        col = "Restart"
        ref = get_filepath(gcc_vs_gcc_refrst, col, end_date)
        dev = get_filepath(gcc_vs_gcc_devrst, col, end_date)

        # Plot mass tables
        bmk.make_benchmark_mass_tables(
            ref,
            gcc_ref_version,
            dev,
            gcc_dev_version,
            dst=gcc_vs_gcc_tablesdir,
            overwrite=True,
            spcdb_dir=spcdb_dir
        )

    #---------------------------------------------------------------
    # GCC vs GCC budgets tables
    #---------------------------------------------------------------
    if ops_budget_table:
        print("\n%%% Creating GCC vs. GCC operations budget tables %%%")

        # Diagnostic collection files to read
        col = "Budget"
        ref = get_filepath(gcc_vs_gcc_refdir, col, gcc_date)
        dev = get_filepath(gcc_vs_gcc_devdir, col, gcc_date)

        # Make budget table. Include calculation of Strat and Accumulation
        bmk.make_benchmark_operations_budget(
            gcc_ref_version,
            ref,
            gcc_dev_version,
            dev,
            sec_in_bmk_month,
            benchmark_type=bmk_type,
            label=mon_yr_str,
            dst=gcc_vs_gcc_tablesdir
        )

    #---------------------------------------------------------------
    # GCC vs. GCC global mean OH, MCF Lifetime, CH4 Lifetime
    #---------------------------------------------------------------
    if OH_metrics:
        print("\n%%% Creating GCC vs. GCC OH metrics %%%")

        # Diagnostic collection files to read
        col  = "ConcAfterChem"
        ref = get_filepath(gcc_vs_gcc_refdir, col, gcc_date)
        dev = get_filepath(gcc_vs_gcc_devdir, col, gcc_date)

        # Meteorology data needed for calculations
        col = "StateMet"
        refmet = get_filepath(gcc_vs_gcc_refdir, col, gcc_date)
        devmet = get_filepath(gcc_vs_gcc_devdir, col, gcc_date)

        # Print OH metrics
        bmk.make_benchmark_oh_metrics(
            ref,
            refmet,
            gcc_ref_version,
            dev,
            devmet,
            gcc_dev_version,
            dst=gcc_vs_gcc_tablesdir,
            overwrite=True
        )

    # --------------------------------------------------------------
    # GCC dev Strat-Trop Exchange
    # --------------------------------------------------------------
    if ste_table:
        print("\n%%% Creating GCC dev Strat-Trop Exchange table %%%")

        # Diagnostic collection files to read
        col = "AdvFluxVert"
        dev = get_filepath(gcc_vs_gcc_devdir, col, gcc_date)

        # Compute monthly and annual average strat-trop exchange of O3
        ste.make_benchmark_ste_table(
            gcc_dev_version,
            dev,
            b_start[0],
            bmk_type=bmk_type,
            dst=gcc_vs_gcc_tablesdir,
            species=['O3'],
            overwrite=True,
            month=b_start[1]
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

        # Diagnostic collection files to read
        col = "SpeciesConc"
        ref = get_filepath(gchp_vs_gcc_refdir, col, gcc_date)
        dev = get_filepath(gchp_vs_gcc_devdir, col, gchp_date, is_gchp=True)

        # Meteorology data needed for GCC calculations
        col = "StateMet"
        refmet = get_filepath(gchp_vs_gcc_refdir, col, gcc_date)

        # Meteorology data needed for GCHP calculations
        col = "StateMet_avg"
        devmet = get_filepath(gchp_vs_gcc_devdir, col, gchp_date, is_gchp=True)

        # Make concentration plots
        bmk.make_benchmark_conc_plots(
            ref,
            gchp_vs_gcc_refstr,
            dev,
            gchp_vs_gcc_devstr,
            refmet=refmet,
            devmet=devmet,
            dst=gchp_vs_gcc_resultsdir,
            weightsdir=weightsdir,
            plot_by_spc_cat=plot_by_spc_cat,
            overwrite=True,
            sigdiff_files=gchp_vs_gcc_sigdiff,
            spcdb_dir=spcdb_dir
        )

    #---------------------------------------------------------------
    # GCHP vs. GCC Emissions plots
    #---------------------------------------------------------------
    if plot_emis:
        print("\n%%% Creating GCHP vs. GCC emissions plots %%%")

        # Diagnostic collection files to read
        col = "Emissions"
        ref = get_filepath(gchp_vs_gcc_refdir, col, gcc_date)
        dev = get_filepath(gchp_vs_gcc_devdir, col, gchp_date, is_gchp=True)

        # Create emissions plots
        bmk.make_benchmark_emis_plots(
            ref,
            gchp_vs_gcc_refstr,
            dev,
            gchp_vs_gcc_devstr,
            dst=gchp_vs_gcc_resultsdir,
            weightsdir=weightsdir,
            plot_by_spc_cat=plot_by_spc_cat,
            plot_by_hco_cat=plot_by_hco_cat,
            overwrite=True,
            sigdiff_files=gchp_vs_gcc_sigdiff,
            spcdb_dir=spcdb_dir
        )

    #---------------------------------------------------------------
    # GCHP vs. GCC tables of emission and inventory totals
    #---------------------------------------------------------------
    if emis_table:
        print("\n%%% Creating GCHP vs. GCC emissions/inventory tables %%%")

        # Diagnostic collection files to read
        col = "Emissions"
        ref = get_filepath(gchp_vs_gcc_refdir, col, gcc_date)
        dev = get_filepath(gchp_vs_gcc_devdir, col, gchp_date, is_gchp=True)

        # Meteorology needed for GCHP
        col = "StateMet_avg"
        devmet = get_filepath(gchp_vs_gcc_devdir, col, gchp_date, is_gchp=True)

        # Print emisisons and inventory tables
        bmk.make_benchmark_emis_tables(
            ref,
            gchp_vs_gcc_refstr,
            dev,
            gchp_vs_gcc_devstr,
            dst=gchp_vs_gcc_resultsdir,
            interval=[sec_in_bmk_month],
            overwrite=True,
            devmet=devmet,
            spcdb_dir=spcdb_dir
        )

    #---------------------------------------------------------------
    # GCHP vs. GCC J-values plots
    #---------------------------------------------------------------
    if plot_jvalues:
        print("\n%%% Creating GCHP vs. GCC J-value plots %%%")

        # Diagnostic collection files to read
        col = "JValues"
        ref = get_filepath(gchp_vs_gcc_refdir, col, gcc_date)
        dev = get_filepath(gchp_vs_gcc_devdir, col, gchp_date, is_gchp=True)

        # Plot J-values
        bmk.make_benchmark_jvalue_plots(
            ref,
            gchp_vs_gcc_refstr,
            dev,
            gchp_vs_gcc_devstr,
            dst=gchp_vs_gcc_resultsdir,
            weightsdir=weightsdir,
            overwrite=True,
            sigdiff_files=gchp_vs_gcc_sigdiff,
            spcdb_dir=spcdb_dir
        )

    #---------------------------------------------------------------
    # GCHP vs GCC column AOD plots
    #---------------------------------------------------------------
    if plot_aod:
        print("\n%%% Creating GCHP vs. GCC column AOD plots %%%")

        # Diagnostic collection files to read
        col = "Aerosols"
        ref = get_filepath(gchp_vs_gcc_refdir, col, gcc_date)
        dev = get_filepath(gchp_vs_gcc_devdir, col, gchp_date, is_gchp=True)

        # Plot AODs
        bmk.make_benchmark_aod_plots(
            ref,
            gchp_vs_gcc_refstr,
            dev,
            gchp_vs_gcc_devstr,
            dst=gchp_vs_gcc_resultsdir,
            weightsdir=weightsdir,
            overwrite=True,
            sigdiff_files=gchp_vs_gcc_sigdiff,
            spcdb_dir=spcdb_dir
        )

    #---------------------------------------------------------------
    # GCHP vs GCC global mass tables
    #---------------------------------------------------------------
    if mass_table:
        print("\n%%% Creating GCHP vs. GCC global mass tables %%%")

        # Diagnostic collection files to read
        col = "Restart"
        ref = get_filepath(gchp_vs_gcc_refrst, col, end_date)
        dev = get_filepath(gchp_vs_gcc_devrst, col, end_date, is_gchp=True)

        # Plot mass tables
        bmk.make_benchmark_mass_tables(
            ref,
            gchp_vs_gcc_refstr,
            dev,
            gchp_vs_gcc_devstr,
            dst=gchp_vs_gcc_tablesdir,
            overwrite=True,
            spcdb_dir=spcdb_dir
        )

    #---------------------------------------------------------------
    # GCHP vs GCC operations budgets tables
    #---------------------------------------------------------------
    if ops_budget_table:
        print("\n%%% Creating GCHP vs GCC operations budget tables %%%")

        # Diagnostic collections to read
        col = "Budget"
        ref = get_filepath(gchp_vs_gcc_refdir, col, gcc_date)
        dev = get_filepath(gchp_vs_gcc_devdir, col, gchp_date, is_gchp=True)

        # Make budget table. Include calculation of Strat. Exclude Transport
        # and Accumulation.
        bmk.make_benchmark_operations_budget(
            gcc_ref_version,
            ref,
            gchp_dev_version,
            dev,
            sec_in_bmk_month,
            benchmark_type=bmk_type,
            label=mon_yr_str,
            operations=["Chemistry","Convection","EmisDryDep","Mixing",
                        "WetDep"],
            compute_accum=False,
            dst=gchp_vs_gcc_tablesdir
        )

    #---------------------------------------------------------------
    # GCHP vs. GCC global mean OH, MCF Lifetime, CH4 Lifetime
    #---------------------------------------------------------------
    if OH_metrics:
        print("\n%%% Skipping GCHP vs GCC OH metrics %%%")

    # --------------------------------------------------------------
    # GCHP vs GCC Strat-Trop Exchange
    # --------------------------------------------------------------
    if ste_table:
        title = "\n%%% Skipping GCHP vs. GCC Strat-Trop Exchange table %%%"
        print(title)

# ======================================================================
# Create GCHP vs GCHP benchmark plots and tables
# ======================================================================
if gchp_vs_gchp:

    #---------------------------------------------------------------
    # GCHP vs GCHP Concentration plots
    #---------------------------------------------------------------
    if plot_conc:
        print("\n%%% Creating GCHP vs. GCHP concentration plots %%%")

        # Diagnostic collection files to read
        col = "SpeciesConc"
        ref = get_filepath(gchp_vs_gchp_refdir, col, gchp_date, is_gchp=True)
        dev = get_filepath(gchp_vs_gchp_devdir, col, gchp_date, is_gchp=True)

        # Meteorology data needed for GCHP calculations
        col = "StateMet_avg"
        refmet = get_filepath(gchp_vs_gchp_devdir, col, gchp_date, is_gchp=True)
        devmet = get_filepath(gchp_vs_gchp_devdir, col, gchp_date, is_gchp=True)

        # Make concentration plots
        bmk.make_benchmark_conc_plots(
            ref,
            gchp_vs_gchp_refstr,
            dev,
            gchp_vs_gchp_devstr,
            refmet=refmet,
            devmet=devmet,
            dst=gchp_vs_gchp_resultsdir,
            weightsdir=weightsdir,
            plot_by_spc_cat=plot_by_spc_cat,
            overwrite=True,
            sigdiff_files=gchp_vs_gchp_sigdiff,
            spcdb_dir=spcdb_dir
        )

    #---------------------------------------------------------------
    # GCHP vs. GCHP Emissions plots
    #---------------------------------------------------------------
    if plot_emis:
        print("\n%%% Creating GCHP vs. GCHP emissions plots %%%")

        # Diagnostic collection files to read
        col = "Emissions"
        ref = get_filepath(gchp_vs_gchp_refdir, col, gchp_date, is_gchp=True)
        dev = get_filepath(gchp_vs_gchp_devdir, col, gchp_date, is_gchp=True)

        # Create emissions plots
        bmk.make_benchmark_emis_plots(
            ref,
            gchp_vs_gchp_refstr,
            dev,
            gchp_vs_gchp_devstr,
            dst=gchp_vs_gchp_resultsdir,
            weightsdir=weightsdir,
            plot_by_spc_cat=plot_by_spc_cat,
            plot_by_hco_cat=plot_by_hco_cat,
            overwrite=True,
            sigdiff_files=gchp_vs_gchp_sigdiff,
            spcdb_dir=spcdb_dir
        )

    #---------------------------------------------------------------
    # GCHP vs. GCHP tables of emission and inventory totals
    #---------------------------------------------------------------
    if emis_table:
        print("\n%%% Creating GCHP vs. GCHP emissions/inventory tables %%%")

        # Diagnostic collection files to read
        col = "Emissions"
        ref = get_filepath(gchp_vs_gchp_refdir,col, gchp_date, is_gchp=True)
        dev = get_filepath(gchp_vs_gchp_devdir, col, gchp_date, is_gchp=True)

        # Meteorology needed for GCHP
        col = "StateMet_avg"
        refmet = get_filepath(gchp_vs_gchp_refdir,col, gchp_date, is_gchp=True)
        devmet = get_filepath(gchp_vs_gchp_devdir, col, gchp_date, is_gchp=True)

        # Print emisisons and inventory tables
        bmk.make_benchmark_emis_tables(
            ref,
            gchp_vs_gchp_refstr,
            dev,
            gchp_vs_gchp_devstr,
            dst=gchp_vs_gchp_resultsdir,
            interval=[sec_in_bmk_month],
            overwrite=True,
            refmet=refmet,
            devmet=devmet,
            spcdb_dir=spcdb_dir
        )

    #---------------------------------------------------------------
    # GCHP vs. GCHP J-values plots
    #---------------------------------------------------------------
    if plot_jvalues:
        print("\n%%% Creating GCHP vs. GCHP J-value plots %%%")

        # Diagnostic collection files to read
        col = "JValues"
        ref = get_filepath(gchp_vs_gchp_refdir, col, gchp_date, is_gchp=True)
        dev = get_filepath(gchp_vs_gchp_devdir, col, gchp_date, is_gchp=True)

        # Plot J-values
        bmk.make_benchmark_jvalue_plots(
            ref,
            gchp_vs_gchp_refstr,
            dev,
            gchp_vs_gchp_devstr,
            dst=gchp_vs_gchp_resultsdir,
            weightsdir=weightsdir,
            overwrite=True,
            sigdiff_files=gchp_vs_gchp_sigdiff,
            spcdb_dir=spcdb_dir
        )

    #---------------------------------------------------------------
    # GCHP vs GCHP column AOD plots
    #---------------------------------------------------------------
    if plot_aod:
        print("\n%%% Creating GCHP vs. GCHP column AOD plots %%%")

        # Diagnostic collection files to read
        col = "Aerosols"
        ref = get_filepath(gchp_vs_gchp_refdir, col, gchp_date, is_gchp=True)
        dev = get_filepath(gchp_vs_gchp_devdir, col, gchp_date, is_gchp=True)

        # Plot AODs
        bmk.make_benchmark_aod_plots(
            ref,
            gchp_vs_gchp_refstr,
            dev,
            gchp_vs_gchp_devstr,
            dst=gchp_vs_gchp_resultsdir,
            weightsdir=weightsdir,
            overwrite=True,
            sigdiff_files=gchp_vs_gchp_sigdiff,
            spcdb_dir=spcdb_dir
        )

    #---------------------------------------------------------------
    # GCHP vs GCHP global mass tables
    #---------------------------------------------------------------
    if mass_table:
        print("\n%%% Creating GCHP vs. GCHP global mass tables %%%")

        # Diagnostic collection files to read
        col = "Restart"
        ref = get_filepath(gchp_vs_gchp_refrst, col, end_date, is_gchp=True)
        dev = get_filepath(gchp_vs_gchp_devrst, col, end_date, is_gchp=True)

        # Plot mass tables
        bmk.make_benchmark_mass_tables(
            ref,
            gchp_vs_gchp_refstr,
            dev,
            gchp_vs_gchp_devstr,
            dst=gchp_vs_gchp_tablesdir,
            overwrite=True,
            spcdb_dir=spcdb_dir
        )

    #---------------------------------------------------------------
    # GCHP vs GCHP operations budgets tables
    #---------------------------------------------------------------
    if ops_budget_table:
        print("\n%%% Creating GCHP vs GCHP operations budget tables %%%")

        # Diagnostic collections to read
        col = "Budget"
        ref = get_filepath(gchp_vs_gchp_refdir, col, gchp_date, is_gchp=True)
        dev = get_filepath(gchp_vs_gchp_devdir, col, gchp_date, is_gchp=True)

        # Make budget table. Include calculation of Strat. Exclude Transport
        # and Accumulation.
        bmk.make_benchmark_operations_budget(
            gchp_ref_version,
            ref,
            gchp_dev_version,
            dev,
            sec_in_bmk_month,
            benchmark_type=bmk_type,
            label=mon_yr_str,
            operations=["Chemistry","Convection","EmisDryDep","Mixing",
                        "WetDep"],
            compute_accum=False,
            dst=gchp_vs_gchp_tablesdir
        )


    #---------------------------------------------------------------
    # GCHP vs. GCHP global mean OH, MCF Lifetime, CH4 Lifetime
    #---------------------------------------------------------------
    if OH_metrics:
        print("\n%%% Skipping GCHP vs GCHP OH metrics %%%")

    # --------------------------------------------------------------
    # GCHP vs GCHP Strat-Trop Exchange
    # --------------------------------------------------------------
    if ste_table:
        print("\n%%% Skipping GCHP vs. GCHP Strat-Trop Exchange table %%%")

# =====================================================================
# Create GCHP vs GCC difference of differences benchmark plots
# =====================================================================
if gchp_vs_gcc_diff_of_diffs:

    if plot_conc:
        print("\n%%% Creating GCHP vs. GCC diff-of-diffs conc plots %%%")

        # Diagnostic collection files to read
        col = "SpeciesConc"
        gcc_ref = get_filepath(gcc_vs_gcc_refdir, col, gcc_date)
        gcc_dev = get_filepath(gcc_vs_gcc_devdir, col, gcc_date)
        gchp_ref = get_filepath(gchp_vs_gchp_refdir, col, gchp_date,
                                is_gchp=True)
        gchp_dev = get_filepath(gchp_vs_gchp_devdir, col, gchp_date,
                                is_gchp=True)

        # Create diff-of-diff plots for species concentrations
        # NOTE: for simplicity, do not convert aerosols to ug/m3.
        bmk.make_benchmark_conc_plots(
            gcc_ref,
            diff_of_diffs_refstr,
            gchp_ref,
            diff_of_diffs_devstr,
            dst=diff_of_diffs_resultsdir,
            weightsdir=weightsdir,
            overwrite=True,
            use_cmap_RdBu=True,
            second_ref=gcc_dev,
            second_dev=gchp_dev,
            cats_in_ugm3=None,
            spcdb_dir=spcdb_dir
        )


