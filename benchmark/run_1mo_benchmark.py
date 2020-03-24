#!/usr/bin/env python
"""
run_1mo_benchmark.py: Driver script for creating benchmark plots.

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

import calendar
import os
from os.path import join
from gcpy import benchmark as bmk
from gcpy.core import get_filepaths
from gcpy.constants import skip_these_vars
import gcpy.budget_ops as opbdg
import gcpy.ste_flux as ste
import numpy as np
import xarray as xr
import warnings

# Tell matplotlib not to look for an X-window
os.environ["QT_QPA_PLATFORM"]="offscreen"

# Suppress harmless run-time warnings (mostly about underflow in division)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)

########################################################################
###           CONFIGURABLE SETTINGS: EDIT THESE ACCORDINGLY          ###
########################################################################

# =====================================================================
# Benchmark information (**EDIT AS NEEDED**)
# Note: When doing GCHP vs GCC comparisions gchp_dev will be compared
# to gcc_dev (not gcc_ref!). This ensures consistency in version names
# when doing GCHP vs GCC diff-of-diffs (mps, 6/27/19)
# =====================================================================

maindir  = "/path/to/main/benchmark/directory"
gcc_ref_version = "gcc_ref_version"
gcc_dev_version = "gcc_dev_version"
gchp_ref_version = "gchp_ref_version"
gchp_dev_version = "gchp_dev_version"

# Path to regridding weights
weightsdir = "/n/holylfs/EXTERNAL_REPOS/GEOS-CHEM/gcgrid/gcdata/ExtData/GCHP/RegriddingWeights"

# =====================================================================
# Comparisons to run (**EDIT AS NEEDED**)
# =====================================================================
bmk_type     = "FullChemBenchmark"
gcc_vs_gcc   = True
gchp_vs_gcc  = False
gchp_vs_gchp = False
gchp_vs_gcc_diff_of_diffs = False

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
OH_metrics   = True
ste_table    = True

# Plot concentrations and emissions by category?
plot_by_spc_cat = True
plot_by_hco_cat = True

# =====================================================================
# Data directories (**EDIT AS NEEDED**)
#
# For gchp_vs_gcc_refdir use gcc_dev_version, not ref (mps, 6/27/19)
# =====================================================================

# Diagnostics
gcc_vs_gcc_refdir   = join(maindir, gcc_ref_version,  "OutputDir")
gcc_vs_gcc_devdir   = join(maindir, gcc_dev_version,  "OutputDir")
gchp_vs_gcc_refdir  = join(maindir, gcc_dev_version,  "OutputDir")
gchp_vs_gcc_devdir  = join(maindir, gchp_dev_version, "OutputDir")
gchp_vs_gchp_refdir = join(maindir, gchp_ref_version, "OutputDir")
gchp_vs_gchp_devdir = join(maindir, gchp_dev_version, "OutputDir")

# Restart files
gcc_vs_gcc_refrst   = join(maindir, gcc_ref_version )
gcc_vs_gcc_devrst   = join(maindir, gcc_dev_version )
gchp_vs_gcc_refrst  = join(maindir, gcc_dev_version )
gchp_vs_gcc_devrst  = join(maindir, gchp_dev_version)
gchp_vs_gchp_refrst = join(maindir, gchp_ref_version)
gchp_vs_gchp_devrst = join(maindir, gchp_dev_version)

# Plots directories (edit as needed)
gcc_vs_gcc_plotsdir    = join(maindir, gcc_dev_version, "Plots")
gchp_vs_gchp_plotsdir  = join(maindir, gchp_dev_version,
                              "Plots", "GCHP_version_comparison")
gchp_vs_gcc_plotsdir   = join(maindir, gchp_dev_version,
                              "Plots", "GCHP_GCC_comparison")
diff_of_diffs_plotsdir = join(maindir, gchp_dev_version,
                              "Plots", "GCHP_GCC_diff_of_diffs")

# Plot title strings (edit as needed)
# For gchp_vs_gcc_refstr use gcc_dev_version, not ref (mps, 6/27/19)
gcc_vs_gcc_refstr    = "{}".format(gcc_ref_version)
gcc_vs_gcc_devstr    = "{}".format(gcc_dev_version)
gchp_vs_gcc_refstr   = "{}".format(gcc_dev_version)
gchp_vs_gcc_devstr   = "{}".format(gchp_dev_version)
gchp_vs_gchp_refstr  = "{}".format(gchp_ref_version)
gchp_vs_gchp_devstr  = "{}".format(gchp_dev_version)
diff_of_diffs_refstr = "{} - {}".format(gcc_dev_version,
                                        gcc_ref_version)
diff_of_diffs_devstr = "{} - {}".format(gchp_dev_version,
                                        gchp_ref_version)

# =====================================================================
# Starting and ending dates (**EDIT AS NEEDED**)
# =====================================================================

# Start and end months of the benchmark
b_start   = (2016, 7)
b_stop    = (2016, 8)

# Convert to strings
s_start   = (str(b_start[0]), str(b_start[1]).zfill(2))
s_stop    = (str(b_stop[0]),  str(b_stop[1]).zfill(2))

# Timestamps for files
gcc_date  = np.datetime64( "{}-{}-01T00:00:00".format(s_start[0], s_start[1]))
gchp_date = np.datetime64("{}-{}-16T12:00:00".format(s_start[0], s_start[1]))
end_date  = np.datetime64("{}-{}-01T00:00:00".format(s_stop[0], s_stop[1]))

# Seconds per month
s_per_mon = [(end_date - gcc_date).astype("float64")]

# String for month and year (e.g. "Jul2016")
mon_yr_str = calendar.month_abbr[b_start[1]] + s_start[0]

########################################################################
###    THE REST OF THESE SETTINGS SHOULD NOT NEED TO BE CHANGED      ###
########################################################################

# ======================================================================
# Files that will contain lists of quantities that have significant
# differences -- we need these for the benchmark approval forms.
# ======================================================================

vstr = "{}_vs_{}".format(gcc_vs_gcc_refstr, gcc_vs_gcc_devstr)
gcc_vs_gcc_sigdiff = [
    join(gcc_vs_gcc_plotsdir, "{}_sig_diffs_sfc.txt".format(vstr)),
    join(gcc_vs_gcc_plotsdir, "{}_sig_diffs_500hpa.txt".format(vstr)),
    join(gcc_vs_gcc_plotsdir, "{}_sig_diffs_zonalmean.txt".format(vstr)),
    join(gcc_vs_gcc_plotsdir, "{}_sig_diffs_emissions.txt".format(vstr))]

vstr = "{}_vs_{}".format(gchp_vs_gcc_refstr, gchp_vs_gcc_devstr)
gchp_vs_gcc_sigdiff = [
    join(gchp_vs_gcc_plotsdir, "{}_sig_diffs_sfc.txt".format(vstr)),
    join(gchp_vs_gcc_plotsdir, "{}_sig_diffs_500hpa.txt".format(vstr)),
    join(gchp_vs_gcc_plotsdir, "{}_sig_diffs_zonalmean.txt".format(vstr)),
    join(gchp_vs_gcc_plotsdir, "{}_sig_diffs_emissions.txt".format(vstr))]

vstr = "{}_vs_{}".format(gchp_vs_gchp_refstr, gchp_vs_gchp_devstr)
gchp_vs_gchp_sigdiff = [
    join(gchp_vs_gchp_plotsdir, "{}_sig_diffs_sfc.txt".format(vstr)),
    join(gchp_vs_gchp_plotsdir, "{}_sig_diffs_500hpa.txt".format(vstr)),
    join(gchp_vs_gchp_plotsdir, "{}_sig_diffs_zonalmean.txt".format(vstr)),
    join(gchp_vs_gchp_plotsdir, "{}_sig_diffs_emissions.txt").format(vstr)]

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

# ======================================================================
# Create GCC vs GCC benchmark plots and tables
# ======================================================================
if gcc_vs_gcc:

    if plot_conc:
        #---------------------------------------------------------------
        # GCC vs GCC Concentration plots
        #
        # Includes lumped species and separates by category if
        # plot_by_spc_cat is true; otherwise excludes lumped species
        # and writes to one file
        #---------------------------------------------------------------
        title = "\n%%% Creating GCC vs. GCC {} concentration plots %%%"
        print(title.format(bmk_type))

        # Files to read
        collection = "SpeciesConc"
        gcc_vs_gcc_reflist = get_filepaths(gcc_vs_gcc_refdir, collection,
                                           [gcc_date], is_gcc=True)
        gcc_vs_gcc_devlist = get_filepaths(gcc_vs_gcc_devdir, collection,
                                           [gcc_date], is_gcc=True)

        # Make concentration plots
        bmk.make_benchmark_plots(gcc_vs_gcc_reflist[0],
                                 gcc_vs_gcc_refstr,
                                 gcc_vs_gcc_devlist[0],
                                 gcc_vs_gcc_devstr,
                                 dst=gcc_vs_gcc_plotsdir,
                                 weightsdir=weightsdir,
                                 plot_by_spc_cat=plot_by_spc_cat,
                                 overwrite=True,
                                 sigdiff_files=gcc_vs_gcc_sigdiff)

    if plot_emis:
        #---------------------------------------------------------------
        # GCC vs. GCC Emissions plots
        #---------------------------------------------------------------
        title = "\n%%% Creating GCC vs. GCC {} emissions plots %%%"
        print(title.format(bmk_type))

        # Files to read
        collection = "Emissions"
        gcc_vs_gcc_reflist = get_filepaths(gcc_vs_gcc_refdir, collection,
                                           [gcc_date], is_gcc=True)
        gcc_vs_gcc_devlist = get_filepaths(gcc_vs_gcc_devdir, collection,
                                           [gcc_date], is_gcc=True)

        # Create emissions plots
        bmk.make_benchmark_emis_plots(gcc_vs_gcc_reflist[0],
                                      gcc_vs_gcc_refstr,
                                      gcc_vs_gcc_devlist[0],
                                      gcc_vs_gcc_devstr,
                                      dst=gcc_vs_gcc_plotsdir,
                                      weightsdir=weightsdir,
                                      plot_by_spc_cat=plot_by_spc_cat,
                                      plot_by_hco_cat=plot_by_hco_cat,
                                      overwrite=True,
                                      sigdiff_files=gcc_vs_gcc_sigdiff)

    if emis_table:
        #---------------------------------------------------------------
        # GCC vs. GCC tables of emission and inventory totals
        #---------------------------------------------------------------
        title = "\n%%% Creating GCC vs. GCC {} emissions/inventory tables %%%"
        print(title.format(bmk_type))

        # Files to read
        collection = "Emissions"
        gcc_vs_gcc_reflist = get_filepaths(gcc_vs_gcc_refdir, collection,
                                           [gcc_date], is_gcc=True)
        gcc_vs_gcc_devlist = get_filepaths(gcc_vs_gcc_devdir, collection,
                                           [gcc_date], is_gcc=True)

        # Print emisisons and inventory tables
        bmk.make_benchmark_emis_tables(gcc_vs_gcc_reflist,
                                       gcc_vs_gcc_refstr,
                                       gcc_vs_gcc_devlist,
                                       gcc_vs_gcc_devstr,
                                       dst=gcc_vs_gcc_plotsdir,
                                       interval=s_per_mon,
                                       overwrite=True)

    if plot_jvalues:
        #---------------------------------------------------------------
        # GCC vs. GCC J-values plots
        #---------------------------------------------------------------
        title = "\n%%% Creating GCC vs. GCC {} J-value plots %%%"
        print(title.format(bmk_type))

        # Files to read
        collection = "JValues"
        gcc_vs_gcc_reflist = get_filepaths(gcc_vs_gcc_refdir, collection,
                                           [gcc_date], is_gcc=True)
        gcc_vs_gcc_devlist = get_filepaths(gcc_vs_gcc_devdir, collection,
                                           [gcc_date], is_gcc=True)
        
        # Plot J-values
        bmk.make_benchmark_jvalue_plots(gcc_vs_gcc_reflist[0],
                                        gcc_vs_gcc_refstr,
                                        gcc_vs_gcc_devlist[0],
                                        gcc_vs_gcc_devstr,
                                        dst=gcc_vs_gcc_plotsdir,
                                        weightsdir=weightsdir,
                                        overwrite=True,
                                        sigdiff_files=gcc_vs_gcc_sigdiff)

    if plot_aod:
        #---------------------------------------------------------------
        # GCC vs GCC column AOD plots
        #---------------------------------------------------------------
        title = "\n%%% Creating GCC vs. GCC {} column AOD plots %%%"
        print(title.format(bmk_type))

        # Files to read
        collection = "Aerosols"
        gcc_vs_gcc_reflist = get_filepaths(gcc_vs_gcc_refdir, collection,
                                           [gcc_date], is_gcc=True)
        gcc_vs_gcc_devlist = get_filepaths(gcc_vs_gcc_devdir, collection,
                                           [gcc_date], is_gcc=True)

        # Plot AODs
        bmk.make_benchmark_aod_plots(gcc_vs_gcc_reflist[0],
                                     gcc_vs_gcc_refstr,
                                     gcc_vs_gcc_devlist[0],
                                     gcc_vs_gcc_devstr,
                                     dst=gcc_vs_gcc_plotsdir,
                                     weightsdir=weightsdir,
                                     overwrite=True,
                                     sigdiff_files=gcc_vs_gcc_sigdiff)

    if mass_table:
        #---------------------------------------------------------------
        # GCC vs GCC global mass tables
        #---------------------------------------------------------------
        title = "\n%%% Creating GCC vs. GCC {} global mass tables %%%"
        print(title.format(bmk_type))

        # Files to read
        collection = "Restart"
        gcc_vs_gcc_reflist = get_filepaths(gcc_vs_gcc_refrst, collection,
                                           [end_date], is_gcc=True)
        gcc_vs_gcc_devlist = get_filepaths(gcc_vs_gcc_devrst, collection,
                                           [end_date], is_gcc=True)

        # Save to "Tables" subfolder of plot dir
        plot_dir = join(gcc_vs_gcc_plotsdir, "Tables")
        
        # Plot mass tables
        bmk.make_benchmark_mass_tables(gcc_vs_gcc_reflist[0],
                                       gcc_vs_gcc_refstr,
                                       gcc_vs_gcc_devlist[0],
                                       gcc_vs_gcc_devstr,
                                       dst=plot_dir,
                                       overwrite=True)
        
    if budget_table:
        #---------------------------------------------------------------
        # GCC vs GCC budgets tables 
        #---------------------------------------------------------------
        title = "\n%%% Creating GCC vs. GCC {} budget tables %%%"
        print(title.format(bmk_type))

        # Files to read
        collection = "Budget"
        gcc_vs_gcc_devlist = get_filepaths(gcc_vs_gcc_devdir, collection,
                                           [gcc_date], is_gcc=True)

        # Save to "Tables" subfolder of plot dir
        plot_dir = join(gcc_vs_gcc_plotsdir, "Budget")
        
        # Print budget tables
        bmk.make_benchmark_budget_tables(gcc_vs_gcc_devlist[0],
                                         gcc_vs_gcc_devstr,
                                         dst=gcc_vs_gcc_plotsdir,
                                         interval=s_per_mon,
                                         overwrite=True)

        # Print "operations budget" for strat, trop, strat+trop 
        opbdg.make_operations_budget_table(gcc_dev_version,
                                           gcc_vs_gcc_devlist[0],
                                           bmk_type,
                                           dst=plot_dir,
                                           label=mon_yr_str,
                                           interval=s_per_mon[0],
                                           overwrite=True)

        
    if OH_metrics:
        #---------------------------------------------------------------
        # GCC vs. GCC global mean OH, MCF Lifetime, CH4 Lifetime
        #---------------------------------------------------------------
        title = "\n%%% Creating GCC vs. GCC {} OH metrics %%%"
        print(title.format(bmk_type))
        
        # Files to read
        collection = ["ConcAfterChem", "StateMet"]
        gcc_vs_gcc_reflist = get_filepaths(gcc_vs_gcc_refdir, collection,
                                           [gcc_date], is_gcc=True)
        gcc_vs_gcc_devlist = get_filepaths(gcc_vs_gcc_devdir, collection,
                                           [gcc_date], is_gcc=True)

        # Save to "Tables" subfolder of plot dir
        plot_dir = join(gcc_vs_gcc_plotsdir, "Tables")
        
        # Print OH metrics
        bmk.make_benchmark_oh_metrics(gcc_vs_gcc_reflist,
                                      gcc_vs_gcc_refstr,
                                      gcc_vs_gcc_devlist,
                                      gcc_vs_gcc_devstr,
                                      dst=plot_dir,
                                      overwrite=True)
        
    if ste_table:
        # --------------------------------------------------------------
        # GCC vs GCC Strat-Trop Exchange
        # --------------------------------------------------------------
        title = "\n%%% Creating GCC vs. GCC {} Strat-Trop Exchange table %%%"
        print(title.format(bmk_type))

        # Files to read
        # Compute monthly and annual average strat-trop exchange of O3
        gcc_vs_gcc_devlist = get_filepaths(gcc_vs_gcc_devdir, "AdvFluxVert",
                                           [gcc_date], is_gcc=True)
        
        # Save to "Tables" subfolder of plot dir
        plot_dir = join(gcc_vs_gcc_plotsdir, "Tables")

        # Create STE table
        ste.make_benchmark_ste_table(gcc_dev_version,
                                     gcc_vs_gcc_devlist[0],
                                     b_start[0],
                                     bmk_type=bmk_type,
                                     dst=plot_dir,
                                     species=['O3'],
                                     overwrite=True,
                                     month=b_start[1])

# ======================================================================
# Create GCHP vs GCC benchmark plots and tables
# ======================================================================
if gchp_vs_gcc:

    if plot_conc:
        #---------------------------------------------------------------
        # GCHP vs GCC Concentration plots
        #
        # Includes lumped species and separates by category if
        # plot_by_spc_cat is true; otherwise excludes lumped species
        # and writes to one file
        #---------------------------------------------------------------
        title = "\n%%% Creating GCHP vs. GCC {} concentration plots %%%"
        print(title.format(bmk_type))

        # Files to read
        collection = "SpeciesConc"
        gchp_vs_gcc_reflist = get_filepaths(gchp_vs_gcc_refdir, collection,
                                            [gcc_date], is_gcc=True)
        gchp_vs_gcc_devlist = get_filepaths(gchp_vs_gcc_devdir, collection,
                                            [gchp_date], is_gchp=True)

        # Make concentration plots
        bmk.make_benchmark_plots(gchp_vs_gcc_reflist[0],
                                 gchp_vs_gcc_refstr,
                                 gchp_vs_gcc_devlist[0],
                                 gchp_vs_gcc_devstr,
                                 dst=gchp_vs_gcc_plotsdir,
                                 weightsdir=weightsdir,
                                 plot_by_spc_cat=plot_by_spc_cat,
                                 overwrite=True,
                                 sigdiff_files=gchp_vs_gcc_sigdiff)

    if plot_emis:
        #---------------------------------------------------------------
        # GCHP vs. GCC Emissions plots
        #---------------------------------------------------------------
        title = "\n%%% Creating GCHP vs. GCC {} emissions plots %%%"
        print(title.format(bmk_type))

        # Files to read
        collection = "Emissions"
        gchp_vs_gcc_reflist = get_filepaths(gchp_vs_gcc_refdir, collection,
                                            [gcc_date], is_gcc=True)
        gchp_vs_gcc_devlist = get_filepaths(gchp_vs_gcc_devdir, collection,
                                            [gchp_date], is_gchp=True)

        # Create emissions plots
        bmk.make_benchmark_emis_plots(gchp_vs_gcc_reflist[0],
                                      gchp_vs_gcc_refstr,
                                      gchp_vs_gcc_devlist[0],
                                      gchp_vs_gcc_devstr,
                                      dst=gchp_vs_gcc_plotsdir,
                                      weightsdir=weightsdir,
                                      plot_by_spc_cat=plot_by_spc_cat,
                                      plot_by_hco_cat=plot_by_hco_cat,
                                      overwrite=True,
                                      sigdiff_files=gchp_vs_gcc_sigdiff)

    if emis_table:
        #---------------------------------------------------------------
        # GCHP vs. GCC tables of emission and inventory totals
        #---------------------------------------------------------------
        title = "\n%%% Creating GCHP vs. GCC {} emissions/inventory tables %%%"
        print(title.format(bmk_type))

        # Files to read
        collection = "Emissions"
        gchp_vs_gcc_reflist = get_filepaths(gchp_vs_gcc_refdir, collection,
                                            [gcc_date], is_gcc=True)
        collection = ["Emissions", "StateMet_avg"]
        gchp_vs_gcc_devlist = get_filepaths(gchp_vs_gcc_devdir, collection,
                                            [gchp_date], is_gchp=True)

        # Print emisisons and inventory tables
        bmk.make_benchmark_emis_tables(gchp_vs_gcc_reflist,
                                       gchp_vs_gcc_refstr,
                                       gchp_vs_gcc_devlist,
                                       gchp_vs_gcc_devstr,
                                       dst=gchp_vs_gcc_plotsdir,
                                       interval=s_per_mon,
                                       overwrite=True)

    if plot_jvalues:
        #---------------------------------------------------------------
        # GCHP vs. GCC J-values plots
        #---------------------------------------------------------------
        title = "\n%%% Creating GCHP vs. GCC {} J-value plots %%%"
        print(title.format(bmk_type))

        # Files to read
        collection = "JValues"
        gchp_vs_gcc_reflist = get_filepaths(gchp_vs_gcc_refdir, collection,
                                            [gcc_date], is_gcc=True)
        gchp_vs_gcc_devlist = get_filepaths(gchp_vs_gcc_devdir, collection,
                                            [gchp_date], is_gchp=True)
        
        # Plot J-values
        bmk.make_benchmark_jvalue_plots(gchp_vs_gcc_reflist[0],
                                        gchp_vs_gcc_refstr,
                                        gchp_vs_gcc_devlist[0],
                                        gchp_vs_gcc_devstr,
                                        dst=gchp_vs_gcc_plotsdir,
                                        weightsdir=weightsdir,
                                        overwrite=True,
                                        sigdiff_files=gchp_vs_gcc_sigdiff)

    if plot_aod:
        #---------------------------------------------------------------
        # GCHP vs GCC column AOD plots
        #---------------------------------------------------------------
        title = "\n%%% Creating GCHP vs. GCC {} column AOD plots %%%"
        print(title.format(bmk_type))

        # Files to read
        collection = "Aerosols"
        gchp_vs_gcc_reflist = get_filepaths(gchp_vs_gcc_refdir, collection,
                                            [gcc_date], is_gcc=True)
        gchp_vs_gcc_devlist = get_filepaths(gchp_vs_gcc_devdir, collection,
                                            [gchp_date], is_gchp=True)

        # Plot AODs
        bmk.make_benchmark_aod_plots(gchp_vs_gcc_reflist[0],
                                     gchp_vs_gcc_refstr,
                                     gchp_vs_gcc_devlist[0],
                                     gchp_vs_gcc_devstr,
                                     dst=gchp_vs_gcc_plotsdir,
                                     weightsdir=weightsdir,
                                     overwrite=True,
                                     sigdiff_files=gchp_vs_gcc_sigdiff)

    if mass_table:
        #---------------------------------------------------------------
        # GCHP vs GCC global mass tables
        #---------------------------------------------------------------
        title = "\n%%% Creating GCHP vs. GCC {} global mass tables %%%"
        print(title.format(bmk_type))

        # Files to read
        collection = "Restart"
        gchp_vs_gcc_reflist = get_filepaths(gchp_vs_gcc_refrst, collection,
                                            [end_date], is_gcc=True)
        collection = ["Restart", "StateMet_inst"]
        gchp_vs_gcc_devlist = get_filepaths(gchp_vs_gcc_devrst, collection,
                                            [end_date], is_gchp=True)

        # Save to "Tables" subfolder of plot dir
        plot_dir = join(gchp_vs_gcc_plotsdir, "Tables")
        
        # Plot mass tables
        bmk.make_benchmark_mass_tables(gchp_vs_gcc_reflist,
                                       gchp_vs_gcc_refstr,
                                       gchp_vs_gcc_devlist,
                                       gchp_vs_gcc_devstr,
                                       dst=plot_dir,
                                       overwrite=True)
        
    if budget_table:
        #---------------------------------------------------------------
        # GCHP vs GCC budgets tables 
        #---------------------------------------------------------------
        title = "\n%%% Creating GCHP vs. GCC {} budget tables %%%"
        print(title.format(bmk_type))

        # Files to read
        collection = "Budget"
        gchp_vs_gcc_devlist = get_filepaths(gchp_vs_gcc_devdir, collection,
                                            [gchp_date], is_gchp=True)

        # Save to "Tables" subfolder of plot dir
        plot_dir = join(gchp_vs_gcc_plotsdir, "Budget")
        
        # Print budget tables
        bmk.make_benchmark_budget_tables(gchp_vs_gcc_devlist[0],
                                         gchp_vs_gcc_devstr,
                                         dst=gchp_vs_gcc_plotsdir,
                                         interval=s_per_mon,
                                         overwrite=True)

        # Print "operations budget" for strat, trop, strat+trop 
        opbdg.make_operations_budget_table(gcc_dev_version,
                                           gchp_vs_gcc_devlist[0],
                                           bmk_type,
                                           dst=plot_dir,
                                           label=mon_yr_str,
                                           interval=s_per_mon[0],
                                           overwrite=True)

        
    if OH_metrics:
        #---------------------------------------------------------------
        # GCHP vs. GCC global mean OH, MCF Lifetime, CH4 Lifetime
        #---------------------------------------------------------------
        title = "\n%%% Creating GCHP vs. GCC {} OH metrics %%%"
        print(title.format(bmk_type))
        
        # Files to read
        collection = ["ConcAfterChem", "StateMet"]
        gchp_vs_gcc_reflist = get_filepaths(gchp_vs_gcc_refdir, collection,
                                            [gcc_date], is_gcc=True)
        collection = ["ConcAfterChem", "StateMet_avg"]
        gchp_vs_gcc_devlist = get_filepaths(gchp_vs_gcc_devdir, collection,
                                            [gchp_date], is_gchp=True)

        # Save to "Tables" subfolder of plot dir
        plot_dir = join(gchp_vs_gcc_plotsdir, "Tables")
        
        # Print OH metrics
        bmk.make_benchmark_oh_metrics(gchp_vs_gcc_reflist,
                                      gchp_vs_gcc_refstr,
                                      gchp_vs_gcc_devlist,
                                      gchp_vs_gcc_devstr,
                                      dst=plot_dir,
                                      overwrite=True)
        
    if ste_table:
        # --------------------------------------------------------------
        # GCHP vs GCC Strat-Trop Exchange
        # --------------------------------------------------------------
        title = "\n%%% Cannot create GCHP vs. GCC Strat-Trop Exchange table %%%"
        print(title)
     
# ======================================================================
# Create GCHP vs GCHP benchmark plots and tables
# ======================================================================
if gchp_vs_gchp:

    if plot_conc:
        #---------------------------------------------------------------
        # GCHP vs GCHP Concentration plots
        #
        # Includes lumped species and separates by category if
        # plot_by_spc_cat is true; otherwise excludes lumped species
        # and writes to one file
        #---------------------------------------------------------------
        title = "\n%%% Creating GCHP vs. GCHP {} concentration plots %%%"
        print(title.format(bmk_type))

        # Files to read
        collection = "SpeciesConc"
        gchp_vs_gchp_reflist = get_filepaths(gchp_vs_gchp_refdir, collection,
                                             [gchp_date], is_gchp=True)
        gchp_vs_gchp_devlist = get_filepaths(gchp_vs_gchp_devdir, collection,
                                             [gchp_date], is_gchp=True)

        # Make concentration plots
        bmk.make_benchmark_plots(gchp_vs_gchp_reflist[0],
                                 gchp_vs_gchp_refstr,
                                 gchp_vs_gchp_devlist[0],
                                 gchp_vs_gchp_devstr,
                                 dst=gchp_vs_gchp_plotsdir,
                                 weightsdir=weightsdir,
                                 plot_by_spc_cat=plot_by_spc_cat,
                                 overwrite=True,
                                 sigdiff_files=gchp_vs_gchp_sigdiff)

    if plot_emis:
        #---------------------------------------------------------------
        # GCHP vs. GCHP Emissions plots
        #---------------------------------------------------------------
        title = "\n%%% Creating GCHP vs. GCC {} emissions plots %%%"
        print(title.format(bmk_type))

        # Files to read
        collection = "Emissions"
        gchp_vs_gchp_reflist = get_filepaths(gchp_vs_gchp_refdir, collection,
                                             [gchp_date], is_gchp=True)
        gchp_vs_gchp_devlist = get_filepaths(gchp_vs_gchp_devdir, collection,
                                             [gchp_date], is_gchp=True)

        # Create emissions plots
        bmk.make_benchmark_emis_plots(gchp_vs_gchp_reflist[0],
                                      gchp_vs_gchp_refstr,
                                      gchp_vs_gchp_devlist[0],
                                      gchp_vs_gchp_devstr,
                                      dst=gchp_vs_gchp_plotsdir,
                                      weightsdir=weightsdir,
                                      plot_by_spc_cat=plot_by_spc_cat,
                                      plot_by_hco_cat=plot_by_hco_cat,
                                      overwrite=True,
                                      sigdiff_files=gchp_vs_gchp_sigdiff)

    if emis_table:
        #---------------------------------------------------------------
        # GCHP vs. GCHP tables of emission and inventory totals
        #---------------------------------------------------------------
        title = "\n%%% Creating GCHP vs. GCHP {} emissions/inventory tables %%%"
        print(title.format(bmk_type))

        # Files to read
        collection = ["Emissions", "StateMet_avg"]
        gchp_vs_gchp_reflist = get_filepaths(gchp_vs_gchp_refdir, collection,
                                             [gchp_date], is_gchp=True)
        gchp_vs_gchp_devlist = get_filepaths(gchp_vs_gchp_devdir, collection,
                                             [gchp_date], is_gchp=True)

        # Print emisisons and inventory tables
        bmk.make_benchmark_emis_tables(gchp_vs_gchp_reflist,
                                       gchp_vs_gchp_refstr,
                                       gchp_vs_gchp_devlist,
                                       gchp_vs_gchp_devstr,
                                       dst=gchp_vs_gchp_plotsdir,
                                       interval=s_per_mon,
                                       overwrite=True)

    if plot_jvalues:
        #---------------------------------------------------------------
        # GCHP vs. GCHP J-values plots
        #---------------------------------------------------------------
        title = "\n%%% Creating GCHP vs. GCHP {} J-value plots %%%"
        print(title.format(bmk_type))

        # Files to read
        collection = "JValues"
        gchp_vs_gchp_reflist = get_filepaths(gchp_vs_gchp_refdir, collection,
                                             [gchp_date], is_gchp=True)
        gchp_vs_gchp_devlist = get_filepaths(gchp_vs_gchp_devdir, collection,
                                             [gchp_date], is_gchp=True)
        
        # Plot J-values
        bmk.make_benchmark_jvalue_plots(gchp_vs_gchp_reflist[0],
                                        gchp_vs_gchp_refstr,
                                        gchp_vs_gchp_devlist[0],
                                        gchp_vs_gchp_devstr,
                                        dst=gchp_vs_gchp_plotsdir,
                                        weightsdir=weightsdir,
                                        overwrite=True,
                                        sigdiff_files=gchp_vs_gchp_sigdiff)

    if plot_aod:
        #---------------------------------------------------------------
        # GCHP vs GCHP column AOD plots
        #---------------------------------------------------------------
        title = "\n%%% Creating GCHP vs. GCHP {} column AOD plots %%%"
        print(title.format(bmk_type))

        # Files to read
        collection = "Aerosols"
        gchp_vs_gchp_reflist = get_filepaths(gchp_vs_gchp_refdir, collection,
                                            [gchp_date], is_gchp=True)
        gchp_vs_gchp_devlist = get_filepaths(gchp_vs_gchp_devdir, collection,
                                            [gchp_date], is_gchp=True)

        # Plot AODs
        bmk.make_benchmark_aod_plots(gchp_vs_gchp_reflist[0],
                                     gchp_vs_gchp_refstr,
                                     gchp_vs_gchp_devlist[0],
                                     gchp_vs_gchp_devstr,
                                     dst=gchp_vs_gchp_plotsdir,
                                     weightsdir=weightsdir,
                                     overwrite=True,
                                     sigdiff_files=gchp_vs_gchp_sigdiff)

    if mass_table:
        #---------------------------------------------------------------
        # GCHP vs GCC global mass tables
        #---------------------------------------------------------------
        title = "\n%%% Creating GCHP vs. GCC {} global mass tables %%%"
        print(title.format(bmk_type))

        # Files to read
        collection = ["Restart", "StateMet_inst"]
        gchp_vs_gchp_reflist = get_filepaths(gchp_vs_gchp_refrst, collection,
                                            [end_date], is_gchp=True)
        gchp_vs_gchp_devlist = get_filepaths(gchp_vs_gchp_devrst, collection,
                                            [end_date], is_gchp=True)

        # Save to "Tables" subfolder of plot dir
        plot_dir = join(gchp_vs_gchp_plotsdir, "Tables")
        
        # Plot mass tables
        bmk.make_benchmark_mass_tables(gchp_vs_gchp_reflist[0],
                                       gchp_vs_gchp_refstr,
                                       gchp_vs_gchp_devlist[0],
                                       gchp_vs_gchp_devstr,
                                       dst=plot_dir,
                                       overwrite=True)
        
    if budget_table:
        #---------------------------------------------------------------
        # GCHP vs GCC budgets tables 
        #---------------------------------------------------------------
        title = "\n%%% Creating GCHP vs. GCC {} budget tables %%%"
        print(title.format(bmk_type))

        # Files to read
        collection = "Budget"
        gchp_vs_gchp_devlist = get_filepaths(gchp_vs_gchp_devdir, collection,
                                            [gchp_date], is_gchp=True)

        # Save to "Tables" subfolder of plot dir
        plot_dir = join(gchp_vs_gchp_plotsdir, "Budget")
        
        # Print budget tables
        bmk.make_benchmark_budget_tables(gchp_vs_gchp_devlist[0],
                                         gchp_vs_gchp_devstr,
                                         dst=gchp_vs_gchp_plotsdir,
                                         interval=s_per_mon,
                                         overwrite=True)

        # Print "operations budget" for strat, trop, strat+trop 
        opbdg.make_operations_budget_table(gcc_dev_version,
                                           gchp_vs_gchp_devlist[0],
                                           bmk_type,
                                           dst=plot_dir,
                                           label=mon_yr_str,
                                           interval=s_per_mon[0],
                                           overwrite=True)

        
    if OH_metrics:
        #---------------------------------------------------------------
        # GCHP vs. GCHP global mean OH, MCF Lifetime, CH4 Lifetime
        #---------------------------------------------------------------
        title = "\n%%% Creating GCHP vs. GCHP {} OH metrics %%%"
        print(title.format(bmk_type))
        
        # Files to read
        collection = ["ConcAfterChem", "StateMet_avg"]
        gchp_vs_gchp_reflist = get_filepaths(gchp_vs_gchp_refdir, collection,
                                            [gchp_date], is_gchp=True)
        gchp_vs_gchp_devlist = get_filepaths(gchp_vs_gchp_devdir, collection,
                                            [gchp_date], is_gchp=True)

        # Save to "Tables" subfolder of plot dir
        plot_dir = join(gchp_vs_gchp_plotsdir, "Tables")
        
        # Print OH metrics
        bmk.make_benchmark_oh_metrics(gchp_vs_gchp_reflist,
                                      gchp_vs_gchp_refstr,
                                      gchp_vs_gchp_devlist,
                                      gchp_vs_gchp_devstr,
                                      dst=plot_dir,
                                      overwrite=True)
        
    if ste_table:
        # --------------------------------------------------------------
        # GCHP vs GCHP Strat-Trop Exchange
        # --------------------------------------------------------------
        title = "\n%%% Cannot create GCHP vs. GCHP Strat-Trop Exchange table %%%"
        print(title)
     
# =====================================================================
# Create GCHP vs GCC difference of differences benchmark plots
# =====================================================================
if gchp_vs_gcc_diff_of_diffs:

    # NOTE: This can be expanded to differences beyond species
    # concentrations by following how this is done for conc plots.
    print("%\n%% Creating GCHP vs. GCC diff-of-diffs concentration plots %%%")

    # Get a list of variables that GCPy should not read
    skip_vars = skip_these_vars

    # Target output files
    diff_of_diffs_refspc = "./gcc_diffs_spc.nc4"
    diff_of_diffs_devspc = "./gchp_diffs_spc.nc4"

    # Create a ref file that contains GCC differences
    gcc_ref  = xr.open_dataset(gcc_vs_gcc_refspc, drop_variables=skip_vars)
    gcc_dev  = xr.open_dataset(gcc_vs_gcc_devspc, drop_variables=skip_vars)
    with xr.set_options(keep_attrs=True):
        gcc_diffs = gcc_dev - gcc_ref
        for v in gcc_dev.data_vars.keys():
            # Ensure the gcc_diffs Dataset includes attributes
            gcc_diffs[v].attrs = gcc_dev[v].attrs
    gcc_diffs.to_netcdf(diff_of_diffs_refspc)

    # Create a dev file that contains GCHP differences. Include special
    # handling if cubed sphere grid dimension names are different since they
    # changed in MAPL v1.0.0.
    gchp_ref = xr.open_dataset(gchp_vs_gchp_refspc, drop_variables=skip_vars)
    gchp_dev = xr.open_dataset(gchp_vs_gchp_devspc, drop_variables=skip_vars)
    refdims = gchp_ref.dims
    devdims = gchp_dev.dims
    if "lat" in refdims and "Xdim" in devdims:
        gchp_ref_newdimnames = gchp_dev.copy()
        for v in gchp_dev.data_vars.keys():
            if "Xdim" in gchp_dev[v].dims:
                gchp_ref_newdimnames[v].values = gchp_ref[v].values.reshape(
                    gchp_dev[v].values.shape)
                # NOTE: the reverse conversion is gchp_dev[v].stack(lat=("nf","Ydim")).transpose("time","lev","lat","Xdim").values
        gchp_ref = gchp_ref_newdimnames.copy()
    with xr.set_options(keep_attrs=True):
        gchp_diffs = gchp_dev.copy()
        for v in gchp_dev.data_vars.keys():
            if "Xdim" in gchp_dev[v].dims or "lat" in gchp_dev[v].dims:
                gchp_diffs[v] = gchp_dev[v] - gchp_ref[v]
                # NOTE: The gchp_diffs Dataset is created without variable
                # attributes; we have to reattach them
                gchp_diffs[v].attrs = gchp_dev[v].attrs
    gchp_diffs.to_netcdf(diff_of_diffs_devspc)


    # Create diff-of-diff plots for species concentrations
    # (includes lumped species and separates by category)
    #
    # NOTE: Since at the present time we are only printing out
    # diff-of-diffs for concentration plots, we can take this
    # call out of the "if plot_conc:" block.
    bmk.make_benchmark_plots(diff_of_diffs_refspc,
                             diff_of_diffs_refstr,
                             diff_of_diffs_devspc,
                             diff_of_diffs_devstr,
                             dst=diff_of_diffs_plotsdir,
                             overwrite=True,
                             use_cmap_RdBu=True)
