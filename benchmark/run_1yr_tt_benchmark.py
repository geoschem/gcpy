#!/usr/bin/env python
"""
run_1yr_benchmark.py: Driver script for creating benchmark plots
for the 1-year TransportTracers benchmark.

Run this script to generate benchmark comparisons between:

    (1) GCC (aka GEOS-Chem "Classic") vs. GCC

Under development:

    (2) GCHP vs GCC
    (3) GCHP vs GCHP

Note that we currently do not create diff-of-diff plots for the transport
tracer benchmark.
 
You can customize this script by editing the following settings in the
"Configurables" section below:

    (1) Edit the path variables so that they point to folders w/ model data
    (2) Edit the version strings for each benchmark simulation
    (3) Edit the switches that turn on/off creating of plots and tables
    (4) If necessary, edit labels for the dev and ref versions

Calling sequence:

    ./run_1yr_tt_benchmark.py

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
from gcpy.constants import skip_these_vars
from gcpy.core import get_filepaths
import gcpy.budget_ops as opbdg
import gcpy.budget_tt as ttbdg
import gcpy.ste_flux as ste
import numpy as np
import warnings

# Tell matplotlib not to look for an X-window
os.environ["QT_QPA_PLATFORM"]="offscreen"

# Suppress harmless run-time warnings (mostly about underflow in division)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# This script has a fixed benchmark type
bmk_type     = "TransportTracersBenchmark"

# Path to regridding weights
weightsdir = "/n/holyscratch01/external_repos/GEOS-CHEM/gcgrid/gcdata/ExtData/GCHP/RegriddingWeights"

########################################################################
###           CONFIGURABLE SETTINGS: EDIT THESE ACCORDINGLY          ###
########################################################################

# =====================================================================
# Benchmark information (**EDIT AS NEEDED**)
# Note: When doing GCHP vs GCC comparisions gchp_dev will be compared
# to gcc_dev (not gcc_ref!). This ensures consistency in version names
# when doing GCHP vs GCC diff-of-diffs (mps, 6/27/19)
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

# =====================================================================
# Specify if this is a gcpy test validation run
# =====================================================================
gcpy_test = True

# =====================================================================
# Comparisons to run (**EDIT AS NEEDED**)
# =====================================================================
gcc_vs_gcc   = True
gchp_vs_gcc  = True
gchp_vs_gchp = True
# GCHP vs GCC diff of diffs not included in transport tracer benchmark

# =====================================================================
# Output to generate (**EDIT AS NEEDED**)
# =====================================================================
plot_conc         = True
plot_wetdep       = True
rnpbbe_budget     = True
operations_budget = True
ste_table         = True # GCC only

# =====================================================================
# Data directories (**EDIT AS NEEDED**)
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
    mainplotsdir          = './Plots'
    gcc_vs_gcc_plotsdir    = join(mainplotsdir,'GCC_version_comparison')
    gchp_vs_gcc_plotsdir   = join(mainplotsdir,'GCHP_GCC_comparison')
    gchp_vs_gchp_plotsdir  = join(mainplotsdir,'GCHP_version_comparison')
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

# Seconds and days in the benchmark year
days_per_month = np.zeros(12)
sec_per_month = np.zeros(12)
for m in range(12):
    days_per_month[m] = monthrange(bmk_year, m + 1)[1]
    sec_per_month[m] = days_per_month[m] * 86400.0
days_per_yr = np.sum(days_per_month)
sec_per_yr = np.sum(sec_per_month)

# Timestamps for GCHP (these are in the middle of the month)
gchp_months = np.zeros(12, dtype="datetime64[h]")
for m in range(12):
    if days_per_month[m] == 31:
        delta = np.timedelta64(((15 * 24) + 12), 'h')
    elif days_per_month[m] == 30:
        delta = np.timedelta64((15 * 24), 'h')
    elif days_per_month[m] == 29:
        delta = np.timedelta64(((14 * 24) + 12), 'h')
    else:
        delta = np.timedelta64((14 * 24), 'h')
    gchp_months[m] = bmk_months[m].astype("datetime64[h]") + delta
gchp_seasons = gchp_months[[0, 3, 6, 9]]

# ======================================================================
# Echo the list of plots & tables that will be made to the screen
# ======================================================================

print("The following plots and tables will be created for {}:".format(bmk_type))
if plot_conc:    print(" - Concentration plots")
if plot_wetdep:  print(" - Convective and large-scale wet deposition plots")
if rnpbbe_budget: print(" - Radionuclides budget table")
if operations_budget: print(" - Operations budget table")
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
    # --------------------------------------------------------------
    if plot_conc:
        print("\n%%% Creating GCC vs. GCC concentration plots %%%")

        # File lists for emissions data (seasonal)
        collection = "SpeciesConc"
        gcc_vs_gcc_refspc = get_filepaths(gcc_vs_gcc_refdir, collection,
                                          bmk_seasons, is_gcc=True)
        gcc_vs_gcc_devspc = get_filepaths(gcc_vs_gcc_devdir, collection,
                                          bmk_seasons, is_gcc=True)

        # Only plot concentration categories for TransportTracers
        restrict_cats = ["RnPbBeTracers", "PassiveTracers"]

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
                                     restrict_cats=restrict_cats,
                                     overwrite=True)

    # --------------------------------------------------------------
    # GCC vs GCC wet deposition plots
    # --------------------------------------------------------------
    if plot_wetdep:
        print("\n%%% Creating GCC vs. GCC wet deposition plots %%%")

        # Loop over wet deposition collections
        collection_list = ["WetLossConv", "WetLossLS"]
        for collection in collection_list:
            gcc_vs_gcc_refwd = get_filepaths(gcc_vs_gcc_refdir, collection,
                                             bmk_seasons, is_gcc=True)
            gcc_vs_gcc_devwd = get_filepaths(gcc_vs_gcc_devdir, collection,
                                             bmk_seasons, is_gcc=True)

            # Create seasonal plots for wet scavenging
            for s in range(bmk_nseasons):
                mon_yr_str = bmk_seasons_names[s]
                bmk.make_benchmark_plots(gcc_vs_gcc_refwd[s],
                                         gcc_vs_gcc_refstr,
                                         gcc_vs_gcc_devwd[s],
                                         gcc_vs_gcc_devstr,
                                         dst=gcc_vs_gcc_plotsdir,
                                         subdst=mon_yr_str,
                                         weightsdir=weightsdir,
                                         benchmark_type=bmk_type,
                                         collection=collection,
                                         restrict_cats=[collection],
                                         overwrite=True)

    # --------------------------------------------------------------
    # GCC vs GCC radionuclides budget tables
    # --------------------------------------------------------------
    if rnpbbe_budget:
        print("\n%%% Creating GCC vs. GCC radionuclides budget table %%%")
        ttbdg.transport_tracers_budgets(gcc_dev_version,
                                        gcc_vs_gcc_devdir,
                                        gcc_vs_gcc_devrstdir,
                                        bmk_year,
                                        dst=gcc_vs_gcc_tablesdir,
                                        overwrite=True)

    # --------------------------------------------------------------
    # GCC vs GCC operations budgets tables
    # --------------------------------------------------------------
    if operations_budget:
        print("\n%%% Creating GCC vs. GCC operations budget tables %%%")
        gcc_vs_gcc_reflist = get_filepaths(gcc_vs_gcc_refdir, "Budget",
                                           bmk_months, is_gcc=True)
        gcc_vs_gcc_devlist = get_filepaths(gcc_vs_gcc_devdir, "Budget",
                                           bmk_months, is_gcc=True)
        opbdg.make_operations_budget_table(gcc_ref_version,
                                           gcc_vs_gcc_reflist,
                                           gcc_dev_version,
                                           gcc_vs_gcc_devlist,
                                           bmk_type,
                                           dst=gcc_vs_gcc_tablesdir,
                                           label=str(bmk_year),
                                           interval=sec_per_yr,
                                           overwrite=True,
                                           pd_float_format="{:13.6e}")

    # --------------------------------------------------------------
    # GCC dev strat-trop exchange table
    # --------------------------------------------------------------
    if ste_table:
        print("\n%%% Creating GCC vs. GCC Strat-Trop Exchange table %%%")
        gcc_vs_gcc_devflx = get_filepaths(gcc_vs_gcc_devdir, "AdvFluxVert",
                                          bmk_months, is_gcc=True)
        species = ["Pb210","Be7","Be10"]
        ste.make_benchmark_ste_table(gcc_dev_version,
                                     gcc_vs_gcc_devflx,
                                     bmk_year,
                                     dst=gcc_vs_gcc_tablesdir,
                                     bmk_type=bmk_type,
                                     species=species,
                                     overwrite=True)

# ======================================================================
# Create GCHP vs GCC benchmark plots and tables
# ======================================================================

if gchp_vs_gcc:

    # --------------------------------------------------------------
    # GCHP vs GCC Concentration plots
    # --------------------------------------------------------------
    if plot_conc:
        print("\n%%% Creating GCHP vs. GCC concentration plots %%%")

        # File lists for emissions data (seasonal)
        collection = "SpeciesConc"
        gchp_vs_gcc_refspc = get_filepaths(gchp_vs_gcc_refdir, collection,
                                           bmk_seasons, is_gcc=True)
        gchp_vs_gcc_devspc = get_filepaths(gchp_vs_gcc_devdir, collection,
                                           gchp_seasons, is_gchp=True)

        # Only plot concentration categories for TransportTracers
        restrict_cats = ["RnPbBeTracers", "PassiveTracers"]

        # Create seasonal concentration plots
        for s in range(bmk_nseasons):
            mon_yr_str = bmk_seasons_names[s]
            bmk.make_benchmark_plots(gchp_vs_gcc_refspc[s],
                                     gchp_vs_gcc_refstr,
                                     gchp_vs_gcc_devspc[s],
                                     gchp_vs_gcc_devstr,
                                     dst=gchp_vs_gcc_plotsdir,
                                     subdst=mon_yr_str,
                                     weightsdir=weightsdir,
                                     benchmark_type=bmk_type,
                                     collection=collection,
                                     restrict_cats=restrict_cats,
                                     overwrite=True)

    # --------------------------------------------------------------
    # GCHP vs GCC wet deposition plots
    # --------------------------------------------------------------
    if plot_wetdep:
        print("\n%%% Creating GCHP vs. GCC wet deposition plots %%%")

        # Get GCHP area array from StateMet diagnostic file since not in
        # the wet loss diagnostics file. Must be called 'AREA' and be m2.
        gchpareapath = get_filepaths(gchp_vs_gcc_devdir, "StateMet_avg",
                                     [gchp_months[0]], is_gchp=True)
        ds_gchp = xr.open_mfdataset(gchpareapath)
        ds_gchp = ds_gchp.rename({'Met_AREAM2': 'AREA'})

        # Store area DataArrays as dictionary. The GCC area can be empty
        # since GCC diagnostics include variable 'AREA' in m2. Since area
        # is time invariant, drop the time dimension to avoid merge issues
        # for data files from other seasons.
        gchp_vs_gcc_areas = {'Ref': [], 
                             'Dev': ds_gchp['AREA'].isel(time=0).drop('time')}

        # Loop over wet deposition collections
        collection_list = ["WetLossConv", "WetLossLS"]
        for collection in collection_list:
            gchp_vs_gcc_refwd = get_filepaths(gchp_vs_gcc_refdir, collection,
                                              bmk_seasons, is_gcc=True)
            gchp_vs_gcc_devwd = get_filepaths(gchp_vs_gcc_devdir, collection,
                                              gchp_seasons, is_gchp=True)

            # Create seasonal plots for wet scavenging
            for s in range(bmk_nseasons):
                mon_yr_str = bmk_seasons_names[s]
                bmk.make_benchmark_plots(gchp_vs_gcc_refwd[s],
                                         gchp_vs_gcc_refstr,
                                         gchp_vs_gcc_devwd[s],
                                         gchp_vs_gcc_devstr,
                                         dst=gchp_vs_gcc_plotsdir,
                                         subdst=mon_yr_str,
                                         weightsdir=weightsdir,
                                         overwrite=True,
                                         benchmark_type=bmk_type,
                                         collection=collection,
                                         restrict_cats=[collection],
                                         normalize_by_area=True,
                                         areas=gchp_vs_gcc_areas)

    # --------------------------------------------------------------
    # GCHP vs GCC radionuclides budget tables
    # --------------------------------------------------------------
    if rnpbbe_budget:
        print("\n%%% Creating GCHP vs. GCC radionuclides budget table %%%")
        ttbdg.transport_tracers_budgets(gchp_dev_version,
                                        gchp_vs_gcc_devdir,
                                        gchp_vs_gcc_devrstdir,
                                        bmk_year,
                                        dst=gchp_vs_gcc_tablesdir,
                                        is_gchp=True,
                                        overwrite=True)

    # --------------------------------------------------------------
    # GCHP vs GCC operations budgets tables
    # --------------------------------------------------------------
    if operations_budget:
        print("\n%%% Creating GCHP vs. GCC operations budget tables %%%")
        gchp_vs_gcc_reflist = get_filepaths(gchp_vs_gcc_refdir, "Budget",
                                            bmk_months, is_gcc=True)
        gchp_vs_gcc_devlist = get_filepaths(gchp_vs_gcc_devdir, "Budget",
                                            gchp_months, is_gchp=True)
        opbdg.make_operations_budget_table(gcc_dev_version,
                                           gchp_vs_gcc_reflist,
                                           gchp_dev_version,
                                           gchp_vs_gcc_devlist,
                                           bmk_type,
                                           dst=gchp_vs_gcc_tablesdir,
                                           label=str(bmk_year),
                                           interval=sec_per_yr,
                                           overwrite=True,
                                           pd_float_format="{:13.6e}")

# =====================================================================
# Create GCHP vs GCHP benchmark plots and tables
# =====================================================================

if gchp_vs_gchp:

    # --------------------------------------------------------------
    # GCC vs GCC Concentration plots
    # --------------------------------------------------------------
    if plot_conc:
        print("\n%%% Creating GCHP vs. GCHP concentration plots %%%")

        # File lists for emissions data (seasonal)
        collection = "SpeciesConc"
        gchp_vs_gchp_refspc = get_filepaths(gchp_vs_gchp_refdir, collection,
                                            gchp_seasons, is_gchp=True)
        gchp_vs_gchp_devspc = get_filepaths(gchp_vs_gcc_devdir, collection,
                                            gchp_seasons, is_gchp=True)

        # Only plot concentration categories for TransportTracers
        restrict_cats = ["RnPbBeTracers", "PassiveTracers"]

        # Create seasonal concentration plots
        for s in range(bmk_nseasons):
            mon_yr_str = bmk_seasons_names[s]
            bmk.make_benchmark_plots(gchp_vs_gchp_refspc[s],
                                     gchp_vs_gchp_refstr,
                                     gchp_vs_gchp_devspc[s],
                                     gchp_vs_gchp_devstr,
                                     dst=gchp_vs_gchp_plotsdir,
                                     subdst=mon_yr_str,
                                     weightsdir=weightsdir,
                                     benchmark_type=bmk_type,
                                     collection=collection,
                                     restrict_cats=restrict_cats,
                                     overwrite=True)

    # --------------------------------------------------------------
    # GCHP vs GCHP wet deposition plots
    # --------------------------------------------------------------
    if plot_wetdep:
        print("\n%%% Creating GCHP vs. GCHP wet deposition plots %%%")

        # Loop over wet deposition collections
        collection_list = ["WetLossConv", "WetLossLS"]
        for collection in collection_list:
            gchp_vs_gchp_refwd = get_filepaths(gchp_vs_gchp_refdir, collection,
                                               gchp_seasons, is_gchp=True)
            gchp_vs_gchp_devwd = get_filepaths(gchp_vs_gchp_devdir, collection,
                                               gchp_seasons, is_gchp=True)

            # Create seasonal plots for wet scavenging
            for s in range(bmk_nseasons):
                mon_yr_str = bmk_seasons_names[s]
                bmk.make_benchmark_plots(gchp_vs_gchp_refwd[s],
                                         gchp_vs_gchp_refstr,
                                         gchp_vs_gchp_devwd[s],
                                         gchp_vs_gchp_devstr,
                                         dst=gchp_vs_gchp_plotsdir,
                                         subdst=mon_yr_str,
                                         weightsdir=weightsdir,
                                         overwrite=True,
                                         benchmark_type=bmk_type,
                                         collection=collection,
                                         restrict_cats=[collection])

    # --------------------------------------------------------------
    # GCHP vs GCHP radionuclides budget table
    # --------------------------------------------------------------
    if rnpbbe_budget:
        print("\n%%% Creating GCHP vs. GCHP radionuclides budget table %%%")
        ttbdg.transport_tracers_budgets(gchp_dev_version,
                                        gchp_vs_gchp_devdir,
                                        gchp_vs_gchp_devrstdir,
                                        bmk_year,
                                        dst=gchp_vs_gchp_tablesdir,
                                        is_gchp=True,
                                        overwrite=True)

    # --------------------------------------------------------------
    # GCHP vs GCHP operations budgets tables
    # --------------------------------------------------------------
    if operations_budget:
        print("\n%%% Creating GCHP vs. GCHP operations budget tables %%%")
        gchp_vs_gchp_reflist = get_filepaths(gchp_vs_gchp_refdir, "Budget",
                                             gchp_months, is_gchp=True)
        gchp_vs_gchp_devlist = get_filepaths(gchp_vs_gchp_devdir, "Budget",
                                             gchp_months, is_gchp=True)
        opbdg.make_operations_budget_table(gchp_dev_version,
                                           gchp_vs_gchp_reflist,
                                           gchp_dev_version,
                                           gchp_vs_gchp_devlist,
                                           bmk_type,
                                           dst=gchp_vs_gchp_tablesdir,
                                           label=str(bmk_year),
                                           interval=sec_per_yr,
                                           overwrite=True,
                                           pd_float_format="{:13.6e}")
