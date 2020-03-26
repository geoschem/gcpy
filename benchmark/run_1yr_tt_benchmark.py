#!/usr/bin/env python
"""
run_1yr_benchmark.py: Driver script for creating benchmark plots
for the 1-year TransportTracers benchmark.

Run this script to generate benchmark comparisons between:

    (1) GCC (aka GEOS-Chem "Classic") vs. GCC

Under development:

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

# =====================================================================
# Configurables
# =====================================================================

# Benchmark information (*MUST EDIT*)
# Note: When doing GCHP vs GCC comparisions gchp_dev will be compared to
#  gcc_dev (not gcc_ref!). This ensures consistency in version names when
#  doing GCHP vs GCC diff-of-diffs (mps, 6/27/19)
maindir  = "/path/to/main/directory"
gcc_ref_version = "gcc_ref_version_string"
gcc_dev_version = "gcc_dev_version_string"
gchp_ref_version = "gchp_ref_version_string"
gchp_dev_version = "gchp_dev_version_string"

# Weights directory (**MUST EDIT**)
weightsdir = "/n/holylfs/EXTERNAL_REPOS/GEOS-CHEM/gcgrid/gcdata/ExtData/GCHP/RegriddingWeights"

# Comparisons to run (edit as needed)
gcc_vs_gcc   = True
gchp_vs_gcc  = False
gchp_vs_gchp = False
gchp_vs_gcc_diff_of_diffs = False

########################################################################
#### TransportTracers BENCHMARK OPTIONS                             ####
########################################################################
bmk_type     = "TransportTracersBenchmark"
plot_conc    = True
plot_wetdep  = True
budget_table = True
ste_table    = True

########################################################################
#### Echo back selections to stdout                                 ####
########################################################################
print("The following plots and tables will be created for {}:".format(bmk_type))
if plot_conc:    print(" - Concentration plots")
if plot_wetdep:  print(" - Convective and large-scale wet deposition plots")
if budget_table: print(" - Budget tables")
if ste_table:    print(" - Table of strat-trop exchange")

# Start and end of the benchmark year (as numpy.datetime64 dates)
bmk_year = 2016
bmk_start_str = "{}-01-01".format(str(bmk_year))
bmk_end_str = "{}-01-01".format(str(bmk_year+1))
bmk_start = np.datetime64(bmk_start_str)
bmk_end = np.datetime64(bmk_end_str)

# Data directories (edit as needed)
gcc_vs_gcc_refdir   = join(maindir, gcc_ref_version,  "OutputDir")
gcc_vs_gcc_devdir   = join(maindir, gcc_dev_version,  "OutputDir")
gchp_vs_gcc_refdir  = join(maindir, gcc_dev_version,  "OutputDir")
gchp_vs_gcc_devdir  = join(maindir, gchp_dev_version, "OutputDir")
gchp_vs_gchp_refdir = join(maindir, gchp_ref_version, "OutputDir")
gchp_vs_gchp_devdir = join(maindir, gchp_dev_version, "OutputDir")

# Restart file directories (edit as needed)
gcc_vs_gcc_refrstdir   = join(maindir, gcc_ref_version, "restarts")
gcc_vs_gcc_devrstdir   = join(maindir, gcc_dev_version, "restarts")
gchp_vs_gcc_refrstdir  = join(maindir, gcc_dev_version, "restarts")
gchp_vs_gcc_devrstdir  = join(maindir, gchp_dev_version)
gchp_vs_gchp_refrstdir = join(maindir, gchp_ref_version)
gchp_vs_gchp_devrstdir = join(maindir, gchp_dev_version)

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
# The rest of these settings should not need to be changed
# =====================================================================

# ----------------------------------------------------------------------
# Months, days, years etc.
# ----------------------------------------------------------------------

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
    
# ----------------------------------------------------------------------
# Files that will contain lists of quantities that have significant
# differences -- we need these for the benchmark approval forms.
# ----------------------------------------------------------------------
if gcc_vs_gcc:
    gcc_vs_gcc_sigdir = join(gcc_vs_gcc_plotsdir, "Sig_Diffs")
    if not os.path.isdir(gcc_vs_gcc_sigdir):
        os.makedirs(gcc_vs_gcc_sigdir)
if gchp_vs_gcc:
    gchp_vs_gcc_sigdir = join(gchp_vs_gcc_plotsdir, "Sig_Diffs")
    if not os.path.isdir(gchp_vs_gcc_sigdir):
        os.makedirs(gchp_vs_gcc_sigdir)
if gchp_vs_gchp:
    gchp_vs_gchp_sigdir = join(gchp_vs_gchp_plotsdir, "Sig_Diffs")
    if not os.path.isdir(gchp_vs_gchp_sigdir):
        os.makedirs(gchp_vs_gchp_sigdir)

# Populate the lists of sigdiff files
gcc_vs_gcc_sigdiff = {}
gchp_vs_gcc_sigdiff = {}
gchp_vs_gchp_sigdiff = {}
for mon_yr_str in bmk_seasons_names:
    
    if gcc_vs_gcc:
        vstr = "{}_vs_{}.{}".format(
            gcc_vs_gcc_refstr, gcc_vs_gcc_devstr, mon_yr_str)
        sigdiff_files = [
            join(gcc_vs_gcc_sigdir, "{}.sig_diffs_sfc.txt".format(vstr)),
            join(gcc_vs_gcc_sigdir, "{}.sig_diffs_500hpa.txt".format(vstr)),
            join(gcc_vs_gcc_sigdir, "{}.sig_diffs_zonalmean.txt".format(vstr)),
            join(gcc_vs_gcc_sigdir, "{}.sig_diffs_emissions.txt".format(vstr))
        ]
        gcc_vs_gcc_sigdiff[mon_yr_str] = sigdiff_files

    if gchp_vs_gcc:
        vstr = "{}_vs_{}.{}".format(
            gchp_vs_gcc_refstr, gchp_vs_gcc_devstr, mon_yr_str)
        sigdiff_files = [
            join(gchp_vs_gcc_sigdir, "{}_sig_diffs_sfc.txt".format(vstr)),
            join(gchp_vs_gcc_sigdir, "{}_sig_diffs_500hpa.txt".format(vstr)),
            join(gchp_vs_gcc_sigdir, "{}_sig_diffs_zonalmean.txt".format(vstr)),
            join(gchp_vs_gcc_sigdir, "{}_sig_diffs_emissions.txt".format(vstr))
        ]
        gchp_vs_gcc_sigdiff[mon_yr_str] = sigdiff_files

    if gchp_vs_gchp:
        vstr = "{}_vs_{}.{}".format(
            gchp_vs_gchp_refstr, gchp_vs_gchp_devstr, mon_yr_str)
        sigdiff_files = [
            join(gchp_vs_gchp_sigdir, "{}_sig_diffs_sfc.txt".format(vstr)),
            join(gchp_vs_gchp_sigdir, "{}_sig_diffs_500hpa.txt".format(vstr)),
            join(gchp_vs_gchp_sigdir, "{}_sig_diffs_zonalmean.txt".format(vstr)),
            join(gchp_vs_gchp_sigdir, "{}_sig_diffs_emissions.txt".format(vstr))
        ]
        gchp_vs_gchp_sigdiff[mon_yr_str] = sigdiff_files

# List of species for the operations budgets
ops_budget_species = ["Rn222",
                      "Pb210",
                      "Pb210Strat",
                      "Be7",
                      "Be7Strat",
                      "Be10",
                      "Be10Strat",
                      "PassiveTracer",
                      "SF6Tracer",
                      "CH3ITracer",
                      "COAnthroEmis25dayTracer",
                      "COAnthroEmis50dayTracer",
                      "COUniformEmis25dayTracer",
                      "GlobEmis90dayTracer",
                      "NHEmis90dayTracer",
                      "SHEmis90dayTracer"]

# ======================================================================
# Functions for area normalization
# ======================================================================
def get_gcc_area(data_dir): 
    """
    Returns the area variable for GEOS-Chem "Classic"
    and renames it to "AREAM2".
    """
    if gcc_vs_gcc or gchp_vs_gcc:
        path = get_filepaths(data_dir, "StateMet",
                            [bmk_months[0]], is_gcc=True)
        with xr.set_options(keep_attrs=True):
            ds = xr.open_mfdataset(path)
            ds = ds.rename({"AREA": "AREAM2"})
            ds = ds["AREAM2"]
            return ds
    else:
        return None

    
def get_gchp_area(data_dir):
    """
    """
    if gchp_vs_gcc or gchp_vs_gchp:
        path = get_filepaths(data_dir, "StateMet_avg",
                             [gchp_months[0]], is_gchp=True)
        with xr.set_options(keep_attrs=True):
            ds = xr.open_mfdataset(path)
            ds = ds.rename({"Met_AREAM2": "AREAM2"})
            ds = ds["AREAM2"]
            return ds
    else:
        return None
    
# ======================================================================
# Create GCC vs GCC benchmark plots and tables
# ======================================================================

if gcc_vs_gcc:

    # Get surface area variables on Ref and Dev grids
    gcc_vs_gcc_areas = {"Ref": get_gcc_area(gcc_vs_gcc_refdir),
                        "Dev": get_gcc_area(gcc_vs_gcc_devdir)}
    
    if plot_conc:
        # --------------------------------------------------------------
        # GCC vs GCC Concentration plots
        # --------------------------------------------------------------
        title = "\n%%% Creating GCC vs. GCC {} concentration plots %%%"
        print(title.format(bmk_type))
   
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
            sigdiff_files= gcc_vs_gcc_sigdiff[mon_yr_str]
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
                                     overwrite=True,
                                     sigdiff_files=sigdiff_files)

    if plot_wetdep:
        # --------------------------------------------------------------
        # GCC vs GCC wet deposition
        # --------------------------------------------------------------
        title = "\n%%% Creating GCC vs. GCC {} wet deposition plots %%%"
        print(title.format(bmk_type))
                
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
                sigdiff_files= gcc_vs_gcc_sigdiff[mon_yr_str]
                bmk.make_benchmark_plots(gcc_vs_gcc_refwd[s],
                                         gcc_vs_gcc_refstr,
                                         gcc_vs_gcc_devwd[s],
                                         gcc_vs_gcc_devstr,
                                         dst=gcc_vs_gcc_plotsdir,
                                         subdst=mon_yr_str,
                                         weightsdir=weightsdir,
                                         overwrite=True,
                                         benchmark_type=bmk_type,
                                         collection=collection,
                                         restrict_cats=[collection],
                                         normalize_by_area=True,
                                         areas=gcc_vs_gcc_areas,
                                         sigdiff_files=sigdiff_files)

    if budget_table:
        # --------------------------------------------------------------
        # GCC vs GCC budgets tables
        # --------------------------------------------------------------
        title = "\n%%% Creating GCC vs. GCC {} budgets %%%"
        print(title.format(bmk_type))

        # Budgets of Radionuclide species
        plot_dir = join(gcc_vs_gcc_plotsdir, "Tables")
        ttbdg.transport_tracers_budgets(gcc_dev_version,
                                        gcc_vs_gcc_devdir,
                                        gcc_vs_gcc_devrstdir,
                                        bmk_year,
                                        dst=plot_dir,
                                        overwrite=True)

        # Change in mass of species after each operation
        gcc_vs_gcc_devbgt = get_filepaths(gcc_vs_gcc_devdir, "Budget",
                                          bmk_months, is_gcc=True)
        label = "{}".format(bmk_year)
        opbdg.make_operations_budget_table(gcc_dev_version,
                                           gcc_vs_gcc_devbgt,
                                           bmk_type,
                                           label,
                                           dst=plot_dir,
                                           interval=sec_per_yr,
                                           species=ops_budget_species,
                                           overwrite=True)

    if ste_table:
        # --------------------------------------------------------------
        # GCC Strat-Trop Exchange
        # --------------------------------------------------------------
        title = \
        "\n%%% Creating GCC vs. GCC {} Strat-Trop Exchange table %%%".format(
            bmk_type)
        print(title)
        
        # Strat-trop exchange of radionuclide species
        gcc_vs_gcc_devflx = get_filepaths(gcc_vs_gcc_devdir, "AdvFluxVert",
                                          bmk_months, is_gcc=True)
        plot_dir = join(gcc_vs_gcc_plotsdir, "Tables")
        species = ["Pb210", "Be7", "Be10"]
        ste.make_benchmark_ste_table(gcc_dev_version,
                                     gcc_vs_gcc_devflx,
                                     bmk_year,
                                     dst=plot_dir,
                                     bmk_type=bmk_type,
                                     species=species,
                                     overwrite=True)

# ======================================================================
# Create GCHP vs GCC benchmark plots and tables
# ======================================================================

if gchp_vs_gcc:

    # Get surface area variables on Ref and Dev grids
    gchp_vs_gcc_areas = {"Ref": get_gcc_area(gchp_vs_gcc_refdir),
                         "Dev": get_gchp_area(gchp_vs_gcc_devdir)}
    
    if plot_conc:
        # --------------------------------------------------------------
        # GCHP vs GCC Concentration plots
        # --------------------------------------------------------------
        title = "\n%%% Creating GCHP vs. GCC {} concentration plots %%%"
        print(title.format(bmk_type))
   
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
            sigdiff_files= gchp_vs_gcc_sigdiff[mon_yr_str]
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
                                     overwrite=True,
                                     sigdiff_files=sigdiff_files)

    if plot_wetdep:
        # --------------------------------------------------------------
        # GCHP vs GCC wet deposition
        # --------------------------------------------------------------
        title = "\n%%% Creating GCHP vs. GCC {} wet deposition plots %%%"
        print(title.format(bmk_type))
        
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
                sigdiff_files= gchp_vs_gcc_sigdiff[mon_yr_str]
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
                                         areas=gchp_vs_gcc_areas,
                                         sigdiff_files=sigdiff_files)

    if budget_table:
        # --------------------------------------------------------------
        # GCHP vs GCC budgets tables
        # --------------------------------------------------------------
        title = "\n%%% Creating GCHP vs. GCC {} budgets %%%"
        print(title.format(bmk_type))

        # Budgets of Radionuclide species
        print("-- Budgets of radionuclide species")
        plot_dir = join(gchp_vs_gcc_plotsdir, "Tables")
        ttbdg.transport_tracers_budgets(gchp_dev_version,
                                        gchp_vs_gcc_devdir,
                                        gchp_vs_gcc_devrstdir,
                                        bmk_year,
                                        dst=plot_dir,
                                        is_gchp=True,
                                        overwrite=True)

        # Change in mass of species after each operation
        print("-- Change in mass after each operation")
        gchp_vs_gcc_devbgt = get_filepaths(gchp_vs_gcc_devdir, "Budget",
                                           gchp_months, is_gchp=True)
        label = "{}".format(bmk_year)
        opbdg.make_operations_budget_table(gchp_dev_version,
                                           gchp_vs_gcc_devbgt,
                                           bmk_type,
                                           label,
                                           dst=plot_dir,
                                           interval=sec_per_yr,
                                           species=ops_budget_species,
                                           overwrite=True)

    if ste_table:
        # --------------------------------------------------------------
        # GCHP vs GCC Strat-Trop Exchange
        # --------------------------------------------------------------
        title = "%%% Cannot create Strat-Trop Exchange tables for GCHP %%%"
        print(title)
        
# =====================================================================
# Create GCHP vs GCHP benchmark plots and tables
# =====================================================================

if gchp_vs_gchp:

    # Get surface area variables on Ref and Dev grids
    gchp_vs_gchp_areas = {"Ref": get_gchp_area(gchp_vs_gchp_refdir),
                          "Dev": get_gchp_area(gchp_vs_gchp_devdir)}
 
    if plot_conc:
        # --------------------------------------------------------------
        # GCC vs GCC Concentration plots
        # --------------------------------------------------------------
        title = "\n%%% Creating GCHP vs. GCHP {} concentration plots %%%"
        print(title.format(bmk_type))
   
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
            sigdiff_files= gchp_vs_gchp_sigdiff[mon_yr_str]
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
                                     overwrite=True,
                                     sigdiff_files=sigdiff_files)

    if plot_wetdep:
        # --------------------------------------------------------------
        # GCHP vs GCHP wet deposition
        # --------------------------------------------------------------
        title = "\n%%% Creating GCHP vs. GCHP {} wet deposition plots %%%"
        print(title.format(bmk_type))
                
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
                sigdiff_files= gchp_vs_gchp_sigdiff[mon_yr_str]
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
                                         restrict_cats=[collection],
                                         normalize_by_area=True,
                                         areas=gchp_vs_gchp_areas,
                                         sigdiff_files=sigdiff_files)

    if budget_table:
        # --------------------------------------------------------------
        # GCHP vs GCHP budgets tables
        # --------------------------------------------------------------
        title = "\n%%% Creating GCHP vs. GCHP {} budgets %%%"
        print(title.format(bmk_type))

        # Budgets of Radionuclide species
        print("-- Budgets of radionuclide species")
        plot_dir = join(gchp_vs_gchp_plotsdir, "Tables")
        ttbdg.transport_tracers_budgets(gchp_dev_version,
                                        gchp_vs_gchp_devdir,
                                        gchp_vs_gchp_devrstdir,
                                        bmk_year,
                                        dst=plot_dir,
                                        is_gchp=True,
                                        overwrite=True)

        # Change in mass of species after each operation
        print("-- Change in mass after each operation")
        gchp_vs_gcc_devbgt = get_filepaths(gchp_vs_gchp_devdir, "Budget",
                                           gchp_months, is_gchp=True)
        label = "{}".format(bmk_year)
        opbdg.make_operations_budget_table(gchp_dev_version,
                                           gchp_vs_gchp_devbgt,
                                           bmk_type,
                                           label,
                                           dst=plot_dir,
                                           interval=sec_per_yr,
                                           species=species,
                                           overwrite=True)

    if ste_table:
        # --------------------------------------------------------------
        # GCHP vs GCC Strat-Trop Exchange
        # --------------------------------------------------------------
        title = "%%% Cannot create Strat-Trop Exchange tables for GCHP %%%"
        print(title)

## =====================================================================
## Create GCHP vs GCC difference of differences benchmark plots
## =====================================================================
#
#if gchp_vs_gcc_diff_of_diffs:
#
#    # NOTE: This can be expanded to differences beyond species
#    # concentrations by following how this is done for conc plots.
#    print("%\n%% Creating GCHP vs. GCC diff-of-diffs concentration plots %%%")
#
#    # Get a list of variables that GCPy should not read
#    skip_vars = skip_these_vars
#
#    # Target output files
#    diff_of_diffs_refspc = "./gcc_diffs_spc.nc4"
#    diff_of_diffs_devspc = "./gchp_diffs_spc.nc4"
#
#    # Create a ref file that contains GCC differences
#    gcc_ref  = xr.open_dataset(gcc_vs_gcc_refspc, drop_variables=skip_vars)
#    gcc_dev  = xr.open_dataset(gcc_vs_gcc_devspc, drop_variables=skip_vars)
#    with xr.set_options(keep_attrs=True):
#        gcc_diffs = gcc_dev - gcc_ref
#        for v in gcc_dev.data_vars.keys():
#            # Ensure the gcc_diffs Dataset includes attributes
#            gcc_diffs[v].attrs = gcc_dev[v].attrs
#    gcc_diffs.to_netcdf(diff_of_diffs_refspc)
#
#    # Create a dev file that contains GCHP differences. Include special
#    # handling if cubed sphere grid dimension names are different since they
#    # changed in MAPL v1.0.0.
#    gchp_ref = xr.open_dataset(gchp_vs_gchp_refspc, drop_variables=skip_vars)
#    gchp_dev = xr.open_dataset(gchp_vs_gchp_devspc, drop_variables=skip_vars)
#    refdims = gchp_ref.dims
#    devdims = gchp_dev.dims
#    if "lat" in refdims and "Xdim" in devdims:
#        gchp_ref_newdimnames = gchp_dev.copy()
#        for v in gchp_dev.data_vars.keys():
#            if "Xdim" in gchp_dev[v].dims:
#                gchp_ref_newdimnames[v].values = gchp_ref[v].values.reshape(
#                    gchp_dev[v].values.shape)
#                # NOTE: the reverse conversion is gchp_dev[v].stack(lat=("nf","Ydim")).transpose("time","lev","lat","Xdim").values
#        gchp_ref = gchp_ref_newdimnames.copy()
#    with xr.set_options(keep_attrs=True):
#        gchp_diffs = gchp_dev.copy()
#        for v in gchp_dev.data_vars.keys():
#            if "Xdim" in gchp_dev[v].dims or "lat" in gchp_dev[v].dims:
#                gchp_diffs[v] = gchp_dev[v] - gchp_ref[v]
#                # NOTE: The gchp_diffs Dataset is created without variable
#                # attributes; we have to reattach them
#                gchp_diffs[v].attrs = gchp_dev[v].attrs
#    gchp_diffs.to_netcdf(diff_of_diffs_devspc)
#
#
#    # Create diff-of-diff plots for species concentrations
#    # (includes lumped species and separates by category)
#    #
#    # NOTE: Since at the present time we are only printing out
#    # diff-of-diffs for concentration plots, we can take this
#    # call out of the "if plot_conc:" block.
#    bmk.make_benchmark_conc_plots(diff_of_diffs_refspc,
#                                  diff_of_diffs_refstr,
#                                  diff_of_diffs_devspc,
#                                  diff_of_diffs_devstr,
#                                  dst=diff_of_diffs_plotsdir,
#                                  overwrite=True,
#                                  use_cmap_RdBu=True)
#
###############################################################################
