#!/usr/bin/env python
'''
run_1yr_benchmark.py: Driver script for creating benchmark plots.

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

    ./run_1yr_benchmark.py

Remarks:

    By default, matplotlib will try to open an X window for plotting.
    If you are running this script in an environment where you do not have
    an active X display (such as in a computational queue), then you will
    need to use these commands to disable the X-window functionality.

        import os
        os.environ['QT_QPA_PLATFORM']='offscreen'

    For more information, please see this issue posted at the ipython site:

        https://github.com/ipython/ipython/issues/10627

    This issue might be fixed in matplotlib 3.0.
'''

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
import gcpy.budget_aer as aerbdg
import gcpy.budget_ops as opbdg
import gcpy.budget_tt as ttbdg
import gcpy.ste_flux as ste
import pandas as pd
import numpy as np
import warnings

# Tell matplotlib not to look for an X-window
os.environ['QT_QPA_PLATFORM']='offscreen'

# Suppress harmless run-time warnings (mostly about underflow in division)
warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', category=UserWarning)

# =====================================================================
# Configurables
# =====================================================================

# Benchmark information (*MUST EDIT*)
# Note: When doing GCHP vs GCC comparisions gchp_dev will be compared to
#  gcc_dev (not gcc_ref!). This ensures consistency in version names when
#  doing GCHP vs GCC diff-of-diffs (mps, 6/27/19)
maindir  = '/path/to/main/directory'
gcc_ref_version = 'gcc_ref_version_string'
gcc_dev_version = 'gcc_dev_version_string'
gchp_ref_version = 'gchp_ref_version_string'
gchp_dev_version = 'gchp_dev_version_string'

# Comparisons to run (edit as needed)
# NOTE: For now, just develop GCC vs GCC and add others later (bmy, 10/11/19)
gcc_vs_gcc   = True
gchp_vs_gcc  = False
gchp_vs_gchp = False
gchp_vs_gcc_diff_of_diffs = False

########################################################################
##### FULL-CHEMISTRY BENCHMARK OPTIONS (comment out if not needed) #####
########################################################################
# Output to generate (edit as needed)
# Plots/tables will be created in this order:
bmk_type     = 'FullChemBenchmark'
plot_conc    = True
plot_emis    = True
emis_table   = True
plot_jvalues = True
plot_aod     = True
mass_table   = True
budget_table = True
ste_table    = True
OH_metrics   = True
plot_wetdep  = False

########################################################################
#### TransportTracers BENCHMARK OPTIONS (comment out if not needed) ####
########################################################################
#bmk_type     = 'TransportTracersBenchmark'
#plot_conc    = False
#plot_wetdep  = False
#budget_table = True
#ste_table    = True
#plot_emis    = False
#emis_table   = False
#plot_jvalues = False
#plot_aod     = False
#mass_table   = False
#OH_metrics   = False

########################################################################
#### Echo back selections to stdout                                 ####
########################################################################
print("The following plots and tables will be created for {}:".format(bmk_type))
if plot_conc:    print(" - Concentration plots")
if plot_wetdep:  print(" - Convective and large-scale wet deposition plots")
if plot_emis:    print(" - Emissions plots")
if plot_jvalues: print(" - J-values (photolysis rates) plots")
if plot_aod:     print(" - Aerosol optical depth plots")
if budget_table: print(" - Budget tables")
if emis_table:   print(" - Table of emissions totals by species and inventory")
if mass_table:   print(" - Table of species mass")
if ste_table:    print(" - Table of strat-trop exchange")

# Start and end of the benchmark year (as numpy.datetime64 dates)
bmk_year = 2016
bmk_start_str = "{}-01-01".format(str(bmk_year))
bmk_end_str = "{}-01-01".format(str(bmk_year+1))
bmk_start = np.datetime64(bmk_start_str)
bmk_end = np.datetime64(bmk_end_str)

# Data directories (edit as needed)
gcc_vs_gcc_refdir   = join(maindir, gcc_ref_version,  'OutputDir')
gcc_vs_gcc_devdir   = join(maindir, gcc_dev_version,  'OutputDir')
gchp_vs_gcc_refdir  = join(maindir, gcc_dev_version,  'OutputDir')
gchp_vs_gcc_devdir  = join(maindir, gchp_dev_version, 'OutputDir')
gchp_vs_gchp_refdir = join(maindir, gchp_ref_version, 'OutputDir')
gchp_vs_gchp_devdir = join(maindir, gchp_dev_version, 'OutputDir')

# Restart file directories (edit as needed)
gcc_vs_gcc_refrstdir   = join(maindir, gcc_ref_version,  'restarts')
gcc_vs_gcc_devrstdir   = join(maindir, gcc_dev_version,  'restarts')
gchp_vs_gcc_refrstdir  = join(maindir, gcc_dev_version,  'restarts')
gchp_vs_gcc_devrstdir  = join(maindir, gchp_dev_version, 'restarts')
gchp_vs_gchp_refrstdir = join(maindir, gchp_ref_version, 'restarts')
gchp_vs_gchp_devrstdir = join(maindir, gchp_dev_version, 'restarts')

# Plots directories (edit as needed)
gcc_vs_gcc_plotsdir    = join(maindir, gcc_dev_version, 'Plots')
gchp_vs_gchp_plotsdir  = join(maindir, gchp_dev_version,
                              'Plots/GCHP_version_comparison')
gchp_vs_gcc_plotsdir   = join(maindir, gchp_dev_version,
                              'Plots/GCHP_GCC_comparison')
diff_of_diffs_plotsdir = join(maindir, gchp_dev_version,
                              'Plots/GCHP_GCC_diff_of_diffs')

# Plot title strings (edit as needed)
# For gchp_vs_gcc_refstr use gcc_dev_version, not ref (mps, 6/27/19)
gcc_vs_gcc_refstr    = '{}'.format(gcc_ref_version)
gcc_vs_gcc_devstr    = '{}'.format(gcc_dev_version)
gchp_vs_gcc_refstr   = '{}'.format(gcc_dev_version)
gchp_vs_gcc_devstr   = '{}'.format(gchp_dev_version)
gchp_vs_gchp_refstr  = '{}'.format(gchp_ref_version)
gchp_vs_gchp_devstr  = '{}'.format(gchp_dev_version)
diff_of_diffs_refstr = '{} - {}'.format(gcc_dev_version,
                                        gcc_ref_version)
diff_of_diffs_devstr = '{} - {}'.format(gchp_dev_version,
                                        gchp_ref_version)

# =====================================================================
# The rest of these settings should not need to be changed
# =====================================================================

###############################################################################
# Under development, comment these out for now
## Paths to concentration after chemistry files
#gcc_vs_gcc_refcac = core.get_gcc_filepath(gcc_vs_gcc_refdir, 'ConcAfterChem',
#                                          year=bmk_year, month=bmk_months)
#
#gcc_vs_gcc_devcac = core.get_gcc_filepath(gcc_vs_gcc_devdir, 'ConcAfterChem',
#                                          year=bmk_year, month=bmk_months)
###############################################################################

# Monthly array of dates
bmk_delta_1m = np.timedelta64(1, 'M')
bmk_months = np.arange(bmk_start, bmk_end,
                       step=bmk_delta_1m, dtype='datetime64[M]')

# Get the benchmark year from the datetime64 object
bmk_year = bmk_months[0].astype('datetime64[Y]').astype(int) + 1970

# Seasonal array of dates
bmk_delta_3m = np.timedelta64(3, 'M')
bmk_seasons = np.arange(bmk_start, bmk_end,
                        step=bmk_delta_3m, dtype='datetime64[M]')
bmk_nseasons = len(bmk_seasons)

# Names for each season (e.g. Jan2016, Apr2016, Jul2016, Oct2016)
bmk_seasons_names = ['Jan', 'Apr', 'Jul', 'Oct']
bmk_seasons_names = [v + str(bmk_year) for v in bmk_seasons_names]

# Seconds in each month of the benchmark year
sec_per_month = np.zeros(12)
for m in range(1, 13):
    month_info = monthrange(bmk_year, m)
    sec_per_month[m-1] = month_info[1] * 86400.0

# Files that will contain lists of quantities that have significant
# differences -- we need these for the benchmark approval forms.
sigdiff_dir = join(gcc_vs_gcc_plotsdir, "Sig_Diffs")
if not os.path.isdir(sigdiff_dir):
    os.mkdir(sigdiff_dir)
gcc_vs_gcc_sigdiff = {}
for mon_yr_str in bmk_seasons_names:
    vstr = '{}_vs_{}.{}'.format(
        gcc_vs_gcc_refstr, gcc_vs_gcc_devstr, mon_yr_str)

    sigdiff_files = [
        join(sigdiff_dir, "{}/{}.sig_diffs_sfc.txt".format(sigdiff_dir, vstr)),
        join(sigdiff_dir, "{}/{}.sig_diffs_500hpa.txt".format(sigdiff_dir, vstr)),
        join(sigdiff_dir, "{}/{}.sig_diffs_zonalmean.txt".format(sigdiff_dir, vstr)),
        join(sigdiff_dir, "{}/{}.sig_diffs_emissions.txt".format(sigdiff_dir, vstr))
    ]
    
    gcc_vs_gcc_sigdiff[mon_yr_str] = sigdiff_files

##############################################################################
# For now we are only developing the GCC vs. GCC 1-yr benchmark code
# so leave this commented out (bmy, 10/15/19)
#vstr = '{}_vs_{}'{}.format(gchp_vs_gcc_refstr, gchp_vs_gcc_devstr)
#gchp_vs_gcc_sigdiff = [
#    join(gchp_vs_gcc_plotsdir, '{}_sig_diffs_sfc.txt'.format(vstr)),
#    join(gchp_vs_gcc_plotsdir, '{}_sig_diffs_500hpa.txt'.format(vstr)),
#    join(gchp_vs_gcc_plotsdir, '{}_sig_diffs_zonalmean.txt'.format(vstr)),
#    join(gchp_vs_gcc_plotsdir, '{}_sig_diffs_emissions.txt'.format(vstr))]
#
#vstr = '{}_vs_{}'.format(gchp_vs_gchp_refstr, gchp_vs_gchp_devstr)
#gchp_vs_gchp_sigdiff = [
#    join(gchp_vs_gchp_plotsdir, '{}_sig_diffs_sfc.txt'.format(vstr)),
#    join(gchp_vs_gchp_plotsdir, '{}_sig_diffs_500hpa.txt'.format(vstr)),
#    join(gchp_vs_gchp_plotsdir, '{}_sig_diffs_zonalmean.txt'.format(vstr)),
#    join(gchp_vs_gchp_plotsdir, '{}_sig_diffs_emissions.txt').format(vstr)]
##############################################################################

# ======================================================================
# Create GCC vs GCC benchmark plots and tables
# ======================================================================

if gcc_vs_gcc:

    if plot_conc:
        # --------------------------------------------------------------
        # GCC vs GCC Concentration plots
        # (FullChemBenchmark or TransportTracersBenchmark)
        #
        # FullChemBenchmark includes lumped species and
        # and also separates plots by category
        # --------------------------------------------------------------
        title = '\n%%% Creating GCC vs. GCC {} concentration plots %%%'.format(
            bmk_type)
        print(title)
   
        # File lists for emissions data (seasonal)
        collection = 'SpeciesConc'
        gcc_vs_gcc_refspc = get_filepaths(gcc_vs_gcc_refdir, collection,
                                          bmk_seasons, is_gcc=True)
        gcc_vs_gcc_devspc = get_filepaths(gcc_vs_gcc_devdir, collection,
                                          bmk_seasons, is_gcc=True)

        # Only plot concentration categories for TransportTracers
        if 'TransportTracers' in bmk_type:
            restrict_cats = ['RnPbBeTracers', 'PassiveTracers']
        else:
            restrict_cats = []  # Add FullChemBenchmark restrictions here

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
                                     benchmark_type=bmk_type,
                                     collection=collection,
                                     restrict_cats=restrict_cats,
                                     overwrite=True,
                                     sigdiff_files=sigdiff_files)

    if plot_emis and "FullChem" in bmk_type:
        # --------------------------------------------------------------
        # GCC vs GCC emissions plots
        # (FullChemBenchmark only)
        # --------------------------------------------------------------
        title = '\n%%% Creating GCC vs. GCC {} emissions plots %%%'.format(
            bmk_type)
        print(title)

        # File lists for emissions data (seasonal)
        gcc_vs_gcc_refhco = get_filepaths(gcc_vs_gcc_refdir, 'Emissions',
                                          bmk_seasons, is_gcc=True)
        gcc_vs_gcc_devhco = get_filepaths(gcc_vs_gcc_devdir, 'Emissions',
                                          bmk_seasons, is_gcc=True)

        # Create seasonal emissions plots
        for s in range(bmk_nseasons):
            mon_yr_str = bmk_seasons_names[s]
            sigdiff_files= gcc_vs_gcc_sigdiff[mon_yr_str]
            bmk.make_benchmark_emis_plots(gcc_vs_gcc_refhco[s],
                                          gcc_vs_gcc_refstr,
                                          gcc_vs_gcc_devhco[s],
                                          gcc_vs_gcc_devstr,
                                          dst=gcc_vs_gcc_plotsdir,
                                          subdst=mon_yr_str,
                                          plot_by_benchmark_cat=True,
                                          plot_by_hco_cat=True,
                                          overwrite=True,
                                          sigdiff_files=sigdiff_files)

    if emis_table and "FullChem" in bmk_type:
        # --------------------------------------------------------------
        # GCC vs GCC tables of emission and inventory totals
        # (FullChemBenchmark only)
        # --------------------------------------------------------------
        title = '\n%%% Creating GCC vs. GCC {} emissions and inventory totals %%%'.format(bmk_type)
        print(title)

        # File lists for J-values data (monthly)
        gcc_vs_gcc_refhco = get_filepaths(gcc_vs_gcc_refdir, 'Emissions',
                                          bmk_months, is_gcc=True)
        gcc_vs_gcc_devhco = get_filepaths(gcc_vs_gcc_devdir, 'Emissions',
                                          bmk_months, is_gcc=True)

        # Create emission tables
        bmk.make_benchmark_emis_tables(gcc_vs_gcc_refhco,
                                       gcc_vs_gcc_refstr,
                                       gcc_vs_gcc_devhco,
                                       gcc_vs_gcc_devstr,
                                       dst=gcc_vs_gcc_plotsdir,
                                       interval=sec_per_month,
                                       overwrite=True)

    if plot_jvalues and "FullChem" in bmk_type:
        # --------------------------------------------------------------
        # GCC vs GCC J-value plots
        # (FullChemBenchmark only)
        # --------------------------------------------------------------
        title = '\n%%% Creating GCC vs. GCC {} J-value plots %%%'.format(
            bmk_type)
        print(title)

        # Paths to J-value data (seasonal)
        gcc_vs_gcc_refjv = get_filepaths(gcc_vs_gcc_refdir, 'JValues',
                                          bmk_seasons, is_gcc=True)
        gcc_vs_gcc_devjv = get_filepaths(gcc_vs_gcc_devdir, 'JValues',
                                          bmk_seasons, is_gcc=True)

        # Create seasonal J-values plots
        for s in range(bmk_nseasons):
            mon_yr_str = bmk_seasons_names[s]
            sigdiff_files= gcc_vs_gcc_sigdiff[mon_yr_str]
            bmk.make_benchmark_jvalue_plots(gcc_vs_gcc_refjv[s],
                                            gcc_vs_gcc_refstr,
                                            gcc_vs_gcc_devjv[s],
                                            gcc_vs_gcc_devstr,
                                            dst=gcc_vs_gcc_plotsdir,
                                            subdst=mon_yr_str,
                                            overwrite=True,
                                            sigdiff_files=sigdiff_files)

    if plot_aod and "FullChem" in bmk_type:
        # --------------------------------------------------------------
        # GCC vs. GCC column AOD plots
        # (FullChemBenchmark only)
        # --------------------------------------------------------------
        title = '\n%%% Creating GCC vs. GCC {} column AOD plots %%%'.format(
            bmk_type)
        print(title)

        ## Paths to aerosol optical depth data
        gcc_vs_gcc_refaod = get_filepaths(gcc_vs_gcc_refdir, 'Aerosols',
                                          bmk_seasons, is_gcc=True)
        gcc_vs_gcc_devaod = get_filepaths(gcc_vs_gcc_devdir, 'Aerosols',
                                          bmk_seasons, is_gcc=True)

        # Create seasonal column AOD plots
        for s in range(bmk_nseasons):
            mon_yr_str = bmk_seasons_names[s]
            sigdiff_files= gcc_vs_gcc_sigdiff[mon_yr_str]
            bmk.make_benchmark_aod_plots(gcc_vs_gcc_refaod[s],
                                         gcc_vs_gcc_refstr,
                                         gcc_vs_gcc_devaod[s],
                                         gcc_vs_gcc_devstr,
                                         dst=gcc_vs_gcc_plotsdir,
                                         subdst=mon_yr_str,
                                         overwrite=True,
                                         sigdiff_files=sigdiff_files)

    if mass_table and "FullChem" in bmk_type:
        # --------------------------------------------------------------
        # GCC vs GCC budgets tables
        # (FullChemBenchmark)
        # --------------------------------------------------------------
        title = '\n%%% Creating GCC vs. GCC {} mass tables %%%'.format(
            bmk_type)
        print(title)

        ## Paths to restart files
        gcc_vs_gcc_refrst = get_filepaths(gcc_vs_gcc_refrstdir, 'Restart',
                                          bmk_seasons, is_gcc=True)
        gcc_vs_gcc_devrst = get_filepaths(gcc_vs_gcc_devrstdir, 'Restart',
                                          bmk_seasons, is_gcc=True)

        # Create seasonal budget tables (mass at end of each season month)
        seasons = ['Jan', 'Apr', 'Jul', 'Oct']
        table_dir = join(gcc_vs_gcc_plotsdir, 'Tables')
        for s in range(bmk_nseasons):
            mon_yr_str = seasons[s]
            bmk.make_benchmark_mass_tables(gcc_vs_gcc_refrst[s],
                                           gcc_vs_gcc_refstr,
                                           gcc_vs_gcc_devrst[s],
                                           gcc_vs_gcc_devstr,
                                           dst=table_dir,
                                           overwrite=True,
                                           subdst=mon_yr_str)

    if budget_table and "FullChem" in bmk_type:
        # --------------------------------------------------------------
        # GCC vs GCC budgets tables
        # (FullChemBenchmark)
        # --------------------------------------------------------------
        title = '\n%%% Creating GCC vs. GCC {} budgets %%%'.format(bmk_type)
        print(title)

        ## Paths to budget data
        gcc_vs_gcc_devbgt = get_filepaths(gcc_vs_gcc_devdir, 'Budget',
                                          bmk_seasons, is_gcc=True)
        seasons = ['Jan', 'Apr', 'Jul', 'Oct']
        #Create seasonal budget tables (mass at end of each season month)
        for s in range(bmk_nseasons):
            mon_yr_str = seasons[s]
            bmk.make_benchmark_budget_tables(gcc_vs_gcc_devbgt[s],
                                             gcc_vs_gcc_devstr,
                                             dst=gcc_vs_gcc_plotsdir,
                                             overwrite=True,
                                             interval=sec_per_month[s*3],
                                             subdst=mon_yr_str)

        # Compute annual mean AOD budgets and aerosol burdens
        aerbdg.aerosol_budgets_and_burdens(gcc_dev_version,
                                           gcc_vs_gcc_devdir,
                                           gcc_vs_gcc_plotsdir,
                                           bmk_year,
                                           overwrite=True)

    if budget_table and "TransportTracers" in bmk_type:
        # --------------------------------------------------------------
        # GCC vs GCC budgets tables
        # (TransportTracers)
        # --------------------------------------------------------------
        title = '\n%%% Creating GCC vs. GCC {} budgets %%%'.format(bmk_type)
        print(title)

        # Budgets of Radionuclide species
        ttbdg.transport_tracers_budgets(maindir,
                                        gcc_dev_version,
                                        gcc_vs_gcc_plotsdir,
                                        bmk_year,
                                        overwrite=True)

        # Change in mass of species after each operation
        species = ["Rn222", "Pb210", "Pb210Strat", "Be7", "Be7Strat",
                   "Be10", "Be10Strat", "PassiveTracer", "SF6Tracer",
                   "CH3ITracer", "COAnthroEmis25dayTracer",
                   "COAnthroEmis50dayTracer", "COUniformEmis25dayTracer",
                   "GlobEmis90dayTracer", "NHEmis90dayTracer",
                   "SHEmis90dayTracer"]
        opbdg.make_operations_budget_table(gcc_dev_version,
                                           gcc_vs_gcc_devdir,
                                           gcc_vs_gcc_plotsdir,
                                           bmk_type,
                                           bmk_year,
                                           species,
                                           overwrite=True)

    if OH_metrics and "FullChem" in bmk_type:
        # --------------------------------------------------------------
        # GCC vs GCC Global mean OH, MCF Lifetime, CH4 Lifetime
        # (FullChemBenchmark only)
        # --------------------------------------------------------------
        title = '\n%%% Creating GCC vs. GCC {} OH metrics %%%'.format(bmk_type)
        print(title)
        gcc_vs_gcc_reflist = [gcc_vs_gcc_refcac, gcc_vs_gcc_refmet]
        gcc_vs_gcc_devlist = [gcc_vs_gcc_devcac, gcc_vs_gcc_devmet]
        bmk.make_benchmark_oh_metrics(gcc_vs_gcc_reflist,
                                      gcc_vs_gcc_refstr,
                                      gcc_vs_gcc_devlist,
                                      gcc_vs_gcc_devstr,
                                      dst=gcc_vs_gcc_plotsdir,
                                      overwrite=True)

    if plot_wetdep and "TransportTracers" in bmk_type:
        # --------------------------------------------------------------
        # GCC vs GCC wet deposition
        # (TransportTracers only)
        # --------------------------------------------------------------
        print('\n%%% Creating GCC vs. GCC wet deposition plots %%%')
                
        # Loop over wet deposition collections
        collection_list = ['WetLossConv', 'WetLossLS']
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
                                         overwrite=True,
                                         benchmark_type=bmk_type,
                                         collection=collection,
                                         restrict_cats=[collection],
                                         sigdiff_files=sigdiff_files)

    if ste_table:
        # --------------------------------------------------------------
        # GCC Strat-Trop Exchange
        # --------------------------------------------------------------
        title = '\n%%% Creating GCC vs. GCC {} Strat-Trop Exchange table %%%'.format(
            bmk_type)
        print(title)
        
        # Pick the species list depending on the benchmark type
        if "TransportTracers" in bmk_type:
            species = ["Pb210", "Be7", "Be10"]
        elif "FullChem" in bmk_type:
            species = ["O3"]
        
        # Strat-trop exchange
        ste.make_benchmark_ste_table(gcc_dev_version,
                                     gcc_vs_gcc_devdir,
                                     gcc_vs_gcc_plotsdir,
                                     bmk_year,
                                     species=species,
                                     overwrite=True)
    
###############################################################################
# NOTE: For now, we are developing the GCC vs. GCC 1-yr benchmark plots.
# Leave this commented out for now.
## =====================================================================
## Create GCHP vs GCC benchmark plots and tables
## =====================================================================
#
#if gchp_vs_gcc:
#    if plot_conc:
#        # Concentration plots
#        # (includes lumped species and separates by category)
#        print('\n%%% Creating GCHP vs. GCC concentration plots %%%')
#        bmk.make_benchmark_conc_plots(gchp_vs_gcc_refspc,
#                                      gchp_vs_gcc_refstr,
#                                      gchp_vs_gcc_devspc,
#                                      gchp_vs_gcc_devstr,
#                                      dst=gchp_vs_gcc_plotsdir,
#                                      overwrite=True,
#                                      sigdiff_files=gchp_vs_gcc_sigdiff)
#
#    if plot_emis:
#        # Emissions plots
#        print('\n%%% Creating GCHP vs. GCC emissions plots %%%')
#        bmk.make_benchmark_emis_plots(gchp_vs_gcc_refhco,
#                                      gchp_vs_gcc_refstr,
#                                      gchp_vs_gcc_devhco,
#                                      gchp_vs_gcc_devstr,
#                                      dst=gchp_vs_gcc_plotsdir,
#                                      plot_by_benchmark_cat=True,
#                                      plot_by_hco_cat=True,
#                                      overwrite=True,
#                                      flip_dev=True,
#                                      sigdiff_files=gchp_vs_gcc_sigdiff)
#
#    if emis_table:
#        # Tables of emissions and inventory totals
#        print('\n%%% Creating GCHP vs. GCC emissions and inventory tables %%%')
#        gchp_vs_gcc_reflist = [gchp_vs_gcc_refhco]
#        gchp_vs_gcc_devlist = [gchp_vs_gcc_devhco, gchp_vs_gcc_devmet]
#        bmk.make_benchmark_emis_tables(gchp_vs_gcc_reflist,
#                                       gchp_vs_gcc_refstr,
#                                       gchp_vs_gcc_devlist,
#                                       gchp_vs_gcc_devstr,
#                                       dst=gchp_vs_gcc_plotsdir,
#                                       overwrite=True)
#
#    if plot_jvalues:
#        # Local noon J-values plots
#        print('\n%%% Creating GCHP vs. GCC J-value plots %%%')
#        bmk.make_benchmark_jvalue_plots(gchp_vs_gcc_refjv,
#                                        gchp_vs_gcc_refstr,
#                                        gchp_vs_gcc_devjv,
#                                        gchp_vs_gcc_devstr,
#                                        dst=gchp_vs_gcc_plotsdir,
#                                        overwrite=True,
#                                        sigdiff_files=gchp_vs_gcc_sigdiff)
#
#    if plot_aod:
#        # Column AOD plots
#        print('\n%%% Creating GCHP vs. GCC column AOD plots %%%')
#        bmk.make_benchmark_aod_plots(gchp_vs_gcc_refaod,
#                                     gchp_vs_gcc_refstr,
#                                     gchp_vs_gcc_devaod,
#                                     gchp_vs_gcc_devstr,
#                                     dst=gchp_vs_gcc_plotsdir,
#                                     overwrite=True,
#                                     sigdiff_files=gchp_vs_gcc_sigdiff)
#
## Under development, leave commented out for now (bmy, 9/5/19)
##    if mass_table:
##        # Global mass tables
##        print('\n%%% Creating GCHP vs. GCC global mass tables %%%')
##        gchp_vs_gcc_reflist = [gchp_vs_gcc_refrst]
##        gchp_vs_gcc_devlist = [gchp_vs_gcc_devrst, gchp_vs_gcc_devmetinst]
##        bmk.make_benchmark_mass_tables(gchp_vs_gcc_reflist,
##                                       gchp_vs_gcc_refstr,
##                                       gchp_vs_gcc_devlist,
##                                       gchp_vs_gcc_devstr,
##                                       dst=gchp_vs_gcc_plotsdir,
##                                       overwrite=True)
#        
#    if OH_metrics:
#        # Global mean OH, MCF Lifetime, CH4 Lifetime
#        print('\n%%% Creating GCHP vs. GCC OH metrics %%%')
#        gchp_vs_gcc_reflist = [gchp_vs_gcc_refcac, gchp_vs_gcc_refmet]
#        gchp_vs_gcc_devlist = [gchp_vs_gcc_devcac, gchp_vs_gcc_devmet]
#        bmk.make_benchmark_oh_metrics(gchp_vs_gcc_reflist,
#                                      gchp_vs_gcc_refstr,
#                                      gchp_vs_gcc_devlist,
#                                      gchp_vs_gcc_devstr,
#                                      dst=gchp_vs_gcc_plotsdir,
#                                      overwrite=True)
#        
## =====================================================================
## Create GCHP vs GCHP benchmark plots and tables
## =====================================================================
#
#if gchp_vs_gchp:
#
#    if plot_conc:
#        # Concentration plots
#        # (includes lumped species and separates by category)
#        print('\n%%% Creating GCHP vs. GCHP concentration plots %%%')
#        bmk.make_benchmark_conc_plots(gchp_vs_gchp_refspc,
#                                      gchp_vs_gchp_refstr,
#                                      gchp_vs_gchp_devspc,
#                                      gchp_vs_gchp_devstr,
#                                      dst=gchp_vs_gchp_plotsdir,
#                                      overwrite=True,
#                                      sigdiff_files=gchp_vs_gchp_sigdiff)
#
#    if plot_emis:
#        # Emissions plots
#        print('\n%%% Creating GCHP vs. GCHP emissions plots %%%')
#        bmk.make_benchmark_emis_plots(gchp_vs_gchp_refhco,
#                                      gchp_vs_gchp_refstr,
#                                      gchp_vs_gchp_devhco,
#                                      gchp_vs_gchp_devstr,
#                                      dst=gchp_vs_gchp_plotsdir,
#                                      plot_by_benchmark_cat=True,
#                                      plot_by_hco_cat=True,
#                                      overwrite=True,
#                                      flip_ref=True,
#                                      flip_dev=True,
#                                      sigdiff_files=gchp_vs_gchp_sigdiff)
#
#    if emis_table:
#        # Tables of emissions and inventory totals
#        print('\n%%% Creating GCHP vs. GCHP emissions and inventory tables %%%')
#        gchp_vs_gchp_reflist = [gchp_vs_gchp_refhco, gchp_vs_gchp_refmet]
#        gchp_vs_gchp_devlist = [gchp_vs_gchp_devhco, gchp_vs_gchp_devmet]
#        bmk.make_benchmark_emis_tables(gchp_vs_gchp_reflist,
#                                       gchp_vs_gchp_refstr,
#                                       gchp_vs_gchp_devlist,
#                                       gchp_vs_gchp_devstr,
#                                       dst=gchp_vs_gchp_plotsdir,
#                                       overwrite=True)
#
#    if plot_jvalues:
#        # Local noon J-values plots
#        print('\n%%% Creating GCHP vs. GCHP J-value plots %%%')
#        bmk.make_benchmark_jvalue_plots(gchp_vs_gchp_refjv,
#                                        gchp_vs_gchp_refstr,
#                                        gchp_vs_gchp_devjv,
#                                        gchp_vs_gchp_devstr,
#                                        dst=gchp_vs_gchp_plotsdir,
#                                        overwrite=True,
#                                        sigdiff_files=gchp_vs_gchp_sigdiff)
#
#    if plot_aod:
#        # Column AOD plots
#        print('\n%%% Creating GCHP vs. GCHP column AOD plots %%%')
#        bmk.make_benchmark_aod_plots(gchp_vs_gchp_refaod,
#                                     gchp_vs_gchp_refstr,
#                                     gchp_vs_gchp_devaod,
#                                     gchp_vs_gchp_devstr,
#                                     dst=gchp_vs_gchp_plotsdir,
#                                     overwrite=True,
#                                     sigdiff_files=gchp_vs_gchp_sigdiff)
#
## Under development, leave commented out for now (bmy, 9/5/19)
##    if mass_table:
##        # Global mass tables
##        print('\n%%% Creating GCHP vs. GCHP global mass tables %%%')
##        gchp_vs_gchp_reflist = [gchp_vs_gchp_refrst, gchp_vs_gchp_refmetinst]
##        gchp_vs_gchp_devlist = [gchp_vs_gchp_devrst, gchp_vs_gchp_devmetinst]
##        bmk.make_benchmark_mass_tables(gchp_vs_gchp_reflist,
##                                       gchp_vs_gchp_refstr,
##                                       gchp_vs_gchp_devlist,
##                                       gchp_vs_gchp_devstr,
##                                       dst=gchp_vs_gchp_plotsdir,
##                                       overwrite=True)
#
#    if OH_metrics:
#        # Global mean OH, MCF Lifetime, CH4 Lifetime
#        print('\n%%% Creating GCHP vs. GCHP OH metrics %%%')
#        gchp_vs_gcc_reflist = [gchp_vs_gchp_refcac, gchp_vs_gchp_refmet]
#        gchp_vs_gcc_devlist = [gchp_vs_gchp_devcac, gchp_vs_gchp_devmet]
#        bmk.make_benchmark_oh_metrics(gchp_vs_gchp_reflist,
#                                      gchp_vs_gchp_refstr,
#                                      gchp_vs_gchp_devlist,
#                                      gchp_vs_gchp_devstr,
#                                      dst=gchp_vs_gchp_plotsdir,
#                                      overwrite=True)
#
## =====================================================================
## Create GCHP vs GCC difference of differences benchmark plots
## =====================================================================
#
#if gchp_vs_gcc_diff_of_diffs:
#
#    # NOTE: This can be expanded to differences beyond species
#    # concentrations by following how this is done for conc plots.
#    print('%\n%% Creating GCHP vs. GCC diff-of-diffs concentration plots %%%')
#
#    # Get a list of variables that GCPy should not read
#    skip_vars = skip_these_vars
#
#    # Target output files
#    diff_of_diffs_refspc = './gcc_diffs_spc.nc4'
#    diff_of_diffs_devspc = './gchp_diffs_spc.nc4'
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
#    if 'lat' in refdims and 'Xdim' in devdims:
#        gchp_ref_newdimnames = gchp_dev.copy()
#        for v in gchp_dev.data_vars.keys():
#            if 'Xdim' in gchp_dev[v].dims:
#                gchp_ref_newdimnames[v].values = gchp_ref[v].values.reshape(
#                    gchp_dev[v].values.shape)
#                # NOTE: the reverse conversion is gchp_dev[v].stack(lat=('nf','Ydim')).transpose('time','lev','lat','Xdim').values
#        gchp_ref = gchp_ref_newdimnames.copy()
#    with xr.set_options(keep_attrs=True):
#        gchp_diffs = gchp_dev.copy()
#        for v in gchp_dev.data_vars.keys():
#            if 'Xdim' in gchp_dev[v].dims or 'lat' in gchp_dev[v].dims:
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
