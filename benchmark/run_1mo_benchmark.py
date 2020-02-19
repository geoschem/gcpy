#!/usr/bin/env python
'''
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
        os.environ['QT_QPA_PLATFORM']='offscreen'

    For more information, please see this issue posted at the ipython site:

        https://github.com/ipython/ipython/issues/10627

    This issue might be fixed in matplotlib 3.0.
'''

# =====================================================================
# Imports and global settings (you should not need to edit these)
# =====================================================================

import os
from os.path import join
import xarray as xr
from gcpy import benchmark as bmk
from gcpy.constants import skip_these_vars
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
gcc_vs_gcc   = True
gchp_vs_gcc  = False
gchp_vs_gchp = False
gchp_vs_gcc_diff_of_diffs = False

# Output to generate (edit as needed)
# Plots/tables will be created in this order:
plot_conc    = True
plot_emis    = True
emis_table   = True
plot_jvalues = True
plot_aod     = True
mass_table   = True
budget_table = True
OH_metrics   = True

# Filename date strings (edit as needed)
gcc_datestr  = '20160701'
gchp_datestr = '20160716'
end_datestr  = '20160801'

# Filename hour strings (edit as needed)
gcc_hourstr  = '0000'
gchp_hourstr = '1200'

# Seconds in the month (July 2016)
sec_per_month = [86400.0 * 31.0]

# Data directories (edit as needed)
# For gchp_vs_gcc_refdir use gcc_dev_version, not ref (mps, 6/27/19)
gcc_vs_gcc_refdir   = join(maindir, gcc_ref_version,  'OutputDir')
gcc_vs_gcc_devdir   = join(maindir, gcc_dev_version,  'OutputDir')
gchp_vs_gcc_refdir  = join(maindir, gcc_dev_version,  'OutputDir')
gchp_vs_gcc_devdir  = join(maindir, gchp_dev_version, 'OutputDir')
gchp_vs_gchp_refdir = join(maindir, gchp_ref_version, 'OutputDir')
gchp_vs_gchp_devdir = join(maindir, gchp_dev_version, 'OutputDir')

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

# Files that will contain lists of quantities that have significant
# differences -- we need these for the benchmark approval forms.
vstr = '{}_vs_{}'.format(gcc_vs_gcc_refstr, gcc_vs_gcc_devstr)
gcc_vs_gcc_sigdiff = [
    join(gcc_vs_gcc_plotsdir, '{}_sig_diffs_sfc.txt'.format(vstr)),
    join(gcc_vs_gcc_plotsdir, '{}_sig_diffs_500hpa.txt'.format(vstr)),
    join(gcc_vs_gcc_plotsdir, '{}_sig_diffs_zonalmean.txt'.format(vstr)),
    join(gcc_vs_gcc_plotsdir, '{}_sig_diffs_emissions.txt'.format(vstr))]

vstr = '{}_vs_{}'.format(gchp_vs_gcc_refstr, gchp_vs_gcc_devstr)
gchp_vs_gcc_sigdiff = [
    join(gchp_vs_gcc_plotsdir, '{}_sig_diffs_sfc.txt'.format(vstr)),
    join(gchp_vs_gcc_plotsdir, '{}_sig_diffs_500hpa.txt'.format(vstr)),
    join(gchp_vs_gcc_plotsdir, '{}_sig_diffs_zonalmean.txt'.format(vstr)),
    join(gchp_vs_gcc_plotsdir, '{}_sig_diffs_emissions.txt'.format(vstr))]

vstr = '{}_vs_{}'.format(gchp_vs_gchp_refstr, gchp_vs_gchp_devstr)
gchp_vs_gchp_sigdiff = [
    join(gchp_vs_gchp_plotsdir, '{}_sig_diffs_sfc.txt'.format(vstr)),
    join(gchp_vs_gchp_plotsdir, '{}_sig_diffs_500hpa.txt'.format(vstr)),
    join(gchp_vs_gchp_plotsdir, '{}_sig_diffs_zonalmean.txt'.format(vstr)),
    join(gchp_vs_gchp_plotsdir, '{}_sig_diffs_emissions.txt').format(vstr)]

# =====================================================================
# The rest of these settings should not need to be changed
# =====================================================================

# Species concentration filenames
gcc_spcfile  = 'GEOSChem.SpeciesConc.{}_{}z.nc4'.format(gcc_datestr,
                                                        gcc_hourstr)   
gchp_spcfile = 'GCHP.SpeciesConc.{}_{}z.nc4'.format(gchp_datestr,
                                                    gchp_hourstr) 

# HEMCO diagnostic filenames
gcc_hcofile  = 'HEMCO_diagnostics.{}{}.nc'.format(gcc_datestr, gcc_hourstr)
gchp_hcofile = 'GCHP.Emissions.{}_{}z.nc4'.format(gchp_datestr, gchp_hourstr)

# 24-hr avg J-value diagnostic filenames
gcc_jvfile  = 'GEOSChem.JValues.{}_{}z.nc4'.format(gcc_datestr, gcc_hourstr)
gchp_jvfile = 'GCHP.JValues.{}_{}z.nc4'.format(gchp_datestr, gchp_hourstr)

# Aerosol optical depth diagnostic filenames
gcc_aodfile  = 'GEOSChem.Aerosols.{}_{}z.nc4'.format(gcc_datestr, gcc_hourstr)
gchp_aodfile = 'GCHP.Aerosols.{}_{}z.nc4'.format(gchp_datestr, gchp_hourstr)

# StateMet diagnostic filenames
gcc_metfile  = 'GEOSChem.StateMet.{}_{}z.nc4'.format(gcc_datestr,
                                                     gcc_hourstr)
gchp_metfile = 'GCHP.StateMet_avg.{}_{}z.nc4'.format(gchp_datestr,
                                                     gchp_hourstr)
gchp_metfileinst = 'GCHP.StateMet_inst.{}_{}z.nc4'.format(end_datestr,
                                                          gcc_hourstr)

# Restart filenames
gcc_rstfile  = 'GEOSChem.Restart.{}_{}z.nc4'.format(end_datestr,
                                                    gcc_hourstr)
gchp_rstfile = 'gcchem_internal_checkpoint.restart.{}_000000.nc4'.format(end_datestr)

# Budget diagnostic filenames
gcc_bgtfile = 'GEOSChem.Budget.{}_{}z.nc4'.format(gcc_datestr,
                                                  gcc_hourstr)
gchp_bgtfile = 'GCHP.Budget.{}_{}z.nc4'.format(gchp_datestr,
                                               gchp_hourstr)

# Concentration after chemistry filename (for OH metrics)
gcc_cacfile = 'GEOSChem.ConcAfterChem.{}_{}z.nc4'.format(gcc_datestr,
                                                         gcc_hourstr)
gchp_cacfile = 'GCHP.ConcAfterChem.{}_{}z.nc4'.format(gchp_datestr,
                                                      gchp_hourstr)

# Paths to species concentration data
gcc_vs_gcc_refspc   = join(maindir, gcc_vs_gcc_refdir,   gcc_spcfile)
gcc_vs_gcc_devspc   = join(maindir, gcc_vs_gcc_devdir,   gcc_spcfile)
gchp_vs_gcc_refspc  = join(maindir, gchp_vs_gcc_refdir,  gcc_spcfile)
gchp_vs_gcc_devspc  = join(maindir, gchp_vs_gcc_devdir,  gchp_spcfile)
gchp_vs_gchp_refspc = join(maindir, gchp_vs_gchp_refdir, gchp_spcfile)
gchp_vs_gchp_devspc = join(maindir, gchp_vs_gchp_devdir, gchp_spcfile)

# Paths to HEMCO diagnostics data
gcc_vs_gcc_refhco   = join(maindir, gcc_vs_gcc_refdir,   gcc_hcofile)
gcc_vs_gcc_devhco   = join(maindir, gcc_vs_gcc_devdir,   gcc_hcofile)
gchp_vs_gcc_refhco  = join(maindir, gchp_vs_gcc_refdir,  gcc_hcofile)
gchp_vs_gcc_devhco  = join(maindir, gchp_vs_gcc_devdir,  gchp_hcofile)
gchp_vs_gchp_refhco = join(maindir, gchp_vs_gchp_refdir, gchp_hcofile)
gchp_vs_gchp_devhco = join(maindir, gchp_vs_gchp_devdir, gchp_hcofile)

# Paths to local noon J-value data
gcc_vs_gcc_refjv    = join(maindir, gcc_vs_gcc_refdir,   gcc_jvfile)
gcc_vs_gcc_devjv    = join(maindir, gcc_vs_gcc_devdir,   gcc_jvfile)
gchp_vs_gcc_refjv   = join(maindir, gchp_vs_gcc_refdir,  gcc_jvfile)
gchp_vs_gcc_devjv   = join(maindir, gchp_vs_gcc_devdir,  gchp_jvfile)
gchp_vs_gchp_refjv  = join(maindir, gchp_vs_gchp_refdir, gchp_jvfile)
gchp_vs_gchp_devjv  = join(maindir, gchp_vs_gchp_devdir, gchp_jvfile)

# Paths to aerosol optical depth data
gcc_vs_gcc_refaod   = join(maindir, gcc_vs_gcc_refdir,   gcc_aodfile)
gcc_vs_gcc_devaod   = join(maindir, gcc_vs_gcc_devdir,   gcc_aodfile)
gchp_vs_gcc_refaod  = join(maindir, gchp_vs_gcc_refdir,  gcc_aodfile)
gchp_vs_gcc_devaod  = join(maindir, gchp_vs_gcc_devdir,  gchp_aodfile)
gchp_vs_gchp_refaod = join(maindir, gchp_vs_gchp_refdir, gchp_aodfile)
gchp_vs_gchp_devaod = join(maindir, gchp_vs_gchp_devdir, gchp_aodfile)

# Paths to StateMet data
gcc_vs_gcc_refmet       = join(maindir, gcc_vs_gcc_refdir,   gcc_metfile)
gcc_vs_gcc_devmet       = join(maindir, gcc_vs_gcc_devdir,   gcc_metfile)
gchp_vs_gcc_refmet      = join(maindir, gchp_vs_gcc_refdir,  gcc_metfile)
gchp_vs_gcc_devmet      = join(maindir, gchp_vs_gcc_devdir,  gchp_metfile)
gchp_vs_gcc_devmetinst  = join(maindir, gchp_vs_gcc_devdir,  gchp_metfileinst)
gchp_vs_gchp_refmet     = join(maindir, gchp_vs_gchp_refdir, gchp_metfile)
gchp_vs_gchp_refmetinst = join(maindir, gchp_vs_gchp_refdir, gchp_metfileinst)
gchp_vs_gchp_devmet     = join(maindir, gchp_vs_gchp_devdir, gchp_metfile)
gchp_vs_gchp_devmetinst = join(maindir, gchp_vs_gchp_devdir, gchp_metfileinst)

# Paths to restart files
gcc_vs_gcc_refrst   = join(maindir, gcc_ref_version,  gcc_rstfile)
gcc_vs_gcc_devrst   = join(maindir, gcc_dev_version,  gcc_rstfile)
gchp_vs_gcc_refrst  = join(maindir, gcc_dev_version,  gcc_rstfile)
gchp_vs_gcc_devrst  = join(maindir, gchp_dev_version, gchp_rstfile)
gchp_vs_gchp_refrst = join(maindir, gchp_ref_version, gchp_rstfile)
gchp_vs_gchp_devrst = join(maindir, gchp_dev_version, gchp_rstfile)

# Paths to budget data
gcc_vs_gcc_devbgt   = join(maindir, gcc_vs_gcc_devdir,   gcc_bgtfile)
gchp_vs_gcc_devbgt  = join(maindir, gchp_vs_gcc_devdir,  gchp_bgtfile)
gchp_vs_gchp_devbgt = join(maindir, gchp_vs_gchp_devdir, gchp_bgtfile)

# Paths to concentration after chemistry files
gcc_vs_gcc_refcac   = join(maindir, gcc_vs_gcc_refdir,   gcc_cacfile)
gcc_vs_gcc_devcac   = join(maindir, gcc_vs_gcc_devdir,   gcc_cacfile)
gchp_vs_gcc_refcac  = join(maindir, gchp_vs_gcc_refdir,  gcc_cacfile)
gchp_vs_gcc_devcac  = join(maindir, gchp_vs_gcc_devdir,  gchp_cacfile)
gchp_vs_gchp_refcac = join(maindir, gchp_vs_gchp_refdir, gchp_cacfile)
gchp_vs_gchp_devcac = join(maindir, gchp_vs_gchp_devdir, gchp_cacfile)

# =====================================================================
# Create GCC vs GCC benchmark plots and tables
# =====================================================================

if gcc_vs_gcc:

    if plot_conc:
        # Concentration plots
        # (includes lumped species and separates by category)
        print('\n%%% Creating GCC vs. GCC concentration plots %%%')
        bmk.make_benchmark_plots(gcc_vs_gcc_refspc,
                                 gcc_vs_gcc_refstr,
                                 gcc_vs_gcc_devspc,
                                 gcc_vs_gcc_devstr,
                                 dst=gcc_vs_gcc_plotsdir,
                                 overwrite=True,
                                 sigdiff_files=gcc_vs_gcc_sigdiff)

    if plot_emis:
        # Emissions plots
        print('\n%%% Creating GCC vs. GCC emissions plots %%%')
        bmk.make_benchmark_emis_plots(gcc_vs_gcc_refhco,
                                      gcc_vs_gcc_refstr,
                                      gcc_vs_gcc_devhco,
                                      gcc_vs_gcc_devstr,
                                      dst=gcc_vs_gcc_plotsdir,
                                      plot_by_benchmark_cat=True,
                                      plot_by_hco_cat=True,
                                      overwrite=True,
                                      sigdiff_files=gcc_vs_gcc_sigdiff)

    if emis_table:
        # Table of emission and inventory totals
        print('\n%%% Creating GCC vs. GCC emissions and inventory tables %%%')
        gcc_vs_gcc_reflist = [gcc_vs_gcc_refhco]
        gcc_vs_gcc_devlist = [gcc_vs_gcc_devhco]
        bmk.make_benchmark_emis_tables(gcc_vs_gcc_reflist,
                                       gcc_vs_gcc_refstr,
                                       gcc_vs_gcc_devlist,
                                       gcc_vs_gcc_devstr,
                                       dst=gcc_vs_gcc_plotsdir,
                                       interval=sec_per_month,
                                       overwrite=True)

    if plot_jvalues:
        # Local noon J-values plots
        print('\n%%% Creating GCC vs. GCC J-value plots %%%')
        bmk.make_benchmark_jvalue_plots(gcc_vs_gcc_refjv,
                                        gcc_vs_gcc_refstr,
                                        gcc_vs_gcc_devjv,
                                        gcc_vs_gcc_devstr,
                                        dst=gcc_vs_gcc_plotsdir,
                                        overwrite=True,
                                        sigdiff_files=gcc_vs_gcc_sigdiff)

    if plot_aod:
        # Column AOD plots
        print('\n%%% Creating GCC vs. GCC column AOD plots %%%')
        bmk.make_benchmark_aod_plots(gcc_vs_gcc_refaod,
                                     gcc_vs_gcc_refstr,
                                     gcc_vs_gcc_devaod,
                                     gcc_vs_gcc_devstr,
                                     dst=gcc_vs_gcc_plotsdir,
                                     overwrite=True,
                                     sigdiff_files=gcc_vs_gcc_sigdiff)

    if mass_table:
        # Global mass tables
        print('\n%%% Creating GCC vs. GCC global mass tables %%%')
        bmk.make_benchmark_mass_tables(gcc_vs_gcc_refrst,
                                       gcc_vs_gcc_refstr,
                                       gcc_vs_gcc_devrst,
                                       gcc_vs_gcc_devstr,
                                       dst=gcc_vs_gcc_plotsdir,
                                       overwrite=True)
        
    if budget_table:
        # Budgets tables
        print('\n%%% Creating GCC vs. GCC budget tables %%%')
        bmk.make_benchmark_budget_tables(gcc_vs_gcc_devbgt,
                                         gcc_vs_gcc_devstr,
                                         dst=gcc_vs_gcc_plotsdir,
                                         interval=sec_per_month,
                                         overwrite=True)

    if OH_metrics:
        # Global mean OH, MCF Lifetime, CH4 Lifetime
        print('\n%%% Creating GCC vs. GCC OH metrics %%%')
        gcc_vs_gcc_reflist = [gcc_vs_gcc_refcac, gcc_vs_gcc_refmet]
        gcc_vs_gcc_devlist = [gcc_vs_gcc_devcac, gcc_vs_gcc_devmet]
        bmk.make_benchmark_oh_metrics(gcc_vs_gcc_reflist,
                                      gcc_vs_gcc_refstr,
                                      gcc_vs_gcc_devlist,
                                      gcc_vs_gcc_devstr,
                                      dst=gcc_vs_gcc_plotsdir,
                                      overwrite=True)
        
# =====================================================================
# Create GCHP vs GCC benchmark plots and tables
# =====================================================================

if gchp_vs_gcc:
    if plot_conc:
        # Concentration plots
        # (includes lumped species and separates by category)
        print('\n%%% Creating GCHP vs. GCC concentration plots %%%')
        bmk.make_benchmark_plots(gchp_vs_gcc_refspc,
                                 gchp_vs_gcc_refstr,
                                 gchp_vs_gcc_devspc,
                                 gchp_vs_gcc_devstr,
                                 dst=gchp_vs_gcc_plotsdir,
                                 overwrite=True,
                                 sigdiff_files=gchp_vs_gcc_sigdiff)

    if plot_emis:
        # Emissions plots
        print('\n%%% Creating GCHP vs. GCC emissions plots %%%')
        bmk.make_benchmark_emis_plots(gchp_vs_gcc_refhco,
                                      gchp_vs_gcc_refstr,
                                      gchp_vs_gcc_devhco,
                                      gchp_vs_gcc_devstr,
                                      dst=gchp_vs_gcc_plotsdir,
                                      plot_by_benchmark_cat=True,
                                      plot_by_hco_cat=True,
                                      overwrite=True,
                                      flip_dev=True,
                                      sigdiff_files=gchp_vs_gcc_sigdiff)

    if emis_table:
        # Tables of emissions and inventory totals
        print('\n%%% Creating GCHP vs. GCC emissions and inventory tables %%%')
        gchp_vs_gcc_reflist = [gchp_vs_gcc_refhco]
        gchp_vs_gcc_devlist = [gchp_vs_gcc_devhco, gchp_vs_gcc_devmet]
        bmk.make_benchmark_emis_tables(gchp_vs_gcc_reflist,
                                       gchp_vs_gcc_refstr,
                                       gchp_vs_gcc_devlist,
                                       gchp_vs_gcc_devstr,
                                       dst=gchp_vs_gcc_plotsdir,
                                       overwrite=True)

    if plot_jvalues:
        # Local noon J-values plots
        print('\n%%% Creating GCHP vs. GCC J-value plots %%%')
        bmk.make_benchmark_jvalue_plots(gchp_vs_gcc_refjv,
                                        gchp_vs_gcc_refstr,
                                        gchp_vs_gcc_devjv,
                                        gchp_vs_gcc_devstr,
                                        dst=gchp_vs_gcc_plotsdir,
                                        overwrite=True,
                                        sigdiff_files=gchp_vs_gcc_sigdiff)

    if plot_aod:
        # Column AOD plots
        print('\n%%% Creating GCHP vs. GCC column AOD plots %%%')
        bmk.make_benchmark_aod_plots(gchp_vs_gcc_refaod,
                                     gchp_vs_gcc_refstr,
                                     gchp_vs_gcc_devaod,
                                     gchp_vs_gcc_devstr,
                                     dst=gchp_vs_gcc_plotsdir,
                                     overwrite=True,
                                     sigdiff_files=gchp_vs_gcc_sigdiff)

# Under development, leave commented out for now (bmy, 9/5/19)
#    if mass_table:
#        # Global mass tables
#        print('\n%%% Creating GCHP vs. GCC global mass tables %%%')
#        gchp_vs_gcc_reflist = [gchp_vs_gcc_refrst]
#        gchp_vs_gcc_devlist = [gchp_vs_gcc_devrst, gchp_vs_gcc_devmetinst]
#        bmk.make_benchmark_mass_tables(gchp_vs_gcc_reflist,
#                                       gchp_vs_gcc_refstr,
#                                       gchp_vs_gcc_devlist,
#                                       gchp_vs_gcc_devstr,
#                                       dst=gchp_vs_gcc_plotsdir,
#                                       overwrite=True)
        
    if OH_metrics:
        # Global mean OH, MCF Lifetime, CH4 Lifetime
        print('\n%%% Creating GCHP vs. GCC OH metrics %%%')
        gchp_vs_gcc_reflist = [gchp_vs_gcc_refcac, gchp_vs_gcc_refmet]
        gchp_vs_gcc_devlist = [gchp_vs_gcc_devcac, gchp_vs_gcc_devmet]
        bmk.make_benchmark_oh_metrics(gchp_vs_gcc_reflist,
                                      gchp_vs_gcc_refstr,
                                      gchp_vs_gcc_devlist,
                                      gchp_vs_gcc_devstr,
                                      dst=gchp_vs_gcc_plotsdir,
                                      overwrite=True)
        
# =====================================================================
# Create GCHP vs GCHP benchmark plots and tables
# =====================================================================

if gchp_vs_gchp:

    if plot_conc:
        # Concentration plots
        # (includes lumped species and separates by category)
        print('\n%%% Creating GCHP vs. GCHP concentration plots %%%')
        bmk.make_benchmark_plots(gchp_vs_gchp_refspc,
                                 gchp_vs_gchp_refstr,
                                 gchp_vs_gchp_devspc,
                                 gchp_vs_gchp_devstr,
                                 dst=gchp_vs_gchp_plotsdir,
                                 overwrite=True,
                                 sigdiff_files=gchp_vs_gchp_sigdiff)

    if plot_emis:
        # Emissions plots
        print('\n%%% Creating GCHP vs. GCHP emissions plots %%%')
        bmk.make_benchmark_emis_plots(gchp_vs_gchp_refhco,
                                      gchp_vs_gchp_refstr,
                                      gchp_vs_gchp_devhco,
                                      gchp_vs_gchp_devstr,
                                      dst=gchp_vs_gchp_plotsdir,
                                      plot_by_benchmark_cat=True,
                                      plot_by_hco_cat=True,
                                      overwrite=True,
                                      flip_ref=True,
                                      flip_dev=True,
                                      sigdiff_files=gchp_vs_gchp_sigdiff)

    if emis_table:
        # Tables of emissions and inventory totals
        print('\n%%% Creating GCHP vs. GCHP emissions and inventory tables %%%')
        gchp_vs_gchp_reflist = [gchp_vs_gchp_refhco, gchp_vs_gchp_refmet]
        gchp_vs_gchp_devlist = [gchp_vs_gchp_devhco, gchp_vs_gchp_devmet]
        bmk.make_benchmark_emis_tables(gchp_vs_gchp_reflist,
                                       gchp_vs_gchp_refstr,
                                       gchp_vs_gchp_devlist,
                                       gchp_vs_gchp_devstr,
                                       dst=gchp_vs_gchp_plotsdir,
                                       overwrite=True)

    if plot_jvalues:
        # Local noon J-values plots
        print('\n%%% Creating GCHP vs. GCHP J-value plots %%%')
        bmk.make_benchmark_jvalue_plots(gchp_vs_gchp_refjv,
                                        gchp_vs_gchp_refstr,
                                        gchp_vs_gchp_devjv,
                                        gchp_vs_gchp_devstr,
                                        dst=gchp_vs_gchp_plotsdir,
                                        overwrite=True,
                                        sigdiff_files=gchp_vs_gchp_sigdiff)

    if plot_aod:
        # Column AOD plots
        print('\n%%% Creating GCHP vs. GCHP column AOD plots %%%')
        bmk.make_benchmark_aod_plots(gchp_vs_gchp_refaod,
                                     gchp_vs_gchp_refstr,
                                     gchp_vs_gchp_devaod,
                                     gchp_vs_gchp_devstr,
                                     dst=gchp_vs_gchp_plotsdir,
                                     overwrite=True,
                                     sigdiff_files=gchp_vs_gchp_sigdiff)

# Under development, leave commented out for now (bmy, 9/5/19)
#    if mass_table:
#        # Global mass tables
#        print('\n%%% Creating GCHP vs. GCHP global mass tables %%%')
#        gchp_vs_gchp_reflist = [gchp_vs_gchp_refrst, gchp_vs_gchp_refmetinst]
#        gchp_vs_gchp_devlist = [gchp_vs_gchp_devrst, gchp_vs_gchp_devmetinst]
#        bmk.make_benchmark_mass_tables(gchp_vs_gchp_reflist,
#                                       gchp_vs_gchp_refstr,
#                                       gchp_vs_gchp_devlist,
#                                       gchp_vs_gchp_devstr,
#                                       dst=gchp_vs_gchp_plotsdir,
#                                       overwrite=True)

    if OH_metrics:
        # Global mean OH, MCF Lifetime, CH4 Lifetime
        print('\n%%% Creating GCHP vs. GCHP OH metrics %%%')
        gchp_vs_gchp_reflist = [gchp_vs_gchp_refcac, gchp_vs_gchp_refmet]
        gchp_vs_gchp_devlist = [gchp_vs_gchp_devcac, gchp_vs_gchp_devmet]
        bmk.make_benchmark_oh_metrics(gchp_vs_gchp_reflist,
                                      gchp_vs_gchp_refstr,
                                      gchp_vs_gchp_devlist,
                                      gchp_vs_gchp_devstr,
                                      dst=gchp_vs_gchp_plotsdir,
                                      overwrite=True)

# =====================================================================
# Create GCHP vs GCC difference of differences benchmark plots
# =====================================================================

if gchp_vs_gcc_diff_of_diffs:

    # NOTE: This can be expanded to differences beyond species
    # concentrations by following how this is done for conc plots.
    print('%\n%% Creating GCHP vs. GCC diff-of-diffs concentration plots %%%')

    # Get a list of variables that GCPy should not read
    skip_vars = skip_these_vars

    # Target output files
    diff_of_diffs_refspc = './gcc_diffs_spc.nc4'
    diff_of_diffs_devspc = './gchp_diffs_spc.nc4'

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
    if 'lat' in refdims and 'Xdim' in devdims:
        gchp_ref_newdimnames = gchp_dev.copy()
        for v in gchp_dev.data_vars.keys():
            if 'Xdim' in gchp_dev[v].dims:
                gchp_ref_newdimnames[v].values = gchp_ref[v].values.reshape(
                    gchp_dev[v].values.shape)
                # NOTE: the reverse conversion is gchp_dev[v].stack(lat=('nf','Ydim')).transpose('time','lev','lat','Xdim').values
        gchp_ref = gchp_ref_newdimnames.copy()
    with xr.set_options(keep_attrs=True):
        gchp_diffs = gchp_dev.copy()
        for v in gchp_dev.data_vars.keys():
            if 'Xdim' in gchp_dev[v].dims or 'lat' in gchp_dev[v].dims:
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
    bmk.make_benchmark_lots(diff_of_diffs_refspc,
                            diff_of_diffs_refstr,
                            diff_of_diffs_devspc,
                            diff_of_diffs_devstr,
                            dst=diff_of_diffs_plotsdir,
                            overwrite=True,
                            use_cmap_RdBu=True)
