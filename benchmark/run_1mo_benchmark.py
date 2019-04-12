#!/usr/bin/env python
'''
Driver script for creating benchmark plots.  The options are:
(1) GEOS-Chem "Classic" vs. GEOS-Chem "Classic"
(2) GCHP vs GEOS-Chem "Classic"
(3) GCHP vs GCHP

Set the path variables below to the folders containing model data
files from benchmark simulations.  Then set the switches according
to the plotting options that you want.
'''

import os
from gcpy import benchmark
import warnings

# Suppress harmless run-time warnings (mostly about underflow in division)
warnings.filterwarnings("ignore", category=RuntimeWarning)

# =====================================================================
# Configurables
# =====================================================================

# Benchmark information (*MUST EDIT*)
maindir  = '/path/to/main/directory'
ref_version = 'ref_version_string'
dev_version = 'dev_version_string'

# Comparisons to run (edit as needed)
gcc_vs_gcc   = True
gchp_vs_gcc  = False
gchp_vs_gchp = False

# Output to generate (edit as needed)
plot_aod     = True
plot_conc    = False
plot_emis    = False
plot_jvalues = False
emis_table   = False

# Filename date strings (edit as needed)
gcc_datestr  = '20160701'
gchp_datestr = '20160716'

# Filename hour strings (edit as needed)
gcc_hourstr  = '0000'
gchp_hourstr = '1200'

# Data directories (edit as needed)
gcc_vs_gcc_refdir   = os.path.join(maindir, ref_version)
gcc_vs_gcc_devdir   = os.path.join(maindir, dev_version)
gchp_vs_gcc_refdir  = os.path.join(maindir, dev_version)
gchp_vs_gcc_devdir  = os.path.join(maindir, dev_version, 'gchp/OutputDir')
gchp_vs_gchp_refdir = os.path.join(maindir, ref_version, 'gchp/OutputDir')
gchp_vs_gchp_devdir = os.path.join(maindir, dev_version, 'gchp/OutputDir')

# Plots directories (edit as needed)
gcc_vs_gcc_plotsdir   = os.path.join(maindir, dev_version, 'output')
gchp_vs_gchp_plotsdir = os.path.join(maindir, dev_version, 'output/GCHP_version_comparison')
gchp_vs_gcc_plotsdir  = os.path.join(maindir, dev_version, 'output/GCHP_GCC_comparison')

# Plot title strings (edit as needed)
gcc_vs_gcc_refstr   = '{}'.format(ref_version)
gcc_vs_gcc_devstr   = '{}'.format(dev_version)
gchp_vs_gcc_refstr  = 'GCC {}'.format(dev_version)
gchp_vs_gcc_devstr  = 'GCHP {}'.format(dev_version)
gchp_vs_gchp_refstr = 'GCHP {}'.format(ref_version)
gchp_vs_gchp_devstr = 'GCHP {}'.format(dev_version)

# =====================================================================
# The rest of these settings should not need to be changed
# =====================================================================

# Species concentration filenames
gcc_spcfile  = 'GEOSChem.SpeciesConc.{}_{}z.nc4'.format(gcc_datestr,      \
                                                        gcc_hourstr)   
gchp_spcfile = 'GCHP.SpeciesConc.{}_{}z.nc4'.format(gchp_datestr,         \
                                                    gchp_hourstr) 

# HEMCO diagnostic filenames
gcc_hcofile  = 'HEMCO_diagnostics.{}{}.nc'.format(gcc_datestr, gcc_hourstr)
gchp_hcofile = 'GCHP.Emissions.{}_{}z.nc4'.format(gchp_datestr, gchp_hourstr)

# Local noon J-value diagnostic filenames
gcc_jvfile  = 'GEOSChem.JValuesLocalNoon.{}_{}z.nc4'.format(gcc_datestr,  \
                                                            gcc_hourstr)
gchp_jvfile = 'GCHP.JValuesLocalNoon.{}_{}z.nc4'.format(gchp_datestr,     \
                                                        gchp_hourstr)

# Aerosol optical depth diagnostic filenames
gcc_aodfile  = 'GEOSChem.Aerosols.{}_{}z.nc4'.format(gcc_datestr, gcc_hourstr)
gchp_aodfile = 'GCHP.Aerosols.{}_{}z.nc4'.format(gchp_datestr, gchp_hourstr)

# Paths to species concentration data
gcc_vs_gcc_refspc   = os.path.join(maindir, gcc_vs_gcc_refdir,   gcc_spcfile)
gcc_vs_gcc_devspc   = os.path.join(maindir, gcc_vs_gcc_devdir,   gcc_spcfile)
gchp_vs_gcc_refspc  = os.path.join(maindir, gchp_vs_gcc_refdir,  gcc_spcfile)
gchp_vs_gcc_devspc  = os.path.join(maindir, gchp_vs_gcc_devdir,  gchp_spcfile)
gchp_vs_gchp_refspc = os.path.join(maindir, gchp_vs_gchp_refdir, gchp_spcfile)
gchp_vs_gchp_devspc = os.path.join(maindir, gchp_vs_gchp_devdir, gchp_spcfile)

# Paths to HEMCO diagnostics data
gcc_vs_gcc_refhco   = os.path.join(maindir, gcc_vs_gcc_refdir,   gcc_hcofile)
gcc_vs_gcc_devhco   = os.path.join(maindir, gcc_vs_gcc_devdir,   gcc_hcofile)
gchp_vs_gcc_refhco  = os.path.join(maindir, gchp_vs_gcc_refdir,  gcc_hcofile)
gchp_vs_gcc_devhco  = os.path.join(maindir, gchp_vs_gcc_devdir,  gchp_hcofile)
gchp_vs_gchp_refhco = os.path.join(maindir, gchp_vs_gchp_refdir, gchp_hcofile)
gchp_vs_gchp_devhco = os.path.join(maindir, gchp_vs_gchp_devdir, gchp_hcofile)

# Paths to local noon J-value data
gcc_vs_gcc_refjv    = os.path.join(maindir, gcc_vs_gcc_refdir,   gcc_jvfile)
gcc_vs_gcc_devjv    = os.path.join(maindir, gcc_vs_gcc_devdir,   gcc_jvfile)
gchp_vs_gcc_refjv   = os.path.join(maindir, gchp_vs_gcc_refdir,  gcc_jvfile)
gchp_vs_gcc_devjv   = os.path.join(maindir, gchp_vs_gcc_devdir,  gchp_jvfile)
gchp_vs_gchp_refjv  = os.path.join(maindir, gchp_vs_gchp_refdir, gchp_jvfile)
gchp_vs_gchp_devjv  = os.path.join(maindir, gchp_vs_gchp_devdir, gchp_jvfile)

# Paths to aerosol optical depth data
gcc_vs_gcc_refaod   = os.path.join(maindir, gcc_vs_gcc_refdir,   gcc_aodfile)
gcc_vs_gcc_devaod   = os.path.join(maindir, gcc_vs_gcc_devdir,   gcc_aodfile)
gchp_vs_gcc_refaod  = os.path.join(maindir, gchp_vs_gcc_refdir,  gcc_aodfile)
gchp_vs_gcc_devaod  = os.path.join(maindir, gchp_vs_gcc_devdir,  gchp_aodfile)
gchp_vs_gchp_refaod = os.path.join(maindir, gchp_vs_gchp_refdir, gchp_aodfile)
gchp_vs_gchp_devaod = os.path.join(maindir, gchp_vs_gchp_devdir, gchp_aodfile)

# =====================================================================
# Create GCC vs GCC benchmark plots and tables
# =====================================================================

if gcc_vs_gcc:

    if plot_conc:
        # Concentration plots
        # (includes lumped species and separates by category)
        benchmark.make_benchmark_conc_plots(gcc_vs_gcc_refspc,           \
                                            gcc_vs_gcc_refstr,           \
                                            gcc_vs_gcc_devspc,           \
                                            gcc_vs_gcc_devstr,           \
                                            dst=gcc_vs_gcc_plotsdir,     \
                                            overwrite=True)

    if plot_emis:
        # Emissions plots
        benchmark.make_benchmark_emis_plots(gcc_vs_gcc_refhco,           \
                                            gcc_vs_gcc_refstr,           \
                                            gcc_vs_gcc_devhco,           \
                                            gcc_vs_gcc_devstr,           \
                                            dst=gcc_vs_gcc_plotsdir,     \
                                            plot_by_benchmark_cat=True,  \
                                            plot_by_hco_cat=True,        \
                                            overwrite=True)

    if emis_table:
        # Emissions tables
        benchmark.make_benchmark_emis_tables(gcc_vs_gcc_refhco,          \
                                             gcc_vs_gcc_refstr,          \
                                             gcc_vs_gcc_devhco,          \
                                             gcc_vs_gcc_devstr,          \
                                             dst=gcc_vs_gcc_plotsdir,    \
                                             overwrite=True)

    if plot_jvalues:
        # Local noon J-values plots
        benchmark.make_benchmark_jvalue_plots(gcc_vs_gcc_refjv,          \
                                              gcc_vs_gcc_refstr,         \
                                              gcc_vs_gcc_devjv,          \
                                              gcc_vs_gcc_devstr,         \
                                              dst=gcc_vs_gcc_plotsdir,   \
                                              local_noon_jvalues=True,   \
                                              overwrite=True)
    if plot_aod:
        # Column AOD plots
        benchmark.make_benchmark_aod_plots(gcc_vs_gcc_refaod,            \
                                           gcc_vs_gcc_refstr,            \
                                           gcc_vs_gcc_devaod,            \
                                           gcc_vs_gcc_devstr,            \
                                           dst=gcc_vs_gcc_plotsdir,      \
                                           overwrite=True)

# =====================================================================
# Create GCHP vs GCC benchmark plots and tables
# =====================================================================

if gchp_vs_gcc:
    if plot_conc:
        # Concentration plots
        # (includes lumped species and separates by category)
        benchmark.make_benchmark_conc_plots(gchp_vs_gcc_refspc,          \
                                            gchp_vs_gcc_refstr,          \
                                            gchp_vs_gcc_devspc,          \
                                            gchp_vs_gcc_devstr,          \
                                            dst=gchp_vs_gcc_plotsdir,    \
                                            overwrite=True)

    if plot_emis:
        # Emissions plots
        benchmark.make_benchmark_emis_plots(gchp_vs_gcc_refhco,          \
                                            gchp_vs_gcc_refstr,          \
                                            gchp_vs_gcc_devhco,          \
                                            gchp_vs_gcc_devstr,          \
                                            dst=gchp_vs_gcc_plotsdir,    \
                                            plot_by_benchmark_cat=True,  \
                                            plot_by_hco_cat=True,        \
                                            overwrite=True,              \
                                            flip_dev=True)

    if emis_table:
        # Table of emissions totals
        pass

    if plot_jvalues:
        # Local noon J-values plots
        benchmark.make_benchmark_jvalue_plots(gchp_vs_gcc_refjv,         \
                                              gchp_vs_gcc_refstr,        \
                                              gchp_vs_gcc_devjv,         \
                                              gchp_vs_gcc_devstr,        \
                                              dst=gchp_vs_gcc_plotsdir,  \
                                              local_noon_jvalues=True,   \
                                              overwrite=True)
    if plot_aod:
        # Column AOD plots
        benchmark.make_benchmark_aod_plots(gchp_vs_gcc_refaod,           \
                                           gchp_vs_gcc_refstr,           \
                                           gchp_vs_gcc_devaod,           \
                                           gchp_vs_gcc_devstr,           \
                                           dst=gchp_vs_gcc_plotsdir,     \
                                           overwrite=True)

# =====================================================================
# Create GCHP vs GCHP benchmark plots and tables
# =====================================================================

if gchp_vs_gchp:
    if plot_conc:
        # Concentration plots
        # (includes lumped species and separates by category)
        benchmark.make_benchmark_conc_plots(gchp_vs_gchp_refspc,         \
                                            gchp_vs_gchp_refstr,         \
                                            gchp_vs_gchp_devspc,         \
                                            gchp_vs_gchp_devstr,         \
                                            dst=gchp_vs_gchp_plotsdir,   \
                                            overwrite=True)

    if plot_emis:
        # Emissions plots
        benchmark.make_benchmark_emis_plots(gchp_vs_gchp_refhco,         \
                                            gchp_vs_gchp_refstr,         \
                                            gchp_vs_gchp_devhco,         \
                                            gchp_vs_gchp_devstr,         \
                                            dst=gchp_vs_gchp_plotsdir,   \
                                            plot_by_benchmark_cat=True,  \
                                            plot_by_hco_cat=True,        \
                                            overwrite=True,              \
                                            flip_ref=True,               \
                                            flip_dev=True)

    if emis_table:
        # Table of emissions totals
        pass

    if plot_jvalues:
        # Local noon J-values plots
        benchmark.make_benchmark_jvalue_plots(gchp_vs_gchp_refjv,        \
                                              gchp_vs_gchp_refstr,       \
                                              gchp_vs_gchp_devjv,        \
                                              gchp_vs_gchp_devstr,       \
                                              dst=gchp_vs_gchp_plotsdir, \
                                              local_noon_jvalues=True,   \
                                              overwrite=True)

if plot_aod:
        # Column AOD plots
        benchmark.make_benchmark_aod_plots(gchp_vs_gchp_refaod,          \
                                           gchp_vs_gchp_refstr,          \
                                           gchp_vs_gchp_devaod,          \
                                           gchp_vs_gchp_devstr,          \
                                           dst=gchp_vs_gchp_plotsdir,    \
                                           overwrite=True)
