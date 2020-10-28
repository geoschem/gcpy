Benchmark Plotting / Tabling
============================

This script is an abbreviated version of `run_1mo_benchmark.py`, the comprehensive benchmark comparison
script that also covers GCHP. `run_1mo_benchmark.py` is available in the `benchmarks` folder of the GCPy repository.


.. code-block:: python


   #!/usr/bin/env python
   """

   Run this script to generate benchmark comparisons between
   two runs of GEOS-Chem Classic. It is intended for 1-month diagnostics,
   but some parts work for other time periods.

   You can customize this by editing the following settings in the
   "Configurables" section below:

      (1) Edit the path variables so that they point to folders w/ model data
      (2) Edit the version strings for each benchmark simulation
      (3) Edit the switches that turn on/off creating of plots and tables
      (4) If necessary, edit labels for the dev and ref versions

   To test gcpy, copy this script anywhere you want to run the test and
   set gcpy_test to True at the top of the script. Benchmark artifacts will
   be created locally in new folder called Plots.

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
   import gcpy.oh_metrics as oh

   # Tell matplotlib not to look for an X-window
   os.environ["QT_QPA_PLATFORM"]="offscreen"

   # Suppress harmless run-time warnings (mostly about underflow in division)
   warnings.filterwarnings("ignore", category=RuntimeWarning)
   warnings.filterwarnings("ignore", category=UserWarning)

   # This script has a fixed benchmark type
   bmk_type     = "FullChemBenchmark"

   ########################################################################
   ###           CONFIGURABLE SETTINGS: ***EDIT AS NEEDED***            ###
   ########################################################################

   # =====================================================================
   # Benchmark information (**EDIT AS NEEDED**)
   # =====================================================================

   # High-level directory containing subdirectories with data
   maindir  = "/path/to/high/level/directory"

   # Version strings that denote your file names
   # NOTE: these will be used in some filenames and so should not have spaces
   # or other characters not appropriate for a filename.
   gcc_ref_version = "GCC_ref"
   gcc_dev_version = "GCC_dev"

   #dates of ref and dev runs
   ref_bmk_year     = '2019'
   ref_bmk_mon      = '7'
   dev_bmk_year     = '2019'
   dev_bmk_mon      = '7'

   # Name to be used for directory of output from this script
   results_dir = "BenchmarkResults"

   # Path to regridding weights (may create new weights in the given directory)
   weightsdir = "/path/to/regridding/weights"

   # =====================================================================
   # Specify if this is a gcpy test validation run
   # =====================================================================
   gcpy_test = True

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
   OH_metrics       = True
   ste_table        = True

   # Plot concentrations and emissions by category?
   plot_by_spc_cat = True
   plot_by_hco_cat = True

   # =====================================================================
   # Data directories
   # =====================================================================

   # Directory names (edit if not same as version strings)
   gcc_ref_dir = gcc_ref_version
   gcc_dev_dir = gcc_dev_version

   # Diagnostic file directory paths
   gcc_vs_gcc_refdir   = join(maindir, gcc_ref_dir,  "OutputDir")
   gcc_vs_gcc_devdir   = join(maindir, gcc_dev_dir,  "OutputDir")

   # Restart file directory paths
   gcc_vs_gcc_refrst   = join(maindir, gcc_ref_dir )
   gcc_vs_gcc_devrst   = join(maindir, gcc_dev_dir )

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
      if not os.path.exists(mainresultsdir): os.mkdir(mainresultsdir)
   else:
      gcc_vs_gcc_resultsdir    = join(maindir, gcc_dev_dir, results_dir)
      if not os.path.exists(gcc_vs_gcc_resultsdir): os.mkdir(gcc_vs_gcc_resultsdir)

   gcc_vs_gcc_tablesdir    = join(gcc_vs_gcc_resultsdir, "Tables")

   # =====================================================================
   # Plot title strings
   # =====================================================================
   gcc_vs_gcc_refstr    = gcc_ref_version
   gcc_vs_gcc_devstr    = gcc_dev_version

   ########################################################################
   ###    THE REST OF THESE SETTINGS SHOULD NOT NEED TO BE CHANGED      ###
   ########################################################################

   # =====================================================================
   # Dates and times
   # =====================================================================

   # Start and end months of the benchmark
   ref_b_start   = (int(ref_bmk_year), int(ref_bmk_mon))
   ref_b_stop    = (int(ref_bmk_year), int(ref_bmk_mon) + 1)

   dev_b_start   = (int(dev_bmk_year), int(dev_bmk_mon))
   dev_b_stop    = (int(dev_bmk_year), int(dev_bmk_mon) + 1)

   # Convert to strings
   ref_s_start   = (str(ref_b_start[0]), str(ref_b_start[1]).zfill(2))
   ref_s_stop    = (str(ref_b_stop[0]),  str(ref_b_stop[1]).zfill(2))

   dev_s_start   = (str(dev_b_start[0]), str(dev_b_start[1]).zfill(2))
   dev_s_stop    = (str(dev_b_stop[0]),  str(dev_b_stop[1]).zfill(2))

   # Timestamps for files
   gcc_ref_date  = np.datetime64( "{}-{}-01T00:00:00".format(ref_s_start[0], ref_s_start[1]))
   end_ref_date  = np.datetime64("{}-{}-01T00:00:00".format(ref_s_stop[0], ref_s_stop[1]))

   gcc_dev_date  = np.datetime64( "{}-{}-01T00:00:00".format(dev_s_start[0], dev_s_start[1]))
   end_dev_date  = np.datetime64("{}-{}-01T00:00:00".format(dev_s_stop[0], dev_s_stop[1]))

   # Seconds per month
   ref_sec_in_bmk_month = (end_ref_date - gcc_ref_date).astype("float64")
   dev_sec_in_bmk_month = (end_dev_date - gcc_dev_date).astype("float64")

   if not np.equal(ref_sec_in_bmk_month, dev_sec_in_bmk_month):
      print('Skipping emissions tables and operations budget tables because months are' + \
           'different lengths')
      emis_table=False
      ops_budget_table=False

   # String for month and year (e.g. "Jul2016")
   if np.equal(gcc_ref_date, gcc_dev_date):
      mon_yr_str = calendar.month_abbr[ref_b_start[1]] + ref_s_start[0]
   else:
      mon_yr_str = calendar.month_abbr[ref_b_start[1]] + ref_s_start[0] + 'Vs' + \
                calendar.month_abbr[dev_b_start[1]] + dev_s_start[0]

   # ======================================================================
   # Significant difference filenames
   # ======================================================================

   vstr = "{}_vs_{}".format(gcc_ref_version, gcc_dev_version)
   gcc_vs_gcc_sigdiff = [
      join(gcc_vs_gcc_resultsdir, "{}_sig_diffs_sfc.txt".format(vstr)),
      join(gcc_vs_gcc_resultsdir, "{}_sig_diffs_500hpa.txt".format(vstr)),
      join(gcc_vs_gcc_resultsdir, "{}_sig_diffs_zonalmean.txt".format(vstr)),
      join(gcc_vs_gcc_resultsdir, "{}_sig_diffs_emissions.txt".format(vstr))]

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

   # ======================================================================
   # Create GCC vs GCC benchmark plots and tables
   # ======================================================================

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
        ref = get_filepath(gcc_vs_gcc_refdir, col, gcc_ref_date)
        dev = get_filepath(gcc_vs_gcc_devdir, col, gcc_dev_date)

        # Meteorology data needed for calculations
        colmet = "StateMet"
        refmet = get_filepath(gcc_vs_gcc_refdir, colmet, gcc_ref_date)
        devmet = get_filepath(gcc_vs_gcc_devdir, colmet, gcc_dev_date)

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
        ref = get_filepath(gcc_vs_gcc_refdir, col, gcc_ref_date)
        dev = get_filepath(gcc_vs_gcc_devdir, col, gcc_dev_date)

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
        ref = get_filepath(gcc_vs_gcc_refdir, col, gcc_ref_date)
        dev = get_filepath(gcc_vs_gcc_devdir, col, gcc_dev_date)

        # Print emisisons and inventory tables
        bmk.make_benchmark_emis_tables(
            ref,
            gcc_vs_gcc_refstr,
            dev,
            gcc_vs_gcc_devstr,
            dst=gcc_vs_gcc_resultsdir,
            ref_interval=[ref_sec_in_bmk_month],
            dev_interval=[dev_sec_in_bmk_month],
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
        ref = get_filepath(gcc_vs_gcc_refdir, col, gcc_ref_date)
        dev = get_filepath(gcc_vs_gcc_devdir, col, gcc_dev_date)

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
        ref = get_filepath(gcc_vs_gcc_refdir, col, gcc_ref_date)
        dev = get_filepath(gcc_vs_gcc_devdir, col, gcc_dev_date)

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
        ref = get_filepath(gcc_vs_gcc_refrst, col, end_ref_date)
        dev = get_filepath(gcc_vs_gcc_devrst, col, end_dev_date)

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
        ref = get_filepath(gcc_vs_gcc_refdir, col, gcc_ref_date)
        dev = get_filepath(gcc_vs_gcc_devdir, col, gcc_dev_date)

        # Make budget table. Include calculation of Strat and Accumulation
        bmk.make_benchmark_operations_budget(
            gcc_ref_version,
            ref,
            gcc_dev_version,
            dev,
            ref_sec_in_bmk_month,
            dev_sec_in_bmk_month,
            benchmark_type=bmk_type,
            label=mon_yr_str,
            dst=gcc_vs_gcc_tablesdir
        )

    #---------------------------------------------------------------
    # GCC vs. GCC global mean OH, MCF Lifetime, CH4 Lifetime
    #---------------------------------------------------------------
    if OH_metrics:
        print("\n%%% Creating GCC vs. GCC OH metrics table %%%")

      #Use this for benchmarks prior to GEOS-Chem 13.0.0
      '''
        # Diagnostic collection files to read
        col  = "ConcAfterChem"
        ref = get_filepath(gcc_vs_gcc_refdir, col, gcc_ref_date)
        dev = get_filepath(gcc_vs_gcc_devdir, col, gcc_dev_date)

        # Meteorology data needed for calculations
        col = "StateMet"
        refmet = get_filepath(gcc_vs_gcc_refdir, col, gcc_ref_date)
        devmet = get_filepath(gcc_vs_gcc_devdir, col, gcc_dev_date)

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
      '''
      
        # Diagnostic collection files to read
        col = "Metrics"
        ref = get_filepath(gcc_vs_gcc_refdir, col, gcc_ref_date)
        dev = get_filepath(gcc_vs_gcc_devdir, col, gcc_dev_date)

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

    # --------------------------------------------------------------
    # GCC dev Strat-Trop Exchange
    # --------------------------------------------------------------
    if ste_table:
        print("\n%%% Creating GCC dev Strat-Trop Exchange table %%%")

        # Diagnostic collection files to read
        col = "AdvFluxVert"
        dev = get_filepath(gcc_vs_gcc_devdir, col, gcc_dev_date)

        # Compute monthly and annual average strat-trop exchange of O3
        ste.make_benchmark_ste_table(
            gcc_dev_version,
            dev,
            dev_b_start[0],
            bmk_type=bmk_type,
            dst=gcc_vs_gcc_tablesdir,
            species=['O3'],
            overwrite=True,
            month=dev_b_start[1]
        )
