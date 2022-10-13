############################
Benchmark Plotting / Tabling
############################

Below is an example configuration file used to input the desired
options for the comprehensive benchmark comparison script
:code:`run_benchmark.py`. Additional configuration file examples can
be found in the :file:`benchmarks` directory of GCpy.

The :file:`run_benchmark.py` script allows one to perform benchmark
comparisons between any simulation duration supplied in the
configuration file provided the ref and dev simulations time periods
match.  Additionally, if the durations specified are exactly one year,
then the corresponding :literal:`bmk_type` specialty comparison script
will be run (either :file:`run_1yr_fullchem_benchmark.py` or
:file:`run_1yr_tt_benchmark.py`). Any other duration will run the standard
suite of benchmark comparisons.

To generate plots from a 1-month benchmark simulation, you would call
:file:`run_benchmark.py` as follows:

.. code-block:: console

   (gcpy_env) $ run_benchmark.py 1mo_benchmark.yml

Where :file:`1mo_benchmark.yml` contains the following inputs:

.. code-block:: yaml

   &---
   # =====================================================================
   # Benchmark configuration file (**EDIT AS NEEDED**)
   # customize in the following manner:
   # (1) Edit the path variables so that they point to folders w/ model data
   # (2) Edit the version strings for each benchmark simulation
   # (3) Edit the switches that turn on/off creating of plots and tables
   # (4) If necessary, edit labels for the dev and ref versions
   # Note: When doing GCHP vs GCC comparisions gchp_dev will be compared
   # to gcc_dev (not gcc_ref!). This ensures consistency in version names
   # when doing GCHP vs GCC diff-of-diffs (mps, 6/27/19)
   # =====================================================================
   #
   # Configuration for 1 month FullChemBenchmark
   #
   # paths:
   #   main_dir:    High-level directory containing ref & dev rundirs
   #   results_dir: Directory where plots/tables will be created
   #   weights_dir: Path to regridding weights
   #   spcdb_dir:   Folder in which the species_database.yml file is
   #                located.  If set to "default", then will look for
   #                species_database.yml in one of the Dev rundirs.
   #
   paths:
     main_dir: /n/holyscratch01/external_repos/GEOS-CHEM/gcgrid/geos-chem/validation/gcpy_test_data/1mon
     results_dir: /path/to/BenchmarkResults
     weights_dir: /n/holyscratch01/external_repos/GEOS-CHEM/gcgrid/gcdata/ExtData/GCHP/RegriddingWeights
     spcdb_dir: default
   #
   # data: Contains configurations for ref and dev runs
   #   version:         Version string (must not contain spaces)
   #   dir:             Path to run directory
   #   outputs_subdir:  Subdirectory w/ GEOS-Chem diagnostic files
   #   restarts_subdir: Subdirectory w/ GEOS-Chem restarts
   #   bmk_start:       Simulation start date (YYYY-MM-DDThh:mm:ss)
   #   bmk_end:         Simulation end date (YYYY-MM-DDThh:mm:ss)
   #   resolution:      GCHP resolution string
   #
   data:
     ref:
       gcc:
         version: GCC_ref
         dir: GCC_ref
         outputs_subdir: OutputDir
         restarts_subdir: Restarts
         bmk_start: "2019-07-01T00:00:00"
         bmk_end: "2019-08-01T00:00:00"
       gchp:
         version: GCHP_ref
         dir: GCHP_ref
         outputs_subdir: OutputDir
         restarts_subdir: Restarts
         bmk_start: "2019-07-01T00:00:00"
         bmk_end: "2019-08-01T00:00:00"
         is_pre_13.1: False
         is_pre_14.0: False
         resolution: c24
     dev:
       gcc:
         version: GCC_dev
         dir: GCC_dev
         outputs_subdir: OutputDir
         restarts_subdir: Restarts
         bmk_start: "2019-07-01T00:00:00"
         bmk_end: "2019-08-01T00:00:00"
       gchp:
         version: GCHP_dev
         dir: GCHP_dev
         outputs_subdir: OutputDir
         restarts_subdir: Restarts
         bmk_start: "2019-07-01T00:00:00"
         bmk_end: "2019-08-01T00:00:00"
         is_pre_13.1: False
         is_pre_14.0: False
         resolution: c24
   #
   # options: Specify the types of comparisons to perform
   #
   options:
     bmk_type: FullChemBenchmark
     gcpy_test: True # Specify if this is a gcpy test validation run
     comparisons:
       gcc_vs_gcc:
         run: True # True to run this comparison
         dir: GCC_version_comparison
         tables_subdir: Tables
       gchp_vs_gcc:
         run: True
         dir: GCHP_GCC_comparison
         tables_subdir: Tables
       gchp_vs_gchp:
         run: True
         dir: GCHP_version_comparison
         tables_subdir: Tables
       gchp_vs_gcc_diff_of_diffs:
         run: True
         dir: GCHP_GCC_diff_of_diffs
   #
   # outputs: Types of output to generate (plots/tables)
   #
     outputs:
       plot_conc: True
       plot_emis: True
       emis_table: True
       plot_jvalues: True
       plot_aod: True
       mass_table: True
       ops_budget_table: False
       OH_metrics: True
       ste_table: True # GCC only
       plot_options: # Plot concentrations and emissions by category?
         by_spc_cat: True
         by_hco_cat: True

YAML configuration files for 1-year benchmarks
(:file:`1yr_fullchem_benchmark.yml`, :file:`1yr_tt_benchmark.yml`) are
also provided in the :file:`benchmarks` folder.
