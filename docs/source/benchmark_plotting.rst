Benchmark Plotting / Tabling
============================

Below is an example configuration file used to input the desired options for the comprehensive benchmark comparison script `run_benchmark.py`. Additional configuration file examples can be found in the `benchmarks` directory of GCpy.

The `run_benchmark.py` script allows one to perform benchmark comparisons between any simulation duration supplied in the configuration file provided the ref and dev simulations time periods match. 
Additionally, if the durations specified are exactly one year, then the corresponding `bmk_type` specialty comparison script will be run (either `run_1yr_fullchem_benchmark.py` or `run_1yr_tt_benchmark.py`). Any other duration will run the standard suite of benchmark comparisons. 


.. code-block:: yaml


  ---
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
# configuration for 1 month benchmark
paths:
  # High-level directory containing subdirectories with data
  main_dir: /n/holyscratch01/external_repos/GEOS-CHEM/gcgrid/geos-chem/validation/gcpy_test_data/1mon
  results_dir: BenchmarkResults # Name to be used for directory of output from this script
  # Path to regridding weights
  weights_dir: /n/holyscratch01/external_repos/GEOS-CHEM/gcgrid/gcdata/ExtData/GCHP/RegriddingWeights
data:
  # contains configurations for ref and dev runs (version, dates, prior to v13.0.0)
  # timestamp format YYYY-MM-DDThh:mm:ss
  ref:
    gcc:
      version: GCC_ref # NOTE: these will be used in some filenames and so should not have spaces
      dir: GCC_ref # Directory name
      subdir: OutputDir
      bmk_start: "2019-07-01T00:00:00" 
      bmk_end: "2019-08-01T00:00:00"
    gchp:
      version: GCHP_ref
      dir: GCHP_ref
      subdir: OutputDir
      bmk_start: "2019-07-01T00:00:00"
      bmk_end: "2019-08-01T00:00:00"
      is_pre_13.1: False # Whether GCHP files are format before version 13.1
      is_pre_14.0: True # Whether GCHP files are format before version 14.0
  dev:
    gcc:
      version: GCC_dev 
      dir: GCC_dev
      subdir: OutputDir
      bmk_start: "2019-07-01T00:00:00" 
      bmk_end: "2019-08-01T00:00:00"
    gchp:
      version: GCHP_dev
      dir: GCHP_dev
      subdir: OutputDir
      bmk_start: "2019-07-01T00:00:00" 
      bmk_end: "2019-08-01T00:00:00"
      is_legacy: False
      is_pre_13.1: False # Whether GCHP files are format before version 13.1
      is_pre_14.0: True # Whether GCHP files are format before version 14.0
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

  outputs: # Output to generate (plots/tables will be created in this order):
    plot_conc: True
    plot_emis: True
    emis_table: True
    plot_jvalues: True
    plot_aod: True
    mass_table: True
    ops_budget_table: False
    OH_metrics: True # GCC only
    ste_table: True # GCC only
    plot_options: # Plot concentrations and emissions by category?
      by_spc_cat: True
      by_hco_cat: True

