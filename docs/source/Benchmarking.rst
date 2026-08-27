.. |br| raw:: html

   <br/>

.. _bmk:

############
Benchmarking
############

The `GEOS-Chem Support Team
<https://geoschem.github.io/support-team>`_ uses GCPy to produce
comparison plots and summary tables from GEOS-Chem benchmark
simulations.  In this chapter we will describe this capability of
GCPy.

.. _bmk-scripts:

======================================
Location of benchmark plotting scripts
======================================

The source code for creating benchmark plots is located in the
:file:`gcpy/benchmark` directory tree.

.. list-table:: **Contents of the gcpy/benchmark directory**
   :header-rows: 1
   :widths: 25 75

   * - File or folder
     - Description
   * - :file:`run_benchmark.py`
     - Benchmark driver script :mod:`gcpy.benchmark.run_benchmark`
   * - :file:`benchmark_slurm.sh`
     - Bash script to submit :file:`run_benchmark.py` as a SLURM batch job
   * - :file:`cloud/`
     - Directory containing template config files (in YAML format) for
       1-hour and 1-month benchmark plot jobs on the AWS cloud.
   * - :file:`config/`
     - Directory containing editable config files (in YAML format) for
       1-month and 1-year benchmark plot jobs.
   * - :file:`__init__.py`
     - Python import script
   * - :file:`modules/` 
     - Contains Python modules imported into the
       :file:`run_benchmark.py` script.  See
       :mod:`gcpy.benchmark.modules` for a detailed listing.
   * - :file:`README.md`
     - Readme file in Markdown format

.. _bmk-steps:

===============================
Steps to create benchmark plots
===============================

Follow these instructions to create comparison plots and summary
tables from GEOS-Chem benchmark simulations.

#. Copy a configuration file from the :file:`gcpy/benchmark/config`
   directory.

   In this example we will use the configuration file that will create
   plots from 1-year full-chemistry benchmark
   simulations. (Configuration files for other benchmark types have a
   similar layout.)

   .. code-block:: console

      $ cp /path/to/GCPy/gcpy/benchmark/config/1yr_fullchem_benchmark.yml .

   |br|

#. Edit the :literal:`paths` section of the configuration file to
   specify the proper directory paths for your system.

   .. code-block:: yaml

      # Configuration for 1-year FullChemBenchmark
      #
      # paths:
      #   main_dir:     High-level directory containing ref & dev rundirs
      #   results_dir:  Directory where plots/tables will be created
      #   weights_dir:  Path to regridding weights
      #   obs_data:     Paths to observations (for models vs. obs plots)
      #
      paths:
        main_dir: /path/to/benchmark/main/dir
        results_dir: /path/to/BenchmarkResults
        weights_dir: /n/holylfs06/LABS/jacob_lab/Shared/GEOS-CHEM/gcgrid/gcdata/ExtData/GCHP/RegriddingWeights
        #
        # Observational data dirs are on Harvard Cannon, edit if necessary
        #
        obs_data:
          ebas_o3:
            data_dir: /n/lab_storage/jacob_lab/Lab/obs_data_for_bmk/sondes_2010-2019
            data_label: "O3 (EBAS, 2019)"
          sondes:
            data_dir: /n/lab_storage/jacob_lab/Lab/obs_data_for_bmk/sondes_2010-2019
            data_file: allozonesondes_2010-2019.csv
            site_file: allozonesondes_site_elev.csv

   |br|

#. Edit the :literal:`data` section to specify the directories (and
   labels) for the Ref and Dev versions for GEOS-Chem Classic and GCHP.

   .. code-block:: yaml

      #
      # data: Contains configurations for ref and dev runs
      #   version:          Version string (must not contain spaces)
      #   dir:              Path to run directory
      #   outputs_subdir:   Subdirectory w/ GEOS-Chem diagnostic files
      #   restarts_subdir:  Subdirectory w/ GEOS-Chem restarts
      #   logs_subdir:      Subdirectory w/ GEOS-Chem log files
      #   logs_template:    Template for log file names (may include tokens)
      #   bmk_start:        Simulation start date (YYYY-MM-DDThh:mm:ss)
      #   bmk_end:          Simulation end date (YYYY-MM-DDThh:mm:ss)
      #   resolution:       GCHP resolution string
      #
      #
      data:
        ref:
          gcc:
            version: GCC_ref
            dir: GCC_ref
            outputs_subdir: OutputDir
            restarts_subdir: Restarts
            logs_subdir: Logs
            logs_template: "log.%Y%m%d"
            bmk_start: "2019-01-01T00:00:00"
            bmk_end: "2020-01-01T00:00:00"
          gchp:
            version: GCHP_ref
            dir: GCHP_ref
            outputs_subdir: OutputDir
            restarts_subdir: Restarts
            logs_subdir: Logs
            logs_template: "gchp.%Y%m%d_0000z.log"
            bmk_start: "2019-01-01T00:00:00"
            bmk_end: "2020-01-01T00:00:00"
            is_pre_14.0: False
            resolution: c24
        dev:
          gcc:
            version: GCC_dev
            dir: GCC_dev
            outputs_subdir: OutputDir
            restarts_subdir: Restarts
            logs_subdir: Logs
            logs_template: "log.%Y%m%d"
            bmk_start: "2019-01-01T00:00:00"
            bmk_end: "2020-01-01T00:00:00"
          gchp:
            version: GCHP_dev
            dir: GCHP_dev
            outputs_subdir: OutputDir
            restarts_subdir: Restarts
            logs_subdir: Logs
            logs_template: "gchp.%Y%m%d_0000z.log"
            bmk_start: "2019-01-01T00:00:00"
            bmk_end: "2020-01-01T00:00:00"
            is_pre_14.0: False
            resolution: c24

   |br|

#. Edit the :literal:`comparisons` section to specify the types of
   comparisons you would like to perform.

   .. code-block:: yaml

      #
      # comparisons: Specifies the comparisons to perform.
      #
      comparisons:
        gcc_vs_gcc:
          run: True
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

   |br|

#. Edit the :literal:`outputs` section to select the plots and tables
   that you would like to generate.

   .. code-block:: yaml

      #
      # outputs: Specifies the plots and tables to generate
      #
      outputs:
        #
        # Benchmark plots
        #
        plot_aod: True
        plot_conc: True
        plot_drydep: True
        plot_emis: True
        plot_jvalues: True
        plot_models_vs_obs: True
        plot_options:
          by_spc_cat: True
          by_hco_cat: True
          # yaxis_units: "pressure" or "level" for the zonal-mean plot Y-axis
          yaxis_units: "pressure"
        #
        # Benchmark tables
        #
        aer_budget_table: True
        emis_table: True
        mass_accum_table: False
        mass_table: True
        Ox_budget_table: True
        OH_metrics: True
        ops_budget_table: False
        sanity_check_table: True
        ste_table: True # GCC only
        summary_table: False
        timing_table: False
        #
        # Comparison plots for selected collections
        # (not normally used for benchmarks)
        #
        plot_budget: False
        plot_2d_met: False
        plot_3d_met: False
        plot_uvflux: False

   |br|

#. Edit the :literal:`n_cores` setting if you wish to change the
   number of computational cores to use.  If not, leave
   :literal:`n_cores` set to :literal:`-1`, which will use as many
   cores as possible.

   .. code-block:: yaml

      #
      # n_cores: Specify the number of cores to use.
      # -1: Use $OMP_NUM_THREADS         cores
      # -2: Use $OMP_NUM_THREADS - 1     cores
      # -N: Use $OMP_NUM_THREADS - (N-1) cores
      #  1: Disable parallelization (use a single core)
      #
      n_cores: -1

   |br|

#. Run :mod:`gcpy.benchmark.run_benchmark`.  You may do this in 2
   ways:

   #. Direct execution from the command line:

      .. code-block:: console

         (gcpy_env) $ python -m gcpy.benchmark.run_benchmark 1yr_fullchem_benchmark.yml

   #. Batch execution with the SLURM scheduler.  First, copy the
      :file:`benchmark_slurm.sh` script to your current directory:

      .. code-block:: console

         (gcpy_env) $ cp /path/to/GCPy/gcpy/benchmark/benchmark_slurm.sh .

      Next, edit your local copy of :file:`benchmark_slurm.sh` to
      specify your SLURM partition name, number of cores, the name of
      your Python environment and the configuration file to use.

      .. code-block:: bash

         #!/bin/bash
         
         #SBATCH -c 8
         #SBATCH -N 1
         #SBATCH -t 0-6:00
         #SBATCH -p sapphire,huce_cascade,seas_compute,shared
         #SBATCH --mem=180000
         #SBATCH --mail-type=END
         
         #============================================================================
         # This us a sample SLURM script that you can use to run the GCPy
         # benchmark plotting code as a SLURM batch job.
         #
         # You can modify the SLURM parameters above for your setup.
         #
         # Tips:
         # -----
         # (1) Use fewer cores to reduce the memory footprint. This may prevent
         #     your job from running out of memory.  Python under Linux seems
         #     to have an issue where not all memory is released back to the OS.
         #
         # (2) We recommend that you generate only one benchmark comparison
         #     (GCC vs GCC, GCHP vs GCC, GCHP vs GCC, or diff of diffs)
         #     at a time.  Otherwise your job will probaly run out of memory.
         #
         # (3) For diff-of-diffs plots, we recommend using 6 cores.
         #============================================================================
         
         # Apply all bash initialization settings
         . ~/.bashrc
         
         # Make sure to set multiple threads; Joblib will use multiple
         # cores to parallelize certain plotting operations.
         export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
         export OMP_STACKSIZE=500m
         
         # Use a non-interactive backend for matplotlib (we're printing to file)
         export MPLBACKEND=agg
         
         # Turn on Python environment (edit for your setup)
         conda activate gcpy_env
         
         # Specify a YAML file with benchmark options
         # Uncomment the file that you wish:
         config="1mo_benchmark.yml"
         #config="1yr_fullchem_benchmark.yml"
         #config="1yr_tt_benchmark.yml"
         
         # Call the run_benchmark script to make the plots
         python -m gcpy.benchmark.run_benchmark "${config}" > "${config/.yml/.log}" 2>&1
         
         # Turn off python environment
         conda deactivate
         
         exit 0

      Lastly, start the SLURM batch execution with this command:

      .. code-block:: console

         $ sbatch benchmark_slurm.sh

.. _bmk-funcs-plot:

============================
Benchmark plotting functions
============================

Module :code:`gcpy.benchmark_funcs` contains several functions for
creating plots and tables from GEOS-Chem benchmark simulations. The
specific outputs generated have been requested by the `GEOS-Chem
Steering Committee <https://geoschem.github.io/steering-cmte>`_  in
order to facilitate comparing benchmark output from different model
versions.

In this section, we will describe functions that create comparison
plots from GEOS-Chem benchmark simulation output.  The functions to
create summary tables will be described :ref:`in a separate section
<bmk-funcs-table>`.

.. list-table:: **Functions creating six-panel comparison plots**
   :align: center
   :header-rows: 1
   :widths: 70 30

   * - Function
     - Plot that it creates
   * - :func:`gcpy.benchmark.modules.benchmark_funcs.make_benchmark_aod_plots`
     - Aerosol optical depth
   * - :func:`gcpy.benchmark.modules.benchmark_funcs.make_benchmark_conc_plots` [#A]_ [#B]_ [#C]_
     - Species concentrations
   * - :func:`gcpy.benchmark.modules.benchmark_drydep.make_benchmark_drydep_plots`
     - Dry deposition velocities
   * - :func:`gcpy.benchmark.modules.benchmark_funcs.make_benchmark_emis_plots`
     - Emissions (by species and category)
   * - :func:`gcpy.benchmark.modules.benchmark_funcs.make_benchmark_jvalue_plots`
     - J-values (photolysis)
   * - :func:`gcpy.benchmark.modules.benchmark_funcs.make_benchmark_wetdep_plots`
     - Wet deposition of soluble species

.. list-table:: **Functions creating model vs. observation plots**
   :align: center
   :header-rows: 1
   :widths: 40 60

   * - Function
     - Plot that it creates
   * - :mod:`gcpy.benchmark.modules.benchmark_models_vs_obs`
     - Modeled ozone vs. surface observations [#D]_ [#E]_
   * - :mod:`gcpy.benchmark.modules.benchmark_models_vs_sondes`
     - Vertical profiles of modeled ozone vs. ozonesondes [#D]_

.. rubric:: Notes:

.. [#A] This function is the only benchmark plotting function that
	supports diff-of-diffs plotting, in which 4 datasets are
	passed and the differences between two groups of
	:literal:`Ref` datasets vs. two groups of :literal:`Dev`
        datasets is plotted (typically used for   comparing changes in
	GCHP vs. changes in GEOS-Chem Classic across model versions).

.. [#B] This is also the only benchmark plotting function that sends
	plots to separate folders based on category (as denoted by the
	plot_by_spc_cat flag). The full list of species categories is
	denoted in `benchmark_categories.yml
	<https://github.com/geoschem/gcpy/blob/dev/gcpy/benchmark_categories.yml>`_

.. [#C] In this function, parallelization occurs at the species
	category level. In all other functions, parallelization occurs
	within calls to :mod:`gcpy.plot.compare_single_level`  and
	:mod:`gcpy.plot.compare_zonal_mean`.

.. [#D] Only available in 1-year fullchem benchmarks.

.. [#E] GEOS-Chem :literal:`SpeciesConc` O\ :sub:`3` diagnostic
	outputs vs. the `EBAS 2019 <https://ebas-data.nilu.no/>`_
	observations.

The functions listed above create comparison plots of most GEOS-Chem
output variables divided into specific categories, e.g. species
categories such as :literal:`Aerosols` or :literal:`Bromine` for the
:literal:`SpeciesConcVV` diagnostic. In eachcategory, these function
create single level PDFs for the surface and 500hPa and zonal
mean PDFs for the entire atmosphere and only the stratosphere (defined
a 1-100hPa). For
:func:`gcpy.benchmark.modules.benchmark_funcs.make_benchmark_emis_plots`,
only single level plots at the surface are produced. All of these
plotting functions include bookmarks within the generated PDFs that
point to the pages containing each plotted quantity. Thus these
functions serve as tools for quickly creating comprehensive plots
comparing two GEOS-Chem runs. These functions are used to create the
publicly available plots for 1-month and 1-year benchmarks of new
versions of GEOS-Chem.

Many of the plotting functions listed above use pre-defined lists of
variables in YAML files. If one dataset includes a variable but the
other dataset does not, the data for that variable in the latter
dataset will be considered to be NaN and will be plotted as such.

.. _bmk-funcs-table:

===========================
Benchmark tabling functions
===========================

The following functions generate summary tables from GEOS-Chem benchmark output:


.. table:: **Functions creating benchmark tables**

   +--------------------------------------------------------------------------------------------------------------+
   | Function and table that it creates                                                                           |
   +==============================================================================================================+
   | :func:`gcpy.benchmark.modules.budget_ox.global_ox_budget` |br|                                               |
   | |br|                                                                                                         |
   | Ox budget (1yr benchmarks only)                                                                              |
   +--------------------------------------------------------------------------------------------------------------+
   | :func:`gcpy.benchmark.modules.benchmark_funcs.make_benchmark_aerosol_tables` |br|                            |
   | |br|                                                                                                         |
   | Global aerosol burdens (1yr benchmarks only)                                                                 |
   +--------------------------------------------------------------------------------------------------------------+
   | :func:`gcpy.benchmark.modules.benchmark_funcs.make_benchmark_emis_tables` |br|                               |
   | |br|                                                                                                         |
   | Emissions (by species & inventory)                                                                           |
   +--------------------------------------------------------------------------------------------------------------+
   | :func:`gcpy.benchmark.modules.benchmark_scrape_gcclassic_timers.make_benchmark_gcclassic_timing_table` |br|  |
   | |br|                                                                                                         |
   | GEOS-Chem Classic timers output                                                                              |
   +--------------------------------------------------------------------------------------------------------------+
   | :func:`gcpy.benchmark.modules.benchmark_scrape_gchp_timers.make_benchmark_gchp_timing_table` |br|            |
   | |br|                                                                                                         |
   | GCHP timers output                                                                                           |
   +--------------------------------------------------------------------------------------------------------------+
   | :func:`gcpy.benchmark.modules.benchmark_funcs.make_benchmark_mass_tables` |br|                               |
   | |br|                                                                                                         |
   | Total mass of each species                                                                                   |
   +--------------------------------------------------------------------------------------------------------------+
   | :func:`gcpy.benchmark.modules.benchmark_funcs.make_benchmark_mass_accumulation_tables` |br|                  |
   | |br|                                                                                                         |
   | Mass accumulation for each species                                                                           |
   +--------------------------------------------------------------------------------------------------------------+
   | :func:`gcpy.benchmark.modules.benchmark_mass_cons_table.make_benchmark_mass_conservation_table` |br|         |
   | |br|                                                                                                         |
   | Timeseries of the PassiveTracer species                                                                      |
   +--------------------------------------------------------------------------------------------------------------+
   | :func:`gcpy.benchmark.modules.oh_metrics.make_benchmark_oh_metrics` |br|                                     |
   | |br|                                                                                                         |
   | Global OH metrics                                                                                            |
   +--------------------------------------------------------------------------------------------------------------+
   | :func:`gcpy.benchmark.modules.benchmark_funcs.make_benchmark_operations_budget` |br|                         |
   | |br|                                                                                                         |
   | Species mass after each operation                                                                            |
   +--------------------------------------------------------------------------------------------------------------+
   | :func:`gcpy.benchmark.modules.ste_flux.make_benchmark_ste_table` |br|                                        |
   | |br|                                                                                                         |
   | Stratosphere-troposphere flux of O\ :sub:`3`                                                                 |
   +--------------------------------------------------------------------------------------------------------------+
   | :func:`gcpy.benchmark.modules.budget_tt.transport_tracers_budgets` |br|                                      |
   | |br|                                                                                                         |
   | Rn222, Pb210, Be7 budgets (1yr benchmarks only)                                                              |
   +--------------------------------------------------------------------------------------------------------------+

The functions listed above create summary tables for quantities such as
total mass of species, total mass of emissions, and OH metrics.

Many of these functions use pre-defined lists of variables in YAML
files. If one dataset includes a variable but the other dataset does
not, the data for that variable in the latter dataset will be
considered to be NaN and will be plotted as such.

.. _bmk-standalone-scripts:

====================================
Standalone benchmark utility scripts
====================================

Unlike the plotting and tabling functions above, the scripts below are
run directly from the command line rather than being called from
:file:`run_benchmark.py`.

benchmark_gcclassic_stats.py
----------------------------

:mod:`gcpy.benchmark.modules.benchmark_gcclassic_stats` scrapes wall
clock time, peak memory usage, OH metrics, and timer statistics from
the public S3 benchmark artifacts of a 1-month GEOS-Chem Classic
benchmark run, and prints them in a format that can be pasted into the
"GEOS-Chem 1-month Benchmark Stats" Google spreadsheet.

Example:

.. code-block:: console

   $ conda activate gcpy_env
   $ python -m gcpy.benchmark.modules.benchmark_gchp_stats 14.8.0-alpha.5 14.8.0-alpha.6
      
benchmark_gchp_stats.py
-----------------------

:mod:`gcpy.benchmark.modules.benchmark_gchp_stats` scrapes wall clock
time, peak memory usage, OH metrics, and timer statistics from the
public S3 benchmark artifacts of a 1-month GCHP benchmark run, and
prints them in a format that can be pasted into the "GEOS-Chem 1-month
Benchmark Stats" Google spreadsheet.

Unlike GEOS-Chem Classic, GCHP does not wrap its executable in
:literal:`/usr/bin/time -v`, so this script instead parses the
:literal:`Mem/Swap Used (MB)` lines that MAPL prints to the run log,
and reads the :file:`Benchmark_Timers_<ref>_vs_<dev>.txt` table
already produced by
:mod:`gcpy.benchmark.modules.benchmark_scrape_gchp_timers`.

.. code-block:: console

   $ conda activate gcpy_env
   $ python -m gcpy.benchmark.modules.benchmark_gchp_stats 14.8.0-alpha.5 14.8.0-alpha.6

benchmark_species_changes.py
-----------------------------

:file:`gcpy/benchmark/modules/benchmark_species_changes.py` generates
GEOS-Chem wiki-formatted tables listing the species that were added
and removed between two versions (by comparing the species metadata
read from each version's log file), plus a summary table of species
counts by category (total species, dry-deposited, wet-deposited, and
photolyzed) for the Ref and Dev versions, with the change and percent
change between them.

Because :literal:`spcdb_files` takes a list of paths, call
:func:`gcpy.benchmark.modules.benchmark_species_changes.make_benchmark_species_changes_wiki_tables`
directly from Python rather than via the command line:

.. code-block:: python

   from gcpy.benchmark.modules.benchmark_species_changes import (
       make_benchmark_species_changes_wiki_tables
   )

   make_benchmark_species_changes_wiki_tables(
       ref_label="14.4.0",
       ref_log="gcc_14.4.0/14.4.0.log",
       dev_label="14.5.0",
       dev_log="gcc_14.5.0/14.5.0.log",
       spcdb_files=[
           "gcc_14.4.0/species_database.yml",
           "gcc_14.5.0/species_database.yml",
       ],
       output_file="wiki_tables.txt",
   )
