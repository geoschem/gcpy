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
   * - ``run_benchmark.py``
     - Benchmark driver script
   * - ``benchmark_slurm.sh``
     - Bash script to submit ``run_benchmark.py`` as a SLURM batch job
   * - ``cloud/``
     - Directory containing template config files (in YAML format) for
       1-hour and 1-month benchmark plot jobs on the AWS cloud.
   * - ``config/``
     - Directory containing editable config files (in YAML format) for
       1-month and 1-year benchmark plot jobs.
   * - ``__init__.py``
     - Python import script
   * - ``modules/``
     - Contains Python modules imported into the ``run_benchmark.py`` script.
   * - ``README.md``
     - Readme file in Markdown format

.. note::

   As of this writing, the benchmarking scripts still use several
   :ref:`plotting <bmk-funcs-plot>` and :ref:`tabling
   <bmk-funcs-table>` functions from module
   :file:`gcpy.benchmark_funcs`.  We are currently in the process of
   moving the functions contained in  :file:`gcpy.benchmark_funcs` to
   the :file:`gcpy/benchmark/modules` directory.

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
     main_dir: /n/home09/ryantosca/T/BM/1yr/
     results_dir: /n/home09/ryantosca/T/BM/1yr/14.7.0-rc.0/GCClassic/FullChem/BenchmarkResults
     weights_dir: /n/holylfs06/LABS/jacob_lab/Shared/GEOS-CHEM/gcgrid/gcdata/ExtData/GCHP/RegriddingWeights
     #
     # Observational data dirs are on Harvard Cannon, edit if necessary
     #
     obs_data:
       ebas_o3:
         data_dir: /n/jacob_lab/Lab/obs_data_for_bmk/ebas_sfc_o3_2019
         data_label: "O3 (EBAS, 2019)"
       sondes:
         data_dir: /n/jacob_lab/Lab/obs_data_for_bmk/sondes_2010-2019
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
   #   species_metadata: Path to species_database.yml file (in dir)
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
         species_metadata: species_database.yml
         outputs_subdir: OutputDir
         restarts_subdir: Restarts
         logs_subdir: Logs
         logs_template: "log.%Y%m%d"
         bmk_start: "2019-01-01T00:00:00"
         bmk_end: "2020-01-01T00:00:00"
       gchp:
         version: GCHP_ref
         dir: GCHP_ref
         species_metadata: species_database.yml
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
         species_metadata: species_database.yml
         outputs_subdir: OutputDir
         restarts_subdir: Restarts
         logs_subdir: Logs
         logs_template: "log.%Y%m%d"
         bmk_start: "2019-01-01T00:00:00"
         bmk_end: "2020-01-01T00:00:00"
       gchp:
         version: GCHP_dev
         dir: GCHP_dev
         species_metadata: species_database.yml
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
        plot_options:
          by_spc_cat: True
          by_hco_cat: True
        #
        # Benchmark tables
        #
        emis_table: True
        mass_accum_table: False
        mass_table: True
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

#. Run the :file:`run.benchmark.py` script.  You may do this in 2
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
         #SBATCH -t 0-4:00
         #SBATCH -p sapphire,huce_cascade,seas_compute,shared
         #SBATCH --mem=100000
         #SBATCH --mail-type=END

         #============================================================================
         # This us a sample SLURM script that you can use to run the GCPy
         # benchmark plotting code as a SLURM batch job.
         #
         # You can modify the SLURM parameters above for your setup.
         #
         # Tip: Using less cores can reduce the amount of memory required.
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
         #config="1mo_benchmark.yml"
         config="1yr_fullchem_benchmark.yml"
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
   :widths: 40 60

   * - Function
     - Plot that it creates
   * - :ref:`bmk-funcs-plot-aod`
     - Aerosol optical depth
   * - :ref:`bmk-funcs-plot-conc`
     - Species concentrations
   * - :ref:`bmk-funcs-plot-dryd`
     - Dry deposition velocities
   * - :ref:`bmk-funcs-plot-emis`
     - Emissions (by species and category)
   * - :ref:`bmk-funcs-plot-jvalue`
     - J-values (photolysis)
   * - :ref:`bmk-funcs-plot-wetdep`
     - Wet deposition of soluble species

.. list-table:: **Functions creating six-panel comparison plots**
   :align: center
   :header-rows: 1
   :widths: 40 60

   * - Function
     - Plot that it creates
   * - :ref:`bmk-funcs-plot-aod`
     - Aerosol optical depth
   * - :ref:`bmk-funcs-plot-conc`
     - Species concentrations
   * - :ref:`bmk-funcs-plot-dryd`
     - Dry deposition velocities
   * - :ref:`bmk-funcs-plot-emis`
     - Emissions (by species and category)
   * - :ref:`bmk-funcs-plot-jvalue`
     - J-values (photolysis)
   * - :ref:`bmk-funcs-plot-wetdep`
     - Wet deposition of soluble species

.. list-table:: **Functions creating model vs. observation plots**
   :align: center
   :header-rows: 1
   :widths: 40 60

   * - Function
     - Plot that it creates
   * - :ref:`bmk-funcs-plot-mvo`
     - Modeled ozone vs. surface observations
   * - :ref:`bmk-funcs-plot-mvs`
     - Vertical profiles of modeled ozone vs. ozonesondes

The functions listed above create comparison plots of most GEOS-Chem
output variables divided into specific categories, e.g. species
categories such as :literal:`Aerosols` or :literal:`Bromine` for the
:literal:`SpeciesConcVV` diagnostic. In eachcategory, these function
create single level PDFs for the surface and 500hPa and zonal
mean PDFs for the entire atmosphere and only the stratosphere (defined
a 1-100hPa). For :code:`make_benchmark_emis_plots()`, only single
level plots at the surface are produced. All of these plotting
functions include bookmarks within the generated PDFs that point to
the pages containing each plotted quantity. Thus these functions serve
as tools for quickly creating comprehensive plots comparing two
GEOS-Chem runs. These functions are used to create the publicly
available plots for 1-month and 1-year benchmarks of new versions of
GEOS-Chem.

Many of the plotting functions listed above use pre-defined lists of
variables in YAML files. If one dataset includes a variable but the
other dataset does not, the data for that variable in the latter
dataset will be considered to be NaN and will be plotted as such.

.. _bmk-funcs-plot-aod:

make_benchmark_aod_plots
------------------------

**Located in module:** :file:`gcpy.benchmark.modules.benchmark_funcs`

This function creates column optical depth plots using the Aerosols
diagnostic output.

.. code-block:: python

   def make_benchmark_aod_plots(
           ref,
           refstr,
           dev,
           devstr,
           spcdb_files,
           varlist=None,
           dst="./benchmark",
           subdst=None,
           cmpres=None,
           overwrite=False,
           verbose=False,
           log_color_scale=False,
           sigdiff_files=None,
           weightsdir='.',
           n_job=-1,
           time_mean=False,
   ):
       """
       Creates PDF files containing plots of column aerosol optical
       depths (AODs) for model benchmarking purposes.

       Args:
           ref: str
               Path name for the "Ref" (aka "Reference") data set.
           refstr: str
               A string to describe ref (e.g. version number)
           dev: str
               Path name for the "Dev" (aka "Development") data set.
               This data set will be compared against the "Reference"
               data set.
           devstr: str
               A string to describe dev (e.g. version number)
           spcdb_files : list
               Paths to species_database.yml files in Ref & Dev rundirs

       Keyword Args (optional):
           varlist: list of str
               List of AOD variables to plot.  If not passed, then all
               AOD variables common to both Dev and Ref will be plotted.
               Use the varlist argument to restrict the number of
               variables plotted to the pdf file when debugging.
               Default value: None
           dst: str
               A string denoting the destination folder where a
               PDF file  containing plots will be written.
               Default value: ./benchmark.
           subdst: str
               A string denoting the sub-directory of dst where PDF
               files containing plots will be written.  In practice,
               subdst is only needed for the 1-year benchmark output,
               and denotes a date string (such as "Jan2016") that
               corresponds to the month that is being plotted.
               Default value: None
           cmpres: string
               Grid resolution at which to compare ref and dev data, e.g. '1x1.25'
           overwrite: bool
               Set this flag to True to overwrite files in the
               destination folder (specified by the dst argument).
               Default value: False.
           verbose: bool
               Set this flag to True to print extra informational output.
               Default value: False
           log_color_scale: bool
               Set this flag to True to enable plotting data (not diffs)
               on a log color scale.
               Default value: False
           sigdiff_files: list of str
               Filenames that will contain the list of quantities having
               having significant differences in the column AOD plots.
               These lists are needed in order to fill out the benchmark
               approval forms.
               Default value: None
           weightsdir: str
               Directory in which to place (and possibly reuse) xESMF regridder
               netCDF files.
               Default value: '.'
           n_job: int
               Defines the number of simultaneous workers for parallel plotting.
               Set to 1 to disable parallel plotting. Value of -1 allows the
               application to decide.
               Default value: -1
           time_mean : bool
               Determines if we should average the datasets over time
               Default value: False
       """

.. _bmk-funcs-plot-conc:

make_benchmark_conc_plots
-------------------------

**Located in module:** :file:`gcpy.benchmark.modules.benchmark_funcs`

Generates species concentration plots using the :literal:`SpeciesConc`
diagnostic output by default.  In particular:

- This function is the only benchmark plotting function that supports
  diff-of-diffs plotting, in which 4 datasets are passed and the
  differences between two groups of :literal:`Ref` datasets vs. two
  groups of :literal:`Dev` datasets is plotted (typically used for
  comparing changes in GCHP vs. changes in GEOS-Chem Classic across
  model versions). |br|
  |br|

- This is also the only benchmark plotting function that sends plots
  to separate folders based on category (as denoted by the
  plot_by_spc_cat flag). The full list of species categories is
  denoted in `benchmark_categories.yml
  <https://github.com/geoschem/gcpy/blob/dev/gcpy/benchmark_categories.yml>`_
  (included in GCPy). |br|
  |br|

- In this function, parallelization occurs at the species category
  level. In all other functions, parallelization occurs within calls
  to :code:`compare_single_level()`  and :code:`compare_zonal_mean()`.=

.. code-block:: python

   def make_benchmark_conc_plots(
           ref,
           refstr,
           dev,
           devstr,
           spcdb_files,
           dst="./benchmark",
           subdst=None,
           overwrite=False,
           verbose=False,
           collection="SpeciesConc",
           benchmark_type="FullChemBenchmark",
           cmpres=None,
           plot_by_spc_cat=True,
           restrict_cats=[],
           plots=["sfc", "500hpa", "zonalmean"],
           use_cmap_RdBu=False,
           log_color_scale=False,
           sigdiff_files=None,
           normalize_by_area=False,
           cats_in_ugm3=["Aerosols", "Secondary_Organic_Aerosols"],
           areas=None,
           refmet=None,
           devmet=None,
           weightsdir='.',
           n_job=-1,
           second_ref=None,
           second_dev=None,
           time_mean=False,
   ):
       """
       Creates PDF files containing plots of species concentration
       for model benchmarking purposes.

       Args:
           ref: str
               Path name for the "Ref" (aka "Reference") data set.
           refstr: str
               A string to describe ref (e.g. version number)
           dev: str
               Path name for the "Dev" (aka "Development") data set.
               This data set will be compared against the "Reference"
               data set.
           devstr: str
               A string to describe dev (e.g. version number)
           spcdb_files : list
               Paths to species_database.yml files in Ref & Dev rundirs

       Keyword Args (optional):
           dst: str
               A string denoting the destination folder where a PDF
               file containing plots will be written.
               Default value: ./benchmark
           subdst: str
               A string denoting the sub-directory of dst where PDF
               files containing plots will be written.  In practice,
               subdst is only needed for the 1-year benchmark output,
               and denotes a date string (such as "Jan2016") that
               corresponds to the month that is being plotted.
               Default value: None
           overwrite: bool
               Set this flag to True to overwrite files in the
               destination folder (specified by the dst argument).
               Default value: False
           verbose: bool
               Set this flag to True to print extra informational output.
               Default value: False
           collection: str
               Name of collection to use for plotting.
               Default value: "SpeciesConc"
           benchmark_type: str
               A string denoting the type of benchmark output to plot, options are
               FullChemBenchmark, TransportTracersBenchmark, or CH4Benchmark.
               Default value: "FullChemBenchmark"
           cmpres: string
               Grid resolution at which to compare ref and dev data, e.g. '1x1.25'
           plot_by_spc_cat: logical
               Set this flag to False to send plots to one file rather
               than separate file per category.
               Default value: True
           restrict_cats: list of strings
               List of benchmark categories in benchmark_categories.yml to make
               plots for. If empty, plots are made for all categories.
               Default value: empty
           plots: list of strings
               List of plot types to create.
               Default value: ['sfc', '500hpa', 'zonalmean']
           log_color_scale: bool
               Set this flag to True to enable plotting data (not diffs)
               on a log color scale.
               Default value: False
           normalize_by_area: bool
               Set this flag to true to enable normalization of data
               by surfacea area (i.e. kg s-1 --> kg s-1 m-2).
               Default value: False
           cats_in_ugm3: list of str
               List of benchmark categories to to convert to ug/m3
               Default value: ["Aerosols", "Secondary_Organic_Aerosols"]
           areas: dict of xarray DataArray:
               Grid box surface areas in m2 on Ref and Dev grids.
               Default value: None
           refmet: str
               Path name for ref meteorology
               Default value: None
           devmet: str
               Path name for dev meteorology
               Default value: None
           sigdiff_files: list of str
               Filenames that will contain the lists of species having
               significant differences in the 'sfc', '500hpa', and
               'zonalmean' plots.  These lists are needed in order to
               fill out the benchmark approval forms.
               Default value: None
           weightsdir: str
               Directory in which to place (and possibly reuse) xESMF regridder
               netCDF files.
               Default value: '.'
           n_job: int
               Defines the number of simultaneous workers for parallel plotting.
               Set to 1 to disable parallel plotting. Value of -1 allows the
               application to decide.
               Default value: -1
           second_ref: str
               Path name for a second "Ref" (aka "Reference") data set for
               diff-of-diffs plotting. This dataset should have the same model
               type and grid as ref.
               Default value: None
           second_dev: str
               Path name for a second "Ref" (aka "Reference") data set for
               diff-of-diffs plotting. This dataset should have the same model
               type and grid as ref.
               Default value: None
           time_mean : bool
               Determines if we should average the datasets over time
               Default value: False
       """

.. _bmk-funcs-plot-dryd:

make_benchmark_drydep_plots
---------------------------

**Located in module:** :file:`gcpy.benchmark.modules.benchmark_drydep`

Generates plots of dry deposition velocities using the GEOS-Chem
:literal:`DryDep` diagnostic output.

.. code-block:: python

   def make_benchmark_drydep_plots(
           ref,
           refstr,
           dev,
           devstr,
           spcdb_files,
           collection="DryDep",
           dst="./benchmark",
           subdst=None,
           cmpres=None,
           overwrite=False,
           verbose=False,
           log_color_scale=False,
           weightsdir=".",
           sigdiff_files=None,
           n_job=-1,
           time_mean=False,
           varlist=None,
   ):
       """
       Creates six-panel comparison plots (PDF format) from GEOS-Chem
       benchmark simualtion output.  Can be used with data collections
       that do not require special handling (e.g. concentrations).

       Args:
           ref: str
               Path name for the "Ref" (aka "Reference") data set.
           refstr: str
               A string to describe ref (e.g. version number)
           dev: str
               Path name for the "Dev" (aka "Development") data set.
               This data set will be compared against the "Reference"
               data set.
           devstr: str
               A string to describe dev (e.g. version number)
           spcdb_files : list
               Paths to species_database.yml files in Ref & Dev rundirs

       Keyword Args (optional):
           collection : str
               Name of the diagnostic collection (e.g. "DryDep")
           dst: str
               A string denoting the destination folder where a PDF
               file containing plots will be written.
               Default value: ./benchmark
           subdst: str
               A string denoting the sub-directory of dst where PDF
               files containing plots will be written.  In practice,
               subdst is only needed for the 1-year benchmark output,
               and denotes a date string (such as "Jan2016") that
               corresponds to the month that is being plotted.
               Default value: None
           benchmark_type: str
               A string denoting the type of benchmark output to plot, options are
               FullChemBenchmark, TransportTracersBenchmark, or CH4Benchmark.
               Default value: "FullChemBenchmark"
           overwrite: bool
               Set this flag to True to overwrite files in the
               destination folder (specified by the dst argument).
               Default value: False.
           verbose: bool
               Set this flag to True to print extra informational output.
               Default value: False.
           n_job: int
               Defines the number of simultaneous workers for parallel plotting.
               Set to 1 to disable parallel plotting. Value of -1 allows the
               application to decide.
               Default value: -1
           time_mean : bool
               Determines if we should average the datasets over time
               Default value: False
           varlist: list of str
               List of variables to plot.  If varlist is None, then
               all common variables in Ref & Dev will be plotted.
       """

.. _bmk-funcs-plot-emis:

make_benchmark_emis_plots
-------------------------

**Located in module:** :file:`gcpy.benchmark.modules.benchmark_funcs`

Generates plots of total emissions using output from
:file:`HEMCO_diagnostics.*` (for GEOS-Chem Classic) and/or
:file:`GCHP.Emissions.*` output files.

.. code-block:: python

   def make_benchmark_emis_plots(
           ref,
           refstr,
           dev,
           devstr,
           spcdb_files,
           dst="./benchmark",
           subdst=None,
           plot_by_spc_cat=False,
           plot_by_hco_cat=False,
           benchmark_type="FullChemBenchmark",
           cmpres=None,
           overwrite=False,
           verbose=False,
           flip_ref=False,
           flip_dev=False,
           log_color_scale=False,
           sigdiff_files=None,
           weightsdir='.',
           n_job=-1,
           time_mean=False,
   ):
       """
       Creates PDF files containing plots of emissions for model
       benchmarking purposes. This function is compatible with benchmark
       simulation output only. It is not compatible with transport tracers
       emissions diagnostics.

       Args:
           ref: str
               Path name for the "Ref" (aka "Reference") data set.
           refstr: str
               A string to describe ref (e.g. version number)
           dev: str
               Path name for the "Dev" (aka "Development") data set.
               This data set will be compared against the "Reference"
               data set.
           devstr: str
               A string to describe dev (e.g. version number)
           spcdb_files : list
               Paths to species_database.yml files in Ref & Dev rundirs

       Keyword Args (optional):
           dst: str
               A string denoting the destination folder where
               PDF files containing plots will be written.
               Default value: './benchmark
           subdst: str
               A string denoting the sub-directory of dst where PDF
               files containing plots will be written.  In practice,
               and denotes a date string (such as "Jan2016") that
               corresponds to the month that is being plotted.
               Default value: None
           plot_by_spc_cat: bool
               Set this flag to True to separate plots into PDF files
               according to the benchmark species categories (e.g. Oxidants,
               Aerosols, Nitrogen, etc.)  These categories are specified
               in the YAML file benchmark_categories.yml.
               Default value: False
           plot_by_hco_cat: bool
               Set this flag to True to separate plots into PDF files
               according to HEMCO emissions categories (e.g. Anthro,
               Aircraft, Bioburn, etc.)
               Default value: False
           benchmark_type: str
               A string denoting the type of benchmark output to plot, options are
               FullChemBenchmark, TransportTracersBenchmark, or CH4Benchmark.
               Default value: "FullChemBenchmark"
           cmpres: string
               Grid resolution at which to compare ref and dev data, e.g. '1x1.25'
           overwrite: bool
               Set this flag to True to overwrite files in the
               destination folder (specified by the dst argument).
               Default value: False
           verbose: bool
               Set this flag to True to print extra informational output.
               Default value: False
           flip_ref: bool
               Set this flag to True to reverse the vertical level
               ordering in the "Ref" dataset (in case "Ref" starts
               from the top of atmosphere instead of the surface).
               Default value: False
           flip_dev: bool
               Set this flag to True to reverse the vertical level
               ordering in the "Dev" dataset (in case "Dev" starts
               from the top of atmosphere instead of the surface).
               Default value: False
           log_color_scale: bool
               Set this flag to True to enable plotting data (not diffs)
               on a log color scale.
               Default value: False
            sigdiff_files: list of str
               Filenames that will contain the lists of species having
               significant differences in the 'sfc', '500hpa', and
               'zonalmean' plots.  These lists are needed in order to
               fill out the benchmark approval forms.
               Default value: None
           weightsdir: str
               Directory in which to place (and possibly reuse) xESMF regridder
               netCDF files.
               Default value: '.'
           n_job: int
               Defines the number of simultaneous workers for parallel plotting.
               Set to 1 to disable parallel plotting.
               Value of -1 allows the application to decide.
               Default value: -1
           time_mean : bool
               Determines if we should average the datasets over time
               Default value: False

       Remarks:
           (1) If both plot_by_spc_cat and plot_by_hco_cat are
               False, then all emission plots will be placed into the
               same PDF file.

           (2) Emissions that are 3-dimensional will be plotted as
               column sums.
       """

.. _bmk-funcs-plot-jvalue:

make_benchmark_jvalue_plots
---------------------------

**Located in module:** :file:`gcpy.benchmark.modules.benchmark_funcs`

Generates plots of J-values using the GEOS-Chem :literal:`JValues`
diagnostic output.

.. code-block:: python

   def make_benchmark_jvalue_plots(
           ref,
           refstr,
           dev,
           devstr,
           spcdb_files,
           varlist=None,
           dst="./benchmark",
           subdst=None,
           local_noon_jvalues=False,
           cmpres=None,
           plots=["sfc", "500hpa", "zonalmean"],
           overwrite=False,
           verbose=False,
           flip_ref=False,
           flip_dev=False,
           log_color_scale=False,
           sigdiff_files=None,
           weightsdir='.',
           n_job=-1,
           time_mean=False,
   ):
       """
       Creates PDF files containing plots of J-values for model
       benchmarking purposes.

       Args:
           ref: str
               Path name for the "Ref" (aka "Reference") data set.
           refstr: str
               A string to describe ref (e.g. version number)
           dev: str
               Path name for the "Dev" (aka "Development") data set.
               This data set will be compared against the "Reference"
               data set.
           devstr: str
               A string to describe dev (e.g. version number)
           spcdb_files : list
               Paths to species_database.yml files in Ref & Dev rundirs

       Keyword Args (optional):
           varlist: list of str
               List of J-value variables to plot.  If not passed,
               then all J-value variables common to both dev
               and ref will be plotted.  The varlist argument can be
               a useful way of restricting the number of variables
               plotted to the pdf file when debugging.
               Default value: None
           dst: str
               A string denoting the destination folder where a
               PDF file  containing plots will be written.
               Default value: ./benchmark.
           subdst: str
               A string denoting the sub-directory of dst where PDF
               files containing plots will be written.  In practice,
               subdst is only needed for the 1-year benchmark output,
               and denotes a date string (such as "Jan2016") that
               corresponds to the month that is being plotted.
               Default value: None
           local_noon_jvalues: bool
               Set this flag to plot local noon J-values.  This will
               divide all J-value variables by the JNoonFrac counter,
               which is the fraction of the time that it was local noon
               at each location.
               Default value: False
           cmpres: string
               Grid resolution at which to compare ref and dev data, e.g. '1x1.25'
           plots: list of strings
               List of plot types to create.
               Default value: ['sfc', '500hpa', 'zonalmean']
           overwrite: bool
               Set this flag to True to overwrite files in the
               destination folder (specified by the dst argument).
               Default value: False.
           verbose: bool
               Set this flag to True to print extra informational output.
               Default value: False
           flip_ref: bool
               Set this flag to True to reverse the vertical level
               ordering in the "Ref" dataset (in case "Ref" starts
               from the top of atmosphere instead of the surface).
               Default value: False
           flip_dev: bool
               Set this flag to True to reverse the vertical level
               ordering in the "Dev" dataset (in case "Dev" starts
               from the top of atmosphere instead of the surface).
               Default value: False
           log_color_scale: bool
               Set this flag to True if you wish to enable plotting data
               (not diffs) on a log color scale.
               Default value: False
           sigdiff_files: list of str
               Filenames that will contain the lists of J-values having
               significant differences in the 'sfc', '500hpa', and
               'zonalmean' plots.  These lists are needed in order to
               fill out the benchmark approval forms.
               Default value: None
           weightsdir: str
               Directory in which to place (and possibly reuse) xESMF regridder
               netCDF files.
               Default value: '.'
           n_job: int
               Defines the number of simultaneous workers for parallel plotting.
               Set to 1 to disable parallel plotting. Value of -1 allows the
               application to decide.
               Default value: -1
           time_mean : bool
               Determines if we should average the datasets over time
               Default value: False

       Remarks:
            Will create 4 files containing J-value plots:
               (1 ) Surface values
               (2 ) 500 hPa values
               (3a) Full-column zonal mean values.
               (3b) Stratospheric zonal mean values
            These can be toggled on/off with the plots keyword argument.

            At present, we do not yet have the capability to split the
            plots up into separate files per category (e.g. Oxidants,
            Aerosols, etc.).  This is primarily due to the fact that
            we archive J-values from GEOS-Chem for individual species
            but not family species.  We could attempt to add this
            functionality later if there is sufficient demand.
       """

.. _bmk-funcs-plot-wetdep:

make_benchmark_wetdep_plots
---------------------------

**Located in module:** :file:`gcpy.benchmark.modules.benchmark_funcs`

Generates plots of wet deposition using the GEOS-Chem
:literal:`WetLossConv` and :literal:`WetLossLS` diagnostic outputs.
It is currently primarily used for 1-Year Transport Tracer benchmarks,
plotting values for the following species as defined in
`benchmark_categories.yml
<https://github.com/geoschem/gcpy/blob/dev/gcpy/benchmark/modules/benchmark_categories.yml>`_

.. code-block:: python

   def make_benchmark_wetdep_plots(
           ref,
           refstr,
           dev,
           devstr,
           collection,
           spcdb_files,
           dst="./benchmark",
           cmpres=None,
           datestr=None,
           overwrite=False,
           verbose=False,
           benchmark_type="TransportTracersBenchmark",
           plots=["sfc", "500hpa", "zonalmean"],
           log_color_scale=False,
           normalize_by_area=False,
           areas=None,
           refmet=None,
           devmet=None,
           weightsdir='.',
           n_job=-1,
           time_mean=False,
   ):
       """
       Creates PDF files containing plots of species concentration
       for model benchmarking purposes.

       Args:
           ref: str
               Path name for the "Ref" (aka "Reference") data set.
           refstr: str
               A string to describe ref (e.g. version number)
           dev: str
               Path name for the "Dev" (aka "Development") data set.
               This data set will be compared against the "Reference"
               data set.
           devstr: str
               A string to describe dev (e.g. version number)
           collection: str
               String name of collection to plot comparisons for.
           spcdb_files : list
               Paths to species_database.yml files in Ref & Dev rundirs

       Keyword Args (optional):
           dst: str
               A string denoting the destination folder where a PDF
               file containing plots will be written.
               Default value: ./benchmark
           datestr: str
               A string with date information to be included in both the
               plot pdf filename and as a destination folder subdirectory
               for writing plots
               Default value: None
           benchmark_type: str
               A string denoting the type of benchmark output to plot, options are
               FullChemBenchmark, TransportTracersBenchmark, or CH4Benchmark.
               Default value: "FullChemBenchmark"
           overwrite: bool
               Set this flag to True to overwrite files in the
               destination folder (specified by the dst argument).
               Default value: False.
           verbose: bool
               Set this flag to True to print extra informational output.
               Default value: False.
           plots: list of strings
               List of plot types to create.
               Default value: ['sfc', '500hpa', 'zonalmean']
           normalize_by_area: bool
               Set this flag to true to enable normalization of data
               by surfacea area (i.e. kg s-1 --> kg s-1 m-2).
               Default value: False
           areas: dict of xarray DataArray:
               Grid box surface areas in m2 on Ref and Dev grids.
               Default value: None
           refmet: str
               Path name for ref meteorology
               Default value: None
           devmet: str
               Path name for dev meteorology
               Default value: None
           n_job: int
               Defines the number of simultaneous workers for parallel plotting.
               Set to 1 to disable parallel plotting. Value of -1 allows the
               application to decide.
               Default value: -1
           time_mean : bool
               Determines if we should average the datasets over time
               Default value: False
       """

.. _bmk-funcs-plot-mvo:

make_benchmark_model_vs_obs_plots
---------------------------------

**Located in module:** :file:`gcpy.benchmark.modules.benchmark_models_vs_obs`

Gnerates plots of monthly-averaged modeled surface
ozone concentrations (using the GEOS-Chem :literal:`SpeciesConc`
diagnostic outputs) vs. the `EBAS 2019 <https://ebas-data.nilu.no/>`_
observations.

.. note::

   Model vs. observation plots are only available in 1-year
   full-chemistry benchmarks.

.. code-block:: python

   def make_benchmark_models_vs_obs_plots(
           obs_filepaths,
           obs_label,
           ref_filepaths,
           ref_label,
           dev_filepaths,
           dev_label,
           varname="SpeciesConcVV_O3",
           dst="./benchmark",
           verbose=False,
           overwrite=False
   ):
       """
       Driver routine to create model vs. observation plots.

       Args:
       obs_filepaths : str|list : Path(s) to the observational data.
       obs_label     : str      : Label for the observational data
       ref_filepaths : str      : Paths to the Ref model data
       ref_label     : str      : Label for the Ref model data
       dev_filepaths : str      : Paths to the Dev model data
       dev_label     : str      : Label for the Dev model data
       varname       : str      : Variable name for model data
       dst           : str      : Destination folder for plots
       verbose       : bool     : Toggles verbose output on/off
       overwrite     : bool     : Toggles overwriting contents of dst

.. _bmk-funcs-plot-mvs:

make_benchmark_model_vs_sondes_plots
------------------------------------

**Located in module:** :file:`gcpy.benchmark.modules.benchmark_models_vs_sondes`

.. note::

   Model vs. ozonesonde plots are only available in 1-year
   full-chemistry benchmarks.

Generates vertical profiles of modeled ozone concentrations (using the
GEOS-Chem :literal:`SpeciesConc` diagnostic outputs) vs. ozonesonde
observations.

.. code-block:: python

   def make_benchmark_models_vs_sondes_plots(
           obs_data_file,
           obs_site_file,
           ref_filepaths,
           ref_label,
           dev_filepaths,
           dev_label,
           dst="./benchmark",
           overwrite=False,
           varname="SpeciesConcVV_O3",

       ):
       """
       Creates plots of sonde data vs. GEOS-Chem output.  For use in the
       1-year benchmark plotting workflow.

       Args
       obs_data_file : str      : File containing sonde data
       obs_site_file : str      : File containing sonde site metadata
       ref_filepaths : str|list : Files for the GEOS-Chem Ref version
       ref_label     : str      : GEOS-Chem Ref version label
       dev_filepaths : str|list : Files for the GEOS-Chem Dev version
       dev_label     : str      : GEOS-Chem Dev version label

       Keyword Args
       dst           : str      : Folder where PDF w/ plots will be created
       overwrite     : bool     : Overwrite contents of dst folder?
       varname       : str      : GEOS-Chem diagnostic variable name
       verbose       : bool     : Activate verbose printout?
       """

.. _bmk-funcs-table:

===========================
Benchmark tabling functions
===========================

.. list-table:: **Functions creating summary tables**
   :align: center
   :header-rows: 1
   :widths: 40 60

   * - Function
     - Table that it creates
   * - :ref:`bmk-funcs-table-oxbdg`
     - Ox budget (1yr benchmarks only)
   * - :ref:`bmk-funcs-table-aer`
     - Global aerosol burdens (1yr benchmarks only)
   * - :ref:`bmk-funcs-table-emis`
     - Emissions (by species & inventory)
   * - :ref:`bmk-funcs-table-gcc-timers`
     - GEOS-Chem Classic timers output
   * - :ref:`bmk-funcs-table-gchp-timers`
     - GCHP timers output
   * - :ref:`bmk-funcs-table-mass`
     - Total mass of each species
   * - :ref:`bmk-funcs-table-accum`
     - Mass accumulation for each species
   * - :ref:`bmk-funcs-table-cons`
     - Timeseries of the PassiveTracer species
   * - :ref:`bmk-funcs-table-oh`
     - Global OH metrics
   * - :ref:`bmk-funcs-table-ops`
     - Species mass after each operation
   * - :ref:`bmk-funcs-table-ttbdg`
     - Rn222, Pb210, Be7 budgets (1yr benchmarks only)

The functions listed above create summary tables for quantities such as
total mass of species, total mass of emissions, and OH metrics.

Many of these functions use pre-defined lists of variables in YAML
files. If one dataset includes a variable but the other dataset does
not, the data for that variable in the latter dataset will be
considered to be NaN and will be plotted as such.

.. _bmk-funcs-table-oxbdg:

global_ox_budget
----------------

**Located in module:** :file:`gcpy.benchmark.modules.budget_ox`

Generates a  budget table for the Ox (odd ozone) family from 1-year
full-chemistry benchmark output.

.. code-block:: python

   def global_ox_budget(
           devstr,
           devdir,
           devrstdir,
           year,
           spcdb_file,
           dst='./1yr_benchmark',
           overwrite=True,
           is_gchp=False,
           gchp_res="c24",
           gchp_is_pre_14_0=False
   ):
       """
       Main program to compute Ox budgets from GEOS-Chem Classic or
       GCHP benchmark simulations.

       Args
       devstr           : str  : Label for the Dev version
       devdir           : str  : Path to the Dev data directory
       devrstdir        : str  : Path to the Dev restart file directory
       year             : int  : Year of the benchmark simulation
       spcdb_file       : str  : Path to the Dev species_database.yml file

       Keyword Args
       dst              : str  : Directory where tables will be written
       overwrite        : bool : Should existing tables should be overwritten?
       is_gchp          : bool : Is Dev from a GCHP benchmark simulation?
       gchp_res         : str  : GCHP resolution string (e.g. "c24")
       gchp_is_pre_14_0 : bool : Is Dev from a GCHP version prior to 14.0.0?
       """

.. _bmk-funcs-table-aer:

make_benchmark_aerosol_tables
-----------------------------

**Located in module:** :file:`gcpy.benchmark.modules.benchmark_funcs`

Generates a table of global aerosol budgets and burdens from GEOS-Chem
1-year full-chemistry benchmark simulation output.

.. code-block:: python

   def make_benchmark_aerosol_tables(
           devdir,
           devlist_aero,
           devlist_spc,
           devlist_met,
           devstr,
           year,
           days_per_mon,
           spcdb_files,
           dst='./benchmark',
           overwrite=False,
           is_gchp=False,
   ):
       """
       Compute FullChemBenchmark aerosol budgets & burdens

       Args:
           devdir: str
               Path to development ("Dev") data directory
           devlist_aero: list of str
               List of Aerosols collection files (different months)
           devlist_spc: list of str
               List of SpeciesConc collection files (different months)
           devlist_met: list of str
               List of meteorology collection files (different months)
           devstr: str
               Descriptive string for datasets (e.g. version number)
           year: str
               The year of the benchmark simulation (e.g. '2016').
           days_per_month: list of int
               List of number of days per month for all months
           spcdb_files : list
               Paths to species_database.yml files in Ref & Dev rundirs

       Keyword Args (optional):
           dst: str
               Directory where budget tables will be created.
               Default value: './benchmark'
           overwrite: bool
               Overwrite burden & budget tables? (default=True)
               Default value: False
           is_gchp: bool
               Whether datasets are for GCHP
               Default value: False

       """

.. _bmk-funcs-table-emis:

make_benchmark_emis_tables
--------------------------

**Located in module:** :file:`gcpy.benchmark.modules.benchmark_funcs`

Generates tables of emissions (by species and by inventory) from the
the :literal:`HEMCO_diagnostics*` outputs.

.. code-block:: python

   def make_benchmark_emis_tables(
           reflist,
           refstr,
           devlist,
           devstr,
           spcdb_files,
           dst="./benchmark",
           benchmark_type="FullChemBenchmark",
           refmet=None,
           devmet=None,
           overwrite=False,
           ref_interval=[2678400.0],
           dev_interval=[2678400.0],
   ):
       """
       Creates a text file containing emission totals by species and
       category for benchmarking purposes.

       Args:
           reflist: list of str
                List with the path names of the emissions file or files
                (multiple months) that will constitute the "Ref"
                (aka "Reference") data set.
           refstr: str
               A string to describe ref (e.g. version number)
           devlist: list of str
                List with the path names of the emissions file or files
                (multiple months) that will constitute the "Dev"
                (aka "Development") data set
           devstr: str
               A string to describe dev (e.g. version number)
           spcdb_files : list
               Paths to species_database.yml files in Ref & Dev rundirs

       Keyword Args (optional):
           dst: str
               A string denoting the destination folder where the file
               containing emissions totals will be written.
               Default value: ./benchmark
           benchmark_type: str
               A string denoting the type of benchmark output to plot, options are
               FullChemBenchmark, TransportTracersBenchmark or CH4Benchmark.
               Default value: "FullChemBenchmark"
           refmet: str
               Path name for ref meteorology
               Default value: None
           devmet: str
               Path name for dev meteorology
               Default value: None
           overwrite: bool
               Set this flag to True to overwrite files in the
               destination folder (specified by the dst argument).
               Default value: False
           ref_interval: list of float
               The length of the ref data interval in seconds. By default, interval
               is set to [2678400.0], which is the number of seconds in July
               (our 1-month benchmarking month).
               Default value: [2678400.0]
           dev_interval: list of float
               The length of the dev data interval in seconds. By default, interval
               is set to [2678400.0], which is the number of seconds in July
               (our 1-month benchmarking month).
               Default value: [2678400.0]

       """

.. _bmk-funcs-table-gcc-timers:

make_benchmark_gcclassic_timing_table
-------------------------------------

**Located in module:**
:file:`gcpy.benchmark.modules.benchmark_scrape_gcclassic_timers`

Generates a comparison table of GEOS-Chem Classic timer values.  This
can be used to determine if computational bottlenecks have been
introduced.

.. code-block:: python

   def make_benchmark_gcclassic_timing_table(
           ref_files,
           ref_label,
           dev_files,
           dev_label,
           dst="./benchmark",
           overwrite=False,
   ):
       """
       Creates a table of timing information for GEOS-Chem Classic
       benchmark simulations given one or more JSON and/or text files
       as input.

       Args
       ref_files : str|list : File(s) with timing info from the "Ref" model
       ref_label : str      : Version string for the "Ref" model
       dev_files : str|list : File(s) with timing info from the "Ref" model
       dev_label : str      : Version string for the "Dev" model

       Kwargs
       dst       : str      : Directory where output will be written
       overwrite : bool     : Overwrite existing files? (default: False)
       """

.. _bmk-funcs-table-gchp-timers:

make_benchmark_gchp_timing_table
--------------------------------

**Located in module:**
:file:`gcpy.benchmark.modules.benchmark_scrape_gchp_timers`

Generates a comparison table of GCHP Classic timer values.  This
can be used to determine if computational bottlenecks have been
introduced.

.. code-block:: python

   def make_benchmark_gchp_timing_table(
           ref_files,
           ref_label,
           dev_files,
           dev_label,
           dst="./benchmark",
           overwrite=False,
   ):
       """
       Creates a table of timing information for GCHP benchmark
       simulations given one or more text files as input.

       Args
       ref_files : str|list : File(s) with timing info from the "Ref" model
       ref_label : str      : Version string for the "Ref" model
       dev_files : str|list : File(s) with timing info from the "Ref" model
       dev_label : str      : Version string for the "Dev" model

       Kwargs
       dst       : str      : Directory where output will be written
       overwrite : bool     : Overwrite existing files? (default: False)
       """

.. _bmk-funcs-table-mass:

make_benchmark_mass_tables
--------------------------

**Located in module:** :file:`gcpy.benchmark.modules.benchmark_funcs`

Generates a comparison table of total mass for each GEOS-Chem species,
using the GEOS-Chem restart file output.

.. code-block:: python

   def make_benchmark_mass_tables(
           ref,
           refstr,
           dev,
           devstr,
           spcdb_files,
           varlist=None,
           dst="./benchmark",
           subdst=None,
           overwrite=False,
           verbose=False,
           label="at end of simulation",
           ref_met_extra=None,
           dev_met_extra=None
   ):
       """
       Creates a text file containing global mass totals by species and
       category for benchmarking purposes.

       Args:
           reflist: str
               Pathname that will constitute
               the "Ref" (aka "Reference") data set.
           refstr: str
               A string to describe ref (e.g. version number)
           dev: list of str
               Pathname that will constitute
               the "Dev" (aka "Development") data set.  The "Dev"
               data set will be compared against the "Ref" data set.
           devstr: str
               A string to describe dev (e.g. version number)
           spcdb_files : list
               Paths to species_database.yml files in Ref & Dev rundirs

       Keyword Args (optional):
           varlist: list of str
               List of variables to include in the list of totals.
               If omitted, then all variables that are found in either
               "Ref" or "Dev" will be included.  The varlist argument
               can be a useful way of reducing the number of
               variables during debugging and testing.
               Default value: None
           dst: str
               A string denoting the destination folder where the file
               containing emissions totals will be written.
               Default value: ./benchmark
           subdst: str
               A string denoting the sub-directory of dst where PDF
               files containing plots will be written.  In practice,
               subdst is only needed for the 1-year benchmark output,
               and denotes a date string (such as "Jan2016") that
               corresponds to the month that is being plotted.
               Default value: None
           overwrite: bool
               Set this flag to True to overwrite files in the
               destination folder (specified by the dst argument).
               Default value: False
           verbose: bool
               Set this flag to True to print extra informational output.
               Default value: False.
           ref_met_extra: str
               Path to ref Met file containing area data for use with restart files
               which do not contain the Area variable.
               Default value: ''
           dev_met_extra: str
               Path to dev Met file containing area data for use with restart files
               which do not contain the Area variable.
               Default value: ''
       """

.. _bmk-funcs-table-accum:

make_benchmark_mass_accumulation_tables
---------------------------------------

**Located in module:** :file:`gcpy.benchmark.modules.benchmark_funcs`

Generates a comparison table of mass accumulation (i.e. mass difference
between the start and end of Ref and Dev benchmark simulations), using
GEOS-Chem restart files.

.. code-block:: python

   def make_benchmark_mass_accumulation_tables(
           ref_start,
           ref_end,
           refstr,
           refperiodstr,
           dev_start,
           dev_end,
           devstr,
           devperiodstr,
           spcdb_files,
           varlist=None,
           dst="./benchmark",
           subdst=None,
           overwrite=False,
           verbose=False,
           label="at end of simulation",
   ):
       """
       Creates a text file containing global mass totals by species and
       category for benchmarking purposes.

       Args:
           ref_start: list of str
               Pathname that will constitute
               the "Ref" (aka "Reference") data set.
           ref_end: list of str
               Pathname that will constitute
               the "Ref" (aka "Reference") data set.
           refstr: str
               A string to describe ref (e.g. version number)
           refperiodstr: str
               Ref simulation period start and end
           dev_start: list of str
               Pathname that will constitute
               the "Dev" (aka "Development") data set.  The "Dev"
               data set will be compared against the "Ref" data set.
           dev_end: list of str
               Pathname that will constitute
               the "Dev" (aka "Development") data set.  The "Dev"
               data set will be compared against the "Ref" data set.
           devstr: str
               A string to describe dev (e.g. version number)
           devperiodstr: str
               Dev simulation period start and end
           spcdb_files : list
               Paths to species_database.yml files in Ref & Dev rundirs

       Keyword Args (optional):
           varlist: list of str
               List of variables to include in the list of totals.
               If omitted, then all variables that are found in either
               "Ref" or "Dev" will be included.  The varlist argument
               can be a useful way of reducing the number of
               variables during debugging and testing.
               Default value: None
           dst: str
               A string denoting the destination folder where the file
               containing emissions totals will be written.
               Default value: ./benchmark
           subdst: str
               A string denoting the sub-directory of dst where PDF
               files containing plots will be written.  In practice,
               subdst is only needed for the 1-year benchmark output,
               and denotes a date string (such as "Jan2016") that
               corresponds to the month that is being plotted.
               Default value: None
           overwrite: bool
               Set this flag to True to overwrite files in the
               destination folder (specified by the dst argument).
               Default value: False
           verbose: bool
               Set this flag to True to print extra informational output.
               Default value: False.
       """

.. _bmk-funcs-table-cons:

make_benchmark_mass_conservation_table
--------------------------------------

**Located in module:** :file:`gcpy.benchmark.modules.benchmark_mass_cons_table`

Generates a timeseries table of the global mass of the
:literal:`PassiveTracer` species.  Usually used with output from
1-year TransportTracers benchmark simulations.  This is an important
check for mass conservation in GEOS-Chem Classic and GCHP.

.. code-block:: python

   def make_benchmark_mass_conservation_table(
           ref_files,
           ref_label,
           dev_files,
           dev_label,
           spcdb_files,
           dst="./benchmark",
           overwrite=False,
           ref_areapath=None,
           dev_areapath=None,
   ):
       """
       Creates a text file containing global mass of passive species
       contained in GEOS-Chem Classic and/or GCHP restart files.

       Args
       ref_files    : list|str : List of files from the Ref model
       ref_label    : str      : Ref version label
       dev_files    : list|str : List of files from the Dev model
       dev_label    : str      : Dev version label
       spcdb_files  : list     : Paths to Ref & Dev species_database.yml files
       dst          : str      : Destination folder for file output
       overwrite    : bool     : Overwrite pre-existing files?
       ref_areapath : list|str : Path to file w/ Ref area data (optional)
       dev_areapath : list|str : Path to file w/ Dev area data (optional)
       """

.. _bmk-funcs-table-oh:

make_benchmark_oh_metrics
-------------------------

**Located in module:** :file:`gcpy.benchmark.modules.oh_metrics`

Generates a table of OH metrics (mean OH concentration,
methyl chloroform lifetime, CH4 lifetime) from the GEOS-Chem
:literal:`Metrics` diagnostic outputs.

.. code-block:: python

   def make_benchmark_oh_metrics(
           ref,
           refstr,
           dev,
           devstr,
           spcdb_files,
           dst="./benchmark",
           overwrite=True,
   ):
       """
       Creates a text file containing metrics of global mean OH, MCF lifetime,
       and CH4 lifetime for benchmarking purposes.

       Args
       ref         : str  : Path name of "Ref" (aka "Reference") data file
       refstr      : str  : Label to describe Ref
       dev         : str  : Path name of "Dev" (aka "Development") data file
       devstr      : str  : Label to describe Dev
       spcdb_files : list : Paths to Ref & Dev species_database.yml files

       Keyword Args
       dst         : str  : Folder where OH Metrics output will be written
       overwrite   : bool : Overwrite previously-generated files? (T/F)
       """

.. _bmk-funcs-table-ops:

make_benchmark_operations_budget
--------------------------------

**Located in module:** :file:`gcpy.benchmark.module.benchmark_funcs`

Creates a table with the change in species mass after each GEOS-Chem
operation, using diagnostic output from GEOS-Chem benchmark
simulations.

.. code-block:: python

   def make_benchmark_operations_budget(
           refstr,
           reffiles,
           devstr,
           devfiles,
           ref_interval,
           dev_interval,
           spcdb_files,
           benchmark_type=None,
           label=None,
           col_sections=["Full", "Trop", "PBL", "FixedLevs", "Strat"],
           operations=["Chemistry", "Convection", "EmisDryDep",
                       "Mixing", "Transport", "WetDep"],
           compute_accum=True,
           compute_restart=False,
           require_overlap=False,
           dst='.',
           species=None,
           overwrite=True,
           verbose=False,
   ):
       """
       Prints the "operations budget" (i.e. change in mass after
       each operation) from a GEOS-Chem benchmark simulation.

       Args:
           refstr: str
               Labels denoting the "Ref" versions
           reffiles: list of str
               Lists of files to read from the "Ref" version.
           devstr: str
               Labels denoting the "Dev" versions
           devfiles: list of str
               Lists of files to read from "Dev" version.
           interval: float
               Number of seconds in the diagnostic interval.
           spcdb_files : list
               Paths to species_database.yml files in Ref & Dev rundirs

       Keyword Args (optional):
           benchmark_type: str
               A string denoting the type of benchmark output to plot, options are
               FullChemBenchmark, TransportTracersBenchmark, or CH4Benchmark.
               Default value: None
           label: str
               Contains the date or date range for each dataframe title.
               Default value: None
           col_sections: list of str
               List of column sections to calculate global budgets for. May
               include Strat eventhough not calculated in GEOS-Chem, but Full
               and Trop must also be present to calculate Strat.
               Default value: ["Full", "Trop", "PBL", "FixedLevs", "Strat"]
           operations: list of str
               List of operations to calculate global budgets for. Accumulation
               should not be included. It will automatically be calculated if
               all GEOS-Chem budget operations are passed and optional arg
               compute_accum is True.
               Default value: ["Chemistry","Convection","EmisDryDep",
                               "Mixing","Transport","WetDep"]
           compute_accum: bool
               Optionally turn on/off accumulation calculation. If True, will
               only compute accumulation if all six GEOS-Chem operations budgets
               are computed. Otherwise a message will be printed warning that
               accumulation will not be calculated.
               Default value: True
           compute_accum: bool
               Optionally turn on/off accumulation calculation. If True, will
               only compute accumulation if all six GEOS-Chem operations budgets
               are computed. Otherwise a message will be printed warning that
               accumulation will not be calculated.
               Default value: True
           compute_restart: bool
               Optionally turn on/off calculation of mass change based on restart
               file. Only functional for "Full" column section.
               Default value: False
           require_overlap: bool
               Whether to calculate budgets for only species that are present in
               both Ref or Dev.
               Default value: False
           dst: str
               Directory where plots & tables will be created.
               Default value: '.' (directory in which function is called)
           species: list of str
               List of species for which budgets will be created.
               Default value: None (all species)
           overwrite: bool
               Denotes whether to overwrite existing budget file.
               Default value: True
           verbose: bool
               Set this switch to True if you wish to print out extra
               informational messages.
               Default value: False
       """

.. _bmk-funcs-table-ste:

make_benchmark_ste_table
------------------------

**Located in module:** :file:`gcpy.benchmark.modules.ste_flux`

Generates a table with the stratosphere-troposphere flux of ozone from
GEOS-Chem benchmark simulation output.

.. note::

   This table is only available for GEOS-Chem Classic benchmarks.

.. code-block:: python

   def make_benchmark_ste_table(
           devstr,
           files,
           year,
           dst='./1yr_benchmark',
           bmk_type="FullChemBenchmark",
           species=None,
           overwrite=True,
           month=None
   ):
       """
       Driver program.  Computes and prints strat-trop exchange for
       the selected species and benchmark year.

       Args:
           devstr: str
               Label denoting the "Dev" version.
           files: str
               List of files containing vertical fluxes.
           year: str
               Year of the benchmark simulation.

       Keyword Args (optional):
           dst: str
               Directory where plots & tables will be created.
           bmk_type: str
               FullChemBenchmark or TransportTracersBenchmark.
           species: list of str
               Species for which STE fluxes are desired.
           overwrite: bool
               Denotes whether to ovewrite existing budget tables.
           month: float
               If passed, specifies the month of a 1-month benchmark.
               Default: None (denotes a 1-year benchmark)
       """

.. _bmk-funcs-table-ttbdg:

transport_tracers_budgets
-------------------------

**Located in module:** :file:`gcpy.benchmark.modules.budget_tt`

Generates a budget table for Rn222, Pb210, and Be7 species from 1-year
TransportTracers benchmark output.

.. code-block:: python

   def transport_tracers_budgets(
           devstr,
           devdir,
           devrstdir,
           year,
           spcdb_file,
           dst='./1yr_benchmark',
           is_gchp=False,
           gchp_res="c00",
           gchp_is_pre_14_0=False,
           overwrite=True,
   ):
       """
       Main program to compute TransportTracers budgets from
       GEOS-Chem Classic or GCHP benchmark simulations.

       Args
       devstr           : str  : Label for the Dev version
       devdir           : str  : Path to the Dev data directory
       devrstdir        : str  : Path to the Dev restart file directory
       year             : int  : Year of the benchmark simulation
       spcdb_file       : str  : Path to the Dev species_database.yml file

       Keyword Args
       dst              : str  : Directory where tables will be written
       overwrite        : bool : Should existing tables should be overwritten?
       is_gchp          : bool : Is Dev from a GCHP benchmark simulation?
       gchp_res         : str  : GCHP resolution string (e.g. "c24")
       gchp_is_pre_14_0 : bool : Is Dev from a GCHP version prior to 14.0.0?
       """
