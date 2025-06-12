# GCPy Changelog

All notable changes to GCPy will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.6.2] - 2025-06-12
### Added
- Added `create_benchmark_sanity_check_table` routine to `gcpy/benchmark/benchmark_funcs.py` to test if variables are all zero or NaN
- Added a chapter on using code profiling tools in the ReadTheDocs documentatio
n
- Added code profiling scripts (in `gcpy/gcpy/profile`) to read and display output from gprofng and Intel VTune profilers
- Added `check_gchp_emission_diags.py` example script and documentation
- Added new benchmark functions `make_benchmark_collection_2d_var_plots` and `make_benchmark_collection_3d_var_plots` which can be used with any GEOS-Chem output collection
- Added 1-month benchmark comparison plot options for `Budget`, `UVFlux`, and `StateMet` collections (2D and 3D vars separately) which are off by default
- Added `export MPLBACKEND=agg` to `gcpy/benchmark/modules/benchmark_slurm.sh` to request a non-interactive MatPlotLib backend
- Added `method` keyword argument to `make_regridder_*` routines in `regrid.py`, with default value `conservative`

### Changed
- Updated `gcpy_environment_py313.yml` to use `esmf==8.8.1` and `esmpy==8.8.1` to fix package inconsistency issues

### Fixed
- Fixed logic error in `compare_varnames` that caused 2D data variables to be flagged as 3D (esp. for GCHP vs GCC comparisons)
- Replaced incorrect collection name `AOD`  with `Aerosols` in the GCC vs. GCC 1-month AOD plots
- Fix the logic of ravel_1D that creates orbit file for 1D GCHP diagnostic

## [1.6.2] - 2025-06-12
### Changed
- Bumped `requests` from 2.32.3 to 2.32.4 in environment files to fix a security issue flagged by Dependabot

### Fixed
- Fix the logic of ravel_1D that creates orbit file for 1D GCHP diagnostic

## [1.6.2] - 2025-06-12
### Added
- Added `create_benchmark_sanity_check_table` routine to `gcpy/benchmark/benchmark_funcs.py` to test if variables are all zero or NaN
- Added a chapter on using code profiling tools in the ReadTheDocs documentatio
n
- Added code profiling scripts (in `gcpy/gcpy/profile`) to read and display output from gprofng and Intel VTune profilers
- Added `check_gchp_emission_diags.py` example script and documentation
- Added new benchmark functions `make_benchmark_collection_2d_var_plots` and `make_benchmark_collection_3d_var_plots` which can be used with any GEOS-Chem output collection
- Added 1-month benchmark comparison plot options for `Budget`, `UVFlux`, and `StateMet` collections (2D and 3D vars separately) which are off by default
- Added `export MPLBACKEND=agg` to `gcpy/benchmark/modules/benchmark_slurm.sh` to request a non-interactive MatPlotLib backend
- Added `method` keyword argument to `make_regridder_*` routines in `regrid.py`, with default value `conservative`

### Changed
- Updated `gcpy_environment_py313.yml` to use `esmf==8.8.1` and `esmpy==8.8.1` to fix package inconsistency issues

### Fixed
- Fixed logic error in `compare_varnames` that caused 2D data variables to be flagged as 3D (esp. for GCHP vs GCC comparisons)
- Replaced incorrect collection name `AOD`  with `Aerosols` in the GCC vs. GCC 1-month AOD plots

## [1.6.1] - 2025-03-24
### Added
- Added GCPy environment file `docs/environment_files/gcpy_environment_py313.yml`, which is based on Python 3.13
- Added GitHub Actions config file`.github/workflows/build-gcpy-environment-py313.yml`
- Added `Install-Mamba-Conda.rst` ReadTheDocs supplmental guide for Mamba/Conda installation
- Added Anaconda, PyPI, Platforms, Downloads badges to `README.md` and `docs/source/index.rst`

### Changed
- Updated paths to regridding weights in `benchmark/config/*.yml` and `examples/diagnostics/compare_diags.yml` files
- Renamed `docs/environment_files/gcpy_environment.yml` to `gcpy_environment_py312.yml`
- Updated symbolic link `environment.yml` to point to `docs/source/gcpy_environment_py312.yml`
- Renamed GitHub Actions config file `.github/workflows/build-gcpy-environment.yml` to `build-gcpy-environment-py312.yml`
- Updated package version information in ReadTheDocs documentation
- Bumped dask from version 2024.5.2 to 2025.3.0 to fix a security issue (raised by @dependabot)
- Bumped jinja2 from version 3.1.5 to 3.1.6 to fix a security issue (raised by @dependabot)
- Edited installation instructions via Mamba/Conda to more accurately reflect the needed commands
- Updated the ReadTheDocs "Releasing new versions" guide with info on conda-forge

### Removed
- Removed PyPi configuration file  `docs/environment_files/gcpy_requirements.txt` and symbolic link `./requirements.txt`
- Removed Python 3.8 from the `python-version` list in `.github/workflows/*.yml` files
- Remove Mamba/Conda installation instructions from `Getting-Started-from-GCPy.rst`
- Removed "Latest pre-release" bade from `README.md` and `docs/source/index.rst`

## [1.6.0] - 2025-01-29
### Added
- Added example script `gcpy/examples/hemco/make_hemco_sa_spec.py` (creates the HEMCO standalone configuration file `HEMCO_sa_Spec.rc`)
- Added module `benchmark_gcclassic_stats.py` for scraping statistics from GEOS-Chem Classic cloud benchmarks
- Added dry deposition velocity comparison plots in 1-month cloud benchmarks
- Added `gcpy/benchmark/modules/benchmark_species_changes.py` to compute the table of species changes between versions
- Added `gcpy/kpp/` folder containing scripts to plot output from the KPP-Standalone box model
- Added ReadTheDocs documentation for plotting output from the KPP-Standalone box model

### Changed
- Changed format of `% diff` column from `12.3e` to `12.3f` in benchmark timing tables 
- Updated `gcpy/benchmark/modules/emission_species.yml` file with emission species for GEOS-Chem 14.5.0
- Updated `gcpy/benchmark/modules/benchmark_categories.yml` with the latest categories for GEOS-Chem 14.5.0
- Updated `gcpy/benchmark/modules/lumped_species.yml` with speciations for GEOS-Chem 14.5.0
- Add `DryDep` to list of collections included in benchmark summary table
- Updated `checkout` GitHub action to v4
- Updated `CodeQL` GitHub action to v3
- Updated `publish-python` GitHub action to v5
- In environment files `gcpy_environment.yml` and `gcpy_requirements.txt`:
  - Update `python` to 3.12.0
  - Update `xesmf` to 0.8.5
  - Update `esmf` and `esmpy` to 8.6.1
- In environment files `read_the_docs_environment.yml` and `read_the_docs_requirements.txt`
  - Update `jinja` to 3.1.5 (fixes a security issue)
- Update `gcpy/setup.py` with the new Python package version numbers
- Updated code in `gcpy/benchmark/modules/` to replace whitespace in Ref and Dev labels with underscores

### Fixed
- Fixed formatting error in `.github/workflows/stale.yml` that caused the Mark Stale Issues action not to run
- Now flag differences greater than +/- 10% in benchmark timing table outputs
- Fixed error in computation of dynamic ratio plot min & max values in `plot/six_plot.py`
- Fixed erroneous species classification in `gcpy/benchmark/modules/benchmark_categories.yml`
- Fixed type errors in `calc_rectilinear_lon_edge` and `calc_rectangular_lat_edge` by casting the length of the output   array from `float` to `int`

### Removed
- Removed `gcpy/benchmark/modules/species_database.yml` file and corresponding code pointing to this

## [1.5.0] - 2024-05-29
### Added
- Script `gcpy/benchmark/modules/benchmark_utils.py`, with common benchmark utility functions
- Script `gcpy/benchmark/modules/benchmark_drydep.py`, with code to create drydep velocity plots
- YAML tag `plot_drydep` in `gcpy/benchmark/config/*.yml` files
- Badge for `build-gcpy-environment` GitHub Action in `README.md`
- Badges in `docs/source/index.rst`
- GitHub action to push GCPy releases to PyPi
- Script `./release/changeVersionNumbers.sh`, used to update version numbers in various files before release
- Mamba/Conda enviroment file `docs/environment_files/read_the_docs_environment.yml`, for building ReadTheDocs documentation
- Environment files `docs/environment_files/gcpy_requirements.txt` and `docs/environment_files/read_the_docs_requirements.txt`
- New benchmark script `gcpy/benchmark/modules/benchmark_models_vs_sondes.py`
- Added fixed level budget diagnostic to budget operations table
- Function `rename_speciesconc_to_speciesconcvv` in `gcpy/benchmark/modules/benchmark_utils.py`
- Function `copy_file_to_dir` in `gcpy/util.py`.  This is a wrapper for `shutil.copyfile`.
- GitHub Action config file `.github/workflows/stale.yml`, which replaces StaleBot
- Added `gcpy/community/format_hemco_data.py` from @hannahnesser
- Added `gcpy/examples/hemco/format_hemco_demo.py` from @hannahnesser
- Added HCl to `gcpy/benchmark/modules/emission_species.yml` and GTChlorine to `gcpy/benchmark/modules/emission_inventories.yml` for GEOS-Chem 14.4.0
- Example script `gcpy/examples/working_with_files/make_mask_file.py`
- Convenience function `replace_whitespace` in `gcpy/util.py`
- Benchmark script `gcpy/benchmark/modules/benchmark_scrape_gcclassic_timers.py`
- Benchmark script `gcpy/benchmark/modules/benchmark_scrape_gchp_timers.py`
- 1-year benchmark scripts now produce GCC vs GCC and GCHP vs GCHP timing tables
- Functions `gcc_vs_gcc_dirs`, `gchp_vs_gcc_dirs`, `gchp_vs_gchp_dirs`, and `get_log_filepaths` in `gcpy/benchmark/modules/benchmark_utils.py`
- Script `gcpy/benchmark/modules/benchmark_mass_cons_table.py`, with code to create mass conservation tables
- Expanded statistics output in benchmark mass conservation tables
- Function `get_datetimes_from_filenames` in `gcpy/benchmark/modules/benchmark_utils.py`
- Function `replace_whitespace` in `gcpy/util.py`

### Changed
- Bump pip from 23.2.1 to 23.3 (dependabot suggested this)
- Bump pypdf from 3.16.1 to 3.17.0 (dependabot suggested this)
- YAML tag `operations_budget` is now `ops_budget_table` in `gcpy/benchmark/config/1yr_tt_benchmark.yml`
- Now require `matplotlib=3.8.0` in `docs/environment_files/environment.yml` (with other pegged versions)
- Now run the `stale` GitHub action at 00:00 UTC on the 1st of each month
- Renamed `docs/environment_files/environment.yml` to `gcpy_environment.yml`
- `environment.yml` links to `docs/environment_files/gcpy_environment.yml`
- `requirements.txt` links to `docs/environment_files/requirements.txt`
- Python packages for RTD documenation builds from `docs/environment_files/environment.yml`
- Script `benchmark_model_vs_obs.py` now uses grid inquiry functions from `grid.py` to return data nearest to a (lat,lon) location
- Moved routine `get_geoschem_level_metadata` to `gcpy/benchmark/modules/benchmark_utils.py`
- Refactored `get_vert_grid.py` (in `gcpy/grid.py`) to accept the `p_sfc` argument; Also never-nested the if-block logic.
- `benchmark_slurm.sh` script now saves output to a log file with the same base name as the YAML config file
- `benchmark_models_vs_obs.py` now reads the observational data paths and metadata from  `1yr_fullchem_benchmark.yml`
- Moved benchmark code from `gcpy/` to `gcpy/benchmark/modules/`
- Benchmark functions now call `rename_speciesconc_to_speciesconcvv`
- Create radionuclide, STE flux, and mass conservation tables for Ref and Dev versions in TransportTracers benchmarks
- Use new function `copy_file_to_dir` to copy the benchmark script and configuration file to the benchmark results folders
- Updated GitHub stalebot config file `stale.yml` with new issue/PR labels that should not go stale
- Updated benchmark driver scripts and config files to print GCClassic & GCHP timing information
- Updated documentation for GCPy 1.5.0 release
- Converted Github issue templates to issue forms using YAML definition files

### Fixed
- CS inquiry functions in `gcpy/cstools.py` now work properly for `xr.Dataset` and `xr.DataArray` objects
- Prevent an import error by using `seaborn-v0_8-darkgrid` in`gcpy/benchmark/modules/benchmark_models_vs_obs.py`
- `gcpy/file_regrid.py` now creates GCHP files with `DELP_DRY` instead of `DELPDRY`
- Fixed bugs in GCHP vs GCHP sections of 1-yr benchmark run scripts to allow comparison of GCHP runs with different grid resolutions
- Fixed silent bug in transport tracer benchmark GCC vs GCHP mass
  tables preventi ng them from being generated
- Import error in `gcpy/examples/diagnostics/compare_diags.py`
- Added missing `n_cores` to `gcpy/examples/diagnostics/compare_diags.yml`
- Added missing `plot_drydep` option to `gcpy/gcpy/benchmark/config/1yr_ch4_benchmark.yml`
- Added `docs/requirements.txt` symbolic link to `docs/environment_files/read_the_docs_requirements.txt` for RTD builds
- `gcpy/file_regrid.py` now tests if `lon_bnds`, `lat_bnds` are in the dataset before trying to drop them
- Ensured that the header line in the aerosol burden table is written to file

### Removed
- Example script `gcpy/examples/plotting/mda8_o3_timeseries.py`
- Removed `Pylint` GitHub action
- Environment file `docs/environment_files/environment.yml`
- Environment file `docs/environment_files/requirements.txt`
- Removed `awscli` from the GCPy environment; version 2 is no longer available on conda-forge or PyPi
- GitHub config files `.github/stale.yml` and `.github/no-response.yml`
- Routine `make_benchmark_mass_conservation_table` in `benchmark_funcs.py`; this is now obsolete

## [1.4.2] - 2024-01-26
### Added
- Example script `create_test_plot.py`, which can be used to check that GCPy has been installed properly
- GitHub action `build-gcpy-environment` which tests installation of the mamba environment specified in in `docs/environment_files/environment.yml`
- YAML file`docs/environment_files/testing.yml` for building an environment without pegged package versions (for testing)
- GitHub action `build-test-environment` to test the environment specified in `testing.yml`

### Changed
- `build-gcpy-environment` GitHub action now runs with several Python versions

### Fixed
- Prevent overwriting of the `results` variable when parallel plotting is deactivated (`n_cores: 1`)

## [1.4.1] - 2023-12-08
### Fixed
- Now use the proper default value for the `--weightsdir` argument to `gcpy/file_regrid.py`

## [1.4.0] - 2023-11-20
### Added
- Added C2H2 and C2H4 to `emission_species.yml`
- Updated `species_database.yml` for consistency with GEOS-Chem 14.2.0
- Added `.github/ISSUE_TEMPLATE/config.yml` file w/ Github issue options
- Added `CONTRIBUTING.md` and `SUPPORT.md`, replacing `docs/source/Contributing.rst` and `docs/source/Report_Request.rst`
- Added option to pass the benchmark type to plotting routines
- Updated `AUTHORS.txt` as of Apr 2023 (concurrent w/ GEOS-Chem 14.2.0)
- Added ReadTheDocs badge in `README.md`
- Added `.readthedocs.yaml` to configure ReadTheDocs builds
- Added cloud benchmarking YAML configuration files to `benchmark/cloud` folder
- Added `README.md` files in `gcpy/benchmark` directory structure
- Added `benchmark/modules/benchmark_models_vs_obs.py` script
- Added `benchmark/modules/GC_72_vertical_levels.csv` file
- Added `multi_index_lat` keyword to `reshape_MAPL_CS` function in `gcpy/util.py`
- Added FURA to `emission_species.yml` and `benchmark_categories.yml`
- Added new routine `format_number_for_table` in `gcpy/util.py`
- Added module `gcpy/cstools.py` with utility functions for cubed-sphere grids
- Added new routine `verify_variable_type` function in `gcpy/util.py`
- Added new routine `format_number_for_table` in `util.py`
- Added BrSALA and BrSALC to `emission_species.yml`
- Added `options:n_cores` to all benchmark YAML config files
- Added `__init__.py` files in subfolders of `gcpy/gcpy`
- `gcpy/benchmark/modules/*.py` scripts are now chmod 644
- Added `ENCODING = "UTF-8"` to `gcpy/constants.py`
- Added statement `from dask.array import Array as DaskArray` in `gcpy plot.py`
- Added SLURM run script `gcpy/benchmark/benchmark_slurm.sh`
- Added `gcpy/plot/gcpy_plot_style` style sheet for title and label default settings
- Added `gcpy/gcpy_plot_style` style sheet for title and label default settings
- Added new cubed-sphere grid inquiry functions to `gcpy/cstools.py`
- Added functions `get_ilev_coord` and `get_lev_coord` to `gcpy/grid.py`
- Add `tk` package to `docs/environment_files/environment.yml`

### Changed
- Simplified the Github issues templates into two options: `new-feature-or-discussion.md` and `question-issue.md`
- The GitHub PR template is now named `./github/PULL_REQUEST_TEMPLATE.md`
- Updated badge links in `README.md`
- Construct ops budget table filename without using the `label` argument
- Updated species_database.yml for consistency with GEOS-Chem 14.2.0
- Renamed TransportTracers species in `benchmark_categories.yml`, `run_1yr_tt_benchmark.py`, and in documentation
- YAML files in `benchmark/` have been moved to `benchmark/config`
- Models vs. O3 obs plots are now arranged by site latitude from north to south
- Routine `print_totals` now prints small and/or large numbers in scientific notation
- Truncate names in benchmark & emissions tables to improve readability
- Add TransportTracers species names to `gcpy/emissions_*.yml` files
- Now pass `n_job=config["options"]["n_cores"]` to benchmark plotting routines
- Script `benchmark.py` to `benchmark_funcs.py` to remove a name collision
- Folder `gcpy/benchmark` is now `gcpy/gcpy/benchmark`
- Folder `benchmark/modules` is now `gcpy/gcpy/benchmark/modules`
- Folder `gcpy/examples` is now `gcpy/gcpy/examples`
- Pass `sys.argv` to the `main()` routine of `run_benchmark.py`,` compare_diags.py`
- Updated `docs/environment_files/environment.yml` for MambaForge (also added `gridspec`)
- Now use `pypdf` instead of `PyPDF2` in `plot.py` and `util.py`
- Added coding suggestions made by `pylint` where possible
- Abstracted and never-nested code from `six_plot` into functions (in `plot.py`)
- Added `main()` routine to `gcpy/file_regrid.py`; Also added updates suggested by Pylint
- Fixed broken regridding code in `gcpy/file_regrid.py`; also refactored for clarity
- Rewrote `Regridding.rst` page; Confirmed that regridding examples work properly
- Now allow `plot_val` to be of type `dask.array.Array` in `plot.py` routines `six_plot` and `single_panel`
- Now add `if` statements to turn of `Parallel()` commands when `n_jobs==1`.
- Do not hardwire fontsize in `gcpy/plot.py`; get defaults from `gcpy_plot_style`
- `gcpy/plot.py` has been split up into smaller modules in the `gcpy/plot` folder
- Updated and cleaned up code in `gcpy/regrid.py`
- Example scripts`plot_single_level` and `plot_comparisons` can now accept command-line arguments
- Example scripts `plot_single_level.py`, `plot_comparisons.py`, `compare_diags.py` now handle GCHP restart files properly
- Now specify the X11 backend with by setting the `MPLBACKEND` environment variable

### Fixed
- Generalized test for GCHP or GCClassic restart file in `regrid_restart_file.py`
- Fixed bug in transport tracer benchmark mass conservation table file write
- Routine `create_display_name` now splits on only the first `_` in species & diag names
- Prevent plot panels from overlapping in six-panel plots
- Prevent colorbar tick labels from overlapping in dynamic-range ratio plots
- Updated `seaborn` plot style names to conform to the latest matplotlib
- Set `lev:positive` and/or `ilev:positive` properly in `regrid_restart_file.py` and `file_regrid.py`
- Prevent overwriting of `lev` coord in `file_regrid.py` at netCDF write time
- Fixed bug in option to allow different units when making comparison plots

### Removed
- Removed `gchp_is_pre_13_1` arguments & code from benchmarking routines
- Removed `is_pre_13_1` tags from `*_benchmark.yml` config files
- Removed `benchmark_emission_totals.ipynb`, this is obsolete
- Replaced `gcpy/benchmark/README` with `README.md`
- Removed `gcpy_test_dir` option from `examples/diagnostics/compare_diags.*`
- Removed `docs/environment_files/gchp_regridding.yml` environment file
- Removed `gcpy/gcpy/benchmark/plot_driver.sh`
- Made benchmark configuration files consistent

## [1.3.3] -- 2023-03-09
### Added
- Updated installation documentation, we now recommend users to create
  a conda environment using the `environment.yml` file
- Benchmark summary table output (intended for 1hr & 1mo benchmarks)
- Species/emissions/inventories that differ between Dev & Ref versions are now printed at the top of the benchmark emissions, inventory, and global mass tables.  if there are too many species with diffs, an alternate message is printed.
- New functions in `benchmark.py` and `util.py` to facilitate printing of the species/emissions/inventories that differ between Dev & Ref versions.
- Added new RTD documentation for installing Conda 4.12.0 with Miniconda
- Added GCHP regridding environnment file `docs/environment_files/gchp_regridding.yml`
- Added new benchmark type CH4Benchmark

### Changed
- Applied cleanup susggestions from pylint to `benchmark.py`, `util.py`, `plot.py`, `oh_metrics.py`, `ste_flux.py`
- Replaced format with f-strings in `benchmark.py`, `util.py`, `plot.py`, `oh_metrics.py`, `ste_flux.py`
- Abstract some common in `benchmark.py` into functions
- Replaced direct calls to `yaml.load` with `util.read_config.file`
- Restore tag information to benchmark `refstr` and `devstr` labels
- Add a newline to diff-of-diffs refstr and devstr if the string is too long.
- Updated GCHP regridding documentation
- Restored `ipython` and `jupyter ` to environment file `environment.yml`

## [1.3.2] -- 2022-10-25

### Fixes
- Fixed malformed version declaration for cartopy (use `==`
  instead of `=`) in setup.py.  This was preventing upload to
  conda-forge.
- Vertically flip GCHP emissions when computing transport tracers budget

## [1.3.1] -- 2022-10-25

### Changed
- Bug fix: Remove extraneous character from setup.py

## [1.3.0] -- 2022-10-25

### Added
- New features in benchmarking scripts (@lizziel, @yantosca)
  - Force garbage collection at end benchmarking functions (@yantosca)
  - Extra print statements (@lizziel)
  - Diff-of-diffs plots for 1-year benchmarks (@lizziel)
  - sparselt is now a GCPy requirement (@lizziel)
- Removed obsolete environment.yml files (@yantosca)
- Added requirements.yml to docs folder for Sphinx/RTD documentation (@yantosca)
- New regridding script `regrid_restart_file.py` (@liambindle)

### Changed
- Fixed several issues in benchmarking scripts (@laestrada, @lizziel, @yantosca)
  - Fixed bug in `budget_ox.py`; The drydep loss of Ox for GCHP was 12x too high
  - Add OMP_NUM_THREADS and OMP_STACKSIZE in `plot_driver.sh` (@yantosca)
  - Increase requested memory to 50MB in `plot_driver.sh` (@yantosca)
  - Benchmark scripts print a message upon completion (@yantosca)
  - Linted several benchmarking routines with Pylint (@yantosca)
  - Rewrote algorithm of add_lumped_species_to_dataset for speed (@yantosca)
  - Can now specify the path to species_database.yml for 1yr benchmarks (@yantosca)
  - 1-yr benchmarks now save output in subdirs of the same path (@lizziel)
  - Avoid hardwiring restart file paths in benchmark scripts (@yantosca)
  - Now use outputs_subdir tag from YAML file for paths to diagnostic files (@yantosca)
  - Now use restarts_subdir tag from YAML file for paths to restart files (@yantosca)
  - GCPy now uses proper year for dev in 1-yr benchmarks (@laestrada)
  - Fixed date string issue in benchmarking scripts (@lizziel)
  - Updates for new GCHP restart file format (@lizziel)
- Updated environment.yml with package versions that work together (@yantosca)
- Updated the AUTHORS.txt and LICENSE.txt files (@yantosca)

## [1.2.0] - 2021-09-22
### Added
- Added Parameter for single_panel to support return of all 6 cubedsphere plots
- Added flexible time period for benchmark plotting scripts
### Changed
- Modified single_panel to vmin/vmax parameters with newer versions of matplotlib (>3.5.0)
- Modified run_benchmark script to select correct species database depending on benchmark type
- Modified filename for Ox budget
- Modified readthedocs build to use mamba instead of conda to fix build failures
- Modified benchmark plotting scripts to use a single run_benchmark.py script
- Modified benchmark categories and species database yaml files
- Fixed bug in mass conservation table percent difference

## [1.1.0] - 2021-09-22

- Added date_time.py module to help manage datetime utility functions
- Added GLYC, HAC, and pFe to benchmark categories
- Added gcpy/budget_ox.py to compute Ox budgets from 1-yr benchmarks
- Added capability to use GCHP 13.1.0+ or legacy file names in benchmark scripts
- Added new methods dataset_reader and get_dataset_mean to util.py

### Changed
- Modified benchmarking scripts to use yaml config files.
- Modified dry-run scripts to use yaml config files.
- Updated benchmark/run_1yr_fullchem_benchmark.py to call the budget_ox.py for GCC vs GCC benchmark generation.
  - NOTE: we are waiting to make sure that the GCHP benchmarks output wetdep fields before activating this feature for GCHP.
- Modified plotting methods in benchmark.py to compute the mean of datasets over the time dimension, if the "time_mean" keyword is passed.
  - This feature is used to generate annual mean plots from 1-yr benchmark output.
- Modified run_1yr_tt_benchmark.py and run_1yr_fullchem_benchmark.py to generate both annual mean and seasonal plots
- Fixed formatting and import order issues in benchmark.py, util.py, budget_ox.py, and the run_*benchmark.py scripts as identified by pylint.
- Modified budget_ox.py to use Ox instead of O3 for computing budget terms

## [1.0.3] - 2021-03-26

### Fixed
- Automatic benchmark script copying no longer overwrites existing files
- Color scales for non-global plots are no longer calculated from full global data
- Regional datasets can now be plotted with cubed-sphere datasets in plot.compare_single_level

## [1.0.2] - 2021-03-18

### Added
- Added GCPy version number and automatic script copying to benchmark scripts
- Added line clarifying lack of Windows support in ReadTheDocs

### Fixed
- Fixed benchmark month seconds calculation for GCHP in 1-month benchmark script
- Fixed label typo in benchmark script GCHP vs. GCC emission plots
- Fixed grid creation for non-global grids in plot.single_panel
- Fixed issue in get_grid_extents when maxlon was in Western Hemisphere

## [1.0.1] - 2021-02-09

### Added
- Added MSA to Sulfur benchmark category
- Added weightsdir parameter to single_panel()
- Added temporary file creation to file_regrid() to decrease memory consumption during cubed-sphere regridding

### Changed
- Removed carbon-based units from benchmark emissions tables
- Environment files now request xESMF through conda-forge rather than pip

### Fixed
- Fixed Cubed-Sphere to Lat/Lon regridding for 1-level files.
- Fixed single panel zonal mean axis selection

## [1.0.0] - 2021-01-05

### Added
- Added complete documentation to a new ReadTheDocs site
- Added conda-forge installation support
- Added file regridder for regridding NetCDF restart and output between GEOS-Chem's horizontal grid types
- Plotting now supports automatic regridding between lat/lon, cubed-sphere, and stretched-grid formats
- Added additional 1-year benchmark plotting capabilities for GCHP
- Added oh_metrics.py, which generates output using the new Metrics collection in GEOS-Chem 13.0.0
- Extra keyword arguments not defined in plotting functions are now passed to matplotlib.pyplot
- Added a command line tool for appending grid-box corners to cubed-sphere datasets
- Added support for arbitrary vertical grids in zonal mean plotting
- Added regridding functions for arbitrary vertical grids

### Changed
- Some constants in constants.py have been tweaked to match GEOS-Chem definitions
- docs/environment.yml, setup.py, and requirements.txt now reflect up-to-date GCPy library requirements
- Most docstrings now use the same format
- Various code formatting changes have been made to align with PEP8 guidelines

### Deprecated
- mean_oh_from_logs.py is replaced in functionality by oh_metrics.py for GEOS-Chem versions >=13.0.0

### Fixed
- Installation through pip (from the repositoryand conda now works correctly

### Removed
- Removed several functions and files that are no longer used, including budget_aer.py and create_budget_table()

## [0.3.1] - 2020-08-21

### Added
- Added instructions on setting PYTHONPATH to include GCPy directory when installing manually
- Added cross-dateline regional plotting capability for both lat/lon and cubed-sphere plots
- Added function to get lev dimension index that matches a requested pressure value
- Added basic up-to-date map plotting examples
- Added pip and tabulate dependencies in gcpy environment yaml file
- Added RRTMG netcdf diagnostics names for converting from bpch to nc
- Added unit string conversion for RRTMG binary diagnostics to compare easily with netcdf

### Changed
- Temporary PDFs are now generated in the system's temp directory rather than within the working directory
- environment.yml now includes version numbers to ensure compatability

### Fixed
- Fixed single panel zonal mean plotting for GCHP
- Fixed existing non-deleted examples code
- Fixed imports for out-of-scope variables

### Removed
- Removed several code examples that were out-of-date.

## [0.3.0] - 2020-07-30

### Added

- Add new function to compute budgets and create budget table that incorporates new optional features.
- Require python package tabulate for generating budget tables.
- Added parallel support for mass and budget table creation in 1-year benchmarks.
- Added capability of completely disabling parallel plotting when calling make_benchmark_*_plots functions.
- Added capability of converting concentrations to ug/m3 for benchmark plotting.
- Added new function to make benchmark wet deposition plots, previously done from function to make concentration plots.

### Changed
- Reorganized functions of GCPy into a more logical and streamlined file structure.
- Updated species_database.yml and benchmark_categories.yml for GEOS-Chem 12.9.2.
- Replaced "Plots" with "Results" in benchmark directory structure. This value is customizable in the benchmark scripts.
- Updated example scripts to use reorganized GCPy functions.
- Updated all benchmark run scripts for consistency, readability, reduced lines of code, and compatibility with reorganized and new GCPy functions

### Fixed
- Fixed documentation and rearranged argument order for diff-of-diffs plot strings.
- Fixed accidental regridding to lat/lon in comparison plots where two cubed-sphere datasets share the same resoltuion.
### Removed
- Removed budget_ops.py in deference to new make_benchmark_operations_budget function.

## [0.2.1] - 2020-05-07

### Fixed
- Fixed bugs calculating lumped species with some or all missing constituents.
- Added documentation for newer keyword arguments in benchmark.py

## [0.2.0] - 2020-05-06

### Added
- Added strat/trop exchange fluxes to 1-year benchmark output (gcpy/ste_flux.py)
- Added operations budgets to 1-year benchmark output (gcpy/budget_ops.py)
- Added seasonal mass table output for 1-year FullChemBenchmark.
- Added mean OH from log files for 1-year FullChemBenchmark (GEOS-Chem Classic only).
- Added function gcplot in core.py for creating individual (rather than six-panel) plots of GEOS-Chem data.
- 47-level model output can now be plotted in addition to the standard 72-level output.
- Added Loader=yaml.FullLoader to the yaml.load command to avoid generating excess warnings.
- Add benchmark plotting option to write concentration and emissions plots to one file
- Add gcpy testing mode for all GEOS-Chem benchmark run scripts using test data on gcgrid.
- Add function to flip and rename GCHP restart files to match GCC names and level convention.
- Add capability to generate GCHP vs GCC and GCHP vs GCHP mass tables.
- Add handling in convert_units for datasets without a time dimension.
- Add initial and final mass to GHCP radionuclide budget tables.

### Changed
- Significant difference files are now written out to the Plots/Sig_Diffs folder for the 1-year benchmarks.
- Updated file names for Pb/Be budget tables in gcpy/budgets_tt.py.
- Created separate driver routines for 1-year FullChem and TransportTracers benchmarks
- Useless warnings when creating benchmark output should now be suppressed
- Can now create benchmark plots in a single file instead of by category.
- Can now plot non-global output files from GEOS-Chem Classic.
- Can now limit plot extents using lat/lon parameters for GEOS-Chem Classic and GCHP.
- L2L regridder filenames now include the grid extent when the regridder does not span the whole globe.
- Input GEOS-Chem Classic data can now be any lat/lon resolution rather than only 4x5, 2x2.5, or 1x1.25.
- Replaced fractional difference plots ((Dev-Ref)/Ref) with ratio plots (Dev/Ref).
- Moved diff-of-diffs functionality from standalone code in the benchmark scripts to benchmark.py.
- The bottom row of diff-of-diffs plotting now shows (Dev2/Dev1)-(Ref2/Ref1) values.
- Paths in example scripts now point to /n/holyscratch01/external_repos/GEOS-CHEM instead of /n/holylfs/EXTERNAL_REPOS/GEOS-CHEM.
- Cleaned up run_1mo_benchmark.py driver scripts
- Operations budgets are now printed as Ref, Dev, Dev-Ref, %diff
- Updated examples/compare_diags.py to point to test benchmark data
- Updated benchmark_categories.yml, species_database.yml, lumped_species.yml, and emission_inventories.yml for recent changes in GEOS-Chem 12.8.0
- Update benchmark run scripts to use version strings rather than subtitle strings in tables filenames.


### Fixed
- Latitude ticks again appear in benchmark zonal mean plots.
- Colorbar tick formatting now never uses offset format, which made colorbar ticks difficult to interpret for small value ranges.
- The list of non-plotted emissions species now populates properly.
- Fixed sig diffs file creation for AOD and JValues.
- Missing values in mass tables are now NaN
- Fix area normalization issues in benchmark plotting functions when using GCHP data.

### Removed
- Removed runnable docstring content.

## [0.1.1] - 2020-02-28

### Added

- This CHANGELOG file to track notable changes in GCPy.

### Changed
- Pb210, Be7, and Be10 species are now added to species_database.yml.
- gcpy/budget_aer.py and gcpy/budget_tt.py now get molecular weights from species_database.yml.
- Updated the value of MW_AIR in constants.py to add more precision.
- gcpy/benchmark.py now writes OH metrics output to the Plots/Tables folder.
- Updated CHANGELOG.md for 0.1.1.
- Updated download_data.py to properly obtain 2 x 2.5 and nested data sets.

## [0.1.0] - 2020-02-26

### Summary

This is the first labeled version of GCPy. The primary functionality of GCPy is plotting and tabling diagnostics from GEOS-Chem. Main features include:

- Functions for comparing GEOS-Chem benchmark output for multiple versions of GEOS-Chem, including 2D plots and mass and budget table creation.
- Support for plotting benchmark output for both GEOS-Chem Classic (lat/lon data) and GCHP (cubed-sphere data).

The first official release version of GCPy, v1.0.0, will correspond with the release of GEOS-Chem 13.0.0.
