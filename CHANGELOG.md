# GCPy Changelog

All notable changes to GCPy will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## Unreleased
### Added
### Changed
## [1.1.0] - 2021-09-22
### Added

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


## [Unreleased]

### Added

### Changed

### Deprecated

### Fixed

### Removed
