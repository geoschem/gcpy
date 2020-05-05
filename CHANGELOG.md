# GCPy Changelog

All notable changes to GCPy will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [Unreleased]

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

### Deprecated

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
