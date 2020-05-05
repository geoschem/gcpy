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

### Deprecated

### Fixed
- Latitude ticks again appear in benchmark zonal mean plots.
- Colorbar tick formatting now never uses offset format, which made colorbar ticks difficult to interpret for small value ranges.
- The list of non-plotted emissions species now populates properly.
- Fixed sig diffs file creation for AOD and JValues.

### Removed

## [0.1.1] - 2020-02-28

### Added

- This CHANGELOG file to track notable changes in GCPy.

### Changed
- Pb210, Be7, and Be10 species are now added to species_database.yml.
- gcpy/budget_aer.py and gcpy/budget_tt.py now get molecular weights from species_database.yml.
- Updated the value of MW_AIR in constants.py to add more precision.
- gcpy/benchmark.py now writes OH metrics output to the Plots/Tables folder.
- Updated CHANGELOG.md for 0.1.1.

## [0.1.0] - 2020-02-26

### Summary

This is the first labeled version of GCPy. The primary functionality of GCPy is plotting and tabling diagnostics from GEOS-Chem. Main features include:

- Functions for comparing GEOS-Chem benchmark output for multiple versions of GEOS-Chem, including 2D plots and mass and budget table creation.
- Support for plotting benchmark output for both GEOS-Chem Classic (lat/lon data) and GCHP (cubed-sphere data).

The first official release version of GCPy, v1.0.0, will correspond with the release of GEOS-Chem 13.0.0.
