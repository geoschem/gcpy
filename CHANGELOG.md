# GCPy Changelog

All notable changes to GCPy will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [Unreleased]

### Added
- Added strat/trop exchange fluxes to 1-year benchmark output (gcpy/ste_flux.py)
- Added operations budgets to 1-year benchmark output (gcpy/budget_ops.py)
- Added seasonal mass table output for 1-year benchmak (FullChem).

### Changed
- Significant difference files are now written out to the Plots/Sig_Diffs folder for the 1-year benchmarks.
- Updated file names for Pb/Be budget tables in gcpy/budgets_tt.py.
- Created separate driver routines for 1-year FullChem and TransportTracers benchmarks
- Useless warnings when creating benchmark output should now be suppressed

### Deprecated

### Fixed

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
