# GCPy example scripts

This directory contains several subdirectories with example scripts that demonstrate the capabilities of GCPy.

## bpch_to_nc

NOTE: The binary punch ("bpch") data format has been retired from GEOS-Chem.  We keep these scripts here for those who work with the GEOS-Chem Adjoint code, which still uses bpch format.

`bpch2nc.py`

- Script to convert GEOS-Chem binary punch (aka "bpch") data to netCDF.

`bpch_tagco_prodloss_to_nc.py`

- Converts the prod/loss data files in bpch format for the tagged CO simulation to netCDF format.


## diagnostics

`compare_diags.py`

- Script to compare the contents of files from two different model versions: A reference version (aka "Ref") and a development version (aka "Dev").

`compare_diags.yml`

- Configuration file for use with `compare_diags.py`

## dry_run

`download_data.py`

- Downloads data from a GEOS-Chem Classic "dry-run" simulation.

`download_data.yml`

- Configuration file for `download_data.py`.


## hemco

`format_hemco_demo.py`
# GCPy example scripts

This directory contains several subdirectories with example scripts that demonstrate the capabilities of GCPy.

## bpch_to_nc

NOTE: The binary punch ("bpch") data format has been retired from GEOS-Chem.  We keep these scripts here for those who work with the GEOS-Chem Adjoint code, which still uses bpch format.

`bpch2nc.py`

- Script to convert GEOS-Chem binary punch (aka "bpch") data to netCDF.

`bpch_tagco_prodloss_to_nc.py`

- Converts the prod/loss data files in bpch format for the tagged CO simulation to netCDF format.


## diagnostics

`compare_diags.py`

- Script to compare the contents of files from two different model versions: A reference version (aka "Ref") and a development version (aka "Dev").

`compare_diags.yml`

- Configuration file for use with `compare_diags.py`

## dry_run

`download_data.py`

- Downloads data from a GEOS-Chem Classic "dry-run" simulation.

`download_data.yml`

- Configuration file for `download_data.py`.


## hemco

`format_hemco_demo.py`

- Demonstrates how to fix a non-COARDS-compliant file (needed for HEMCO) using the `gcpy/community/format_hemco_data.py` module from Hannah Nesser (@hannahnesser).

`make_mask_file.py`

- Creates mask files for HEMCO emissions for a given country.


## plotting

`create_test_plot.py`

- Script to create a test pattern plot.  Useful for testing if the Python environment has been installed properly.

`plot_comparisons.py`

- Plots data from two different models side-by-side for comparison purposes, in a "six-panel" plot layout.

`plot_single_panel.py`

- Creates several different types of single-panel plots.

`plot_timeseries.py`

- Reads and plots timeseries data.


## working_with_files

`add_blank_var_to_restart_file.py`

- Adds a "dummy" DataArray containing all zeroes to a GEOS-Chem restart file.

`concatenate_files.py`

- Combines several netCDF data files into a single file using xarray.

`insert_field_into_restart_file.py`

- Adds a DataArray field into a GEOS-Chem restart file.

`regrid_restart_ll_to_cs.py`

- Regrids data from the lat-lon grid to a cubed-sphere grid.