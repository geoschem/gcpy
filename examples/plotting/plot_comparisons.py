#!/usr/bin/env python
"""
Six Panel Comparison Plots
--------------------------------------
This example script demonstrates the comparitive plotting capabilities of GCPy,
including single level plots as well as global zonal mean plots.
These comparison plots are frequently used to evaluate results from different runs / versions
of GEOS-Chem, but can also be used to compare results from different points in one run that
are stored in separate xarray datasets.
The example data described here is in lat/lon format, but the same code works equally
well for cubed-sphere (GCHP) data. 
"""

#xarray allows us to read in any NetCDF file, the format of most GEOS-Chem diagnostics,
#as an xarray Dataset
import xarray as xr
ref_ds = xr.open_dataset('first_run/GEOSChem.Restart.20160801_0000z.nc4')
dev_ds = xr.open_dataset('second_run/GEOSChem.Restart.20160801_0000z.nc4')

import gcpy.plot as gcplot

"""
Single level plots
------------------
"""

#compare_single_level generates sets of six panel plots for data at a specified level in your datasets.
#By default, the level at index 0 (likely the surface) is plotted. Here we will plot data at ~500 hPa,
#which is located at index 21 in the standard 72-level and 47-level GMAO vertical grids.
ilev=21

#You likely want to look at the same variables across both of your datasets. If a variable is in
#one dataset but not the other, the plots will show NaN values for the latter.
#You can pass variable names in a list to these comparison plotting functions (otherwise all variables will plot).
varlist = ['SpeciesRst_O3', 'SpeciesRst_CO2']

#compare_single_level has many arguments which can be optionally specified. The first four arguments are required.
#They specify your first xarray Dataset, the name of your first dataset, your second xarray Dataset, and the name of
#your second dataset. Here we will also pass a specific level and the names of the variables you want to plot.
import matplotlib.pyplot as plt
gcplot.compare_single_level(ref_ds, 'Dataset 1', dev_ds, 'Dataset 2', ilev=ilev, varlist=varlist)
plt.show()

#Using plt.show(), you can view the plots interactively. You can also save out the plots to a PDF.
gcplot.compare_single_level(ref_ds, 'Dataset 1', dev_ds, 'Dataset 2', ilev=ilev, varlist=varlist, pdfname='single_level.pdf')

"""
Zonal Mean Plotting
-------------------
"""
#compare_zonal_mean generates sets of six panel plots containing zonal mean data across your dataset.
#compare_zonal_mean shares many of the same arguments as compare_single_level.
#You can specify pressure ranges in hPa for zonal mean plotting (by default every vertical level is plotted)
gcplot.compare_zonal_mean(ref_ds, 'Dataset 1', dev_ds, 'Dataset 2', pres_range=[0, 100], varlist=varlist, pdfname='zonal_mean.pdf')

