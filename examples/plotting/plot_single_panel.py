#!/usr/bin/env python
"""
Global and Regional Single Panel Plots
--------------------------------------
This example script demonstrates the core single panel plotting capabilities of GCPy,
including global and regional single level plots as well as global zonal mean plots.
The example data described here is in lat/lon format, but the same code works equally
well for cubed-sphere (GCHP) data. 
For full documentation on the plotting capabilities of GCPy (including full argument lists), 
please see the GCPy Wiki at https://github.com/geoschem/gcpy/wiki
"""

#xarray allows us to read in any NetCDF file, the format of most GEOS-Chem diagnostics,
#as an xarray Dataset
import xarray as xr
ds = xr.open_dataset('GEOSChem.Restart.20160701_0000z.nc4')

#You can easily view the variables available for plotting using xarray.
#Each of these variables has its own xarray DataArray within the larger Dataset container.
print(ds.data_vars)

#Most variables have some sort of prefix; in this example all variables are
#prefixed with 'SpeciesRst_'. We'll select the DataArray for ozone.
da = ds.SpeciesRst_O3

#Printing a DataArray gives a summary of the dimensions and attributes of the data.
print(da)
#This Restart file has a time dimension of size 1, with 72 vertical levels,
#46 latitude indicies, and 72 longitude indices.
import gcpy.plot as gcplot

"""
Single level plots
------------------
"""
#gcpy.single_panel is the core plotting function of GCPy, able to create a one panel zonal mean or 
#single level plot. Here we will create a single level plot of ozone at ~500 hPa.
#We must manually index into the level that we want to plot (index 22 in the standard 72-layer
#and 47-layer GMAO vertical grids). 
slice_500 = da.isel(lev=22)

#single_panel has many arguments which can be optionally specified. The only argument you must always
#pass to a call to single_panel is the DataArray that you want to plot.
#By default, the created plot includes a colorbar with units read from the DataArray, an automatic title 
#(the data variable name in the DataArray), and an extent equivalent to the full lat/lon extent of the DataArray
import matplotlib.pyplot as plt
gcplot.single_panel(slice_500)
plt.show()

#You can specify a specific area of the globe you would like plotted using the 'extent' argument,
#which uses the format [min_longitude, max_longitude, min_latitude, max_latitude] with bounds [-180, 180, -90, 90] 
gcplot.single_panel(slice_500, extent=[50, -90, -10, 60])
plt.show()

#Other commonly used arguments include specifying a title and a colormap (defaulting to a White-Green-Yellow-Red colormap)
#You can find more colormaps at https://matplotlib.org/tutorials/colors/colormaps.html
gcplot.single_panel(slice_500, title='500mb Ozone over the North Pacific', comap = plt.cm.viridis, 
                 log_color_scale=True, extent=[80, -90, -10, 60])
plt.show()

"""
Zonal Mean Plotting
-------------------
"""

#Use the plot_type argument to specify zonal_mean plotting
gcplot.single_panel(da, plot_type="zonal_mean")
plt.show()

#You can specify pressure ranges in hPa for zonal mean plot (by default every vertical level is plotted)
gcplot.single_panel(da, pres_range=[0, 100], log_yaxis=True, log_color_scale=True)
plt.show()

