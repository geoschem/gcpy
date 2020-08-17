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
import gcpy.plot 
#gcpy.gcplot is the core plotting function of GCPy, able to create a one panel zonal mean or 
#single level plot. Here we will create a single level plot of ozone at ~500 hPa.
#We must manually index into the level that we want to plot (index 22 in the standard 72-layer
#and 47-layer GMAO vertical grids). 
slice_500 = da.isel(lev=22)

#gcplot has many arguments which can be optionally specified. The only argument you must always
#pass to a call to gcplot is the DataArray that you want to plot.
#By default, the created plot includes a colorbar with units read from the DataArray, an automatic title 
#(the data variable name in the DataArray), and an extent equivalent to the full lat/lon extent of the DataArray
import matplotlib.pyplot as plt
gcpy.plot.gcplot(slice_500)
plt.show()

#You can specify a specific area of the globe you would like plotted using the 'extent' argument,
#which users the format [min_longitude, max_longitude, min_latitude, max_latitude]
gcpy.plot.gcplot(slice_500, extent=[50, -90, -10, 60])


################################################
# To create a cartopy plot, we have to pass a **projection** to our matplotlib
# axis. This contains the utility/logic to transform data from the coordinate
# system we give to matplotlib to the actual cartographic/geographic system.
#
# Now, when we plot on this axis, we just need to tell matplotlib how to transform
# our data to its coordinate system. We do this by passing a **transform** object.
# No matter what projection you use, if you're plotting rectangular lat/lon data
# from GEOS-Chem or another model, you should use the ``PlateCarree()`` projection
# for this transform.
#
# For convenience, we'll use the wrapper that each ``DataArray`` has for
# matplotlib functions. This is just for convenience; we could always get the
# vector of lon/lat values and plot the data in **avg_o3** directly.

def gcplot(plot_vals,
           ax=None,
           plot_type="single_level",
           grid={},
           gridtype="",
           title="fill",
           comap=WhGrYlRd,
           norm=[],
           unit="",
           extent=(None, None, None, None),
           masked_data=None,
           use_cmap_RdBu=False,
           log_color_scale=False,
           add_cb=True,
           pres_range = [0, 2000],
           pedge=np.full((1, 1), -1),
           pedge_ind=np.full((1,1), -1),
           log_yaxis=False,
           xtick_positions=np.arange(-90,91,30),
           xticklabels = []
):


import cartopy.crs as ccrs
fig = plt.figure(figsize=(8, 4))
ax = fig.add_subplot(111, projection=ccrs.PlateCarree())

avg_o3.plot.imshow(x='lon', y='lat', vmin=0, vmax=60., cmap='gist_stern',
                   ax=ax, transform=ccrs.PlateCarree())
ax.coastlines()
plt.show()
