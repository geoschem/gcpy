#!/usr/bin/env python
"""
Global Average Plot
===================

Using cartopy_ we can easily visualize gridded model output on maps with
different cartographic projections, and then configure them with any
aesthetics or features such as continents, geopolitical borders, gridlines,
and more.

.. _cartopy: http://scitools.org.uk/cartopy/docs/latest/
"""
# Author: Daniel Rothenberg
# Version: June 2, 2017

import matplotlib.pyplot as plt
plt.style.use(['seaborn-talk', 'seaborn-ticks'])

import xbpch

# First we read in a sample dataset containing GEOS-Chem output
ds = xbpch.open_bpchdataset(
    "/Users/daniel/workspace/bpch/test_data/ND49_20060101_ref_e2006_m2010.bpch",
    diaginfo_file="/Users/daniel/Desktop/sample_nd49/diaginfo.dat",
    tracerinfo_file="/Users/daniel/Desktop/sample_nd49/tracerinfo.dat",
    dask=False, memmap=True
)

# This dataset contains multiple timesteps of multiple tracers, but we can
# extract just ozone, and average over the timesteps to get a slice to plot.
avg_o3 = ds['IJ_AVG_S_O3'].mean('time')

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

import cartopy.crs as ccrs
fig = plt.figure(figsize=(8, 4))
ax = fig.add_subplot(111, projection=ccrs.PlateCarree())

avg_o3.plot.imshow(x='lon', y='lat', vmin=0, vmax=60., cmap='gist_stern',
                   ax=ax, transform=ccrs.PlateCarree())
ax.coastlines()
plt.show()
