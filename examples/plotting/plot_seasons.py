#!/usr/bin/env python
"""
Seasonal Panel Plot
===================

For a given field timeseries, compute seasonal averages over all data and
plot each average on a four-panel figure.
"""
# Author: Daniel Rothenberg
# Version: June 2, 2017

import matplotlib.pyplot as plt
plt.style.use(['seaborn-talk', 'seaborn-ticks'])

import xbpch

# First we read in a sample dataset containing GEOS-Chem output
ds = xbpch.open_bpchdataset(
    "/Users/daniel/workspace/bpch/test_data/ref_e2006_m2008.bpch",
    diaginfo_file="/Users/daniel/Desktop/sample_nd49/diaginfo.dat",
    tracerinfo_file="/Users/daniel/Desktop/sample_nd49/tracerinfo.dat",
    dask=True, memmap=True
)

# Compute seasonal averages by doing splitting along the "seasons" corresponding
# to each timestep, and taking the average over time in that group
seasonal_o3 = (
    ds['IJ_AVG_S_O3']
    .isel(lev=0)  # select just surface values
    .groupby('time.season').mean('time')
)
print(seasonal_o3)

# Note that we now have a new dimension, "season", corresponding to the groups
# we split the dataset by. We can use this dimension as a 'facet' to assemble
# a collection of plots.
# TODO: Cleanup axis proportions
import cartopy.crs as ccrs
g = seasonal_o3.plot.imshow(
    x='lon', y='lat', # Use lat/lon for the x/y axis coordinates
    vmin=0, vmax=60., cmap='gist_stern', # Colormap settings for all panels
    col='season', col_wrap=2, # Facet over 'season', with 2 columns on the grid
    transform=ccrs.PlateCarree(), # Geographic transform for coordinates
    subplot_kws=dict(projection=ccrs.PlateCarree())
        # Have each subpanel use this cartographic projection
)
for ax in g.axes.ravel():
    ax.coastlines()
plt.show()
