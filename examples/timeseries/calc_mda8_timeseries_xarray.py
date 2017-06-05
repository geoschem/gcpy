#!/usr/bin/env python
"""
MDA8 Timeseries Calculations with xarray
========================================

A common statistic used when constructing standards for air quality
criteria pollutants is to look at the ranked distribution of the
daily maxima of rolling 8-hour averages of a substance, or MDA8 for
short.
"""
# Author: Daniel Rothenberg
# Version: June 1, 2017

import matplotlib.pyplot as plt
plt.style.use(['seaborn-talk', 'seaborn-ticks'])

import pandas as pd
import xarray as xr
import xbpch
from dask.diagnostics import ProgressBar

from os.path import join

# First we need to read in some data. We'll read a multi-file ND49 BPCH
# dataset using the xbpch package.
dates = ["200601{:02d}".format(d) for d in range(1, 22)]
ROOT = "/Users/daniel/workspace/bpch/test_data/"
fns = [
    join(ROOT, "ND49_{}_ref_e2006_m2010.bpch".format(date))
    for date in dates
]
nd49_data = xbpch.open_mfbpchdataset(
    fns, diaginfo_file="/Users/daniel/Desktop/sample_nd49/diaginfo.dat",
    tracerinfo_file="/Users/daniel/Desktop/sample_nd49/tracerinfo.dat",
    dask=True, memmap=True,
)
o3_data = nd49_data['IJ_AVG_S_O3']
with ProgressBar():
    print("Loading data into memory")
    o3_data.load()

# Second, we compute the 8-hour rolling averages for the ozone.
avg_8hr_o3 = (
    o3_data
    .rolling(time=8, min_periods=6)
    .mean()
)

# By default, this takes the last timestamp in a rolling interval; i.e. the
# timestamps correspond to the preceding 8 hours. We want them to refer to
# the proeding 8 hours, so we can adjust them using datetime arithmetic
times_np = avg_8hr_o3.time.values
times_pd = pd.to_datetime(times_np) - pd.Timedelta('8h')
avg_8hr_o3.time.values[:] = times_pd

# Finally, aggregate by calendar day and compute the maxima of the set of
# 8-hour averages for each day
mda8_o3 = avg_8hr_o3.resample('D', dim='time', how='max')

# Select data for one specific location, near Boston
boston_mda8_o3 = mda8_o3.sel(lon=-71., lat=42., method='nearest')
boston_o3 = o3_data.sel(lon=-71., lat=42., method='nearest')

# Plot both the original (hourly) and MDA* timeseries on the same plot.
fig = plt.figure(figsize=(9, 3))
ax = fig.add_subplot(111)
boston_o3.plot(ax=ax, color='k')
ax.stem(boston_mda8_o3.time.values, boston_mda8_o3.data,
        ':r', markerfmt='ro')
ax.set_ylim(0)

import matplotlib.dates as mdates
ax.xaxis.set_major_formatter(mdates.DateFormatter("%h %d"))
for tick in ax.xaxis.get_majorticklabels():
    tick.set_horizontalalignment('center')

ax.set_xlabel("")
ax.set_ylabel("(MDA8) O$_3$ [ppb]")

plt.show()
