#!/usr/bin/env python
"""
MDA8 Timeseries Calculations with Pandas
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

# First we need to read in some data. For simplicity, let's pull out
# some ozone data I've previously processed from the EPA CASTNET data.
# TODO: Simplify this to an easily-reproducible/distributable example
import feather
o3_data = feather.read_dataframe(
    "/Users/daniel/workspace/air_quality/data/cast/processed/ozone.feather"
)

# Select the data at just one station
station_data = (
    o3_data
    .loc[o3_data.SITE_ID == 'ABT147']
    .set_index('DATE_TIME')
)
print(station_data.head())

# Second, we compute the 8-hour rolling averages for the ozone.
avg_8hr_o3 = (
    station_data['OZONE']
    .rolling(8, min_periods=6)
    .mean()
)

# By default, this takes the last timestamp in a rolling interval; i.e. the
# timestamps correspond to the preceding 8 hours. We want them to refer to
# the proeding 8 hours, so we can adjust them using datetime arithmetic
times = avg_8hr_o3.index.values - pd.Timedelta('8h')
avg_8hr_o3.index.values[:] = times

# Finally, aggregate by calendar day and compute the maxima of the set of
# 8-hour averages for each day
mda8_o3 = avg_8hr_o3.resample('D').max()

# Plot just one year of data
ax = mda8_o3.loc['1994'].plot(figsize=(9, 3), color='k')
ax.set_ylabel("MDA8 O$_3$ [ppb]")
plt.show()
