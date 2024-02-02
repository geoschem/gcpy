#!/usr/bin/env python
"""
MDA8 Timeseries Calculations
========================================

A common statistic used when constructing standards for air quality
criteria pollutants is to look at the ranked distribution of the
daily maxima of rolling 8-hour averages of a substance, or MDA8 for
short.
"""
# Author: Daniel Rothenberg
# Version: June 1, 2017
#
# NOTE: Now placed into a function so that we only execute
# these commands when running as a standalone script.
# (Bob Yantosca, 09 Aug 2023)

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import pandas as pd
import xarray as xr
plt.style.use(['seaborn-v0_8-talk', 'seaborn-v0_8-ticks'])


def main():
    """
    MDA8 O3 timeseries example
    """

    # Read hourly data ### EDIT THIS IF NECESSARY
    InFile = 'GEOSChem.Hourly_SfcO3.2017.nc'
    ds = xr.open_dataset(InFile)
    o3_data = ds['SpeciesConc_O3']

    # Compute the 8-hour rolling averages for ozone
    avg_8hr_o3 = (o3_data.rolling(time=8, min_periods=6).mean())

    # By default, this takes the last timestamp in a rolling interval; i.e. the
    # timestamps correspond to the preceding 8 hours. We want them to refer to
    # the proeding 8 hours, so we can adjust them using datetime arithmetic
    times_np = avg_8hr_o3.time.values
    times_pd = pd.to_datetime(times_np) - pd.Timedelta('8h')
    avg_8hr_o3.time.values[:] = times_pd

    # Finally, aggregate by calendar day and compute the maxima of the set of
    # 8-hour averages for each day
    mda8_o3 = avg_8hr_o3.resample(time='D').max(dim='time')
    mda8_o3.name='mda8_o3'

    # Save output to new netCDF file
    mda8_o3.to_netcdf(
        'GEOSChem.MDA8_O3.20170.nc', 'w',
        format='NETCDF4',
        encoding={'lat': {'_FillValue': None},
                  'lon': {'_FillValue': None},
                  'time': {'_FillValue': None},
                  'mda8_o3': {'_FillValue': None}}
    )

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

    ax.xaxis.set_major_formatter(mdates.DateFormatter("%h %d"))
    for tick in ax.xaxis.get_majorticklabels():
        tick.set_horizontalalignment('center')

    ax.set_xlabel("")
    ax.set_ylabel("(MDA8) O$_3$ [ppb]")

    plt.show()


# Only execute when running as a standalone script
if __name__ == '__main__':
    main()
