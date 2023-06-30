#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
gcpy/benchmark/modules/benchmark_model_vs_obs.py

Python functions to plot modeled O3 from 1-year fullchem benchmark
simulations against observations for the year 2019.

Author: Matt Rowlinson <matthew.rowlinson@york.ac.uk>

Linted with PyLint and incorporated into GCPy
by Bob Yantosca <yantosca@seas.harvard.edu>
'''
import os
from datetime import datetime, timedelta
import pandas as pd
import numpy as np
import xarray as xr


def read_nas(input_file):
    '''
    Read nasa ames data files from EBAS (https://ebas-data.nilu.no)
        Creates data frame of O3 values converted to ppb and dictionary
        with key site information (name, lat, lon, altitude)

    Parameters
    ----------
    input_file : str
                 path to data file

    Returns
    ------
    data : pandas dataframe
        O3 in ppbV
    coords : dict
        Dictionary containing formatted site name: lon, lat and altitude
    '''
    with open(input_file, encoding='UTF-8') as the_file:
        header = np.array(
            [next(the_file) for x in range(155) ]
        )
        n_hdr = int(
            header[0].split(' ')[0]
        )
        st_ymd = header[6].split(' ')
        st_ymd = list(
            filter(
                None,
                st_ymd
            )
        )
        start_date = datetime(
            int(st_ymd[0]),
            int(st_ymd[1]),
            int(st_ymd[2])
        )
        for line in header:
            if 'Station name' in line:
                site = line.split(':')[1:]
                site = '_'.join(site).replace('\n','').replace('  ',' ').replace('/','-')
                site = site.replace('Atmospheric Observatory','')
                site = site.replace(' Research Station','')
            elif 'Station longitude:' in line:
                lon = float(line.split(' ')[-1].replace('\n',''))
            elif 'Station latitude:' in line:
                lat = float(line.split(' ')[-1].replace('\n',''))
            elif 'Station altitude:' in line:
                alt = float(line.split(' ')[-2].replace('\n',''))

    file_hdr = np.loadtxt(
        input_file,
        skiprows=n_hdr
    )
    data = pd.DataFrame(
        file_hdr,
        index=file_hdr[:,0]
    )
    data, flag = find_times(
        data,
        start_date
    )
    data = pd.DataFrame(
        {
            'Value': data.values/1.99532748,
            'Flag': flag
        },
        index=data.index
    )
    data = data[data.Flag == 0.000]
    data = data.loc['2019']
    data = data.resample('H').mean()
    data = pd.DataFrame(
        {
            site:data.Value
        },
        index=data.index
    )
    coords= { site:
          {
              'lon': lon,
              'lat': lat,
              'alt': alt
          }
    }
    return data, coords


def find_times(data, start_time):
    '''
    Convert timestamps in nasa ames data files to python datetime objects
        Set DataFrame index to the new datetime array

    Parameters
    ----------
    data : pandas DataFrame
        DataFrame with O3 values from GAW site
    start_time : str
        Reference start time for timestamp taken from nasa ames file

    Returns
    ------
    data : pandas dataFrame
        O3 in ppbV with datetime index
    flag : pandas dataframe
        QC flag with datetime index
    '''
    end_time=data[data.columns[1]]
    time_x=[]
    for i in range(len(end_time)):
        time_x.append(start_time + timedelta(days=end_time.values[i]))
    data.index = time_x
    flag = data[data.columns[-1]]
    data = data[data.columns[2]]

    return data, flag


def find_file_list(path, substrs):
    '''
    Find list of output files in given location containing specified string

    Parameters
    ----------
    path : str
        Path to rundir where output is stored
    substrs : list
        List of strings to search for in destination

    Returns
    ------
    file_list : list
        List of available files instring format
    '''
    file_list =[]
    for root, directory, files in os.walk(path):
        for file_name in files:
            for substr in substrs:
                if substr in file_name:
                    file_list.append(os.path.join(root, file_name))
    file_list.sort()
    return file_list


def get_data_as_xr(path, dates='2019'):
    '''
    Read GEOS-Chem SpeciesConc_O3 in ppbV from desired rundir

    Parameters
    ----------
    path : str
        Path to rundir where output is stored
    dates : list
        List of strings to search for in destination,
        used to select specific dates

    Returns
    ------
    ds : xarray DataSet
        GEOS-Chem output from given rundir as an xarray dataset
    '''
    ds_o3 = xr.open_mfdataset(
        find_file_list(
            path,
            [f'SpeciesConc.{dates}']
        )[:],
        combine='by_coords'
    )
    ds_o3 = ds_o3['SpeciesConc_O3'] * 1e9 ## Read O3 data and convert to ppbV

    return ds_o3


def find_nearest_3d(ds_gc, lon_value, lat_value, alt_value):
    '''
    Find GEOS-Chem gridbox closest to the observational dataset.
        Uses lat, lon and alt from obs to select most appropriate GC data.

    Parameters
    ----------
    ds_gc : xarray dataset
        GEOS-Chem output to be processed
    lon_value : float
        GAW site longitude
    lat_value : float
        GAW site latitude
    alt_value : float
        GAW site altitude
        List of strings to search for in destination,
        used to select specific dates

    Returns
    ------
    x_idx, y_idx, z_idx
        GEOS-Chem for single gridbox closest to GAW site specifications
    '''
    gc_alts = pd.read_csv(
        './GC_72_vertical_levels.csv'
    )['Altitude (km)'] * 1e3

    x_idx=(
        np.abs(
            ds_gc.lon - float(lon_value)
        )
    ).argmin()
    
    y_idx=(
        np.abs(
            ds_gc.lat - float(lat_value)
        )
    ).argmin()

    z_idx=(
        np.abs(
            gc_alts.values - float(alt_value)
        )
    ).argmin()

    return x_idx, y_idx, z_idx
