#!/usr/bin/env python
'''
Example script that illustrates how to create a netCDF file
from an old GEOS-Chem binary punch ("bpch") file.
'''

# Imports
import os
from os.path import join
import gcpy
import xarray as xr
import xbpch as xb
import warnings

# Suppress harmless run-time warnings (mostly about underflow in division)
warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', category=UserWarning)

# ----------------------------------------------------------------------
# Set file names (EDIT THESE ACCORDINGLY)
# ----------------------------------------------------------------------
bpch_file = '/path/to/my/bpch_file.bpch'
tinfo_file = '/path/to/my/tracerinfo.dat'
dinfo_file = '/path_to_my/diaginfo.dat'
ncfile = '/path/to/my/netcdf_file.nc'

# Open the bpch file and save it into an xarray Dataset object
# NOTE: For best results, also specify the corresponding
# tracerinfo.dat diaginfo.dat metadata files.
try:
    ds = xb.open_bpchdataset(filename=bpchfile,
                             tracerinfo_file=tinfo_file,
                             diaginfo_file=dinfo_file)
except FileNotFoundError:
    print('Could not find file {}'.format(bpchfile))
    raise  

# ----------------------------------------------------------------------
# Further manipulate the Dataset
# ----------------------------------------------------------------------

# Transpose the order of the xarray Dataset object read by
# xbpch so that its dimensions will be in the same order as
# Dataset objects read from netCDF files.
ds = ds.transpose()

# Convert the bpch variable names to the same naming
# convention as the netCDF ("History") diagnostics.
ds = gcpy.convert_bpch_names_to_netcdf_names(ds)

# xbpch does not include a time dimension, so we'll add one here
coords = ds.coords
coords['time'] = 0.0

# If you are using xarray 0.12.1, then delete the following
# variable attributes to avoid problems writing to netCDF.
for v in ds.data_vars.keys():

    # Append time to the data array
    ds[v] = xr.concat([ds[v]], 'time')

    # Add long_name attribute for COARDS netCDF compliance
    ds[v].attrs['long_name'] = ds[v].attrs['full_name']

    # Remove some extraneous attributes that xbpch sets
    del ds[v].attrs['name']
    del ds[v].attrs['full_name']
    del ds[v].attrs['scale_factor']
    del ds[v].attrs['hydrocarbon']
    del ds[v].attrs['tracer']
    del ds[v].attrs['category']
    del ds[v].attrs['chemical']
    del ds[v].attrs['original_shape']
    del ds[v].attrs['origin']
    del ds[v].attrs['number']
    del ds[v].attrs['molwt']
    del ds[v].attrs['C']

# ------------------------------------------------------------------
# Edit attributes for coordinate dimensions
# ------------------------------------------------------------------

# Time
ds['time'].attrs['long_name'] = 'time'
ds['time'].attrs['units'] = \
    'hours since 2016-{}-01 00:00:00.00 UTC'.format(mstr)
ds['time'].attrs['calendar'] = 'standard'
ds['time'].attrs['axis'] = 'T'

# "lon", "lat", "lev"
ds['lon'].attrs['axis'] = 'X'
ds['lat'].attrs['axis'] = 'Y'
ds['lev'].attrs['axis'] = 'Z'
ds['lev'].attrs['units'] = 'level'

# Global title
ds.attrs['title'] = 'Created by bpch2nc.py'
ds.attrs['conventions'] = 'COARDS'
ds.attrs['references'] = 'www.geos-chem.org; wiki.geos-chem.org'

# ------------------------------------------------------------------
# Create the netCDF file
# ------------------------------------------------------------------
ds.to_netcdf(ncfile)
