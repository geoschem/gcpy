#!/usr/bin/env python
'''
Example script that illustrates how to create a netCDF file
from an old GEOS-Chem binary punch ("bpch") file.
'''

# Imports
import gcpy
import xarray as xr
import xbpch as xb
import warnings

# Suppress harmless run-time warnings (mostly about underflow in division)
warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', category=UserWarning)

# ----------------------------------------------------------------------
# User configurable settings (EDIT THESE ACCORDINGLY)
# ----------------------------------------------------------------------

# Name of Bpch file
bpchfile = '/path/to/bpch/file'

# tracerinfo.dat and diaginfo,dat fiels
tinfo_file = '/path/to/tracerinfo.dat'
dinfo_file = '/path/to/diaginfo.dat'

# Name of netCDF file
ncfile = '/path/to/netcdf/file'

# Date string for the time:units attribute
datestr = 'YYYY-MM-DD'

# Number of seconds in the diagnostic interval (assume 1-month)
interval = 86400.0 * 31.0

# ----------------------------------------------------------------------
# Open the bpch file and save it into an xarray Dataset object
# NOTE: For best results, also specify the corresponding
# tracerinfo.dat diaginfo.dat metadata files.
# ----------------------------------------------------------------------
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

# ------------------------------------------------------------------
# Further edit variable attributes
# ------------------------------------------------------------------
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

    # Make the units attribute consistent with the units
    # attribute from the GEOS-Chem History diagnostics
    # NOTE: There probably is a more Pythonic way to code
    # this, but this will work for sure.
    if 'ug/m3' in ds[v].units:
        ds[v].attrs['units'] = 'ug m-3'
    if 'ug Celsius/m3' in ds[v].units:
        ds[v].attrs['units'] = 'ug C m-3'
    if 'count/cm3' in ds[v].units:
        ds[v].attrs['units'] = 'molec m-3'               
    if 'cm/s' in ds[v].units:
        ds[v].attrs['units'] = 'cm s-1'
    if 'count/cm2/s' in ds[v].units:
        ds[v].attrs['units'] = 'molec cm-2 s-1'
    if 'kg/m2s' in ds[v].units:
        ds[v].attrs['units'] = 'kg m-2 s-1'
    if 'kg/m2/s' in ds[v].units:
        ds[v].attrs['units'] = 'kg m-2 s-1'
    if 'kg/s' in ds[v].units:
        ds[v].attrs['units'] = 'kg s-1'
    if 'W/m2' in ds[v].units:
        ds[v].attrs['units'] = 'W m-2'
    if 'm/s' in ds[v].units:
        ds[v].attrs['units'] = 'm s-1'
    if 'Pa/s' in ds[v].units:
        ds[v].attrs['units'] = 'Pa s-1'
    if 'g/kg' in ds[v].units:
        ds[v].attrs['units'] = 'g kg-1'
    if v.strip() == 'TotalOC':
        ds[v].attrs['units'] = 'ug m-3'
    if v.strip() in [ 'HO2concAfterChem']:
        ds[v].attrs['units'] = 'ppb'
    if v.strip() in ['O1DconcAfterChem', 
                     'O3PconcAfterChem',
                     'OHconcAfterChem']:
        ds[v].attrs['units'] = 'molec cm-3'
    if v.strip() in ['Loss_CO', 'Prod_CO',
                     'Loss_Ox', 'Prod_Ox', 'Prod_SO4']:
        ds[v].attrs['units'] = 'molec/cm3/s'
    if v.strip() in 'Met_CLDTOPS':
        ds[v].attrs['units'] = 'level'
    if v.strip() in 'Met_PHIS':
        ds[v].attrs['units'] = 'm2 s-1'
    if v.strip() in ['Met_PRECCON', 'Met_PRECTOT']:
        ds[v].attrs['units'] = 'kg m-2 s-1'
    if v.strip() in 'Met_AVGW':
        ds[v].attrs['units'] = 'vol vol-1'
    if v.strip() in 'Met_AIRNUMDEN':
        ds[v].attrs['units'] = 'molec cm-3'
    if v.strip() in ['ProdCOfromCH4', 'ProdCOfromNMVOC']:
        ds[v].attrs['units'] = 'molec cm-3 s-1'

    # Convert these prodloss diagnostics from kg (bpch) to kg/s
    # to be consistent with the GEOS-Chem History diagnostics
    # NOTE: Assume a 1-month interval (
    if v.strip() in ['ProdSO4fromH2O2inCloud',    'ProdSO4fromO3inCloud', 
                     'ProdSO4fromO2inCloudMetal', 'ProdSO4fromO3inSeaSalt',
                     'ProdSO4fromHOBrInCloud',    'ProdSO4fromSRO3',
                     'ProdSO4fromSRHObr',         'ProdSO4fromO3s']:
        ds[v].attrs['units'] = 'kg S s-1'
        ds[v] = ds[v] / interval
    if v.strip() in ['LossHNO3onSeaSalt']:
        ds[v].attrs['units'] = 'kg s-1'
        ds[v] = ds[v] / interval

# ------------------------------------------------------------------
# Edit attributes for coordinate dimensions
# ------------------------------------------------------------------

# Time
ds['time'].attrs['long_name'] = 'time'
ds['time'].attrs['units'] = \
    'hours since {} 00:00:00.00 UTC'.format(datestr)
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
