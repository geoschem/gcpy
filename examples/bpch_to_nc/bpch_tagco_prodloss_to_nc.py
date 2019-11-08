#!/usr/bin/env python
'''
Example script that illustrates how to generate netCDF files containing
the prod/loss fields that are needed for the tagged CO simulation from
the ND65 bpch prod/loss diagnostic.
'''

# Imports
from os.path import join
import gcpy
import xarray as xr
import xbpch as xb
import warnings

# Suppress harmless run-time warnings (mostly about underflow in division)
warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', category=UserWarning)

# ----------------------------------------------------------------------
# Set directory paths (EDIT THESE ACCORDINGLY)
# ----------------------------------------------------------------------
maindir = '/path/to/bpch/data/from/benchmark/simulation'
bpchdir = join(maindir, 'bpch')
ncdir = join(maindir, 'nc')

# ----------------------------------------------------------------------
# Loop over months and read data from each monthly bpch file:
# We will create a new netCDF file for each month
# ----------------------------------------------------------------------
for m in range(1,13):

    # Month string, with padded zeroes ("01", "02", ...)
    mstr = str(m).zfill(2)

    # Print month
    print('Processing month: {}'.format(mstr))

    # File names (EDIT THESE ACCORDINGLY)
    bpch_file = join(bpchdir, 'trac_avg.GC_12.4.0.2016{}01'.format(mstr))
    tinfo_file = join(bpchdir, 'tracerinfo.dat')
    dinfo_file = join(bpchdir, 'diaginfo.dat')
    ncfile = join(ncdir, 'GEOSChem.TagCOInputs.2016{}01.nc'.format(mstr))    

    # Open the bpch file and save it into an xarray Dataset object
    # NOTE: For best results, also specify the corresponding
    # tracerinfo.dat diaginfo.dat metadata files.
    try:
        ds = xb.open_bpchdataset(filename=bpch_file,
                                 tracerinfo_file=tinfo_file,
                                 diaginfo_file=dinfo_file)
    except FileNotFoundError:
        print('Could not find file {}'.format(bpch_file))
        raise

    # ------------------------------------------------------------------
    # Further manipulate the Dataset
    # ------------------------------------------------------------------
    
    # Only extract the bpch fields we need for the tagged CO
    varlist = ['PORL_L_S_PCO', 'PORL_L_S_LCO', 
               'PORL_L_S_PCO_CH4', 'PORL_L_S_PCO_NMVO']
    ds = ds[varlist]

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
    # For each variable in the Dataset:
    # (1) Append the "time" dimension;
    # (2) Edit the variable attributes, removing non-necessary attrs.
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
    ds.attrs['title'] = 'P(CO) and L(CO) for tagged CO simulation'
    ds.attrs['conventions'] = 'COARDS'
    ds.attrs['references'] = 'www.geos-chem.org; wiki.geos-chem.org'

    # ------------------------------------------------------------------
    # Create the netCDF file
    # ------------------------------------------------------------------
    ds.to_netcdf(ncfile)
