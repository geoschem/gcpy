#!/usr/bin/env python

'''
This Python script concatenates several individual netCDF files
into a single netCDF file using xarray.

Calling sequence:
    ./concatentate_files.py

Remarks:
    If you have several individual files with one variable per file,
    you should consider concatenating them into a single file.
    This is often more efficient, as opening each netCDF file incurs
    computational overhead.  It is usually faster to read data from
    a file with multiple variables than to having to open several
    files with one variable each.
'''

# Imports
from os.path import join
import xarray as xr
from xarray.coding.variables import SerializationWarning
import numpy as np
import warnings

# Suppress harmless run-time warnings (mostly about underflow or NaNs)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=SerializationWarning)

# Use  open_mfdataset to open all files into a single Dataset
# (Edit file path accordingly for your setup)
indir = '/path/to/my/netcdf/files/'
infiles = join(indir, '.nc*')
ds = xr.open_mfdataset(infiles) 

# Keep all netCDF attributes
with xr.set_options(keep_attrs=True):
    
    # Loop over all variables in the Dataset
    for v in ds.data_vars.keys():
        
        # OPTIONAL STEP:
        # Xarray will try convert missing values to NaN's,
        # so you may need to replace these with zeroes.
        #
        # If your netCDF files represent e.g. emissions, 
        # or other physical quantities, you may want to 
        # replace these with zeros, so that NaNs won't 
        # get read into atmospheric models, etc.
        #
        # NOTE: ds[v].values converts to a numpy ndarray,
        # so that you can use numpy functions.
        ds[v].where(np.isnan(ds[v].values), other=0.0, drop=False)
        
        # OPTIONAL: Print min & max for each variable
        # Comment out if you wish
        print('{} : {} {}'.format(
            v, np.min(ds[v].values), np.max(ds[v].values)))

# Write to the output file to the main path
# (Edit file path and file name accordingly for your setup)
outdir = '/path/to/my/output/file'
outfile = join(outdir, 'my_concatenated_output_file.nc')
ds.to_netcdf(outfile)
