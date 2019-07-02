#!/usr/bin/env python
'''
Example script that illustrates how to create a netCDF file
from an old GEOS-Chem binary punch ("bpch") file.
'''
import gcpy
import xarray as xr
import xbpch as xb

# File names (edit these accordingly)
bpch_file = '/path/to/my/bpch_file.bpch'
tinfo_file = '/path/to/my/tracerinfo.dat'
dinfo_file = '/path_to_my/diaginfo.dat'
ncfile = '/path/to/my/netcdf_file.nc'

# Open the bpch file and save it into an xarray Dataset object
# NOTE: For best results, also specify the corresponding
# tracerinfo.dat diaginfo.dat metadata files.
ds = xb.open_bpchdataset(filename=bpchfile,
                         tracerinfo_file=tinfo_file,
                         diaginfo_file=dinfo_file)

# Transpose the order of the xarray Dataset object read by
# xbpch so that its dimensions will be in the same order as
# Dataset objects read from netCDF files.
ds = ds.transpose()

# Convert the bpch variable names to the same naming
# convention as the netCDF ("History") diagnostics.
ds = gcpy.convert_bpch_names_to_netcdf_names(ds)

# If you are using xarray 0.12.1, then delete the following
# variable attributes to avoid problems writing to netCDF.
for v in ds.data_vars.keys():
    del ds[v].attrs['scale_factor']
    del ds[v].attrs['hydrocarbon']
    del ds[v].attrs['tracer']
    del ds[v].attrs['category']
    del ds[v].attrs['chemical']
    del ds[v].attrs['original_shape']
    del ds[v].attrs['origin']

# Create the netCDF file
ds.to_netcdf(ncfile)
