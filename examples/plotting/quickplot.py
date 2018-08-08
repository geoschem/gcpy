#!/usr/bin/env python

'''
quickplot.py: Can be used to generate a quick-look plot of a variable
in a netCDF file.  The plotting is bare-bones, without a gridded-map
or choice of color scale.  This is useful, for example, when you need to
make a "sanity-check" plot to make sure that a given model output
was generated properly.

NOTE: Python indexing starts at 0!  Make sure to subtract 1 from the time
and lev values that you supply!
'''

# Imports
import sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt


def do_the_plot(netcdf_file, variable, time, lev):

    '''
    Method do_the_plot produces a quick-look lon-lat plot of a given
    VARIABLE in a NETCDF_FILE, for a given TIME slice and level (LEV) slice.
    '''

    # Read the entire file into an xarray Dataset object
    ds = xr.open_dataset(netcdf_file)

    # List of data variables (excludes index variables)
    varlist = ds.data_vars.keys()
    varlist = [v for v in varlist if ds[v].ndim > 2]

    # Make sure the variable is located in the file name
    if variable not in varlist:
        raise Exception(
            'Variable {} is not in {}!'.format(variable, netcdf_file))

    # Make sure the requested time slice is in range
    if time < 0 or time > ds.dims['time']-1:
        max_val = ds.dims['time']-1
        err_msg = 'time={}! Must be between 0 and {}!'.format(time, max_val)
        raise Exception(err_msg)

    # Make sure the requrested level slice (lev) is in range
    if lev < 0 or lev > ds.dims['lev']-1:
        max_val = ds.dims['lev']-1
        err_msg = 'time={}! Must be between 0 and {}!'.format(lev, max_val)
        raise Exception(err_msg)

    # Extract data for the given species into an xarray DataArray object
    dr = ds[variable].isel(time=time, lev=lev)

    # Get the time and level values for printout
    time_val = dr['time'].data
    lev_val = dr['lev'].data

    # Informational printout
    print('netCDF file name:      {}'.format(netcdf_file))
    print('Requested variable:    {}'.format(variable))
    print('Requested time slice:  {}; value = {}'.format(time, time_val))
    print('Requested level slice  {}; value = {}'.format(lev, lev_val))
    print('Min value of data      {} {}'.format(np.min(dr.values), dr.units))
    print('Max value of data      {} {}'.format(np.max(dr.values), dr.units))

    # Define the colorbar label to show the name of the variable and the units.
    # This needs to be defined as a Python dict passed via xr.plot.imshow.
    cb_label = { "label": "{}   ({})".format(variable, dr.units) }

    # Plot the data using the matplotlib default "rainbow" color table
    xr.plot.imshow(dr, cmap="rainbow", cbar_kwargs=cb_label)
    plt.show()
    
if __name__ == "__main__":

    '''
    Script driver.  This section only gets executed if we run this script
    at the Unix command prompt.
    '''

    # Number of arguments
    # IMPORTANT NOTE: Similar to C/C++, the first argument is always the name
    # of this program.  Positional arguments start with sys.argv[1].
    n_args = len(sys.argv)

    # 4 arguments passed
    if n_args == 5:
        netcdf_file = sys.argv[1]
        variable = sys.argv[2]
        time = int(sys.argv[3])  # Need to convert string to int
        lev = int(sys.argv[4])  # or else isel will throw an error!

    # 3 arguments passed (set lev to zero by default)
    elif n_args == 4:
        netcdf_file = sys.argv[1]
        variable = sys.argv[2]
        time = int(sys.argv[3])
        lev = 0

    # 2 arguments passed (set time & lev to zero by default)
    elif n_args == 3:
        netcdf_file = sys.argv[1]
        variable = sys.argv[2]
        time = 0
        lev = 0

    # Otherwise throw an error
    else:
        raise Exception("Usage: quickplot.py FILENAME SPECIES")

    # Make the plot
    do_the_plot(netcdf_file, variable, time, lev)
