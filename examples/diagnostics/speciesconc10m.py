#!/usr/bin/env python

'''
Example script showing how to read, process, and plot SpeciesConc10m
diagnostic data from GEOS-Chem (from the ConcAboveSfc collection).
'''

# Imports
import xarray as xr
from gcpy import core
import matplotlib.pyplot as plt

# ----------------------------
# Read and process the data!
# ----------------------------

# Path to a ConcAboveSpc collection file
# (YOU MUST EDIT THIS FOR YOUR SETUP!!)
filename = '/path/to/GEOSChem.ConcAboveSfc.YYYYMMDD_hhmmz.nc4'

# Open dataset
try:
    ds = xr.open_dataset(filename)
except FileNotFoundError:
    raise('Could not find file {}'.format(filename))

# The time-averaged SpeciesConc10m collection contains the cumulative
# sum of species concentration at 10m altitude (the actual altitude
# can be changed in the GEOS-Chem input.geos file).
#
# Because the computation is only valid where the Monin-Obhukov
# similarity (z/L <= 1) holds, we must divide the SpeciesConc10m
# variables by the DryDepZLfrac variable.  DryDepZLfrac is the
# fraction of the time z/L <= 1 at each grid box.
#
# We will use the convenience routine "divide_dataset_by_dataarray".
# This will divide all of the selected variables in a Dataset object
# by a single variable.
varlist = [v for v in ds.data_vars.keys() if 'SpeciesConc10m' in v]
ds = core.divide_dataset_by_dataarray(ds, ds['DryDepZLfrac'], varlist=varlist)

# ----------------------------
# Plot the data!
# ----------------------------

# We would also like to plot the DryDepZLfrac variable,
# so let's append it to the list of variables.
varlist.append('DryDepZLfrac')

# Plot variables
for v in varlist:
    print('Now plotting {}'.format(v))

    # Take the first time slice of this variable.
    dr = ds[v].isel(time=0)

    # Convert SpeciesConc10m units to ppbv.
    if 'SpeciesConc10m' in v:
        dr = dr * 1e9
        dr.attrs['units'] = 'ppbv'

    # Define the colorbar label to show the name of the variable and the units.
    # This needs to be defined as a Python dict passed via xr.plot.imshow.
    cb_label = { "label": "{}   ({})".format(v, dr.units) }

    # Create the plot using the matplotlib default "rainbow" color table.
    xr.plot.imshow(dr, cmap="rainbow", cbar_kwargs=cb_label)

    # Display the plot!
    plt.show()
