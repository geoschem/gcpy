#!/usr/bin/env python

"""
Adds an extra DataArray for into restart files.

Calling sequence:
    ./insert_species_into_restart.py
"""

# Imports
import warnings
import xarray as xr
from xarray.coding.variables import SerializationWarning
from gcpy import constants

# Suppress harmless run-time warnings (mostly about underflow or NaNs)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=SerializationWarning)

def main():
    """
    Appends extra species to restart files.
    """

    # Data vars to skip
    skip_vars = constants.skip_these_vars

    # List of dates (EDIT accordingly)
    file_list = [
        'GEOSChem.Restart.fullchem.20190101_0000z.nc4',
        'GEOSChem.Restart.fullchem.20190701_0000z.nc4',
        'GEOSChem.Restart.TOMAS15.20190701_0000z.nc4',
        'GEOSChem.Restart.TOMAS40.20190701_0000z.nc4',
        'GCHP.Restart.fullchem.20190101_0000z.c180.nc4',
        'GCHP.Restart.fullchem.20190101_0000z.c24.nc4',
        'GCHP.Restart.fullchem.20190101_0000z.c360.nc4',
        'GCHP.Restart.fullchem.20190101_0000z.c48.nc4',
        'GCHP.Restart.fullchem.20190101_0000z.c90.nc4',
        'GCHP.Restart.fullchem.20190701_0000z.c180.nc4',
        'GCHP.Restart.fullchem.20190701_0000z.c24.nc4',
        'GCHP.Restart.fullchem.20190701_0000z.c360.nc4',
        'GCHP.Restart.fullchem.20190701_0000z.c48.nc4',
        'GCHP.Restart.fullchem.20190701_0000z.c90.nc4'
    ]

    # Keep all netCDF attributes
    with xr.set_options(keep_attrs=True):

        # Loop over dates
        for file_name in file_list:

            # Input and output files
            infile = '../' + file_name
            outfile = file_name

            print("Creating " + outfile)

            # Open input file
            dset = xr.open_dataset(infile, drop_variables=skip_vars)

            # Create a new DataArray from a given species (EDIT ACCORDINGLY)
            if "GCHP" in infile:
                darr = dset["SPC_ETO"]
                darr.name = "SPC_ETOO"
            else:
                darr = dset["SpeciesRst_ETO"]
                darr.name = "SpeciesRst_ETOO"

            # Update attributes (EDIT ACCORDINGLY)
            darr.attrs["FullName"] = "peroxy radical from ethene"
            darr.attrs["Is_Gas"] = "true"
            darr.attrs["long_name"] = "Dry mixing ratio of species ETOO"
            darr.attrs["MW_g"] = 77.06

            # Merge the new DataArray into the Dataset
            dset = xr.merge([dset, darr], compat="override")

            # Create a new file
            dset.to_netcdf(outfile)

            # Free memory by setting dset to a null dataset
            dset = xr.Dataset()


# Only execute when we run as a standalone script
if __name__ == "__main__":
    main()
