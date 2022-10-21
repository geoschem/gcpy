#!/usr/bin/env python

"""
Adds an extra DataArray for into restart files.

Calling sequence:
    ./insert_species_into_restart.py
"""

# Imports
import gcpy.constants as gcon
import xarray as xr
from xarray.coding.variables import SerializationWarning
import warnings

# Suppress harmless run-time warnings (mostly about underflow or NaNs)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=SerializationWarning)

def main():
    """
    Appends extra species to restart files.
    """

    # Data vars to skip
    skip_vars = gcon.skip_these_vars

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
        for f in file_list:

            # Input and output files
            infile = '../' + f
            outfile = f

            print("Creating " + outfile)

            # Open input file
            ds = xr.open_dataset(infile, drop_variables=skip_vars)

            # Create a new DataArray from a given species (EDIT ACCORDINGLY)
            if "GCHP" in infile:
                dr = ds["SPC_ETO"]
                dr.name = "SPC_ETOO"
            else:
                dr = ds["SpeciesRst_ETO"]
                dr.name = "SpeciesRst_ETOO"

            # Update attributes (EDIT ACCORDINGLY)
            dr.attrs["FullName"] = "peroxy radical from ethene"
            dr.attrs["Is_Gas"] = "true"
            dr.attrs["long_name"] = "Dry mixing ratio of species ETOO"
            dr.attrs["MW_g"] = 77.06

            # Merge the new DataArray into the Dataset
            ds = xr.merge([ds, dr], compat="override")

            # Create a new file
            ds.to_netcdf(outfile)

            # Free memory by setting ds to a null dataset
            ds = xr.Dataset()

if __name__ == "__main__":
    main()
