#!/usr/bin/env python

"""
Example program to add a variable containing all zeroes to an
xarray Dataset object.  Useful when you have to insert new
species into GEOS-Chem diagnostic or restart files.

Calling sequence:
-----------------
add_blank_var.py varname infile outfile
"""

from sys import argv
import numpy as np
import xarray as xr
from gcpy.util import create_blank_dataarray


def add_blank_var_to_ncfile(
        varname,
        infile,
        outfile,
        varattrs=None
):
    """
    Adds a variable containing all zeroes to a netCDF file

    Args:
    -----
    varname : str
        Name of the variable to add.
    infile: str
        Name of the input netCDF file (will not be overwritten).
    outfile: str
        Name of the output netCDF file containing varname.

    Keyword Args:
    -------------
    varattrs: dict
        A dictionary containing netCDF variable attributes to be written
        to outfile.
    """
    with xr.set_options(keep_attrs=True):

        ds = xr.open_dataset(infile)

        if varattrs is None:
            varattrs = ds.attrs

        da = create_blank_dataarray(
            varname,
            ds.sizes,
            ds.coords,
            varattrs,
            fill_value=0.0,
            fill_type=np.float32
        )

        ds = xr.merge([ds, da])

        ds.to_netcdf(outfile)


if __name__ == '__main__':

    # Name of the blank varible to add (EDIT AS NEEDED)
    # NOTE: For GCHP, the prefix must be "SPC_" instead of "SpeciesRst_"
    varname = "SpeciesRst_PRO2"

    # Variable attributes (EDIT AS NEEDED)
    varattrs = {
        "MW_g"            : "146.98",
        "long_name"       : "Dummy species to track production rate of RO2",
        "units"           : "mol mol-1 dry",
        "Is_Gas"          : "true",
        "averaging_method": "instantaneous",
        "_FillValue"      : np.nan
    }

    # Add blank variable to restart file (EDIT FILENAMES AS NEEDED)
    add_blank_var_to_ncfile(
        varname,
        'GEOSChem.Restart.20190701_0000z.nc4',
        'new.GEOSChem.Restart.20190701_0000z.nc4',
        varattrs=varattrs
    )
