#!/usr/bin/env python

"""
Example program to add a variable containing all zeroes to an
xarray Dataset object.  Useful when you have to insert new
species into GEOS-Chem diagnostic or restart files.

Calling sequence:
-----------------
add_blank_var.py varname infile outfile
"""

import numpy as np
import xarray as xr
from gcpy.util import create_blank_dataarray
from gcpy.constants import skip_these_vars


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

        dset = xr.open_dataset(
            infile,
            drop_variables=skip_these_vars
        )

        if varattrs is None:
            varattrs = dset.attrs

        darr = create_blank_dataarray(
            varname,
            dset.sizes,
            dset.coords,
            varattrs,
            fill_value=0.0,
            fill_type=np.float32
        )

        dset = xr.merge([dset, darr])

        dset.to_netcdf(outfile)


if __name__ == '__main__':

    # Name of the blank varible to add (EDIT AS NEEDED)
    # NOTE: For GCHP, the prefix must be "SPC_" instead of "SpeciesRst_"
    VAR_NAME = "SpeciesRst_PRO2"

    # Variable attributes (EDIT AS NEEDED)
    VAR_ATTRS = {
        "MW_g"            : "146.98",
        "long_name"       : "Dummy species to track production rate of RO2",
        "units"           : "mol mol-1 dry",
        "Is_Gas"          : "true",
        "averaging_method": "instantaneous",
        "_FillValue"      : np.nan
    }

    # Add blank variable to restart file (EDIT FILENAMES AS NEEDED)
    add_blank_var_to_ncfile(
        VAR_NAME,
        'GEOSChem.Restart.20190701_0000z.nc4',
        'new.GEOSChem.Restart.20190701_0000z.nc4',
        varattrs=VAR_ATTRS
    )
