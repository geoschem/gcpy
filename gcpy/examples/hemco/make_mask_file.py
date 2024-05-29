#!/usr/bin/env python3
"""
Create a mask file (for emissions) from a netCDF file of country IDs:
Download this file before using:
https://gcgrid.s3.amazonaws.com/HEMCO/MASKS/v2014-07/countrymask_0.1x0.1.nc

Usage:
------
./make_mask_file.py [-i filein] -o fileout -c country_id -m true|false

where

-i filein     : File of country IDs (download from link above).
                Default value: "countrymask_0.1x0.1.nc"

-o fileout    : Output file for the mask

-c country_id : ID of your desired country.
                Use a netcdf file viewer to determine this value.

-m true|false : Create a mirrored (i.e. inverted) mask.
                Default value: False/

Examples:
---------

# Create a mask for Canada
./make_mask_file.py -o Canada_Mask.01x01.nc -c 124

# Create a mirrored mask for Mexico
./make_mask_file.py -o Mexico_Mask_Mirror.01x01.nc -c 484 -m true

"""
import argparse
import numpy as np
import xarray as xr

def make_mask(
        filein,
        country_id,
        fileout,
        mirror=False,
):
    """
    Creates a netCDF mask file for a given country.

    Args
    filein     : str  : File with country ID values
    country_id : int  : ID of the country that you want masked
    fileout    : str  : Output mask file
    mirror     : bool : Return a mirrored (i.e. inverted) mask?
    """
    with xr.set_options(keep_attrs=True):

        # Define zero and one values (for normal + mirror masks)
        one = np.float32(1)
        zero = np.float32(0)
        if mirror:
            one = np.float32(0)
            zero = np.float32(1)

        # Open file and rename mask variable to MASK
        dset = xr.open_dataset(filein)
        dset = dset.rename({"CountryID": "MASK"})

        # Mask out the country
        array = np.where(
            dset["MASK"].values == country_id,
            one,
            zero
        )

        # Cast to float to avoid issues w/ GCHP input
        dset["MASK"].values = array
        dset["MASK"] = dset["MASK"].astype(np.float32)

        # Write to disk
        dset.to_netcdf(fileout)


if __name__ == '__main__':

    # Tell parser which arguments to expect
    parser = argparse.ArgumentParser(
        description="General cubed-sphere to cubed-sphere regridder."
    )
    parser.add_argument(
        "-i", "--filein",
        metavar="FILEIN",
        type=str,
        required=False,
        default="countrymask_0.1x0.1.nc",
        help="netCDF file with country IDs"
    )
    parser.add_argument(
        "-o", "--fileout",
        metavar="FILEOUT",
        type=str,
        required=True,
        help="name of output file"
    )
    parser.add_argument(
        "-c", "--country-id",
        metavar="COUNTRY-ID",
        required=True,
        type=int,
        help="Country ID value to match in input file",
    )
    parser.add_argument(
        "-m", "--mirror",
        metavar="MIRROR",
        type=bool,
        required=False,
        default=False,
        help="Create mirrored (reversed) mask"
    )
    args = parser.parse_args()
    make_mask(
        args.filein,
        args.country_id,
        args.fileout,
        args.mirror,
    )
