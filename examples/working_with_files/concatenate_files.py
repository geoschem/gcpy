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
import gcpy.constants as gcon
import os
import numpy as np
import xarray as xr
from xarray.coding.variables import SerializationWarning
import warnings

# Suppress harmless run-time warnings (mostly about underflow or NaNs)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=SerializationWarning)


def find_files_in_dir(path, substrs):
    '''
    Returns a list of all files in a directory that match one or more
    substrings.

    Args:
    -----
        path : str
            Path to the directory in which to search for files.

        substrs : list of str
            List of substrings used in the search for files.

    Returns:
    --------
        file_list : list of str
            List of files in the directory (specified by path)
            that match all substrings (specified in substrs).
    '''

    # Initialize
    file_list = []

    # Walk through the given data directory.  Then for each file found,
    # add it to file_list if it matches text in search_list.
    for root, directory, files in os.walk(path):
        for f in files:
            for s in substrs:
                if s in f:
                    file_list.append(os.path.join(root, f))

    # Return an alphabetically sorted list of files
    file_list.sort()
    return file_list


def replace_nans_with_zeroes(ds, verbose=True):
    '''
    Replaces NaN values with zeroes for each variable
    within an an xarray Dataset.

    Args:
    ----
        ds : xarray Dataset
            The input dataset, containing one or more data variables.

    Keyword Args (optional):
    ------------------------
        verbose : boolean
            Set this switch to print out the variable name, as well
            as the min and max of the variable.  This will illustrate
            the replacement of NaNs with zeroes.
    '''

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
            if verbose:
                print('{} : {} {}'.format(
                    v, np.min(ds[v].values), np.max(ds[v].values)))

    # Return the modified Datast
    return ds


def main():
    '''
    Main program.
    '''

    # File path containing data files
    # (YOU CAN EDIT THIS)
    path_to_dir = '/path/to/my/netcdf/files/'

    # List of search strings that each file must contain
    # (YOU CAN EDIT THIS)
    substrs = ['SpeciesConc']

    # Get a list of variables that GCPy should not read.
    # These are mostly variables introduced into GCHP with the MAPL v1.0.0
    # update.  These variables contain either repeated or non-standard
    # dimensions that can cause problems in xarray when combining datasets.
    skip_vars = gcon.skip_these_vars
    
    # Look for all the netCDF files in the path
    file_list = find_files_in_dir(path_to_dir, substrs)

    # Return a single xarray Dataset containing data from all files
    # NOTE: Need to add combine="nested" for xarray 0.15 and higher
    v = xr.__version__.split(".")
    if int(v[0]) == 0 and int(v[1]) >= 15: 
        ds = xr.open_mfdataset(file_list,
                                 drop_variables=skip_vars,
                                 combine="nested")
    else:
        ds = xr.open_mfdataset(file_list,
                                 drop_variables=skip_vars)

    # Replace NaN values with zeroes
    ds = replace_nans_with_zeroes(ds, verbose=True)

    # Specify the path and filename for the concatenated data
    # (YOU CAN EDIT THIS)
    outdir = '/path/to/my/output/file'
    outfile = os.path.join(outdir, 'my_concatenated_output_file.nc')

    # Write concatenated data to a netCDF file
    ds.to_netcdf(outfile)

if __name__ == "__main__":
    main()
