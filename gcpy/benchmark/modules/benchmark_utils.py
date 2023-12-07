"""
Specific utilities for creating plots from GEOS-Chem benchmark simulations.

TODO: Migrate benchmark-specific utilities from gcpy/benchmark.py to here.
"""
import os
import gc
import numpy as np
from gcpy import util
from gcpy.constants import skip_these_vars
from gcpy.plot.compare_single_level import compare_single_level
from gcpy.plot.compare_zonal_mean import compare_zonal_mean

# Suppress numpy divide by zero warnings to prevent output spam
np.seterr(divide="ignore", invalid="ignore")


def make_collection_subdir(
        dst,
        collection,
        datestr
):
    """
    Creates a subdirectory for the given collection type
    in the destination directory.

    Args:
    ----
    dst        (str) : Destination directory
    collection (str) : Collection name
    datestr    (str) : Date string

    Returns:
    --------
    target_dst (str) : Directory that was created

    """

    # Make a collection subdirectory
    target_dst = os.path.join(dst, collection)
    if not os.path.isdir(target_dst):
        os.mkdir(target_dst)

    # If datestr is passed, create a further subdirectory
    if datestr is not None:
        target_dst = os.path.join(target_dst, datestr)
        if not os.path.isdir(target_dst):
            os.mkdir(target_dst)

    return target_dst


def read_ref_and_dev(
        ref,
        dev,
        time_mean=False,
        verbose=False
):
    """
    Reads files from the Ref and Dev models into xarray.Dataset objects.


    Args:
    -----
    ref (str or list) : Ref data file(s)
    def (str or list) : Ref data file(s)
    time_mean (bool)  : Return the average over the time dimension

    Returns:
    --------
    refds (xarray.Dataset) : Data from the Ref model
    devds (xarray.Dataset) : Data from the Dev model
    """

    # Get the function that will read file(s) into a dataset
    reader = util.dataset_reader(time_mean, verbose=verbose)

    # Open datasets
    refds = reader(ref, drop_variables=skip_these_vars)
    devds = reader(dev, drop_variables=skip_these_vars)

    # Take
    if time_mean:
        refds = util.dataset_mean(refds)
        devds = util.dataset_mean(devds)

    return refds, devds


def get_common_varnames(
        refds,
        devds,
        var_prefix=None,
        verbose=False,
    """
    """
    
    # Get list of variables in collection
    vardict = util.compare_varnames(refds, devds, quiet=not verbose)
    varlist = [v for v in vardict["commonvars3D"] if collection + "_" in v]

    # Select variables having a common prefix
    if var_prefix is not None:
        varlist = [v for v in varlist if var_prefix in v]

    return varlist.sort()
