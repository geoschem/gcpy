#!/usr/bin/env python
"""
Example script that can compare diagnostics from two different netCDF
collections.  Similar to compute_diagnostics.ipynb, but can be used
without having to open a Jupyter notebook.
"""

# Imports
import os
from os.path import join
import warnings
import numpy as np
import xarray as xr
import gcpy.benchmark as bmk
import gcpy.constants as gcon
import gcpy.util as util

# Tell matplotlib not to look for an X-window
os.environ["QT_QPA_PLATFORM"] = "offscreen"

# Suppress harmless run-time warnings (mostly about underflow in division)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)

########################################################################
###           CONFIGURABLE SETTINGS: FILES AND FOLDERS               ###
########################################################################

# Main data directory
maindir = "/n/holyscratch01/external_repos/GEOS-CHEM/gcgrid/geos-chem/validation/gcpy_test_data/1mon/"
refstr = "GCC_ref"
devstr = "GCC_dev"

# Regridding weights directory
weightsdir = "/n/holyscratch01/external_repos/GEOS-CHEM/gcgrid/gcdata/ExtData/GCHP/RegriddingWeights"

# Directory where PDF files will be sent
plotsdir = join(".", "Plots")

# Ref version
refdir = join(maindir, refstr, "OutputDir")
reffile = join(refdir, "GEOSChem.SpeciesConc.20160701_0000z.nc4")

# Dev version
devdir = join(maindir, devstr, "OutputDir")
devfile = join(devdir, "GEOSChem.SpeciesConc.20160701_0000z.nc4")

# PDF names
pdfname_level = join(plotsdir, "single_level_comparison.pdf")
pdfname_zonal = join(plotsdir, "zonal_mean_level_comparison.pdf")

########################################################################
###           CONFIGURABLE SETTINGS: PLOTTING OPTIONS                ###
########################################################################

# Plot options
create_single_level_plot = True
create_zonal_mean_plot = True
print_totals_and_diffs = True

# Specify the level that you wish to plot (starting from 0)
# NOTE: For single level plots only
level_to_plot = 0

# Turn on extra debug warnings
verbose = False

# Only process variables containing this substring,
# or set to None to process all common variables in Ref & Dev.
restrict_vars = None

# ======================================================================
# Methods
# ======================================================================

def compare_data(refdata, devdata):
    """
    Compares data from two different xarray datasets.

    Args:
    -----
        refdata : xarray Dataset
            The first Dataset to be compared.

        devdata : xarray Dataset
            The Dataset to be compared against refdata.
    """

    # For each variable in refdata, but not non devdata, add an
    # array of NaN values to refdata. Ditto for devdata.  This will
    # allow us to show that the variable is missing in either
    # refdata or devdata.
    [refdata, devdata] = util.add_missing_variables(refdata, devdata)

    # Get the list of common variable names
    quiet = not verbose
    vardict = util.compare_varnames(refdata, devdata, quiet=quiet)
    varlist_level = vardict["commonvars2D"] + vardict["commonvars3D"]
    varlist_zonal = vardict["commonvars3D"]

    # Restrict variables to those containing a given substring
    if restrict_vars is not None:
        varlist_level = [v for v in varlist_level if restrict_vars in v]
        varlist_zonal = [v for v in varlist_zonal if restrict_vars in v]

    # ==================================================================
    # Generate the single level comparison plot
    # ==================================================================
    if create_single_level_plot:
        bmk.compare_single_level(
            refdata,
            refstr,
            devdata,
            devstr,
            ilev=level_to_plot,
            varlist=varlist_level,
            pdfname=pdfname_level,
            weightsdir=weightsdir,
            verbose=verbose
        )

    # ==================================================================
    # Generate the zonal mean comparison plot
    # ==================================================================
    if create_zonal_mean_plot:
        bmk.compare_zonal_mean(
            refdata,
            refstr,
            devdata,
            devstr,
            varlist=varlist_zonal,
            pdfname=pdfname_zonal,
            weightsdir=weightsdir,
            verbose=verbose
        )

    # ==================================================================
    # Print totals for each quantity
    # ==================================================================
    if print_totals_and_diffs:

        # Header
        print("{} Ref={} Dev={} {}".format(
            "Variable".ljust(22),
            refstr.ljust(20),
            devstr.ljust(20),
            "Dev-Ref"
        ))

        # Data
        for v in varlist_level:
            refsum = np.sum(refdata[v].values)
            devsum = np.sum(devdata[v].values)
            diff = devsum - refsum
            print("{} : {} | {} | {} ".format(
                v.ljust(20),
                str(refsum).ljust(22),
                str(devsum).ljust(22),
                diff
            ))

def main():
    """
    Main program, reads data and calls compare_data to make plots.
    """

    # Create directories for plots and weights if they do not exist
    if create_single_level_plot or create_zonal_mean_plot:
        if not os.path.isdir(plotsdir):
            os.mkdir(plotsdir)
        if not os.path.isdir(weightsdir):
            os.mkdir(weightsdir)

    # Read the Ref abd Dev data into xarray Dataset objects
    skip_vars = gcon.skip_these_vars
    refdata = xr.open_dataset(reffile, drop_variables=skip_vars)
    devdata = xr.open_dataset(devfile, drop_variables=skip_vars)

    # Create the comparison plots and sums
    # NOTE: all other variables are global and thus
    # are seen by all functions above
    compare_data(refdata, devdata)


if __name__ == "__main__":
    main()
