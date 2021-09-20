#!/usr/bin/env python
"""
Example script that can compare diagnostics from two different netCDF
collections.  Similar to compute_diagnostics.ipynb, but can be used
without having to open a Jupyter notebook.
"""

# Imports
import os
import sys
import warnings
import numpy as np
import xarray as xr
from yaml import load as yaml_load_file
import gcpy.benchmark as bmk
import gcpy.constants as constants
import gcpy.util as util

# Tell matplotlib not to look for an X-window
os.environ["QT_QPA_PLATFORM"] = "offscreen"

# Suppress harmless run-time warnings (mostly about underflow in division)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# ======================================================================
# Methods
# ======================================================================


def create_dirs(config):
    """
    Create directories for plots and weights if they do not exist.

    Arguments:
        config : dict
            Configuration information read from a YAML file
    """

    # Extract fields from the config object
    do_level_plot = config["options"]["level_plot"]["create_plot"]
    do_zonal_mean = config["options"]["zonal_mean"]["create_plot"]
    plotsdir = config["paths"]["plots_dir"]
    weightsdir = config["paths"]["weights_dir"]

    # Create plotsdir and weightsdir if necessary
    if do_level_plot or do_zonal_mean:
        if not os.path.isdir(plotsdir):
            os.mkdir(plotsdir)
            if not os.path.isdir(weightsdir):
                os.mkdir(weightsdir)


def read_data(config):
    """
    Read data from the Ref and Dev datasets

    Arguments:
        config : dict
            Configuration information read from a YAML file

    Returns:
        data : dict
            Contains Ref and Dev data as xarray Dataset fields.
    """

    # If we are using the gcpy_test data, use the gcpy_test_dir
    if config["options"]["gcpy_test"]:
        rootdir = config["paths"]["test_data_dir"]
    else:
        rootdir = config["paths"]["main_dir"]

    # Define paths to Ref & Dev files
    ref_file = os.path.join(
        rootdir,
        config["data"]["ref"]["dir"],
        config["data"]["ref"]["subdir"],
        config["data"]["ref"]["file"]
    )
    dev_file = os.path.join(
        rootdir,
        config["data"]["dev"]["dir"],
        config["data"]["dev"]["subdir"],
        config["data"]["dev"]["file"]
    )

    # Read Ref data
    try:
        refdata = xr.open_dataset(
            ref_file,
            drop_variables=constants.skip_these_vars
        )
    except Exception:
        msg = "Error reading " + ref_file
        raise Exception(msg)

    # Read Dev data
    try:
        devdata = xr.open_dataset(
            dev_file,
            drop_variables=constants.skip_these_vars
        )
    except Exception:
        msg = "Error reading " + dev_file
        raise Exception(msg)

    # Define dictionary for return
    data = {
        "ref": refdata,
        "dev": devdata
    }

    return data


def print_totals_and_diffs(config, refdata, devdata, varlist):
    """
    Prints totals and differences of data.  This is a quick way
    to see if there are nonzero differences between Ref and Dev
    datasets.

    Arguments:
        config : dict
            Configuration information read from a YAML file
        refdata : xarray Dataset
            Contains data from the Ref model run
        devdata : xarray Dataset
            Contains data from the Dev model run
        varlist : list of str
            Contains a list of data variables.
    """

    # Get quantities from the config object
    filename = config["options"]["totals_and_diffs"]["filename"]
    do_screen = config["options"]["totals_and_diffs"]["print_to_screen"]

    # Determine if we will print to a file
    do_file = len(filename) > 0
    if not do_file and not do_screen:
        print("...Printing to screen and file are both disabled!")
        return

    # Open a text file for output if a filename is specified
    if do_file:
        pathname = os.path.join(
            config["paths"]["plots_dir"],
            filename
        )
        f = open(pathname, 'w')

    # If we are printing only variables with zero diffs,
    # then alert the user (print to screen & file)
    if config["options"]["totals_and_diffs"]["skip_zero_diffs"]:
        line = "... Only showing variables with nonzero differences"
        if do_file:
            print(line, file=f)
        else:
            print(line)

    # Print top header
    line = "{} Ref={} Dev={} {}".format(
        "Variable".ljust(22),
        config["data"]["ref"]["label"].ljust(20),
        config["data"]["dev"]["label"].ljust(20),
        "Dev-Ref"
    )
    if do_file:
        print(line, file=f)
    if do_screen:
        print(line)

    # Always print nonzero differences, but only print zero differences
    # if the configuration option "skip_zero_diffs" is False.
    for v in varlist:
        refsum = np.sum(refdata[v].values)
        devsum = np.sum(devdata[v].values)
        diff = devsum - refsum

        if np.abs(diff) > 0.0:
            line = "{} : {} | {} | {} ".format(
                v.ljust(20),
                str(refsum).ljust(22),
                str(devsum).ljust(22),
                diff
            )
            if do_file:
                print(line, file=f)
            if do_screen:
                print(line)
        else:
            if not config["options"]["totals_and_diffs"]["skip_zero_diffs"]:
                line = "{} : {} | {} | {} ".format(
                    v.ljust(20),
                    str(refsum).ljust(22),
                    str(devsum).ljust(22),
                    diff
                )
                if do_file:
                    print(line, file=f)
                if do_screen:
                    print(line)

    # Close file
    if do_file:
        f.close()


def compare_data(config, data):
    """
    Compares data from two different xarray datasets.

    Args:
        data : dict
            Contains Ref and Dev data as xarray Dataset fields.
    """

    # Get xarray datasets from the data object
    refdata = data["ref"]
    devdata = data["dev"]

    # For each variable in refdata, but not non devdata, add an
    # array of NaN values to refdata. Ditto for devdata.  This will
    # allow us to show that the variable is missing in either
    # refdata or devdata.
    [refdata, devdata] = util.add_missing_variables(refdata, devdata)

    # Get the list of common variable names
    verbose = config["options"]["verbose"]
    quiet = not verbose
    vardict = util.compare_varnames(refdata, devdata, quiet=quiet)
    varlist_level = vardict["commonvars2D"] + vardict["commonvars3D"]
    varlist_zonal = vardict["commonvars3D"]

    # Restrict variables to those containing a given substring
    restrict_vars = config["options"]["restrict_vars"]
    if len(restrict_vars) > 0:
        varlist_level = [v for v in varlist_level if v in restrict_vars]
        varlist_zonal = [v for v in varlist_zonal if v in restrict_vars]

    # ==================================================================
    # Generate the single level comparison plot
    # ==================================================================
    if config["options"]["level_plot"]["create_plot"]:
        print('... Creating single-level plots')
        pdfname = os.path.join(
            config["paths"]["plots_dir"],
            config["options"]["level_plot"]["pdfname"]
        )
        bmk.compare_single_level(
            refdata,
            config["data"]["ref"]["label"],
            devdata,
            config["data"]["dev"]["label"],
            ilev=config["options"]["level_plot"]["level_to_plot"],
            varlist=varlist_level,
            pdfname=pdfname,
            weightsdir=config["paths"]["weights_dir"],
            verbose=verbose
        )

    # ==================================================================
    # Generate the zonal mean comparison plot
    # ==================================================================
    if config["options"]["zonal_mean"]["create_plot"]:
        print('... Creating zonal mean plots')
        pdfname = os.path.join(
            config["paths"]["plots_dir"],
            config["options"]["zonal_mean"]["pdfname"]
        )
        bmk.compare_zonal_mean(
            refdata,
            config["data"]["ref"]["label"],
            devdata,
            config["data"]["dev"]["label"],
            varlist=varlist_zonal,
            pdfname=pdfname,
            weightsdir=config["paths"]["weights_dir"],
            verbose=verbose
        )

    # ==================================================================
    # Print totals for each quantity
    # ==================================================================
    if config["options"]["totals_and_diffs"]["create_table"] or \
       config["options"]["totals_and_diffs"]["print_to_screen"]:
        print('... Printing totals and differences')
        print_totals_and_diffs(
            config,
            refdata,
            devdata,
            varlist_level
        )


def main():
    """
    Main program, reads data and calls compare_data to make plots.
    """
    
    # Take the config file as the 2nd argument (or use a default)
    # NOTE: sys.argv[0] is always the program name!
    config_filename = sys.argv[1] if len(sys.argv) == 2 else "compare_diags.yml"

    # Get paths and options from the configuration file
    config = util.read_config_file(config_filename)

    # Create dirs for plots & weights (if necessary)
    create_dirs(config)

    # Read and compare data
    compare_data(config, read_data(config))


if __name__ == "__main__":
    main()
