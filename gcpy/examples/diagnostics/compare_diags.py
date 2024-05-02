#!/usr/bin/env python
"""
Example script that can compare diagnostics from two different netCDF
collections.  Similar to compute_diagnostics.ipynb, but can be used
without having to open a Jupyter notebook.  The parameters for the
configuration are specified in a YAML file whose name is passed
as an argument.
"""
import os
import sys
import warnings
import numpy as np
from gcpy.util import add_missing_variables, compare_varnames, \
    dataset_reader, read_config_file, rename_and_flip_gchp_rst_vars
from gcpy.constants import skip_these_vars
from gcpy.plot.compare_single_level import compare_single_level
from gcpy.plot.compare_zonal_mean import compare_zonal_mean

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

    # Root data path
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

    # Function to read the data
    reader = dataset_reader(
        multi_files=False,
        verbose=False
    )

    # Read Ref data
    try:
        refdata = reader(
            ref_file,
            drop_variables=skip_these_vars
        ).load()
    except FileNotFoundError as exc:
        msg = "Error reading " + ref_file
        raise FileNotFoundError(msg) from exc

    # Read Dev data
    try:
        devdata = reader(
            dev_file,
            drop_variables=skip_these_vars
        ).load()
    except FileNotFoundError as exc:
        msg = "Error reading " + dev_file
        raise FileNotFoundError(msg) from exc

    # If the data is from a GCHP restart file, rename variables and
    # flip levels to match the GEOS-Chem Classic naming and level
    # conventions.  Otherwise no changes will be made.
    refdata = rename_and_flip_gchp_rst_vars(refdata)
    devdata = rename_and_flip_gchp_rst_vars(devdata)

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

    # Determine if we will print percent or absolute differences
    do_percent_diff = False
    if any(x in config["options"]["totals_and_diffs"]["diff_type"] \
           for x in ["percent", "pctdiff", "%"]):
        do_percent_diff = True

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
        ofile = open(pathname, 'w', encoding="UTF-8")


    # Percent threshold for reporting differences
    threshold = config["options"]["totals_and_diffs"]["small_diff_threshold"]

    # If we are skipping variables with small differences,
    # then alert the user (print to screen & file)
    if config["options"]["totals_and_diffs"]["skip_small_diffs"]:
        diff_label = f"|absolute difference| > {threshold}"
        if do_percent_diff:
            diff_label = f"|percent difference| > {threshold} %"
        line = f"... Only showing variables with {diff_label}"
        if do_file:
            print(line, file=ofile)
        else:
            print(line)

    # Print top header
    diff_label = "Dev - Ref"
    if do_percent_diff:
        diff_label = "(Dev-Ref)/Ref % diff"
    line = "{} Ref={} Dev={} {}".format(
        "Variable".ljust(22),
        config["data"]["ref"]["label"].ljust(20),
        config["data"]["dev"]["label"].ljust(20),
        diff_label
    )
    if do_file:
        print(line, file=ofile)
    if do_screen:
        print(line)

    # Always print nonzero differences, but only print zero differences
    # if the configuration option "skip_zero_diffs" is False.
    for var in varlist:

        # Absolute difference
        refsum = np.sum(refdata[var].values)
        devsum = np.sum(devdata[var].values)

        # Absolute difference
        absdiff = np.sum(devdata[var].values - refdata[var].values)

        # Compute percent difference if needed
        # otherwise we'll use the absolute difference for plotting
        diff = absdiff
        if do_percent_diff:
            diff = 0.0
            if np.abs(refsum) > 0.0:
                diff = (absdiff / refsum) * 100.0

        # Line to be printed
        line = "{} : {} | {} | {} ".format(
            var.ljust(20),
            str(refsum).ljust(22),
            str(devsum).ljust(22),
            diff
        )

        # Skip small values
        if np.abs(diff) > threshold:
            if do_file:
                print(line, file=ofile)
            if do_screen:
                print(line)

        # Or don't skip any values
        else:
            if not config["options"]["totals_and_diffs"]["skip_small_diffs"]:
                if do_file:
                    print(line, file=ofile)
                if do_screen:
                    print(line)

    # Close file
    if do_file:
        ofile.close()


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
    [refdata, devdata] = add_missing_variables(refdata, devdata)

    # Get the list of common variable names
    verbose = config["options"]["verbose"]
    quiet = not verbose
    vardict = compare_varnames(refdata, devdata, quiet=quiet)
    varlist_level = vardict["commonvars2D"] + vardict["commonvars3D"]
    varlist_zonal = vardict["commonvars3D"]

    # Restrict variables to those containing a given substring
    restrict_vars = config["options"]["restrict_vars"]
    if len(restrict_vars) > 0:
        varlist_level = [v for v in varlist_level if v in restrict_vars]
        varlist_zonal = [v for v in varlist_zonal if v in restrict_vars]

    # Determine if we need to flip levels in the vertical
    flip_ref = False
    flip_dev = False
    if "flip_levels" in config["data"]["ref"]:
        flip_ref = config["data"]["ref"]["flip_levels"]
    if "flip_levels" in config["data"]["dev"]:
        flip_dev = config["data"]["dev"]["flip_levels"]        
        
    # ==================================================================
    # Generate the single level comparison plot
    # ==================================================================
    if config["options"]["level_plot"]["create_plot"]:
        print('... Creating single-level plots')
        pdfname = os.path.join(
            config["paths"]["plots_dir"],
            config["options"]["level_plot"]["pdfname"]
        )
        compare_single_level(
            refdata,
            config["data"]["ref"]["label"],
            devdata,
            config["data"]["dev"]["label"],
            flip_ref=flip_ref,
            flip_dev=flip_dev,
            ilev=config["options"]["level_plot"]["level_to_plot"],
            varlist=varlist_level,
            pdfname=pdfname,
            weightsdir=config["paths"]["weights_dir"],
            n_job=config["options"]["n_cores"],
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
        compare_zonal_mean(
            refdata,
            config["data"]["ref"]["label"],
            devdata,
            config["data"]["dev"]["label"],
            flip_ref=flip_ref,
            flip_dev=flip_dev,
            varlist=varlist_zonal,
            pdfname=pdfname,
            weightsdir=config["paths"]["weights_dir"],
            n_job=config["options"]["n_cores"],
            verbose=verbose
        )

    # ==================================================================
    # Print totals for each quantity
    # ==================================================================
    if config["options"]["totals_and_diffs"]["create_table"]:
        print('... Printing totals and differences')
        print_totals_and_diffs(
            config,
            refdata,
            devdata,
            varlist_level
        )


def main(argv):
    """
    Main program, reads data and calls compare_data to make plots.
    """

    # Take the config file as the 2nd argument (or use a default)
    # NOTE: sys.argv[0] is always the program name!
    if len(argv) == 2:
        config_file = argv[1]
    else:
        config_file = "compare_diags.yml"

    # Get paths and options from the configuration file
    config = read_config_file(config_file)

    # Create dirs for plots & weights (if necessary)
    create_dirs(config)

    # Read and compare data
    compare_data(config, read_data(config))


# Only execute when we run as a standalone script
if __name__ == "__main__":
    main(sys.argv)
