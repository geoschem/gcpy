#!/usr/bin/env python3
"""
Scrapes GEOS-Chem Classic benchmark timing information from one or
more JSON or text files.
"""
import os
import json
import numpy as np
from gcpy.util import make_directory, replace_whitespace, verify_variable_type


def read_gcclassic(input_files):
    """
    Determines whether we should call a function to parse the given
    input file(s) as JSON or plain text.

    Args
    input_files : str|list     : File or list of files to parse

    Returns
    result      : list of dict : List of dicts with timing info
    """
    try:
        result = read_timing_data(input_files, read_one_json_file)
    except ValueError:
        result = read_timing_data(input_files, read_one_text_file)
    return result


def read_timing_data(input_files, reader):
    """
    Parses the GEOS-Chem Classic timing information in JSON format
    and returns a dictionary with the results.

    Args
    input files : str|list     : JSON or text file(s) to parse
    reader      : function     : Function that will parse the file(s)

    Returns
    timing      : list of dict : Dictionary with timing information
    """
    # Return value
    timing = []

    # If more than one file has been provided, read the timing
    # information and return a list of dictionaries with results
    if isinstance(input_files, list):
        for input_file in input_files:
            result = reader(input_file)
            timing.append(result)
        return timing

    # If only one file has been provided, then read it
    # and return the dictionary in a list
    if isinstance(input_files, str):
        result = reader(input_files)
        timing.append(result)
        return timing

    raise ValueError("Argument 'input_files' is not of type str or list!")


def read_one_json_file(json_file):
    """
    Parses a GEOS-Chem JSON file with timing information
    and returns a dictionary with the results.

    Args
    json_file : str  : JSON file with timing information

    Returns
    result    : dict : Dictionary with timing information
    """

    # Make sure file exists
    if not os.path.exists(json_file):
        raise FileNotFoundError(f"Could not find {json_file}!")

    # If the file is not a JSON file, raise a ValueError, as
    # this will prompt read_gcclassic to parse the file as text.
    try:
        with open(json_file, encoding="utf-8") as ifile:
            result = json.load(ifile)
            return result["GEOS-Chem Classic timers"]
    except ValueError as err:
        raise ValueError from err


def read_one_text_file(text_file):
    """
    Parses the GEOS-Chem Classic log file (plain text) with
    timing information and returns a dictionary with the results.

    Args
    text_file : str  : Text file with timing information

    Returns
    result    : dict : Dictionary with timing information
    """
    keep_line = False
    timers = {}

   # Make sure file exists
    if not os.path.exists(text_file):
        raise FileNotFoundError(f"Could not find {text_file}!")

    # Read the line backwards and get just keep the timing information
    with open(text_file, encoding="utf-8") as ifile:

        for line in reversed(list(ifile)):
            line = line.strip("\n")

            # Set a flag to denote the start & end of timing info
            if "Unit conversions" in line:
                keep_line = True
            if "----------" in line:
                keep_line = False
                break

            # Append timing info lines into a list
            if keep_line:
                substr = line.split(":")
                key = substr[0].strip()
                if "THE TIMER DID NOT RUN" in line:
                    val = np.nan
                else:
                    val = substr[3].split()[1].strip()
                timers[key] = {"seconds": val}

    return timers


def sum_timers(timers):
    """
    Sums the time in seconds for each GEOS-Chem timer.  Input may be
    a single dict with timing information or a list of dicts.

    Args
    timers : dict|list : GEOS-Chem timing information from one or more
                         JSON or log files.

    Returns
    result : dict      : Sum of timing information
    """

    # If timers is of type dict, no summing is needed.
    if isinstance(timers, dict):
        return timers

    # If timers is a list of dicts, sum the times
    # in seconds into a new dict, and then return.
    if isinstance(timers, list):

        # Initialize the result dict
        result = {}
        for timer in timers:
            for (key, val) in timer.items():
                result[key] = 0.0

        # Then sum the time in seconds for each timer
        for timer in timers:
            for (key, val) in timer.items():
                result[key] += float(val["seconds"])

        return result

    raise ValueError("Argument 'timers' must be of type str or dict!")


def print_timer(key, ref, dev, ofile):
    """
    Prints timing info for a single timer to a log file.

    Args
    key   : str  : Dictionary key to print
    ref   : dict : Timing information from the "Ref" model
    dev   : dict : Timing information from the "Dev" model
    ofile : file : File object where info will be written
    """
    pctdiff = np.nan
    if np.abs(ref[key] > 0.0):
        pctdiff = ((dev[key] - ref[key]) / ref[key]) * 100.0
    line = f"{key:<22}  {ref[key]:>18.3f}  {dev[key]:>18.3f}   {pctdiff:>12.3f}"
    if np.abs(pctdiff) >= 10.0:  # Flag diffs > +/- 10%
        line += " *"
    print(line, file=ofile)


def display_timers(ref, ref_label, dev, dev_label, table_file):
    """
    Prints the GEOS-Chem timer information to a table.

    Args
    ref        : dict : Timing information from the "Ref" model
    ref_label  : str  : Version string for the "Ref" model
    dev        : dict : Timing information from the "Dev" model
    dev_label  : str  : Version string for the "Dev" model
    table_file : str  : File name for the timing table output
    """
    with open(table_file, "w", encoding="utf-8") as ofile:

        # Print header
        print("%"*79, file=ofile)
        print("%%% GEOS-Chem Classic Benchmark Timing Information",
              file=ofile)
        print("%%%", file=ofile)
        print(f"%%% Ref = {ref_label}", file=ofile)
        print(f"%%% Dev = {dev_label}", file=ofile)
        print("%"*79, file=ofile)
        print("\n", file=ofile)
        print(f"{'Timer':<22}  {'Ref [s]':>18}  {'Dev [s]':>18}   {'% Diff':>12}", file=ofile)
        print("-"*79, file=ofile)

        # Print timers
        print_timer("GEOS-Chem",             ref, dev, ofile)
        print_timer("HEMCO",                 ref, dev, ofile)
        print_timer("All chemistry",         ref, dev, ofile)
        print_timer("=> Gas-phase chem",     ref, dev, ofile)
        print_timer("=> Photolysis",         ref, dev, ofile)
        print_timer("=> Aerosol chem",       ref, dev, ofile)
        print_timer("=> Linearized chem",    ref, dev, ofile)
        print_timer("Transport",             ref, dev, ofile)
        print_timer("Convection",            ref, dev, ofile)
        print_timer("Boundary layer mixing", ref, dev, ofile)
        print_timer("Dry deposition",        ref, dev, ofile)
        print_timer("Wet deposition",        ref, dev, ofile)
        print_timer("Diagnostics",           ref, dev, ofile)
        print_timer("Unit conversions",      ref, dev, ofile)


def make_benchmark_gcclassic_timing_table(
        ref_files,
        ref_label,
        dev_files,
        dev_label,
        dst="./benchmark",
        overwrite=False,
):
    """
    Creates a table of timing information for GEOS-Chem Classic
    benchmark simulations given one or more JSON and/or text files
    as input.

    Args
    ref_files : str|list : File(s) with timing info from the "Ref" model
    ref_label : str      : Version string for the "Ref" model
    dev_files : str|list : File(s) with timing info from the "Ref" model
    dev_label : str      : Version string for the "Dev" model

    Kwargs
    dst       : str      : Directory where output will be written
    overwrite : bool     : Overwrite existing files? (default: False)

    """
    verify_variable_type(ref_files, (str, list))
    verify_variable_type(ref_label, str)
    verify_variable_type(dev_files, (str, list))
    verify_variable_type(dev_label, str)
    verify_variable_type(dst, str)

    # Create the destination folder
    make_directory(dst, overwrite)

    # Replace whitespace in the ref and dev labels
    ref_label = replace_whitespace(ref_label)
    dev_label = replace_whitespace(dev_label)

    # Strip timing info from JSON/text file(s) and sum the them.
    ref_timers = sum_timers(read_gcclassic(ref_files))
    dev_timers = sum_timers(read_gcclassic(dev_files))

    # Filename for output
    timing_table = replace_whitespace(
        os.path.join(
            dst,
            f"Benchmark_Timers_{ref_label}_vs_{dev_label}.txt"
        )
    )

    # Write timing info to a table
    display_timers(
        ref_timers,
        replace_whitespace(ref_label),
        dev_timers,
        replace_whitespace(dev_label),
        timing_table,
    )
