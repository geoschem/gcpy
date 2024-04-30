#!/usr/bin/env python3
"""
"""
import os
from gcpy.util import verify_variable_type
import json


def read_gcclassic(ifile):
    """
    Determines if the input is a valid JSON.

    Args
    ifile  : str  : file name

    Returns
    result : dict : Dictionary with timing information
    """

    # Make sure file exists
    if not os.path.exists(ifile):
        raise FileNotFoundError(f"Could not find {ifile}!")

    # First try to read the file as a JSON,
    # then try to read the file as text.
    try:
        result = read_gcclassic_json(ifile)
    except ValueError as err:
        result = read_gcclassic_log(ifile)
    return result


def read_gcclassic_json(
        ifile
):
    """
    Parses the GEOS-Chem Classic timing information in JSON format
    and returns a dictionary with the results.

    Args
    ifile  : str  : File name

    Returns
    result : dict : Dictionary with timing information
    """
    try:
        with open(ifile, encoding="utf-8") as json_file:
            result = json.load(json_file)
            return result["GEOS-Chem Classic timers"]
    except ValueError as err:
        raise ValueError from err


def read_gcclassic_log(ifile):
    """
    Parses the GEOS-Chem Classic log file with timing information
    and returns a dictionary with the results.

    Args
    ifile  : str  : File name

    Returns
    result : dict : Dictionary with timing information
    """
    keep_line = False
    timers = {}

    # Read the line backwards and get just keep the timing information
    with open(ifile, encoding="utf-8") as log_file:

        for line in reversed(list(log_file)):
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
                val = substr[3].split()[1].strip()
                timers[key] = {"seconds": val}

    return timers


def print_timer(key, ref, dev, ofile):
    """
    Prints timing info for a single timer to a log file.
    """
    line = f"{key:<25}  {ref[key]['seconds']:>20}  {dev[key]['seconds']:>20}"
    print(line, file=ofile)


def display_timers(ref, ref_label, dev, dev_label, table_file):
    """
    Prints the GEOS-Chem timer information to a table.

    Args
    ref : dict : Timer output from the "Ref" model
    ref : dict : Timer output from the "Dev" model
    """
    with open(table_file, "w", encoding="utf-8") as ofile:

        # Print header
        print(f"{'Timer':<25}  {ref_label:>20}  {dev_label:>20}", file=ofile)
        print(f"{'-'*25:<25}  {'-'*20:>20}  {'-'*20:>20}", file=ofile)
        
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

        
def make_benchmark_timing_table(
        ref_file,
        ref_label,
        dev_file,
        dev_label,
        dst,
):
    """
    """
    verify_variable_type(ref_file, (str, list))
    verify_variable_type(ref_label, str)
    verify_variable_type(dev_file, (str, list))
    verify_variable_type(dev_label, str)
    verify_variable_type(dst, str)

    # Strip timing info from JSON or log ifle
    ref_timers = read_gcclassic(ref_file)
    dev_timers = read_gcclassic(dev_file)

    # Write timing info to a table
    display_timers(
        ref_timers,
        ref_label,
        dev_timers,
        dev_label,
        "sample_output.txt",
    )


if __name__ == '__main__':
    make_benchmark_timing_table(
       "./gcclassic_timers.json",
        "GCC 14.4.0",
        "./execute.gc_4x5_merra2_fullchem_benchmark.log",
        "GCHP 14.4.0",
        "./"
    )
#        "./execute.gchp_merra2_fullchem_benchmark.log",
