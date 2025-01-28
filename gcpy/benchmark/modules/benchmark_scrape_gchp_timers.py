#!/usr/bin/env python3
"""
Scrapes GCHP Classic benchmark timing information from one or
more text files.
"""
import os
import numpy as np
from gcpy.util import make_directory, replace_whitespace, verify_variable_type


def read_timing_data(input_files):
    """
    Parses the GEOS-Chem Classic timing information in JSON format
    and returns a dictionary with the results.

    Args
    input files : str|list     : Text file(s) to parse

    Returns
    timing      : list of dict : Dictionary with timing information
    """
    # Return value
    timing = []

    # If more than one file has been provided, read the timing
    # information and return a list of dictionaries with results
    if isinstance(input_files, list):
        for input_file in input_files:
            result = read_one_text_file(input_file)
            timing.append(result)
        return timing

    # If only one file has been provided, then read it
    # and return the dictionary in a list
    if isinstance(input_files, str):
        result = read_one_text_file(input_files)
        timing.append(result)
        return timing

    raise ValueError("Argument 'input_files' is not of type str or list!")


def count_characters(text, char_to_match="."):
    """
    Returns the number of characters in a string of text.

    Args
    text          : str : The text to parse

    Kwargs
    char_to_match : str : The character to look for in "text"

    Returns
    result        : int : Number of underscores in "text"

    Reference
    https://stackoverflow.com/questions/991350/counting-repeated-characters-in-a-string-in-python
    """
    # Create a dictionary where each character of "text"
    # is a key, and all values are set to zero.
    count = dict.fromkeys(text, 0)

    # Increment each time a character is found
    for char in text:
        count[char] += 1

    # Return the count of underscores
    if char_to_match not in count:
        return 0
    return count[char_to_match]


def read_one_text_file(text_file):
    """
    Parses the GCHP log file (plain text) with timing information
    and returns a dictionary with the results.

    Args
    text_file : str  : Text file with timing information

    Returns
    result    : dict : Dictionary with timing information
    """

    # Make sure file exists
    if not os.path.exists(text_file):
        raise FileNotFoundError(f"Could not find {text_file}!")

    # ==================================================================
    # Parse the GCHP log file
    # ==================================================================

    # Initialize local variables
    keep_line = False
    temp_timers = []
    inclusive = 0
    temp_timers = []

    # Open the log file
    with open(text_file, encoding="utf-8") as ifile:

        # Read each line in the file
        for line in ifile:

            # Strip newlines; skip empty lines
            line = line.strip()
            if len(line) == 0:
                continue

            # GCHP timers section (also skip header lines)
            if 'Times for component <GCHPchem>' in line:
                keep_line = True
                inclusive = 3
                continue
            if keep_line and 'Min                            Mean' in line:
                continue
            if keep_line and '============================' in line:
                continue
            if keep_line and 'Name                          %' in line:
                continue
            if keep_line and '------ ---------- ----------' in line:
                continue
            if keep_line and '---------------------------------' in line:
                keep_line = False
                continue

            # Summary section (also skip header lines)
            if 'Report on process:  0' in line:
                keep_line = True
                inclusive = 2
                continue
            if keep_line and 'Inclusive' in line:
                continue
            if keep_line and '================' in line:
                continue
            if keep_line and 'Name' in line:
                continue
            if keep_line and '-------- --------- ------ --------- ------' \
               in line:
                continue

            # NOTE: This line only appears in cloud benchmarks,
            # which signals the end of GCHP output and the start of
            # job statistics.  Exit when we encounter this.
            if keep_line and "Command being timed:" in line:
                break

            # Append timing info lines into a list of dicts
            if keep_line:
                substr = line.split()
                key = substr[0].strip()
                val = float(substr[inclusive].strip())
                temp_timers.append({key: val})

    # ==================================================================
    # Save timing results into a "flattened" dictionary
    # ==================================================================
    hdr = ["", "", "", "", ""]
    timers = {}
    for timer in temp_timers:
        for (key, val) in timer.items():

            # Denote how deep into the dictionary this key goes
            # as determined by the number of prefixing "-" characters
            depth = count_characters(key, "-") / 2

            # Remove any prefixed "-" characters
            new_key = key.strip("-")

            # Add results into the "timers" dictionary as a
            # "flattened" dictionary, for expediency
            # (This is the only way to update a nested dict)
            if depth == 0:
                hdr[0] = new_key
                timers[new_key] = val
            elif depth == 1:
                hdr[1] = new_key
                new_key = f"{hdr[0]}.{new_key}"
                timers[new_key] = val
            elif depth == 2:
                hdr[2] = new_key
                new_key = f"{hdr[0]}.{hdr[1]}.{new_key}"
                timers[new_key] = val
            elif depth == 3:
                hdr[3] = new_key
                new_key = f"{hdr[0]}.{hdr[1]}.{hdr[2]}.{new_key}"
                timers[new_key] = val
            elif depth == 4:
                hdr[4] = new_key
                new_key = f"{hdr[0]}.{hdr[1]}.{hdr[2]}.{hdr[3]}.{new_key}"
                timers[new_key] = val
            else:
                new_key = \
                    f"{hdr[0]}.{hdr[1]}.{hdr[2]}.{hdr[3]}.{hdr[4]}.{new_key}"
                timers[new_key] = val

    return timers


def sum_timers(timers):
    """
    Sums the time in seconds for each GEOS-Chem timer.  Input may be
    a single dict with timing information or a list of dicts.

    Args
    timers : dict|list : GHCP timing information from one or more
                         log files in plain text format

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
                result[key] += float(val)

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
    # Denote the level of the dictionary key by counting "." chars
    depth = count_characters(key, ".")

    # Prefix "--" characters to the end of the key to denote depth
    # to replicate the label style at the end of the GCHP log file
    label = "--"*depth + key.split(".")[-1]

    # Line to print
    pctdiff = np.nan
    if np.abs(ref[key] > 0.0):
        pctdiff = ((dev[key] - ref[key]) / ref[key]) * 100.0
    line = \
        f"{label:<22}  {ref[key]:>18.3f}  {dev[key]:>18.3f}   {pctdiff:>12.3f}"
    if np.abs(pctdiff) >= 10.0:  # Flag diffs > +/- 10%
        line += " *"
    print(line, file=ofile)


def display_timers(ref, ref_label, dev, dev_label, table_file):
    """
    Prints the GCHP timer information to a table.

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
        print("%%% GCHP Benchmark Timing Information", file=ofile)
        print("%%%", file=ofile)
        print(f"%%% Ref = {ref_label}", file=ofile)
        print(f"%%% Dev = {dev_label}", file=ofile)
        print("%"*79, file=ofile)

        # GCHPchem timers
        print("\n", file=ofile)
        print(f"{'GCHPchem Timer':<22}  {'Ref [s]':>18}  {'Dev [s]':>18}   {'% Diff':>12}", file=ofile)
        print("-"*79, file=ofile)
        for key in dev:
            if key.startswith("GCHPchem"):
                print_timer(key, ref, dev, ofile)

        # Summary timers
        print("\n", file=ofile)
        print(f"{'Summary':<22}  {'Ref [s]':>18}  {'Dev [s]':>18}   {'% Diff':>12}", file=ofile)
        print("-"*79, file=ofile)
        for key in dev:
            if key.startswith("All"):
                print_timer(key, ref, dev, ofile)


def make_benchmark_gchp_timing_table(
        ref_files,
        ref_label,
        dev_files,
        dev_label,
        dst="./benchmark",
        overwrite=False,
):
    """
    Creates a table of timing information for GCHP benchmark
    simulations given one or more text files as input.

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
    ref_timers = sum_timers(read_timing_data(ref_files))
    dev_timers = sum_timers(read_timing_data(dev_files))

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
