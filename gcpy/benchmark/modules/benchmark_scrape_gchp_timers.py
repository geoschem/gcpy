#!/usr/bin/env python3
"""
Scrapes GCHP Classic benchmark timing information from one or
more text files.
"""
import os
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


def count_characters(text, char_to_match="-"):
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
    keep_line = True
    temp_timers = []

    # Make sure file exists
    if not os.path.exists(text_file):
        raise FileNotFoundError(f"Could not find {text_file}!")

    # Read the line backwards and get just keep the timing information
    with open(text_file, encoding="utf-8") as ifile:

        for line in reversed(list(ifile)):
            line = line.strip("\n")

            # Set a flag to denote the start & end of timing info
            if "-------- --------- ------ --------- ------" in line:
                keep_line = False
                break

            # Append timing info lines into a list of dicts
            if keep_line:
                substr = line.split()
                key = substr[0].strip()
                val = float(substr[2].strip())
                temp_timers.append({key: val})

        # Because we were reading the end of the file backwards, the
        # entries in temp_timers are reversed.  Now read through them
        # in the forward order.
        hdr = ["", "", ""]
        timers = {}
        for timer in reversed(temp_timers):
            for (key, val) in timer.items():

                # Denote how deep into the dictionary this key goes
                # as determined by the number of prefixing "-" characters
                depth = count_characters(key, "-") / 2

                # Remove any prefixed "-" characters
                new_key = key.strip("-")

                # Add results into the "timers" dictionary as a
                # "flattened" dictionary, for expediency
                if depth == 0:
                    hdr[0] = new_key
                    timers[new_key] = val
                elif depth == 1:
                    hdr[1] = new_key
                    new_key = f"{hdr[0]}_{new_key}"
                    timers[new_key] = val
                elif depth == 2:
                    hdr[2] = new_key
                    new_key = f"{hdr[0]}_{hdr[1]}_{new_key}"
                    timers[new_key] = val
                else:
                    new_key = f"{hdr[0]}_{hdr[1]}_{hdr[2]}_{new_key}"
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
    # Denote the level of the dictionary key by counting "_" chars
    depth = count_characters(key, "_")

    # Prefix "--" characters to the end of the key to denote depth
    # to replicate the label style at the end of the GCHP log file
    label = "--"*depth + key.split("_")[-1]

    # Line to print
    line = f"{label:<25}  {ref[key]:>20.3f}  {dev[key]:>20.3f}"
    print(line, file=ofile)


def display_timers(ref, ref_label, dev, dev_label, table_file):
    """
    Prints the GEOS-Che timer information to a table.

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
        print("%%% GCHP Classic Benchmark Timing Information", file=ofile)
        print("%%%", file=ofile)
        print(f"%%% Ref = {ref_label}", file=ofile)
        print(f"%%% Dev = {dev_label}", file=ofile)
        print("%"*79, file=ofile)
        print("\n", file=ofile)
        print(f"{'Timer':<25}  {'Ref [s]':>20}  {'Dev [s]':>20}", file=ofile)
        print("-"*79, file=ofile)

        # Print timers
        print_timer("All",                            ref, dev, ofile)
        print_timer("All_SetService",                 ref, dev, ofile)
        print_timer("All_SetService_GCHP",            ref, dev, ofile)
        print_timer("All_SetService_GCHP_GCHPctmEnv", ref, dev, ofile)
        print_timer("All_SetService_GCHP_GCHPchem",   ref, dev, ofile)
        print_timer("All_SetService_GCHP_DYNAMICS",   ref, dev, ofile)
        print_timer("All_Initialize",                 ref, dev, ofile)
        print_timer("All_Initialize_GCHP",            ref, dev, ofile)
        print_timer("All_Initialize_GCHP_GCHPctmEnv", ref, dev, ofile)
        print_timer("All_Initialize_GCHP_DYNAMICS",   ref, dev, ofile)
        print_timer("All_Initialize_EXTDATA",         ref, dev, ofile)
        print_timer("All_Initialize_HIST",            ref, dev, ofile)
        print_timer("All_Run",                        ref, dev, ofile)
        print_timer("All_Run_GCHP",                   ref, dev, ofile)
        print_timer("All_Run_GCHP_GCHPctmEnv",        ref, dev, ofile)
        print_timer("All_Run_GCHP_GCHPchem",          ref, dev, ofile)
        print_timer("All_Run_GCHP_DYNAMICS",          ref, dev, ofile)
        print_timer("All_Run_EXTDATA",                ref, dev, ofile)
        print_timer("All_Run_HIST",                   ref, dev, ofile)
        print_timer("All_Finalize",                   ref, dev, ofile)
        print_timer("All_Finalize_GCHP",              ref, dev, ofile)
        print_timer("All_Finalize_GCHP_GCHPctmEnv",   ref, dev, ofile)
        print_timer("All_Finalize_GCHP_GCHPchem",     ref, dev, ofile)
        print_timer("All_Finalize_GCHP_DYNAMICS",     ref, dev, ofile)
        print_timer("All_Finalize_EXTDATA",           ref, dev, ofile)
        print_timer("All_Finalize_HIST",              ref, dev, ofile)


def make_benchmark_timing_table(
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


if __name__ == '__main__':

    REF_FILES = [
        "./execute.gchp_merra2_fullchem_benchmark.log",
        "./execute.gchp_merra2_fullchem_benchmark.log",
    ]
    DEV_FILES = "./execute.gchp_merra2_fullchem_benchmark.log"

    # Debug test
    make_benchmark_timing_table(
        REF_FILES,
        "GCHP 14.4.0 list input",
        DEV_FILES,
        "GCHP 14.4.0 str input",
        dst="./",
        overwrite=True,
)
