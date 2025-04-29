#!/usr/bin/env python3
"""
Compares GCHP emission diagnostic entries in HISTORY.rc
and in HEMCO_Diagn.rc configuration files for consistency.
"""
from os.path import join, realpath
from sys import argv
from re import sub
from gcpy.constants import ENCODING
from gcpy.util import verify_variable_type


def read_history(filename):
    """
    Reads the HISTORY_Diagn.rc file and returns a list of
    diagnostic container names containing "Emis" or "Inv".

    Args
    filename : str  : Path to the HISTORY.rc file

    Returns
    result   : list : List of diagnostic entries
    """
    verify_variable_type(filename, str)

    with open(filename, "r", encoding=ENCODING) as ifile:
        result = []
        for line in ifile:

            # Only take lines with diagnostic fields
            if "GCHPchem" not in line:
                continue
            if "Budget" in line:
                continue
            if "Emis" not in line and "Inv" not in line:
                continue
            if "Emissions.fields:" in line:
                line = line.replace("Emissions.fields:", "")

            # The diagnostic container is the 1st word
            # (also strip quotes)

            first_word = line.strip().split()[0]
            first_word = first_word.replace("'", "").replace('"', "")

            # Replace multiple # chars with a single # char
            first_word = sub(r"#+", "#", first_word)

            # Add to list of diagnosic entries
            result.append(first_word)

    return result


def read_hemco_diagn(filename):
    """
    Reads the HEMCO_Diagn.rc file and returns a list of
    diagnostic container names containing "Emis" or "Inv".

    Args
    filename : str  : Path to the HEMCO_Diagn.rc file

    Returns
    result   : list : List of diagnostic entries
    """
    verify_variable_type(filename, str)

    with open(filename, "r", encoding=ENCODING) as ifile:
        result = []
        for line in ifile:

            # Strip newlines and split into substrings
            line = line.strip().split()

            # Skip if the line is empty
            if len(line) == 0:
                continue

            # Skip if the 1st word doesn't contain "Emis" or "Inv"
            first_word = line[0]
            if "Emis" not in first_word and "Inv" not in first_word:
                continue

            # Replace quotes 
            first_word.replace("'", "").replace('"', "")

            # Replace multiple "#" characters with a single "#"
            first_word = sub(r"#+", "#", first_word)

            # Append to list of diagnosic entries
            result.append(first_word)

    return result


def compare_containers(history_entries, hco_diagn_entries):
    """
    Compares the list of GCHP emission entries in HISTORY.rc 
    and HEMCO_Diagn.rc.  Returns the list of entries common to 
    both, and in one file but not the other.

    Args
    history_containers   : list : Emission entries in HISTORY.rc
    hco_diagn_containers : list : Emission entries in HEMCO_Diagn.rc

    Returns
    result               : dict : Results of the comparison
    """
    verify_variable_type(history_entries, list)
    verify_variable_type(hco_diagn_entries, list)

    # Convert lists to sets
    set_history = set(history_entries)
    set_hemco = set(hco_diagn_entries)

    # Diagnostic entries common to History and HEMCO
    result = {}
    result["common"] = sorted(list(set_history & set_hemco))

    # Diagnostic entries found in one file but not the other
    result["in_history_not_hemco"] = sorted(list(set_history - set_hemco))
    result["in_hemco_not_history"] = sorted(list(set_hemco - set_history))

    return result


def print_results(result):
    """
    Prints the results of the comparison between HISTORY.rc
    and HEMCO_Diagn.rc

    Args
    result : dict : Dictionary with results of the comparison
    """
    verify_variable_type(result, dict)

    print("="*79)
    print("Common to both HISTORY.rc and HEMCO_Diagn.rc")
    print("="*79)
    for var in result["common"]:
        print(f"  {var}")

    print("")
    print("="*79)
    print("In HISTORY.rc but not in HEMCO_Diagn.rc")
    print("="*79)
    for var in result["in_history_not_hemco"]:
         print(f"  {var}")

    print("")
    print("="*79)
    print("In HEMCO_Diagn.rc but not in HISTORY.rc")
    print("="*79)
    for var in result["in_hemco_not_history"]:
         print(f"  {var}")


def main(path_to_rundir, simulation):
    """
    Main program.  Calls routines to read GCHP HISTORY.rc and
    HEMCO_Diagn.rc files, to compare emission diagnostics, and
    to print the results.

    Args
    path_to_rundir : str : Path to the run/GCHP folder
    simulation     : str : Simulation name (e.g. fullchem)
    """
    verify_variable_type(path_to_rundir, str)
    verify_variable_type(simulation, str)

    # Create paths to HISTORY.rc and HEMCO_Diagn.rc
    history_file = realpath(
        join(
            path_to_rundir,
            "HISTORY.rc.templates",
            f"HISTORY.rc.{simulation}"
        )
    )
    hco_diagn_file = realpath(
        join(
            path_to_rundir,
            "HEMCO_Diagn.rc.templates",
            f"HEMCO_Diagn.rc.{simulation}"
        )
    )

    # Read emission entries
    history_entries = read_history(history_file)
    hco_diag_entries = read_hemco_diagn(hco_diagn_file)

    # Compare and print results
    result = compare_containers(history_entries, hco_diag_entries)
    print_results(result)


if __name__ == '__main__':

    if len(argv) != 3:
        MSG = "Usage: python -m gcpy.examples.diagnostics.check_gchp_emission_diags /path/to/run/GCHP SIMULATION-NAME"
        raise ValueError(MSG)

    main(argv[1], argv[2])
