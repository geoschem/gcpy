#!/usr/bin/env python3

import os
import sys
from gcpy.util import read_config_file
"""
Generates a list of species that differ between versions, that
can be printed on the GEOS-Chem wiki.
"""

def read_one_log_file(log_file):
    """
    Parses the GEOS-Chem Classic log file (plain text) with
    timing information and returns a dictionary with the results.

    Args
    text_file : str  : Text file with timing information

    Returns
    result    : dict : Dictionary with timing information
    """
    keep_line = False
    species = {}
    keys = []
    results = []

   # Make sure file exists
    if not os.path.exists(log_file):
        raise FileNotFoundError(f"Could not find {log_file}!")

    # Read the log file and just keep species information
    with open(log_file, encoding="utf-8") as ifile:
        for line in list(ifile):
            line = line.strip("\n")

            # Set a flag to denote the start of species info
            if "SPECIES NAMES AND INDICES" in line:
                keep_line = True
                continue

            # Append species info lines into a list
            if keep_line:

                # Break out of the loop after the species info section
                if "==========" in line:
                    keep_line = False
                    break

                # Skip hard rule and/or empty lines
                if "----------" in line or len(line) == 0:
                    continue

                # Get the list of columns
                if "Name" in line:
                    keys = line.split()
                    continue

                # Add each species
                substr = line.split()
                for idx in range(1, 7):
                    substr[idx] = not "-" in substr[idx]
                species[substr[0]] = dict(zip(keys[1:7], substr[1:7]))

    return species


def append_keys(species, species_database, keys):
    """
    Copies dictionary keys from the species database
    to the existing dictionary
    """

    for key in keys:
        if key in species_database:
            species[key] = bool(species_database[key])
        else:
            species[key] = False

    return species


def get_species_metadata(input_file, species_database):
    """
    """
    species = {}

    # Read log from the Ref version
    try:
        species = read_one_log_file(input_file)
    except FileNotFoundError as exc:
        msg = f"Could not find {input_file}!"
        raise FileNotFoundError(msg) from exc

    # Also append relevant keys from the species database
    for name in species:
        species[name] = append_keys(
            species[name],
            species_database[name],
            ["Is_Advected", "Is_Aerosol", "Is_Gas"]
        )

    return species


def print_species_metadata():
    """
    """
    return

def main(
        ref_dir,
        ref_log,
        dev_dir,
        dev_log,
        spcdb_dir
):
    """
    Main program
    """

    # Read the species database
    input_file = os.path.join(spcdb_dir, "species_database.yml")
    try:
        species_database = read_config_file(input_file)
    except FileNotFoundError as exc:
        msg = f"Could not find {input_file}!"
        raise FileNotFoundError(msg) from exc

    # Get species metadata for the Ref version
    input_file = os.path.join(ref_dir, ref_log)
    ref_spc = get_species_metadata(input_file, species_database)

    # Get species metadata for the Dev version
    input_file = os.path.join(dev_dir, dev_log)
    dev_spc = get_species_metadata(input_file, species_database)


    print(ref_spc["ACET"])


if __name__ == '__main__':

    if len(sys.argv) != 3:
        MSG = "Usage: ./benchmark_check_species.py REF-LOG-FILE DEV-LOG_FILE"
        raise ValueError(MSG)

#    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    main(
        "GCC_ref",
        "14.5.0-alpha.5.log",
        "GCC_dev",
        "14.5.0-alpha.6.log",
        "GCC_dev",
    )
