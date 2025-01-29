#!/usr/bin/env python3
"""
Generates a list of species that differ between versions, that
can be printed on the GEOS-Chem wiki.

Example:

  python -m gcpy.benchmark.modules.benchmark_species_changes \
   --ref-label "14.4.0"                                      \
   --ref-log   "gcc_14.4.0/14.4.0.log                        \
   --dev-label "14.5.0"                                      \
   --dev-log   "gcc_14.5.0/14.5.0.log"                       \
   --spcdb-dir "gcc_14.5.0/14.5.0.log"                       \
   --output-file "wiki_tables.txt"
"""

from os.path import exists, join
import argparse
import pandas as pd
from gcpy.util import read_config_file, verify_variable_type


def read_one_log_file(log_file):
    """
    Parses the GEOS-Chem Classic log file (plain text) with
    timing information and returns a dictionary with the results.

    Args
    text_file : str  : GEOS-Chem log file

    Returns
    species   : dict : Dictionary with species metadata
    """
    keep_line = False
    species = {}
    keys = []

   # Make sure file exists
    if not exists(log_file):
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

    Args
    species          : dict : Dictionary w/ GEOS-Chem metadata
    species_database : dict : GEOS-Chem species database
    keys             : list : Keys in species_database to append
                              to species
    """
    for key in keys:
        if key in species_database:
            if "Formula" in key or "FullName" in key:
                species[key] = species_database[key]
            else:
                species[key] = bool(species_database[key])
        else:
            species[key] = False

    return species


def get_species_metadata(log_file, species_database):
    """
    Returns the relevant metadata for a given species, taken from a
    GEOS-Chem log file as well as from the species database.

    Args
    log_file         : str          : GEOS-Chem log file
    species_database : dict         : GEOS-Chem species database

    Returns
    species_df       : pd.DataFrame : GEOS-Chem species metadata
    """

    # Read the list of species from the log file
    try:
        species = read_one_log_file(log_file)
    except FileNotFoundError as exc:
        msg = f"Could not find {log_file}!"
        raise FileNotFoundError(msg) from exc

    # Also append relevant keys from the species database
    for name in species:
        species[name] = append_keys(
            species[name],
            species_database[name],
            ["Is_Advected", "Is_Aerosol", "Is_Gas", "Formula", "FullName"]
        )

    return pd.DataFrame.from_dict(species).drop("ModelId")


def bool_to_str(val):
    """
    Converts a boolean True value to an "X" for printing
    in the wiki table.

    Args
    val : bool : Boolean value to test.
    """
    string = ""
    if val:
        string = "X"
    return string


def write_wiki_table_header(ofile):
    """
    Writes the header of a wiki table.

    Args
    ofile   : _.io.TextIOWrapper : Output file handle
    """
    line = "{| border=1 cellspacing=0 cellpadding=5\n"
    line += "!width='100px' bgcolor='#CCCCCC'|Name\n"
    line += "!width='100px' bgcolor='#CCCCCC'|Formula\n"
    line += "!width='200px' bgcolor='#CCCCCC'|Fullname\n"
    line += "!width='30px' bgcolor='#CCCCCC'|Advected\n"
    line += "!width='30px' bgcolor='#CCCCCC'|Dry deposited\n"
    line += "!width='30px' bgcolor='#CCCCCC'|Gas\n"
    line += "!width='30px' bgcolor='#CCCCCC'|Photolyzed\n"
    line += "!width='30px' bgcolor='#CCCCCC'|Wet deposited\n"

    print(line, file=ofile)


def write_wiki_row(species, ofile):
    """
    Prints metadata for a given GEOS-Chem species

    Args
    series : pd.Series : GEOS-Chem species metadata
    ofile  : File      : File object for output file
    """
    line = "|-valign='top'\n"
    line += f"|{species.name}\n"
    line += f"|{species['Formula']}\n"
    line += f"|{species['FullName']}\n"
    line += f"|{bool_to_str(species['Is_Advected'])}\n"
    line += f"|{bool_to_str(species['DryDepId'])}\n"
    line += f"|{bool_to_str(species['Is_Gas'])}\n"
    line += f"|{bool_to_str(species['PhotolId'])}\n"
    line += f"|{bool_to_str(species['WetDepId'])}\n"
    line += " "

    print(line, file=ofile)


def write_wiki_table_footer(ofile):
    """
    Writes the footer for a wiki table.

    Args
    ofile   : _.io.TextIOWrapper : Output file handle
    """
    print("|}\n", file=ofile)


def create_table(keys, species, ofile):
    """
    Creates a wiki table containing selected species.

    Args
    keys    : list               : Names of species to include in table
    species : pd.DataFrame       : Species metadata
    ofile   : _.io.TextIOWrapper : Output file handle
    """
    write_wiki_table_header(ofile)

    for key in keys:
        write_wiki_row(species[key], ofile)

    write_wiki_table_footer(ofile)

    print(" ", file=ofile)


def check_for_species_changes(species, ref, dev):
    """
    Prints a list of species with attributes that have changed
    between the Ref and Dev versions.

    species : list         : List of species in both Ref and Dev
    ref     : pd.DataFrame : Species metadata for the Ref version
    ref     : pd.DataFrame : Species metadata for the Dev version
    """
    changed = {}
    keys = dev.index.tolist()

    # Search for species with changed attributes between versions
    for spc in species:
        for key in keys:
            if ref[spc][key] != dev[spc][key]:
                result = {key: {"Ref": ref[spc][key], "Dev": dev[spc][key]}}
                if spc in changed:
                    changed[spc].update(result)
                else:
                    changed[spc] = result

    # Print results
    print("Species with changes (also check advected species)")
    for (key, value) in changed.items():
        print(key, value)


def make_benchmark_species_changes_wiki_tables(
        ref_label,
        ref_log,
        dev_label,
        dev_log,
        spcdb_dir,
        output_file,
):
    """
    Creates tablees of species that have been added and removed
    between Ref and Dev versions.

    Args
    ref_label   : str : Label for the Ref version
    ref_log     : str : Path to log file for the Ref version
    dev_label   : str : Label for the Dev version
    dev_log     : str : Path to log file for the Dev version
    spcdb_dir   : str : Directory containing species_database.yml
    output_file : str : Path to file with generated wiki tables
    """
    verify_variable_type(ref_label, str)
    verify_variable_type(ref_log, str)
    verify_variable_type(dev_label, str)
    verify_variable_type(dev_log, str)
    verify_variable_type(spcdb_dir, str)

    # Read the species database
    input_file = join(spcdb_dir, "species_database.yml")
    try:
        species_database = read_config_file(input_file)
    except FileNotFoundError as exc:
        msg = f"Could not find {input_file}!"
        raise FileNotFoundError(msg) from exc

    # Get species metadata for the Ref and Dev versions
    ref = get_species_metadata(ref_log, species_database)
    dev = get_species_metadata(dev_log, species_database)

    # Get list of species in Ref, Dev and both Ref & Dev
    in_ref = set(list(ref.columns))
    in_dev = set(list(dev.columns))
    in_both = sorted(list(in_dev & in_ref))

    # Lists of species that were added and removed
    added = sorted(list(in_dev - in_ref))
    removed = sorted(list(in_ref - in_dev))

    # Also print list of species in both Ref & Dev
    # to make the table manually
    check_for_species_changes(in_both, ref, dev)

    # Create the wiki tables
    with open(output_file, "w", encoding="utf-8") as ofile:

        print("=== Species added ===\n", file=ofile)
        print(
            f"Species added between versions {ref_label} and {dev_label}:\n",
            file=ofile
        )
        create_table(added, dev, ofile)

        print("=== Species removed===\n", file=ofile)
        print(
            f"Species removed between versions {ref_label} and {dev_label}:\n",
            file=ofile
        )
        create_table(removed, ref, ofile)


def main():
    """
    Parses command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="benchmark_species_changes.py"
    )
    parser.add_argument(
        "--ref-label",
        metavar="REF_LABEL",
        type=str,
        required=True,
        default="GCC_ref",
        help="Label for the Ref version"
    )
    parser.add_argument(
        "--ref-log",
        metavar="REF_LOG",
        type=str,
        required=True,
        help="Log file from the Ref version"
    )
    parser.add_argument(
        "--dev-label",
        metavar="DEV_LABEL",
        type=str,
        required=True,
        default="GCC_dev",
        help="Label for the Dev version"
    )
    parser.add_argument(
        "--dev-log",
        metavar="DEV_LOG",
        type=str,
        required=True,
        help="Log file from the Dev version"
    )
    parser.add_argument(
        "--spcdb-dir",
        metavar="SPCDB_DIR",
        type=str,
        required=True,
        default="./",
        help="Directory where the species_database.yml file resides"
    )
    parser.add_argument(
        "-o", "--output-file",
        metavar="OUTPUT_FILE",
        type=str,
        required=True,
        default="wiki_tables.txt",
        help="File where the wiki table will be written"
    )

    args = parser.parse_args()

    make_benchmark_species_changes_wiki_tables(
        args.ref_label,
        args.ref_log,
        args.dev_label,
        args.dev_log,
        args.spcdb_dir,
        args.output_file
    )


if __name__ == '__main__':
    main()
