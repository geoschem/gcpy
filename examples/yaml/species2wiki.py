#!/usr/bin/env python

"""
species2wiki.py: Creates a MediaWiki table from a GEOS-Chem
species database file in YAML format.  This prints out information
for the table

Calling sequence:
-----------------
./species2wiki.py species   # Creates wiki table w/ general metadata
./species2wiki.py henry     # Creates wiki table w/ Henry's law metadata
./species2wiki.py wetdep    # Creates wiki table w/ wetdep metadata
./species2wiki.py drydep    # Creates wiki table w/ drydep metadata
"""

# Imports
import os
import sys
import yaml

# ======================================================================
# Configurables (MUST EDIT)
# ======================================================================

# YAML file to read
yaml_file = '../../gcpy/species_database.yml'

# ======================================================================
# Methods
# ======================================================================

def print_species_table(metadata):
    """
    Extracts general metadatas from the species database and
    creates a text file in MediaWiki table format.

    Args:
    -----
       metadata : dict
           Dictionary with species metadata
    """

    # Fields to print
    keys_to_print = [
        "Formula", "FullName", "MW_g", "Gas/Aer", "Is_Kpp",
        "Is_Advected", "Is_DryDep", "Is_WetDep", "Is_Photolysis"
    ]

    # Output file
    wiki_table_file = "wiki_species_table.txt"

    # Open file
    with open(wiki_table_file, "w") as f:

        # Loop over species names
        for s in metadata.keys():

            # Skip anchors for other variables
            if "_PROP" in s:
                continue

            # Print the species name as the first column
            spc_db = metadata[s]
            print('\n|-valign="top"', file=f)
            print("|{}".format(s), file=f)

            # Loop over other tags of the species database
            # Special handling for MW_g
            for t in keys_to_print:
                if "EmMW_g" in t or "MolecRatio" in t:
                    pass

                elif "MW_g" in t:
                    if t in spc_db.keys():
                        print("|{}".format(spc_db[t]), end="", file=f)
                        if "EmMW_g" in t and "MolecRatio" in t:
                            print("<br>({}, {}C)".format(
                                spc_db["EmMW_g"],
                                spc_db["MolecRatio"]
                            ), file=f)
                        else:
                            print(file=f)
                    else:
                        print("| -", file=f)

                elif "Gas/Aer" in t:
                    if "Is_Gas" in spc_db.keys():
                        print("|Gas", file=f)
                    elif "Is_Aerosol" in spc_db.keys():
                        print("|Aer", file=f)

                elif "Is_" in t:
                    if t in spc_db.keys():
                        if spc_db[t] is True:
                            print("|X", file=f)
                        else:
                            print("| -", file=f)
                    else:
                        print("| -", file=f)

                else:
                    if t in spc_db.keys():
                        print("|{}".format(spc_db[t]), file=f)
                    else:
                        print("| -", file=f)


        # Close the file
        f.close()


def print_henry_table(metadata):
    """
    Extracts Henry's law metadata from the species database and
    creates a text file in MediaWiki table format.

    Args:
    -----
       metadata : dict
           Dictionary with species metadata
    """

    # Fields to print -- for "GEOS-Chem Species" table on the wiki
    keys_to_print = ["FullName", "DD_Hstar", "Henry_K0", "Henry_CR"]

    # Output file
    wiki_table_file = "wiki_henry_table.txt"

    # Open file
    with open(wiki_table_file, "w") as f:

        # Loop over species names
        for s in metadata.keys():

            # Skip anchors for other variables
            if "_PROP" in s:
                continue

            # Print the species name as the first column
            spc_db = metadata[s]
            print('\n|-valign="top"', file=f)
            print("|{}".format(s), file=f)

            # Loop over other tags of the species database
            for t in keys_to_print:
                if t in spc_db.keys():
                    print("|{}".format(spc_db[t]), file=f)
                else:
                    print("| -", file=f)

        # Close the file
        f.close()


def print_drydep_table(metadata):
    """
    Extracts drydep metadata from the species database and
    creates a text file in MediaWiki table format.

    Args:
    -----
       metadata : dict
           Dictionary with species metadata
    """

    # Fields to print -- for "GEOS-Chem Species" table on the wiki
    keys_to_print = [
        "MW_g", "EmMW_g", "MolecRatio", "Radius", "Density",
        "DD_AeroDryDep", "DD_DustDryDep", "DD_DvzAerSnow",
        "DD_DvzMinVal",  "DD_Hstar", "DD_F0"
    ]

    # Output file
    wiki_table_file = "wiki_drydep_table.txt"

    # Open file
    with open(wiki_table_file, "w") as f:

        # Loop over species names
        for s in metadata.keys():

            # Skip anchors for other variables
            if "_PROP" in s:
                continue

            # Print the species name as the first column
            spc_db = metadata[s]
            print('\n|-valign="top"', file=f)
            print("|{}".format(s), file=f)

            # Loop over other tags of the species database
            for t in keys_to_print:

                if "DD_DvzMinVal" in t:
                    if t in spc_db.keys():
                        print("|{} snow<br>{} land".format(
                            spc_db[t][0],
                            spc_db[t][1]
                        ), file=f)
                    else:
                        print("| -", file=f)

                else:
                    if t in spc_db.keys():
                        print("|{}".format(spc_db[t]), file=f)
                    else:
                        print("| -", file=f)

        # Close the file
        f.close()


def print_wetdep_table(metadata):
    """
    Extracts wetdep metadata from the species database and
    creates a text file in MediaWiki table format.

    Args:
    -----
       metadata : dict
           Dictionary with species metadata
    """

    # Fields to print -- for "GEOS-Chem Species" table on the wiki
    keys_to_print = [
        "MW_g", "EmMW_g", "MolecRatio", "WD_CoarseAer", "WD_AerScafEff",
        "WD_KcScaleFac", "WD_RainoutEff", "WD_RetFactor"
    ]

    # Output file
    wiki_table_file = "wiki_wetdep_table.txt"

    # Open file
    with open(wiki_table_file, "w") as f:

        # Loop over species names
        for s in metadata.keys():

            # Skip anchors for other variables
            if "_PROP" in s:
                continue

            # Print the species name as the first column
            spc_db = metadata[s]
            print('\n|-valign="top"', file=f)
            print("|{}".format(s), file=f)

            # Loop over other tags of the species database
            for t in keys_to_print:

                if "WD_KcScaleFac" in t or "WD_RainoutEff" in t:
                    if t in spc_db.keys():
                        print("|{}".format(spc_db[t][0]), file=f)
                        print("|{}".format(spc_db[t][1]), file=f)
                        print("|{}".format(spc_db[t][2]), file=f)
                    else:
                        print("| -", file=f)
                        print("| -", file=f)
                        print("| -", file=f)

                else:
                    if t in spc_db.keys():
                        print("|{}".format(spc_db[t]), file=f)
                    else:
                        print("| -", file=f)

        # Close the file
        f.close()


def main():
    """
    Main program.  Parses arguments and calls the proper routine
    to create the table in MediaWiki format
    """
    
    # Parse arguments
    n_args = len(sys.argv)
    if n_args == 0 or n_args > 2:
        msg = "Usage: species2wiki.py [species|henry|wetdep|drydep]"
        raise ValueError(msg)
    table = sys.argv[1].upper()

    # Read the YAML file into a dict
    try:
        metadata = yaml.load(open(yaml_file), Loader=yaml.FullLoader)
    except FileNotFoundError:
        msg = "Could not find filename: {}".format(filename)
        raise FileNotFoundError(msg)

    # Print the selected table in MediaWiki format
    if "SPECIES" in table:
        print_species_table(metadata)
    elif "HENRY" in table:
        print_henry_table(metadata)
    elif "WETDEP" in table:
        print_wetdep_table(metadata)
    elif "DRYDEP" in table:
        print_drydep_table(metadata)
    else:
        msg = "Argument must be one of species|henry|wetdep|drydep!"
        raise ValueError(msg)


if __name__ == "__main__":
    main()
