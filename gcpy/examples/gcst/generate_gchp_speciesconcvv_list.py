#!/usr/bin/env python3
"""
Reads a GEOS-Chem classic log file listing species names and writes
out the corresponding SpeciesConcVV_* entries for use in a GCHP
HISTORY.rc collection.
"""
from sys import argv
from gcpy.constants import ENCODING


def read_file(input_file):
    """
    Reads a GEOS-Chem log file and extracts species names, converting
    each into a GCHP SpeciesConcVV_* diagnostic entry.

    Parameters
    ----------
    input_file : str
        Path to the GEOS-Chem classic log file listing species names.

    Returns
    -------
    varlist : list of str
        SpeciesConcVV_* entries, one per species, in reverse order
        (as they should appear in the GCHP HISTORY.rc file).
    """
    # Open file
    with open(input_file, "r", encoding=ENCODING) as ifile:

        varlist = []
        
        # Loop over lines
        for line in ifile:

            # Skip empty lines or hard break lines
            if len(line) == 0:
                continue
            if "=" in line:
                continue

            # Take the species name only and prefix SpeciesConcVV_
            line = f"SpeciesConcVV_{line.strip().split(" ")[0]}"

            # Exclude entries that are not species names
            if "SpeciesConcVV_-" in line:
                continue
            if "SpeciesConcVV_Name" in line:
                continue
            if "SpeciesConcVV_SPECIES" in line:
                continue
            if "SpeciesConcVV_ " in line:
                continue

            # Append into the list
            varlist.append(line)

        # Return the reversed list (as is in the GCHP HISTORY.rc)
        return varlist[::-1]


def write_file(varlist, output_file):
    """
    Writes SpeciesConcVV_* entries to a file, formatted for use in a
    GCHP HISTORY.rc collection's .fields list.

    Parameters
    ----------
    varlist : list of str
        SpeciesConcVV_* entries to write, e.g. as returned by
        read_file().
    output_file : str
        Path to the file where the entries will be written.
    """
    with open(output_file, "w", encoding=ENCODING) as ofile:
        for var in varlist:
            print(f"                             '{var:<29}', \"GCHPchem\",", file=ofile)


def main(input_file, output_file):
    """
    Reads species names from a GEOS-Chem log file and writes the
    corresponding SpeciesConcVV_* entries to an output file for use
    in a GCHP HISTORY.rc collection.

    Parameters
    ----------
    input_file : str
        Path to the GEOS-Chem classic log file listing species names.
    output_file : str
        Path to the file where the SpeciesConcVV_* entries will be
        written.
    """

    # Read the species list from a GEOS-CHem run
    varlist = read_file(input_file)

    write_file(varlist, output_file)

    
if __name__ == '__main__':

    if len(argv) != 3:
        raise ValueError("Usage: python -m generate_gchp_speciesconcvv_list log_file output file")

    main(argv[1], argv[2])
