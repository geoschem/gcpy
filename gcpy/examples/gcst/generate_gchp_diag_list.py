#!/usr/bin/env python3
"""
Reads a GEOS-Chem classic log file listing species names and writes
out the corresponding :literal:`SpeciesConcVV_*` entries for use in a
GCHP :file:`HISTORY.rc` collection.

Examples
--------

Examples
--------
.. code-block:: bash

   $ conda activate gcpy_env
   $ python generate_gchp_speciesconcvv_list.py \
       --diag_prefix <prefix> \
       --input-file </path/to/log/file> \
       --output-file </path/to/output/file>

Parameters
----------
--diag-prefix : str
    Prefix of the diagnostic field (e.g. SpeciesConcVV, JValues, etc.)
--input-file : str
    Path to a GCHP log file.
--output-file : str
    Path of the file to which diagnostic entries will be written.

"""
import argparse
import sys
from gcpy.constants import ENCODING


def read_file(diag_prefix, input_file):
    """
    Reads a GEOS-Chem log file and extracts species names, converting
    each into a GCHP SpeciesConcVV_* diagnostic entry.

    Parameters
    ----------
    diag_prefix : str
        Prefix for the diagnostic field (e.g. "SpeciesConcVV")
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
        found = False

        # Loop over lines
        for line in ifile:

            # Skip until we get to the section w/ names and indices
            if not found:
                if "SPECIES NAMES AND INDICES" in line:
                    found = True
                    continue

            # If we are in the section w/ names and indices ...
            if found:
                
                # Break out of this loop when we encounter a line
                # of "=" characters, that denotes end-of-section
                if "=" in line:
                    found = False
                    break

                # Skip empty lines or hard break lines
                if len(line) == 0:
                    continue

                # Take the species name only and add the diag_prefix
                line = line.strip().split(" ")[0]
                if len(line) == 0:
                    continue
                line = f"{diag_prefix}_{line}"

                # Exclude entries that are not species names
                if f"{diag_prefix}_-" in line:
                    continue
                if f"{diag_prefix}_Name" in line:
                    continue
                if f"{diag_prefix}_SPECIES" in line:
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
            print(f"                            '{var:<29}', \'GCHPchem\',", file=ofile)


def generate_list(diag_prefix, input_file, output_file):
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
    varlist = read_file(diag_prefix, input_file)
    write_file(varlist, output_file)


def main():
    """
    Command-line entry point.
    """
    parser = argparse.ArgumentParser(
        prog='generate_inttest_report.py',
        description='Generate a plain-text GEOS-Chem integration-test report.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            'Examples:\n'
            '\n'
            '  python generate_inttest_report.py \\\n'
            '      --diag_prefix SpeciesConcVV \\\n'
            '      --input-file  gchp.YYYYMMDD_hhmmss.log \\\n'
            '      --output-file HISTORY.rc.snippets.txt \\\n'
        ),
    )
    parser.add_argument(
        '--diag-prefix',
        required=True,
        metavar='DIAG_PREFIX',
        help='Diagnostic field prefix (e.g. "SpeciesConcVV"),',
    )
    parser.add_argument(
        '--input-file',
        required=True,
        metavar='INPUT_FILE',
        help='GCHP log file containing species names and indices.',
    )
    parser.add_argument(
        '--output-file',
        required=True,
        metavar='OUTPUT_FILE',
        help='File where GCHP diagnostic entries will be written.',
    )

    args = parser.parse_args()

    try:
        generate_list(args.diag_prefix, args.input_file, args.output_file)
    except (ValueError, OSError) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
