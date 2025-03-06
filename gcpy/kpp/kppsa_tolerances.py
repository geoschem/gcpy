#!/usr/bin/env python3
"""
Reads KPP-Standalone tolerance loop output and displays results
in sorted order.  Useful for determining the best combination of
absolute and relative tolerances for a particular integrator.
"""

import argparse
from numpy import int32
import pandas as pd
from gcpy.util import verify_variable_type
from gcpy.kpp.kppsa_utils import kppsa_read_tolerance_loop


def kppsa_print_sorted_tolerances(
        dframe,
        sort_by,
        ascending,
        values,
):
    """
    Prints the results from the KPP-Standalone tolerance loop
    in sorted order.

    Args:
    dframe    : str  : DataFrame containing sorted values
    sort-by   : str  : Key by which "dframe" is sorted
    ascending : bool : Is "dframe" sorted in ascending order (T/F)?
    values    : int  : Number of sorted values to display
    """
    verify_variable_type(dframe, pd.DataFrame)
    verify_variable_type(sort_by, str)
    verify_variable_type(ascending, bool)
    verify_variable_type(values, int)

    # Print header
    if ascending:
        order = ("smallest", "ascending")
    else:
        order = ("largest", "descending")
    print("="*79)
    msg = f"Showing the {values} {order[0]} values, "
    msg += f"sorted by key '{sort_by}', in {order[1]} order"
    print(msg)
    print("="*79)
    print()

    # Print values as integers or reals
    key = ""
    count = 0
    for index, row in dframe.iterrows():
        count += 1
        print(f"{count}: Combination {index.strip('combo')}")
        for key in row.keys():
            if key in ["FunCount", "JacCount", "TotSteps",
                       "AccSteps", "RejSteps", "LuDecomps", "Substs"]:
                print(f"  {key.ljust(9)} : {int32(row[key])}")
            else:
                print(f"  {key.ljust(9)} : {row[key]}")
        print()


def kppsa_tolerances(
        filename,
        sort_by,
        ascending=True,
        values=5,
):
    """
    Reads data from a KPP-Standalone loop over absolute and
    relative tolerances and displays

    Args:
    file      : str  : File w/ KPP-Standalone tolerance loop output
    sort-by   : str  : DataFrame key to sort by

    Kwargs
    ascending : bool : Sort in ascending order (Default: true)
    values    : int  : Number of values to display (Default: 5)
    """
    verify_variable_type(filename, str)
    verify_variable_type(sort_by, (str, type(None)))
    verify_variable_type(ascending, bool)
    verify_variable_type(values, int)

    dframe = kppsa_read_tolerance_loop(filename)
    dframe = dframe.sort_values(sort_by, ascending=ascending).head(values)

    kppsa_print_sorted_tolerances(dframe, sort_by, ascending, values)


def main():
    """
    Parses arguments from the command line and calls
    kppsa_plot_species_at_sites.

    Command-line arguments
    --filename  (or -f) : File w/ KPP-Standalone tolerance loop output
    --sort-by   (or -s) : DataFrame key to sort by
    --ascending (or -a) : Sort in ascending order (true/false)?
    --values    (or -v) : Number of values to display
    """
    # Tell the parser which arguments to look for
    parser = argparse.ArgumentParser(
        description="KPP-Standalone tolerance loop query program"
    )
    parser.add_argument(
        "--filename", '-f',
        metavar="FILENAME",
        type=str,
        required=True,
        help="Filename containing KPP-Standalone tolerance loop output"
    )
    parser.add_argument(
        "--sort-by", '-s',
        metavar="SORTBY",
        type=str,
        required=True,
        help="DataFrame key to sort by",
    )
    parser.add_argument(
        "--ascending", "-a",
        required=False,
        action="store_true",
        help="Sort in ascending order",
    )
    parser.add_argument(
        "--descending", "-d",
        required=False,
        dest="ascending",
        action="store_false",
        help="Sort in ascending order",
    )
    parser.add_argument(
        "--values", '-v',
        metavar="VALUES",
        type=int,
        required=False,
        default=5,
        help="Number of values to display",
    )

    # Parse command-line arguments
    args = parser.parse_args()

    # Sort data as requested
    kppsa_tolerances(
        args.filename,
        args.sort_by,
        ascending=args.ascending,
        values=args.values,
    )


if __name__ == '__main__':
    main()
