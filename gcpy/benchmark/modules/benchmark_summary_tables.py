#!/usr/bin/env python3
"""
Creates summary tables fron GEOS-Chem benchmark output.
"""
import os
import numpy as np
from gcpy.util import \
    add_missing_variables, array_equals, compare_varnames, dataset_reader, \
    get_filepath, make_directory, replace_whitespace, \
    unique_values
from gcpy.constants import \
    ENCODING, SKIP_THESE_VARS

# Suppress numpy divide by zero warnings to prevent output spam
np.seterr(divide="ignore", invalid="ignore")


def create_benchmark_summary_table(
        refpath,
        refstr,
        refdate,
        devpath,
        devstr,
        devdate,
        collections,
        dst="./benchmark",
        overwrite=False,
        outfilename="Summary.txt",
        verbose=False,
        ref_gchp=False,
        dev_gchp=False
):
    """
    Creates a benchmark summary table that shows which data collections
    have difference.  Useful for scanning the 1-hr and 1-month benchmark
    outputs.

    Args:
        refpath: str
            Path to the first data set to be compared (aka "Ref").
        refstr: str
            A string that can be used to identify refdata
            (e.g. a model version number or other identifier).
        refdate: np.datetime64
            Date/time stamp used by the "Ref" data files.
        ref_gchp: bool
            Set to True if the "Ref" data comes from a GCHP run.
            Default value: False
        devpath: str
            Path to the second data set to be compared (aka "Dev").
        devstr: str
            A string that can be used to identify the data set specified
            by devfile (e.g. a model version number or other identifier).
        dev_gchp: bool
            Set to True if the "Ref" data comes from a GCHP run.
            Default value: False
        devdate: np.datetime64
            Date/time stamp used by the "Dev" data files.
        collections: list of strings
            List of diagnostic collections to examine.

    Keyword Args (optional):
        dst: str
            A string denoting the destination folder where the file
            containing emissions totals will be written.
            Default value: "./benchmark"
        overwrite: bool
            Set this flag to True to overwrite files in the
            destination folder (specified by the dst argument).
            Default value: False
        outfilename: str
            Name of the text file which will contain the table of
            emissions totals.
            Default value: "Summary.txt"
        verbose: bool
            Set this switch to True if you wish to print out extra
            informational messages.
            Default value: False

    Remarks:
        This method is mainly intended for model benchmarking purposes,
        rather than as a general-purpose tool.

        Species properties (such as molecular weights) are read from a
        YAML file called "species_database.yml".
    """

    # ==================================================================
    # Open file for output
    # ==================================================================

    # Replace whitespace in the ref and dev labels
    refstr = replace_whitespace(refstr)
    devstr = replace_whitespace(devstr)

    # Create the directory for output
    make_directory(dst, overwrite)

    # Create file
    try:
        f = open(os.path.join(dst, outfilename), "w", encoding=ENCODING)
    except (IOError, OSError, FileNotFoundError) as e:
        msg = f"Could not open {outfilename} for writing!"
        raise e(msg) from e

    # Title strings
    title1 = "### Benchmark summary table"
    title2 = f"### Ref = {refstr}"
    title3 = f"### Dev = {devstr}"

    # Print header to file
    print("#" * 80, file=f)
    print(f"{title1 : <77}{'###'}", file=f)
    print(f"{'###'  : <77}{'###'}", file=f)
    print(f"{title2 : <77}{'###'}", file=f)
    print(f"{title3 : <77}{'###'}", file=f)
    print("#" * 80, file=f)
    print(file=f)

    # ==================================================================
    # Read data and look differences btw Ref & Dev versions
    # ==================================================================

    # Variables to skip
    skip_vars = SKIP_THESE_VARS
    skip_vars.append("AREA")

    # Pick the proper function to read the data
    reader = dataset_reader(
        multi_files=False,
        verbose=verbose
    )

    # Loop over diagnostic files
    for col in collections:

        # Read Ref data
        refdata = reader(
            get_filepath(
                refpath,
                col,
                refdate,
                is_gchp=ref_gchp
            ),
            drop_variables=skip_vars
        )

        # Get Dev data
        devdata = reader(
            get_filepath(
                devpath,
                col,
                devdate,
                is_gchp=dev_gchp
            ),
            drop_variables=skip_vars
        )

        # Make sure that Ref and Dev datasets have the same variables.
        # Variables that are in Ref but not in Dev will be added to Dev
        # with all missing values (NaNs). And vice-versa.
        [refdata, devdata] = add_missing_variables(
            refdata,
            devdata
        )

        # Find all common variables between the two datasets
        vardict = compare_varnames(
            refdata,
            devdata,
            quiet=True
        )

        # List of differences for this collection
        diff_list = []

        # Keep track of which variables are different
        # NOTE: Use 32-point float for comparisons since this is
        # the precision used for History diagnostics.
        for v in vardict["commonvarsData"]:
            if not array_equals(
                    refdata[v],
                    devdata[v],
                    dtype=np.float32
            ):
                diff_list.append(v)

        # Drop duplicate values from diff_list
        diff_list = unique_values(diff_list, drop=[None])

        if len(diff_list) == 0:
            print("-" *  79, file=f)
            print(f"{col}: {devstr} is identical to {refstr}", file=f)
            print(file=f)
        else:
            print("-" *  79, file=f)
            print(f"{col}: {devstr} differs from {refstr}", file=f)
            print("\n  Diagnostics that differ", file=f)
            for i, v in enumerate(diff_list):
                print(f"    {v}", file=f)
                if i > 10:
                    print(f"    ... and {len(diff_list) - 10} others", file=f)
                    break
            print(file=f)

    # ==================================================================
    # Close files
    # ==================================================================
    f.close()


def create_benchmark_sanity_check_table(
        devpath,
        devstr,
        devdate,
        collections,
        dst="./benchmark",
        is_gchp=False,
        overwrite=False,
        outfilename="Diagnostic_Sanity_Check.txt",
        verbose=False,
):
    """
    Creates a diagnostic sanity check table that shows which diagnostic
    variables are zero or NaN everywhere.  This can help to identify
    bugs in diagnostic output.

    Args:
        devpath: str
            Path to the data set to be compared (aka "Dev").
        devstr: str
            A string that can be used to identify the data set specified
            by devfile (e.g. a model version number or other identifier).
        devdate: np.datetime64
            Date/time stamp used by the "Dev" data files.
        collections: list of strings
            List of diagnostic collections to examine.

    Keyword Args (optional):
        dst: str
            A string denoting the destination folder where the file
            containing emissions totals will be written.
            Default value: "./benchmark"
        is_gchp : bool
           Set this flag to true to denote if the data is from GCHP.
        overwrite: bool
            Set this flag to True to overwrite files in the
            destination folder (specified by the dst argument).
            Default value: False
        outfilename: str
            Name of the text file which will contain the table of
            emissions totals.
            Default value: "Summary.txt"
        verbose: bool
            Set this switch to True if you wish to print out extra
            informational messages.
            Default value: False

    Remarks:
        This method is mainly intended for model benchmarking purposes,
        rather than as a general-purpose tool.
    """

    # ==================================================================
    # Initial preparations
    # ==================================================================

    # Replace whitespace in the ref and dev labels
    devstr = replace_whitespace(devstr)

    # Create the directory for output (if necessary)
    make_directory(dst, overwrite)
    outfilename = os.path.join(dst, outfilename)

    # Pick the proper function to read the data
    reader = dataset_reader(
        multi_files=False,
        verbose=verbose
    )

    # Variables to skip
    skip_vars = SKIP_THESE_VARS
    skip_vars.append("AREA")

    # ==================================================================
    # Open output file and write header
    # ==================================================================
    with open(outfilename, "w", encoding=ENCODING) as ofile:

        # Title strings
        title1 = "### Benchmark diagnostic sanity check table"
        title2 = f"### Dev = {devstr}"

        # Print header to file
        print("#" * 80, file=ofile)
        print(f"{title1 : <77}{'###'}", file=ofile)
        print(f"{'###'  : <77}{'###'}", file=ofile)
        print(f"{title2 : <77}{'###'}", file=ofile)
        print("#" * 80, file=ofile)

        # ==============================================================
        # Loop over diagnostic collections and scan files
        # ==============================================================
        for col in collections:

            # Read data into an xr.DataSet object
            file_name = get_filepath(
                devpath,
                col,
                devdate,
                is_gchp=is_gchp,
            )
            dset = reader(
                file_name,
                drop_variables=skip_vars
            )

            # Determine which variables are all zeroes or NaN
            all_zeros_or_nans = []
            for var in dset.data_vars:
                data = dset[var].values
                if np.all(data == 0) or np.all(data == np.nan):
                    all_zeros_or_nans.append(var)

            # ===========================================================
            # Print results for each collection
            # ===========================================================
            print("", file=ofile)
            print("="*80, file=ofile)
            print(f"{os.path.basename(file_name)}", file=ofile)
            print("="*80, file=ofile)
            print("", file=ofile)

            if len(all_zeros_or_nans) == 0:
                print("No variables were all zero or all NaN", file=ofile)
            else:
                print("These variables were all zero or all NaN:", file=ofile)
                for var in all_zeros_or_nans:
                    print(f"   {var}", file=ofile)
