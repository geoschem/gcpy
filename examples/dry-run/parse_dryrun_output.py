#!/usr/bin/env python

"""
Description:
------------
This Python script (assumes Python3) reads a GEOS-Chem or
HEMCO-standalone log file containing dry-run output and
creates the following:

(1) A list of unique files that are required for the GEOS-Chem or
    HEMCO-standalone simulation;

(2) A script to download missing files from the AWS s3://gcgrid bucket;

(3) A script to download missing files from Compute Canada.

This script only requires the "os" and "sys" packages, which are
part of core Python.

Remarks:
--------
NOTE: As of 19 Nov 2019, the Compute Canada script has not yet
been implemented.  Still working on it...
"""

import os
import sys


def extract_pathnames_from_log(filename):
    """
    Returns a list of pathnames from a GEOS-Chem log file.

    Args:
    -----
        filename : str
            GEOS-Chem or HEMCO-standalone log file with dry-run output.

    Returns:
    --------
        data paths : dict
            data_paths['comments']: Dry-run comment lines
            data_paths['found'] : List of file paths found on disk
            data_paths['missing']: List of file paths that are missing

    Author:
    -------
        Jiawei Zhuang (jiaweizhuang@g.harvard.edu)
        Modified by Bob Yantosca (yantosca@seas.harvard.edu)
    """

    # Initialization
    comments = ['!'*79,
                '!!! LIST OF (UNIQUE) FILES REQUIRED FOR THE SIMULATION']
    data_found = set()
    data_missing = set()

    # Open file (or die with error)
    try:
        f = open(filename, "r")
    except FileNotFoundError:
        raise FileNotFoundError("Could not find file {}".format(filename))

    # Read data from the file line by line.
    # Add file paths to the data_list set.
    line = f.readline()
    while line:
        upcaseline = line.upper()

        # Search for data paths that have been found
        if (": OPENING" in upcaseline) or (": READING" in upcaseline):
            data_found.add(line.split()[-1])

        # Search for data paths that are missing
        if ("FILE NOT FOUND" in upcaseline):
            data_missing.add(line.split()[-1])

        # Search for certain dry-run comment strings
        # (and make sure to prevent duplicates)
        if ('!!! Sta' in line) or ('!!! End' in line) or \
           ('!!! Sim' in line) or ('!!! Met' in line) or \
           ('!!! Gri' in line):
            if line.rstrip() not in comments:
                comments.append(line.rstrip())

        # Read next line
        line = f.readline()

    # Add another line to the comment list
    comments.append('!'*79)

    # Close file and return
    # The "sorted" command will return unique values
    f.close()
    return {'comments': comments,
            'found': sorted(list(data_found)),
            'missing': sorted(list(data_missing))}


def write_unique_paths(data_paths, filename):
    """
    Writes unique data paths from dry-run output to a file.

    Args:
    -----
        data_paths : dict
            data_paths['comments']: Dry-run comment lines
            data_paths['found'] : List of file paths found on disk
            data_paths['missing']: List of file paths that are missing

        filename : str
            Name of the file that will be created.
    """
    combined_paths = data_paths['found'] + data_paths['missing']
    combined_paths.sort()

    with open(filename, "w") as f:
        for comment in data_paths['comments']:
            print(comment, file=f)
        for path in combined_paths:
            print(path, file=f)
        for comment in data_paths['comments']:
            print(comment, file=f)
        f.close()


def create_script_for_aws_s3(data_paths, script_name):
    """
    Creates a data download script to obtain missing files
    from the GEOS-Chem s3://gcgrid bucket on the AWS cloud.

    Args:
    -----
        data_paths : dict
            data_paths['comments']: Dry-run comment lines
            data_paths['found'] : List of file paths found on disk
            data_paths['missing']: List of file paths that are missing

        script_name : str
            Name of the script that will be created.
    """
    with open(script_name, "w") as f:

        # Write shebang line to script
        print('#!/bin/bash', file=f)
        print(file=f)

        # Write download commands for only the missing data files
        for pathname in data_paths['missing']:

            # Search for
            if 'ExtData' in pathname:

                # Create AWS data download coommand
                index = pathname.find('ExtData')+7
                local_dir = os.path.dirname(pathname)
                remote_path = 's3://gcgrid' + pathname[index:]
                cmd = 'aws s3 cp ' + remote_path + ' ' + local_dir + '/'

                # Write download command to file
                print(cmd, file=f)
                print(file=f)

        # Close file and make it executable
        f.close()
        os.chmod(script_name, 0o755)


def main():
    """
    Main program.
    """

    # Number of arguments
    n_args = len(sys.argv)

    # Test the # of arguments
    if n_args == 0 or n_args > 2:
        raise FileNotFoundError("Usage: parse_dryrun_output.py LOG-FILE-NAME")

    # Name of input file
    log_file_name = sys.argv[1]

    # Get a list of data paths, both found and missing:
    data_paths = extract_pathnames_from_log(log_file_name)

    # Write a list of unique file paths
    write_unique_paths(data_paths, log_file_name + '.unique')

    # Create script to download missing files from AWS S3
    create_script_for_aws_s3(data_paths, 'download_from_aws_s3.sh')


if __name__ == "__main__":
    main()
