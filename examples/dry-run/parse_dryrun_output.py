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
NOTE: As of 26 Nov 2019, the Compute Canada script has not yet
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
            data_paths["comments"]: Dry-run comment lines
            data_paths["found"] : List of file paths found on disk
            data_paths["missing"]: List of file paths that are missing

    Author:
    -------
        Jiawei Zhuang (jiaweizhuang@g.harvard.edu)
        Modified by Bob Yantosca (yantosca@seas.harvard.edu)
    """

    # Initialization
    comments = ["!"*79,
                "!!! LIST OF (UNIQUE) FILES REQUIRED FOR THE SIMULATION"]
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

        # Convert line to uppercase for string match
        upcaseline = line.upper()

        # Search for data paths that have been found
        if (": OPENING" in upcaseline) or (": READING" in upcaseline):
            data_found.add(line.split()[-1])

        # Search for data paths that are missing
        elif ("FILE NOT FOUND" in upcaseline):
            data_missing.add(line.split()[-1])

        # Search for certain dry-run comment strings
        # (and make sure to prevent duplicates)
        elif ("!!! STA" in upcaseline) or ("!!! END" in upcaseline) or \
             ("!!! SIM" in upcaseline) or ("!!! MET" in upcaseline) or \
             ("!!! GRI" in upcaseline):
            if line.rstrip() not in comments:
                comments.append(line.rstrip())

        else:
            pass

        # Read next line
        line = f.readline()

    # Add another line to the comment list
    comments.append("!"*79)

    # Close file and return
    # The "sorted" command will return unique values
    f.close()
    return {"comments": comments,
            "found": sorted(list(data_found)),
            "missing": sorted(list(data_missing))}


def get_run_info(input_geos_file="./input.geos"):
    """
    Searches through the input.geos file for GEOS-Chem run parameters.

    Returns:
    -------
        run_info : dict of str
            Contains the GEOS-Chem run parameters: start_date,
            start_time, end_date, end_time, met, grid, and sim.
    """
    run_info = {}

    try:
        with open(input_geos_file, "r") as f:
            for line in f:
                if "Start YYYYMMDD" in line:
                    substr = line.split(":")[1]
                    run_info["start_date"] = (substr.split(" ")[1]).strip()
                    run_info["start_time"] = (substr.split(" ")[2]).strip()
                elif "End   YYYYMMDD" in line:
                    substr = line.split(":")[1]
                    run_info["end_date"] = (substr.split(" ")[1]).strip()
                    run_info["end_time"] = (substr.split(" ")[2]).strip()
                elif "Met field" in line:
                    run_info["met"] = (line.split(":")[1]).strip()
                elif "Simulation name" in line:
                    run_info["sim"] = (line.split(":")[1]).strip()
                elif "Grid resolution" in line:
                    grid = (line.split(":")[1]).strip()

                    # Adjust grid string to match file names
                    if '4.0x5.0' in grid:
                        run_info['grid'] = '4x5'
                    elif '2.0x2.5' in grid:
                        run_info['grid'] = '2x2.5'
                    elif '0.5x0.625' in grid:
                        run_info['grid'] = '05x0625'
                    elif '0.25x0.3125' in grid:
                        run_info['grid'] = '025x03125'
                    break
            f.close()
    except FileNotFoundError:
        raise FileNotFoundError("Could not open {}".format(input_geos_file))

    return run_info


def expand_restart_file_names(data_paths, run_info):
    """
    Tests if the GEOS-Chem restart file is a symbolic link to
    ExtData.  If so, will append the link to the remote file
    to the line in which the restart file name is found.
    """
    prefix = ""

    # Get the prefix to ExtData
    for path in data_paths["found"] + data_paths["missing"]:
        if 'ExtData' in path:
             index = path.find("ExtData")+8
             prefix = path[0:index] + 'GEOSCHEM_RESTARTS/v2018-11/'
             break

    # Search for the restart file name in the found files
    new_list = []
    for path in data_paths["found"]:
        if 'GEOSChem.Restart' in path:
            realpath = prefix + 'initial_GEOSChem_rst.' + \
                       run_info["grid"] + '_' + run_info["sim"] + '.nc'
            path = path + ' --> ' + realpath
        new_list.append(path)
    data_paths['found'] = sorted(new_list)

    # Search for the restart file name in the missing files
    new_list = []
    for path in data_paths["missing"]:
        if 'GEOSChem.Restart' in path:
            realpath = prefix + 'initial_GEOSChem_rst.' + \
                       run_info["grid"] + '_' + run_info["sim"] + '.nc'
            path = path + ' --> ' + realpath
        new_list.append(path)
    data_paths['missing'] = sorted(new_list)

    # Return the updated data paths
    return data_paths


def write_unique_paths(data_paths, filename):
    """
    Writes unique data paths from dry-run output to a file.

    Args:
    -----
        data_paths : dict
            data_paths["comments"]: Dry-run comment lines
            data_paths["found"] : List of file paths found on disk
            data_paths["missing"]: List of file paths that are missing

        filename : str
            Name of the file that will be created.
    """
    combined_paths = data_paths["found"] + data_paths["missing"]
    combined_paths.sort()

    try:
        with open(filename, "w") as f:
            for comment in data_paths["comments"]:
                print(comment, file=f)
            for path in combined_paths:
                print(path, file=f)
            for comment in data_paths["comments"]:
                print(comment, file=f)
        f.close()
    except FileNotFoundError:
        raise FileNotFoundError("Could not write {}".format(filename))


def create_script_for_aws_s3(data_paths, script_name):
    """
    Creates a data download script to obtain missing files
    from the GEOS-Chem s3://gcgrid bucket on the AWS cloud.

    Args:
    -----
        data_paths : dict
            data_paths["comments"]: Dry-run comment lines
            data_paths["found"] : List of file paths found on disk
            data_paths["missing"]: List of file paths that are missing

        script_name : str
            Name of the script that will be created.
    """
    with open(script_name, "w") as f:

        # Write shebang line to script
        print("#!/bin/bash\n", file=f)

        # Write comment header
        title = "#"*60 + "\n" + \
                "# Script to download data from AWS s3://gcgrid\n" + \
                "# to a local disk, such as an AWS EBS volume.\n" + \
                "#\n" + \
                "# This script is automatically generated with:\n" + \
                "#    ./geos --dryrun > log.dryrun\n" + \
                "#    ./parse_dryrun_output.py log.dryrun\n" + \
                "#\n" + \
                "# To download data files, run this script\n" + \
                "# from your GEOS-Chem run-directory:\n" + \
                "#    ./download_from_aws_s3.sh\n" + \
                "#\n" + \
                "# USE WITH CAUTION!  Downloading data from AWS\n" + \
                "# to non-AWS systems incurs a data egress fee!\n" + \
                "#"*60 + "\n"
        print(title, file=f)

        # Write download commands for only the missing data files
        for path in data_paths["missing"]:

            # For missing GEOS-Chem restart files, create a symbolic
            # link to the coresponding initial restart file in the
            # s3://gcgrid/GEOSCHEM_RESTARTS folder.
            if "-->" in path:

                # First copy the restart file over from
                # s3://gcgrid/GEOSCHEM_RESTARTS to the
                remote_rst = (path.split("-->")[1]).strip()
                local_rst = (path.split("-->")[0]).strip()
                index1 = remote_rst.find("initial")
                index2 = remote_rst.find("ExtData") + 7
                prefix = remote_rst[0:index1]
                remote_rst = "s3://gcgrid" + remote_rst[index2:]
                cmd = "aws s3 cp --request-payer=requester " + \
                    remote_rst + " " + prefix

                # Write download command to script
                print(cmd, file=f)
                print(file=f)

                # Remove the prior link for safety's sake
                cmd = 'unlink ' + local_rst
                print(cmd, file=f)

                # Then create a symbolic link to the copied-over file
                cmd = 'ln -s ' + prefix + remote_rst[index1-9:] + \
                      ' ' + local_rst

                # Write symbolic link command to script
                print(cmd, file=f)
                print(file=f)

            # For files that are in ExtData, create the corresponding
            # command needed to download them from AWS s3://gcgrid.
            elif "ExtData" in path:
                index = path.find("ExtData")+7
                local_dir = os.path.dirname(path)
                remote_path = "s3://gcgrid" + path[index:]
                cmd = "aws s3 cp --request-payer=requester " + \
                    remote_path + " " + local_dir + "/"

                # Write download command to script
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
        raise ValueError("Usage: parse_dryrun_output.py LOG-FILE-NAME")

    # Name of input file
    log_file_name = sys.argv[1]

    # Get information about the run
    run_info = get_run_info()

    # Get a unique list of data paths, both found and missing:
    # Expand the data paths to include links to restart files
    data_paths = extract_pathnames_from_log(log_file_name)
    data_paths = expand_restart_file_names(data_paths, run_info)

    # Write a list of unique file paths
    write_unique_paths(data_paths, log_file_name + ".unique")

    # Create script to download missing files from AWS S3
    create_script_for_aws_s3(data_paths, "download_from_aws_s3.sh")


if __name__ == "__main__":
    main()
