'''Specific utilities for GEOS-Chem and GCPy on the AWS Cloud'''

import os
from . import core


def s3_download_cmds_from_log(log_files,
                              s3_cp_cmd='aws s3 cp --request-payer=requester',
                              s3_root='s3://gcgrid',
                              prefix_filter='/home/ubuntu/ExtData'):
    '''
    Reads a GEOS-Chem log file and creates a list of bash commands to
    download data from the GEOS-Chem S3 bucket on AWS.  This is
    a convenient way of only downloading the data that is actually
    read by GEOS-Chem to your AWS instance.

    Args:
    -----
        log_files : list of str | str
            List of log files to read (e.g. from GEOS-Chem and HEMCO).
            Can also be a single string.

    Keyword Args (optional):
    ------------------------
        s3_cp_cmd : str
            Command to copy files from AWS S3.
            Default value: aws s3 cp --request-payer=requester

        s3_root : str
            Root folder of the GEOS-Chem S3 bucket on AWS.
            Default value: s3://gcgrid

        prefix_filter : str
            Restricts the output to file paths starting with
            this prefix.
            Default value: /home/ubuntu/ExtData

    Returns:
    --------
        s3_cmds : list of str
            AWS S3 download commands for each file that is
            read by GEOS-Chem, as indicated in the log file.
    '''

    # ==================================================================
    # Initialization
    # ==================================================================
    
    # Error check arguments.  If log_files is a string, make it a
    # list.  Then make sure log_files is a list or exit w/ error.
    if isinstance(log_files, str):
        log_files = [log_files]
    if not isinstance(log_files, list):
        raise TypeError('The log_files argument must be a list!')

    # Define empty lists
    paths = []
    cmd_list = []

    # ==================================================================
    # Parse the log files and return a list of path names
    # ==================================================================
    for log_file in log_files:
        paths = paths + core.extract_pathnames_from_log(log_file,
                                                        prefix_filter)


    # ==================================================================
    # Create the list of bash commands to download files from S3
    # ==================================================================
    for path in paths:

        # Local path name and dir name
        local_path = prefix_filter + path
        local_dir = os.path.dirname(local_path)
        local_file = os.path.basename(local_path)

        # S3 path name
        s3_path = '{}{}'.format(s3_root, path)

        # Create the download command.  For some files we have to
        # link to other files, because these don't usually get
        # synced to the AWS gcgrid S3 bucket.
        if 'gmi.clim.IPMN.geos5.2x25.nc' in local_file:

            # ------------------------------
            # PMN -> IPMN and NPMN
            # ------------------------------

            # First download the PMN file
            needed_file = 'gmi.clim.PMN.geos5.2x25.nc'
            needed_path = '{}/{}'.format(local_dir, needed_file)
            needed_s3_path = '{}/{}'.format(s3_root, needed_file)
            cmd_pref = 'if [[ !(-f {}) ]]; then'.format(needed_path)
            cmd = '{} {} {} {}/; fi'.format(
                cmd_pref, s3_cp_cmd, needed_s3_path, local_dir)
            cmd_list.append(cmd)

            # Then create links for the other species
            for species in ['IPMN', 'NPMN']:
                cmd = 'ln -s {}/{} {}/gmi.clim.{}.geos.2x2.5.nc'.format(
                    local_dir, needed_file, local_dir, species)
                cmd_list.append(cmd)

        elif 'gmi.clim.RIPA.geos5.2x25.nc' in local_file:

            # ------------------------------
            # RIP -> RIPA, RIPB, RIPD
            # ------------------------------

            # First download the RIP file
            needed_file = 'gmi.clim.RIP.geos5.2x25.nc'
            needed_path = '{}/{}'.format(local_dir, needed_file)
            needed_s3_path = '{}/{}'.format(s3_root, needed_file)
            cmd_pref = 'if [[ !(-f {}) ]]; then'.format(needed_path)
            cmd = '{} {} {} {}/; fi'.format(
                cmd_pref, s3_cp_cmd, needed_s3_path, local_dir)
            cmd_list.append(cmd)

            # Then create links for the other species
            for species in ['RIPA', 'RIPB', 'RIPD']:
                cmd = 'ln -s {}/{} {}/gmi.clim.{}.geos.2x2.5.nc'.format(
                    local_dir, needed_file, local_dir, species)
                cmd_list.append(cmd)

        else:

            # ------------------------------
            # No special handling needed
            # ------------------------------
            cmd_pref = 'if [[ !(-f {}) ]]; then'.format(local_path)
            cmd = '{} {} {} {}/; fi'.format(
                cmd_pref, s3_cp_cmd, s3_path, local_dir)
            cmd_list.append(cmd)

    return cmd_list
        

def s3_list_cmds_from_log(log_files,
                          s3_ls_cmd='aws s3 ls --request-payer=requester',
                          s3_root='s3://gcgrid',
                          prefix_filter='/home/ubuntu/ExtData'):
    '''
    Reads a GEOS-Chem log file and creates a list of bash commands to
    list data from the GEOS-Chem S3 bucket on AWS.

    Args:
    -----
        log_files : list of str | str
            List of log files to read (e.g. from GEOS-Chem and HEMCO).
            Can also be a single string.

    Keyword Args (optional):
    ------------------------
        s3_ls_cmd : str
            Command to copy files from AWS S3.
            Default value: aws s3 ls --request-payer=requester --recursive

        s3_root : str
            Root folder of the GEOS-Chem S3 bucket on AWS.
            Default value: s3://gcgrid

        prefix_filter : str
            Restricts the output to file paths starting with
            this prefix.
            Default value: /home/ubuntu/ExtData

    Returns:
    --------
        s3_cmds : list of str
            AWS S3 file list commands for each file that is
            read by GEOS-Chem, as indicated in the log file.
    '''

    # ==================================================================
    # Initialization
    # ==================================================================

    # Error check arguments.  If log_files is a string, make it a
    # list.  Then make sure log_files is a list or exit w/ error.
    if isinstance(log_files, str):
        log_files = [log_files]
    if not isinstance(log_files, list):
        raise TypeError('The log_files argument must be a list!')

    # Define empty lists
    paths = []
    cmd_list = []

    # ==================================================================
    # Parse the log files and return a list of path names
    # ==================================================================
    for log_file in log_files:
        paths = paths + core.extract_pathnames_from_log(log_file,
                                                        prefix_filter)

    # ==================================================================
    # Create the list of bash commands to list files in S3
    # ==================================================================
    for path in paths:

        # S3 path name
        s3_path = '{}{}'.format(s3_root, path)

        # Populate a list of bash commands for listing files
        cmd = '{} {}'.format(s3_ls_cmd, s3_path)
        cmd_list.append(cmd)

    return cmd_list
        
    
def s3_script_create(s3_cmds, script_name='aws_cmd_script.sh'):
    '''
    Creates a bash script containing AWS download or list commands.

    Args:
    -----
        s3_cmds : list of str
            List of AWS download (or list) commands.

        script_name : str
            Name for the script that will be created.

    Calling Sequence:
    -----------------
    >>> from gcpy import aws
    >>> s3_cmds = aws.list_cmds_from_log('GC_12.6.0.log')
    >>> aws.s3_script_create(s3_cmds, 's3_download_commands.sh')
    '''

    # ==================================================================
    # Initialization
    # ==================================================================

    # Error check arguments
    if not isinstance(s3_cmds, list):
        raise TypeError('The s3_cmds argument must be a list!')

    # ==================================================================
    # Write the script of bash commansds
    # ==================================================================

    # Open the script file for writing
    try:
        f = open(script_name, 'w+')
    except FileNotFoundError:
        raise FileNotFoundError(
            'Could not open {} for writing!'.format(script_name))

    # Write the shebang line and header comments
    print('#!/bin/bash\n', file=f)
    print('#', file=f, end='')
    print( '='*78, file=f)
    print('# This script was created by gcpy.aws.s3_script_create', file=f)
    print('#', file=f, end='')
    print( '='*78 + '\n', file=f)

    # Write the individual download commands to the file
    for cmd in s3_cmds:
        print('{}\n'.format(cmd), file=f)

    # Close the script file and make it executable
    f.close()
    os.chmod(script_name, 0o775)
    
