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

    # Make sure log_files has at least one element for iterating
    if len(log_files) == 1:
        log_files = [log_files]
    
    # Get the file paths from the log files
    paths = []
    for log_file in log_files:
        paths = paths + core.extract_pathnames_from_log(log_file,
                                                        prefix_filter)

    # Return a list of bash commands to download each file from S3
    cmd_list = []
    for path in paths:

        # Local path name and dir name
        local_path = prefix_filter + path
        local_dir = os.path.dirname(local_path)

        # S3 path name
        s3_path = '{}{}'.format(s3_root, path)

        # Create copy command and append to the list
        cmd = 'if [[ !(-f {}) ]]; then {} {} {}/; fi'.format(
            local_path, s3_cp_cmd, s3_path, local_dir)
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

    # Make sure log_files has at least one element for iterating
    if len(log_files) == 1:
        log_files = [log_files]
    
    # Get the file paths from the log files
    paths = []
    for log_file in log_files:
        paths = paths + core.extract_pathnames_from_log(log_file,
                                                        prefix_filter)

    # Return a list of bash commands to download each file from S3
    cmd_list = []
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

    # Open the script file for writing
    try:
        f = open(script_name, 'w+')
    except FileNotFoundError:
        print('Could not open: {}'.format(script_name))
        raise

    # Write the shebang line
    print('#!/bin/bash', file=f)
    print(file=f)

    # Write the individual download commands to the file
    for cmd in s3_cmds:
        print(cmd, file=f)
        print(file=f)

    # Close the script file and make it executable
    f.close()
    os.chmod(script_name, 0o775)
    
