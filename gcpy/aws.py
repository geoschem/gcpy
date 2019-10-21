'''Specific utilities for GEOS-Chem and GCPy on the AWS Cloud'''

import os
from . import core


def s3_download_cmds_from_log(log_file_name,
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
        filename : str
            Name of the GEOS-Chem log file to read.

    Keyword Args (optional):
    ------------------------
        s3_cp_cmd : str
            Command to copy files from AWS S3.
            Default value: "aws s3 cp --request-payer=requester'

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

    # Get the file names from the log
    filelist = core.extract_pathnames_from_log(log_file_name, prefix_filter)

    # Return a list of bash commands
    cmd_list = []
    for f in filelist:
        cmd = 'if [[ !(-f {}{}) ]]; then {} {}/{} {}{}; fi'.format(
            prefix_filter, f, s3_cp_cmd, s3_root, f, prefix_filter, f)
        cmd_list.append(cmd)

    return cmd_list
        

def s3_list_cmds_from_log(log_file_name,
                          s3_ls_cmd='aws s3 ls --request-payer=requester',
                          src='s3://gcgrid',
                          prefix_filter='/home/ubuntu/ExtData'):
    '''
    Reads a GEOS-Chem log file and creates a list of bash commands to
    list data from the GEOS-Chem S3 bucket on AWS.

    Args:
    -----
        filename : str
            Name of the GEOS-Chem log file to read.

    Keyword Args (optional):
    ------------------------
        s3_ls_cmd : str
            Command to copy files from AWS S3.
            Default value: aws s3 ls --request-payer=requester

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

    # Get the file names from the log
    filelist = core.extract_pathnames_from_log(log_file_name, prefix_filter)

    # Return a list of bash commands
    cmd_list = []
    for f in filelist:
        cmd = 'if [[ !(-f {}{}) ]]; then {} {}/{} {}{}; fi'.format(
            prefix_filter, f, s3_cp_cmd, src, f, prefix_filter, f)
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
    f = open(script_name, 'w+')

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
    
