#!/usr/bin/env python

'''
Python script to create a bash script with commands to download data
from the s3://gcgrid bucket to the ExtData folder in your AWS instance.
Only those files that are read by GEOS-Chem and HEMCO will be included.
This can be useful in creating a new AMI.
'''

# Imports
import gcpy.aws as aws

# Specify a log file from a GEOS-Chem simulation that was done
# elsewhere.  The file paths in this log file have all been
# renamed to /home/ubuntu/ExtData.
log_files = ['/path/to/GEOSChem/log', '/path/to/HEMCO/log']


# Set force_download to False if you would only like to copy those
# files from S3 that aren't already found on your EBS volume
# or local disk.  (This can save some time)
#
# Set force_download to True if you would like to download files
# from S3 even if they are already present on your EBS volume.
# (or local disk).
force_download = False

# Get a list of commands to download s3 data from the log file
s3_cmds = aws.s3_download_cmds_from_log(log_files, force_download)

# Write a bash script with the proper commands
script_name = 's3_download_script.sh'
aws.s3_script_create(s3_cmds, script_name)
