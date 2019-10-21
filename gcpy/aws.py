'''Specific utilities for GEOS-Chem and GCPy on the AWS Cloud'''

import os
from . import core


def s3_download_cmds_from_log(filename, copy_tmpl=None,
                              dst='ExtData', src='s3://gcgrid', **kwargs):
    '''
    '''
    filelist = core.extract_pathnames_from_log(filename, **kwargs)

    if copy_tmpl is None:
        copy_tmpl = 'aws s3 cp --request-payer=requester '

    cmd_list = []
    for f in filelist:

        # Extract the file name from the path
        basename = os.path.basename(f)

        # Create copy commands
        cmd = 'if [[ !(-f {a}/{b})]]; then {c}/{d} {e}/{f}; fi'.format(
            a=dst, b=f, c=copy_tmpl, d=f, e=dst, f=basename)
        cmd_list.append(cmd)

    return cmd_list
        
    '''
    alias s3copy="aws s3 cp --request-payer=requester "
    alias s3ls="aws s3 ls --request-payer=requester "
    '''
