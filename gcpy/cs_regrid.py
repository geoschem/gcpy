import argparse
import os.path

import numpy as np
import xarray as xr
try:
    import xesmf as xe
    from distutils.version import LooseVersion
    if LooseVersion(xe.__version__) < LooseVersion("0.2.1"):
        raise ImportError("cs_regrid.py requires xESMF version 0.2.1 or higher.")
except ImportError as e:
    print('cs_regrid.py requires xESMF version 0.2.1 or higher!\n\nSee the installation instructions here: https://xesmf.readthedocs.io/en/latest/installation.html\n')
import pandas as pd

from gcpy.grid import make_grid_SG
from gcpy.regrid import make_regridder_S2S, sg_hash, reformat_dims

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='General cubed-sphere to cubed-sphere regridder.')
    parser.add_argument('-i', '--filein',
                        metavar='FIN',
                        type=str,
                        required=True,
                        help='input NetCDF file')
    parser.add_argument('-o', '--fileout',
                        metavar='FOUT',
                        type=str,
                        required=True,
                        help='name of output file')
    parser.add_argument('--sg_params_in',
                        metavar='P',
                        type=float,
                        nargs=3,
                        default=[1.0, 170.0, -90.0],
                        help='input grid stretching parameters (stretch-factor, target longitude, target latitude)')
    parser.add_argument('--sg_params_out',
                        metavar='P',
                        type=float,
                        nargs=3,
                        default=[1.0, 170.0, -90.0],
                        help='output grid stretching parameters (stretch-factor, target longitude, target latitude)')
    parser.add_argument('--cs_res_out',
                        metavar='RES',
                        type=int,
                        required=True,
                        help='output grid\'s cubed-sphere resolution')
    parser.add_argument('--dim_format_in',
                        metavar='WHICH',
                        type=str,
                        choices=['checkpoint', 'diagnostic'],
                        required=True,
                        help='format of the input file\'s dimensions (choose from: checkpoint, diagnostic)')
    parser.add_argument('--dim_format_out',
                        metavar='WHICH',
                        type=str,
                        choices=['checkpoint', 'diagnostic'],
                        required=True,
                        help='format of the output file\'s dimensions (choose from: checkpoint, diagnostic)')
    args = parser.parse_args()

    # Load dataset
    ds_in = xr.open_dataset(args.filein, decode_cf=False)
    ds_in = ds_in.load()
    # Reformat dimensions to T, Z, F, Y, X
    ds_in = reformat_dims(ds_in, format=args.dim_format_in, towards_common=True)

    # Drop variables that don't look like fields
    non_fields = [v for v in ds_in.variables.keys() if len(set(ds_in[v].dims) - {'T', 'Z', 'F', 'Y', 'X'}) > 0
                  or len(ds_in[v].dims) == 0]
    ds_in = ds_in.drop(non_fields)

    # Transpose to T, Z, F, Y, X
    ds_in = ds_in.transpose('T', 'Z', 'F', 'Y', 'X')

    assert ds_in.dims['X'] == ds_in.dims['Y']
    cs_res_in = ds_in.dims['X']

    if cs_res_in == args.cs_res_out and all([v1 == v2 for v1, v2 in zip(args.sg_params_in, args.sg_params_out)]):
        print('Skipping regridding since grid parameters are identical')
        ds_out = ds_in
    else:
        # Make regridders
        regridders = make_regridder_S2S(
            cs_res_in, args.cs_res_out,
            sf_in=args.sg_params_in[0], tlon_in=args.sg_params_in[1], tlat_in=args.sg_params_in[2],
            sf_out=args.sg_params_out[0], tlon_out=args.sg_params_out[1], tlat_out=args.sg_params_out[2]
        )
        # For each output face, sum regridded input faces
        oface_datasets = []
        for oface in range(6):
            oface_regridded = []
            for iface, regridder in regridders[oface].items():
                ds_iface = ds_in.isel(F=iface)
                if 'F' in ds_iface.coords:
                    ds_iface = ds_iface.drop('F')
                oface_regridded.append(regridder(ds_iface, keep_attrs=True))
            oface_regridded = xr.concat(oface_regridded, dim='intersecting_ifaces').sum('intersecting_ifaces',
                                                                                        keep_attrs=True)
            oface_datasets.append(oface_regridded)
        ds_out = xr.concat(oface_datasets, dim='F')
        # Put regridded dataset back into a familiar format
        ds_out = ds_out.rename({
            'y': 'Y',
            'x': 'X',
        })
        ds_out = ds_out.drop(['lat', 'lon'])  # lat, lon are from xESMF which we don't want

    # Reformat dimensions to desired output format
    ds_out = reformat_dims(ds_out, format=args.dim_format_out, towards_common=False)
    
    # Store stretched-grid parameters as metadata
    ds_out.attrs['stretch_factor'] = args.sg_params_out[0]
    ds_out.attrs['target_longitude'] = args.sg_params_out[1]
    ds_out.attrs['target_latitude'] = args.sg_params_out[2]

    # Write dataset
    ds_out.to_netcdf(
        args.fileout,
        format='NETCDF4_CLASSIC'
    )
    # Print the resulting dataset
    print(ds_out)
