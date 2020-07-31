import argparse
import hashlib
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

from .grid import make_grid_SG


def sg_hash(cs_res, stretch_factor: float, target_lat: float, target_lon: float):
    return hashlib.sha1('cs={cs_res},sf={stretch_factor:.5f},tx={target_lon:.5f},ty={target_lat:.5f}'.format(
        stretch_factor=stretch_factor,
        target_lat=target_lat,
        target_lon=target_lon,
        cs_res=cs_res
    ).encode()).hexdigest()[:7]


def make_regridder_S2S(csres_in, csres_out, sf_in=1, tlat_in=-90, tlon_in=170, sf_out=1, tlat_out=-90, tlon_out=170, weightsdir='.'):
    igrid, igrid_list = make_grid_SG(csres_in, stretch_factor=sf_in, target_lat=tlat_in, target_lon=tlon_in)
    ogrid, ogrid_list = make_grid_SG(csres_out, stretch_factor=sf_out, target_lat=tlat_out, target_lon=tlon_out)
    regridder_list = []
    for o_face in range(6):
        regridder_list.append({})
        for i_face in range(6):
            weights_fname = f'conservative_sg{sg_hash(csres_in, sf_in, tlat_in, tlon_in)}_F{i_face}_sg{sg_hash(csres_out, sf_out, tlat_out, tlon_out)}_F{o_face}.nc'
            weights_file = os.path.join(weightsdir, weights_fname)
            reuse_weights = os.path.exists(weights_file)
            try:
                regridder = xe.Regridder(igrid_list[i_face],
                                         ogrid_list[o_face],
                                         method='conservative',
                                         filename=weights_file,
                                         reuse_weights=reuse_weights)
                regridder_list[-1][i_face] = regridder
            except ValueError:
                print(f"iface {i_face} doesn't intersect oface {o_face}")
    return regridder_list


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create a restart file for GCHP')
    parser.add_argument('-i', '--filein',
                        metavar='FILEIN',
                        type=str,
                        required=True,
                        help='input file')
    parser.add_argument('-o', '--fileout',
                        metavar='RES',
                        type=str,
                        required=True,
                        help='output file name')
    parser.add_argument('--sg_params_in',
                        metavar='SF X Y',
                        type=float,
                        nargs=3,
                        default=[1.0, 170.0, -90.0],
                        help='input grid stretched-grid definition')
    parser.add_argument('--sg_params_out',
                        metavar='SF X Y',
                        type=float,
                        nargs=3,
                        default=[1.0, 170.0, -90.0],
                        help='output grid stretched-grid definition')
    parser.add_argument('--cs_res_out',
                        metavar='RES',
                        type=int,
                        required=True,
                        help='output grid cubed-sphere resolution')
    parser.add_argument('--stacked_in',
                        action='store_true',
                        help='output file in stacked format')
    parser.add_argument('--stacked_out',
                        action='store_true',
                        help='output file in stacked format')
    args = parser.parse_args()

    ds_in = xr.open_dataset(args.filein, decode_cf=False)

    if args.stacked_in:
        cs_res_in = ds_in.dims['lat'] // 6
        y = np.linspace(1, cs_res_in, cs_res_in)
        nf = np.linspace(1, 6, 6)
        mi = pd.MultiIndex.from_product([nf, y])
        ds_in = ds_in.assign_coords({'lat': mi})
        ds_in = ds_in.unstack('lat').rename({
            'lat_level_0': 'face',
            'lat_level_1': 'Y',
            'lon': 'X',
        })
        ds_in = ds_in.transpose('time', 'lev', 'face', 'Y', 'X')
    else:
        cs_res_in = ds_in.dims['Ydim']
        ds_in = ds_in.unstack('lat').rename({
            'nf': 'face',
            'Ydim': 'Y',
            'Xdim': 'X',
        })
        ds_in = ds_in.transpose('time', 'lev', 'face', 'Y', 'X')

    ds_in['conservative_check'] = xr.DataArray(
        np.ones((ds_in.dims['time'], ds_in.dims['lev'], ds_in.dims['face'], ds_in.dims['Y'], ds_in.dims['X'])),
        dims=['time', 'lev', 'face', 'Y', 'X']
    )

    regridders = make_regridder_S2S(
        cs_res_in, args.cs_res_out,
        sf_in=args.sg_params_in[0], tlon_in=args.sg_params_in[1], tlat_in=args.sg_params_in[2],
        sf_out=args.sg_params_out[0], tlon_out=args.sg_params_out[1], tlat_out=args.sg_params_out[2]
    )

    oface_datasets = []
    for oface in range(6):
        oface_regridded = [regridder(ds_in.isel(face=iface).drop('face'), keep_attrs=True) for iface, regridder in regridders[oface].items()]
        oface_regridded = xr.concat(oface_regridded, dim='intersecting_ifaces').sum('intersecting_ifaces')
        oface_datasets.append(oface_regridded)

    ds_out = xr.concat(oface_datasets, dim='face')

    ds_out = ds_out.rename({
        'face': 'nf',
        'y': 'Ydim',
        'x': 'Xdim',
    })
    ds_out = ds_out.transpose('time', 'lev', 'nf', 'Ydim', 'Xdim')
    ds_out = ds_out.drop(['lat', 'lon'])

    if args.stacked_out:
        cs_res_out = ds_out.dims['Xdim']
        ds_out = ds_out.stack(lat=['nf', 'Ydim'])
        ds_out = ds_out.rename({'Xdim': 'lon'})
        ds_out = ds_out.assign_coords({
            'lat': np.linspace(1, 6*cs_res_out, 6*cs_res_out), 'lon': np.linspace(1, cs_res_out, cs_res_out)
        })
        ds_out = ds_out.transpose('time', 'lev', 'lat', 'lon')
    else:
        cs_res_out = ds_out.dims['Xdim']
        ds_out = ds_out.assign_coords({
            'nf': np.linspace(1, 6, 6),
            'Ydim': np.linspace(1, cs_res_out, cs_res_out),
            'Xdim': np.linspace(1, cs_res_out, cs_res_out),
        })
        ds_out = ds_out.transpose('time', 'lev', 'nf', 'Ydim', 'Xdim')

    # Review conservative check
    score = 1 + (1 - abs(ds_out['conservative_check']))
    print('Conservative check:')
    print(f'    - P95 score: {score.quantile(0.95).item()}')
    print(f'    - P99 score: {score.quantile(0.99).item()}')
    ds_out = ds_out.drop(['conservative_check'])

    # Write dataset
    ds_out.to_netcdf(
        args.fileout,
        format='NETCDF4_CLASSIC'
    )

    print(ds_out)