import argparse
import numpy as np
import xarray as xr
import pandas as pd
from gcpy.grid import make_grid_CS


def create_track_func(args):
    nf = np.linspace(1, 6, 6)
    Ydim = np.linspace(1, args.cs_res, args.cs_res)
    Xdim = np.linspace(1, args.cs_res, args.cs_res)

    grid, _ = make_grid_CS(args.cs_res)
    #grid, _ = make_grid_SG(args.cs_res, *args.sg_params)

    lon = xr.DataArray(
        grid['lon'] %
        360,
        coords={
            'nf': nf,
            'Ydim': Ydim,
            'Xdim': Xdim},
        dims=[
            'nf',
            'Ydim',
            'Xdim'])
    lon.values[lon.values > 180] -= 360
    lat = xr.DataArray(
        grid['lat'],
        coords={
            'nf': nf,
            'Ydim': Ydim,
            'Xdim': Xdim},
        dims=[
            'nf',
            'Ydim',
            'Xdim'])

    ds = xr.Dataset({'longitude': lon, 'latitude': lat})

    #ds['time'] = (((ds['longitude'] / 15) + args.overpass_time) + \
    #(ds['latitude'] / 180) * args.vertical_scan_time) % 24
    # vary overpass time with latitude
    overpass_offset = ds['latitude'] / 90 * 24 / args.orbits_per_day / 4 * 60
    if args.direction == 'ascending':
        # overpass delayed at high northern latitudes if ascending
        overpass_offset = -overpass_offset

    overpass_time_timedelta_min = ds['longitude'] / \
        360 * 24 * 60 + overpass_offset

    overhead_time = pd.to_datetime(args.overpass_time, format='%H:%M').time()
    ds['time'] = (overhead_time.hour + overhead_time.minute /
                  60 + overpass_time_timedelta_min / 60) % 24

    ds = ds.stack(track=['nf', 'Ydim', 'Xdim'])
    ds = ds.sortby('time')

    ds = ds.reset_index('track')
    ds = ds.assign_coords({'track': ds.time}).drop(
        'time').rename({'track': 'time'})

    ds['longitude'].attrs['long_name'] = 'longitude'
    ds['longitude'].attrs['units'] = 'degrees_east'

    ds['latitude'].attrs['long_name'] = 'latitude'
    ds['latitude'].attrs['units'] = 'degrees_north'

    ds2 = ds.drop(['nf', 'Ydim', 'Xdim'])
    ds2['nf'] = xr.DataArray(ds.nf.values, dims=['time'])
    ds2['Ydim'] = xr.DataArray(ds.Ydim.values, dims=['time'])
    ds2['Xdim'] = xr.DataArray(ds.Xdim.values, dims=['time'])

    ds2['time'].attrs['long_name'] = 'time'
    ds2['time'].attrs['units'] = 'hours since 1900-01-01 00:00:00'

    del ds

    encoding = {k: {'dtype': np.float32, 'complevel': 9, 'zlib': True}
                for k in ds2.variables}
    ds2.to_netcdf(args.o, encoding=encoding, format='NETCDF4_CLASSIC')


def unravel_func(args):
    track = xr.open_dataset(args.track)
    track_mi = pd.MultiIndex.from_arrays([track.nf, track.Ydim, track.Xdim])
    track = track.drop(['nf', 'Ydim', 'Xdim', 'latitude', 'longitude'])

    ds = xr.open_dataset(args.i)

    ds = ds.reindex({'time': track.time}, method='nearest')
    ds = xr.merge([ds, track], compat='no_conflicts')

    ds = ds.assign_coords({'time': track_mi})
    ds = ds.unstack('time')
    ds = ds.sortby(['nf', 'Ydim', 'Xdim'])

    ds.to_netcdf(args.o)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Utility for 1D diagnostics')
    subparsers = parser.add_subparsers(dest='command')

    create_track = subparsers.add_parser(
        'create_track', help='Generate satellite track file')
    create_track.add_argument('--cs_res',
                              metavar='RES',
                              type=int,
                              required=True,
                              help='grid\'s cubed-sphere resolution')
    create_track.add_argument(
        '--sg_params', metavar='P', type=float, nargs=3,
        default=[1.0, 170.0, -90.0],
        help='grid stretching parameters (stretch-factor, target longitude, target latitude)')
    create_track.add_argument('--overpass_time',
                              metavar='HH:MM',
                              type=str,
                              required=True,
                              help='Overpass time')
    create_track.add_argument('--orbits_per_day',
                              metavar='N',
                              type=int,
                              required=True,
                              help='number of orbits per day')
    create_track.add_argument('--direction',
                              type=str,
                              choices=['ascending', 'descending'],
                              required=True,
                              help='direction of orbit')
    create_track.add_argument('-o',
                              type=str,
                              required=True,
                              help='output filename')

    unravel = subparsers.add_parser(
        'unravel', help='Unravel 1D diagnostic file')
    unravel.add_argument('--track',
                         metavar='F',
                         type=str,
                         required=True,
                         help='track file')
    unravel.add_argument('-i',
                         metavar='F',
                         type=str,
                         required=True,
                         help='input file (1D diagnostic file)')
    unravel.add_argument('-o',
                         metavar='F',
                         type=str,
                         required=True,
                         help='output filename')

    args = parser.parse_args()
    if args.command is None:
        parser.print_help()
    elif args.command == 'create_track':
        create_track_func(args)
    elif args.command == 'unravel':
        unravel_func(args)
