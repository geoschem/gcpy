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
        grid['lon'] % 360,
        coords={
            'nf': nf,
            'Ydim': Ydim,
            'Xdim': Xdim},
        dims=[
            'nf',
            'Ydim',
            'Xdim'])
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

    longitude = ds.longitude.values
    longitude[longitude > 180] -= 360
    
    # overpass time is early in the east and late in the west
    overpass_time_timedelta_min = -longitude / 360 * 24 * 60 + overpass_offset

    overhead_time = pd.to_datetime(args.overpass_time, format='%H:%M').time()
    ds['time'] = (overhead_time.hour + overhead_time.minute / 60 + overpass_time_timedelta_min / 60) % 24

    ds = ds.stack(track=['nf', 'Ydim', 'Xdim'])
    ds = ds.sortby('time')

    ds = ds.reset_index('track')
    ds = ds.assign_coords({'track': ds.time}).drop(
        'time').rename({'track': 'time'})

    ds['longitude'].attrs['long_name'] = 'longitude'
    ds['longitude'].attrs['units'] = 'degrees_east'
    ds['longitude'].values = ds['longitude'].values % 360

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
    def time_to_time_of_day(ds, numpy_time_unit='ns'):
        # Get mode date
        (unique_dates, date_counts) = np.unique(ds.time.astype('datetime64[D]'), return_counts=True)
        mode_date = unique_dates[np.argmax(date_counts)]
        time_of_day = (ds.time - mode_date).astype(f'timedelta64[{numpy_time_unit}]')
        return ds.assign_coords(time=time_of_day), mode_date

    def find_shift(tracked_output, track_file):
        output_time = tracked_output.time.values.astype(float) / 1e9 / 3600 % 24
        track_time = track_file.time.values.astype(float) / 1e9 / 3600 % 24

        def score(shift):
            return np.sum(np.abs(np.concatenate([output_time[shift:], output_time[:shift]]) - track_time))

        def has_converged(absdiff, curr_min, n=10):
            return len(absdiff) > n and np.all(absdiff[-n:] - curr_min > 0)

        absdiff = []
        shift = 0
        curr_min = score(shift)

        for shift in range(len(track_time)):
            absdiff.append(score(shift))
            curr_min = min(curr_min, absdiff[-1])
            if has_converged(absdiff, curr_min, n=20):
                break

        return np.argmin(absdiff)

    # Load track file and 1D output
    track = xr.open_dataset(args.track)
    tracked_output = xr.open_dataset(args.i)

    # Convert time coordinate to time of day
    track, _ = time_to_time_of_day(track)
    tracked_output, tracked_output_date = time_to_time_of_day(tracked_output)

    # The 1D output is shifted---find how much we need to roll backwards
    shift = find_shift(tracked_output, track)
    tracked_output = tracked_output.roll(dict(time=-shift))

    # Create a multiindex for nf,Ydim,Xdim
    track_mi = pd.MultiIndex.from_arrays([track.nf.values, track.Ydim.values, track.Xdim.values], names=['nf', 'Ydim', 'Xdim'])
    track = track.drop(['nf', 'Ydim', 'Xdim', 'latitude', 'longitude'])

    # Merge datasets
    tracked_output = tracked_output
    tracked_output = tracked_output.assign_coords({'time': track.time})
    tracked_output = xr.merge([tracked_output, track], compat='equals')

    # Unstack time coordinate (1D index)
    tracked_output = tracked_output.assign_coords({'time': track_mi})
    tracked_output = tracked_output.unstack('time')

    # Remake time coordinate
    tracked_output = tracked_output.expand_dims(dim='time')
    tracked_output.assign_coords({'time': [tracked_output_date]})

    tracked_output.to_netcdf(args.o)




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
