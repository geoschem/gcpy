import argparse
import xarray as xr
from .grid import make_grid_SG

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Append grid-box corners to file')
    parser.add_argument('filein')
    parser.add_argument(
        '--sg_params', metavar='P', type=float, nargs=3,
        default=[1.0, 170.0, -90.0],
        help='input grid stretching parameters (stretch-factor, target longitude, target latitude)')
    args = parser.parse_args()

    ds = xr.open_dataset(args.filein)

    csgrid, csgrid_list = make_grid_SG(
        ds.dims['Ydim'],
        stretch_factor=args.sg_params[0],
        target_lon=args.sg_params[1],
        target_lat=args.sg_params[2]
    )

    ds = ds.assign_coords(
        {'corner_lons': xr.DataArray(
            csgrid['lon_b'],
            dims=['nf', 'YCdim', 'XCdim']),
         'corner_lats': xr.DataArray(
             csgrid['lat_b'],
             dims=['nf', 'YCdim', 'XCdim'])})
    ds.load()
    ds.close()
    ds.to_netcdf(args.filein)

    print(ds)
