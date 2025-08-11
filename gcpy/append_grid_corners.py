#!/usr/bin/env python3
"""
Program to append corner longitudes and latitudes to an xarray dataset
contianing cubed-sphere "stretched-grid" data.
"""
import argparse
import xarray as xr
from gcpy.constants import DEFAULT_SG_PARAMS
from gcpy.grid import make_grid_sg

if __name__ == '__main__':

    # Define arguments to be parsed
    parser = argparse.ArgumentParser(
        description='Append grid-box corners to file'
    )
    parser.add_argument('filein')
    parser.add_argument(
        '--sg_params', metavar='P', type=float, nargs=3,
        default=DEFAULT_SG_PARAMS,
        help='input grid stretching parameters (stretch-factor, target longitude, target latitude)'
    )
    args = parser.parse_args()

    # Open netCDF file into an xarray Dataset
    ds = xr.open_dataset(args.filein)

    # Create stretched grid
    csgrid, csgrid_list = make_grid_sg(
        ds.dims['Ydim'],
        stretch_factor=args.sg_params[0],
        target_lon=args.sg_params[1],
        target_lat=args.sg_params[2]
    )

    # Assign the corner lons & lats to the Dataset as coordinates
    ds = ds.assign_coords({
        'corner_lons': xr.DataArray(
            csgrid['lon_b'],
            dims=['nf', 'YCdim', 'XCdim']
        ),
        'corner_lats': xr.DataArray(
            csgrid['lat_b'],
            dims=['nf', 'YCdim', 'XCdim']
        )
    })

    # Write to netCDF
    ds.load()
    ds.close()
    ds.to_netcdf(args.filein)

    print(ds)
