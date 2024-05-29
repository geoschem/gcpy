#!/usr/bin/env python3
"""
Example to create an ObsPack input file with lon, lat, alt, etc coordinates
to specify the locations where GEOS-Chem should archive output.

Authors:
Alli Moon (GitHub: @alli-moon) and Yuk Chun Chan (GitHub: @yc-chan)
"""

import datetime as dt
import numpy as np
import xarray as xr
import pandas as pd

def main():
    """
    Creates an input file for the GEOS-Chem ObsPack diagnostic for a
    given station site location and time range.
    """

    # ==================================================================
    # Configurable settings -- edit for your own use case!
    DATASET_ID = 'GC-MODEL'
    OUTPUT_DIR = './'

    SAMPLING_DATES_SUMMER = pd.date_range(
        start='2022-05-31',
        end='2022-07-01',
        freq='d'
    )
    SAMPLING_DATES_WINTER = pd.date_range(
        start='2023-01-13',
        end='2023-2-13',
        freq='d'
    )

    SITE_LAT = 32.2644
    SITE_LON = -64.8767
    SITE_ALT = 0
    ASSUMED_INLET_HEIGHT_ABOVE_GROUND = 30 # Unit:m
    # ==================================================================

    for idx_date, i_date in enumerate(SAMPLING_DATES_SUMMER):
        i_num_location = 1
        sites_info_dict = {}
        for idx_site in range(i_num_location):
            sites_info_dict[idx_site] = {}
            sites_info_dict[idx_site]['lat'] = SITE_LAT
            sites_info_dict[idx_site]['lon'] = SITE_LON
            sites_info_dict[idx_site]['alt'] = SITE_ALT + ASSUMED_INLET_HEIGHT_ABOVE_GROUND
        i_num_obs = i_num_location*24
        i_coords = {
            'calendar_components': np.arange(6).astype('int8'),
            'obs': np.arange(i_num_obs).astype('int32'),
        }
        i_lat_arr = (np.ones(i_num_obs)*np.nan).astype('float32')
        i_lon_arr = (np.ones(i_num_obs)*np.nan).astype('float32')
        i_alt_arr = (np.ones(i_num_obs)*np.nan).astype('float32')
        i_time_arr = np.zeros([i_num_obs]).astype('int32')

        # Options: 2 = 1-hour avg; 4= instantaneous
        i_samplemethod = np.ones(i_num_obs).astype('int8')*2
        i_obspackid_arr = np.chararray(i_num_obs, itemsize=200)

        idx_obs_tdy = 0
        i_date_str = i_date.strftime('%Y-%m-%d')
        i_date_24hour_index = pd.date_range(
            start=i_date_str+' 00:30',
            end=i_date_str+' 23:30',
            freq='h'
        )
        i_date_in_dt_unit_arr = (
            i_date_24hour_index - dt.datetime(1970,1,1)
        ).total_seconds()
        for i_hour in range(24):
            for idx_site in sites_info_dict.keys():
                i_lat = sites_info_dict[idx_site]['lat']
                i_lon = sites_info_dict[idx_site]['lon']
                i_lat_arr[idx_obs_tdy] = i_lat
                i_lon_arr[idx_obs_tdy] = i_lon
                i_alt_arr[idx_obs_tdy] = sites_info_dict[idx_site]['alt']
                i_time_arr[idx_obs_tdy] = i_date_in_dt_unit_arr[i_hour]
                i_obspackid_arr[idx_obs_tdy] = f"{i_date.strftime('%Y-%m-%d')}_{i_hour:02}30_BermudaTudorHill_{DATASET_ID}".encode("utf-8") + b"\x00"
                idx_obs_tdy +=1

        lat_data_var = xr.DataArray(
            i_lat_arr,
            dims=['obs'],
            coords={'obs': i_coords['obs']},
            name='latitude',
        )
        lat_data_var.attrs['units'] = 'degrees_north'
        lat_data_var.attrs['_FillValue'] = -1.e+34
        lat_data_var.attrs['long_name'] = 'Sample latitude'

        lon_data_var = xr.DataArray(
            i_lon_arr,
            dims=['obs'],
            coords={'obs': i_coords['obs']},
            name='longitude'
        )
        lon_data_var.attrs['units'] = 'degrees_east'
        lon_data_var.attrs['_FillValue'] = -1.e+34
        lon_data_var.attrs['long_name'] = 'Sample latitude'

        alt_data_var = xr.DataArray(
            i_alt_arr,
            dims=['obs'],
            coords={'obs':i_coords['obs']},
            name='altitude'
        )
        alt_data_var.attrs['units'] = 'meters'
        alt_data_var.attrs['_FillValue'] = -1.e+34
        alt_data_var.attrs['long_name'] = 'sample altitude in meters above sea level'
        alt_data_var.attrs['comment'] = 'Altitude is surface elevation plus sample intake height in meters above sea level'

        obspack_id_data_var = xr.DataArray(
            i_obspackid_arr,
            dims=['obs'],
            coords={'obs': i_coords['obs']},
            name='obspack_id',
        )
        obspack_id_data_var.attrs['long_name'] = "Unique ObsPack observation id"
        obspack_id_data_var.attrs['comment'] = "Unique observation id string that includes obs_id, dataset_id and obspack_num."

        time_data_var = xr.DataArray(
            i_time_arr,
            dims=['obs'],
            coords={'obs':i_coords['obs']},
            name='time',
        )
        time_data_var.attrs['units'] = 'Seconds since 1970-01-01 00:00:00 UTC'
        time_data_var.attrs['_FillValue'] = -999999999
        time_data_var.attrs['long_name'] = 'Seconds since 1970-01-01 00:00:00 UTC'
        samplemethod_data_var = xr.DataArray(
            i_samplemethod,
            dims=['obs'],
            coords={"obs": i_coords['obs']},
            name='CT_sampling_strategy',
        )
        samplemethod_data_var.attrs['_FillValue'] = -9
        samplemethod_data_var.attrs['long_name'] = 'model sampling strategy'
        samplemethod_data_var.attrs['values'] = 'How to sample model. 1=4-hour avg; 2=1-hour avg; 3=90-min avg; 4=instantaneous'

        i_config_ds = xr.merge(
            [
                lat_data_var,
                lon_data_var,
                alt_data_var,
                time_data_var,
                obspack_id_data_var,
                samplemethod_data_var
            ]
        )

        i_config_ds = i_config_ds.assign_coords(i_coords)

        i_output_filename = f"ObsPack_config_{i_date.strftime('%Y%m%d')}.nc"
        i_output_filepath = OUTPUT_DIR+i_output_filename

        i_config_ds.to_netcdf(
            i_output_filepath,
            unlimited_dims='obs'
        )


if __name__ == '__main__':
    main()
