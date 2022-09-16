# Example usage:
#
#   regrid_restart_file \
#       GEOSChem.Restart.fullchem.20190701_0000z.nc4 \
#       4x5_to_c24_weights.nc \
#       GCHP.Restart.fullchem.20190701_0000z.c24.nc4
#
#   regrid_restart_file \
#       http://ftp.as.harvard.edu/gcgrid/geos-chem/10yr_benchmarks/13.0.0/GCClassic/restarts/GEOSChem.Restart.20161101_0000z.nc4 \
#       4x5_to_c24_weights.nc \
#       GCHP.Restart.fullchem.20190701_0000z.c24.nc4

import os
import sys
import re
import logging
import tempfile
import xarray as xr
import numpy as np
import sparselt.esmf
import sparselt.xr
from pathlib import Path
import requests


temp_files=[]


def cleanup_tempfile():
    global temp_files
    if len(temp_files) > 0:
        logging.debug(f"Deleting {len(temp_files)} temp files")
    for filepath in temp_files:
        Path(filepath).unlink(missing_ok=True)
    temp_files = []


def is_gchp_restart_file(ds):
    is_gchp_restart = 'SPC_O3' in ds.data_vars
    is_gcclassic = 'SpeciesRst_O3' in ds.data_vars
    if not any((is_gchp_restart, is_gcclassic)):
        raise ValueError("Couldn't determine if the provided file is a GC-Classic or GCHP restart file.")
    return is_gchp_restart


def open_dataset(file_or_url, CHUNK_SIZE=8192):
    global temp_files
    is_url = bool(re.match(r'https?://', file_or_url))
    if is_url:
        logging.debug(f"Downloading {file_or_url}")
        with requests.get(file_or_url, stream=True) as r:
            r.raise_for_status()  # raise HTTPError
            tempfile_fd, tempfile_path = tempfile.mkstemp()
            with open(tempfile_fd, 'wb') as f:
                bytes_downloaded = 0
                for chunk in r.iter_content(chunk_size=CHUNK_SIZE):
                    bytes_downloaded += len(chunk)
                    f.write(chunk)
            temp_files.append(tempfile_path)
    file_path = tempfile_path if is_url else file_or_url
    logging.debug(f"Opening {file_path}")
    return xr.open_dataset(file_path)


def rename_variables(ds, to_gchp=True):
    to_gchp_re_sub = [
        (r'SpeciesRst_(.+)', r'SPC_\1'),
        (r'Met_(.+)', r'\1'),
        (r'Met_DELPDRY', r'DELP_DRY'),
        (r'Chem_(WetDepNitrogen|DryDepNitrogen|H2O2AfterChem|SO2AfterChem|KPPHvalue)', r'\1'),
    ]
    to_gcclassic_re_sub = [
        (r'SPC_(.+)', r'SpeciesRst_\1'),
        (r'(TropLev|BXHEIGHT)', r'Met_\1')
    ]
    re_sub_arg_list = to_gchp_re_sub if to_gchp else to_gcclassic_re_sub

    rename_dict = {}
    for re_sub_args in re_sub_arg_list:
        rename_dict.update({
            name: re.sub(*re_sub_args, name) for name in ds.data_vars if re.match(re_sub_args[0], name)
        })
    logging.info(f"Renaming {len(rename_dict)} variables")
    return ds.rename(rename_dict)


def reverse_lev(ds):
    logging.info(f"Reversing coordinate 'lev'")
    ds = ds.reindex(lev=ds.lev[::-1])
    ds = ds.assign_coords(lev=ds.lev.values[::-1])
    return ds


def drop_variables(ds, output_template):
    input_var_set = set(ds.data_vars)
    output_var_set = set(output_template.data_vars)
    drop_vars = input_var_set - output_var_set
    missing_vars = output_var_set - input_var_set
    if len(drop_vars) > 0:
        logging.info(f"Dropping {len(drop_vars)} variables from the input restart file that dont exist in the output template")
        logging.debug(f"Variables being dropped from the input restart file: {drop_vars}")
        ds = ds.drop(drop_vars)
    if len(missing_vars) > 0:
        logging.warning(f"The input restart file is missing {len(missing_vars)} variables that exist in the output template")
        logging.debug(f"Variables missing in the input restart file: {missing_vars}")
        output_template = output_template.drop(missing_vars)
    return ds, output_template


def regrid(ds, output_template, weights_file):
    weights = open_dataset(weights_file)
    input_dims = [('lat', 'lon'), (ds.dims['lat'], ds.dims['lon'])]

    output_template_shape = (output_template.dims['lat'], output_template.dims['lon'])
    resize_output_template = np.prod(output_template_shape) != weights.dst_grid_dims.item()
    if resize_output_template:
        if is_gchp_restart_file(output_template):
            # This is useful for stretched-grid simulations because they usually don't have a "normal" grid size    
            cs_res = np.sqrt(weights.dst_grid_dims.item() / 6).astype(int)
            logging.info(f"Reshaping the output restart file template to grid size C{cs_res}")
            output_shape = (6 * cs_res, cs_res)
            func = lambda *args, **kwargs: np.ones(output_shape)*np.nan
            vfunc = np.vectorize(func, signature='(lat,lon)->(lat1,lon1)')
            new_output_template = xr.apply_ufunc(
                vfunc, output_template, keep_attrs=True,
                input_core_dims=[['lat', 'lon']], output_core_dims=[['lat1', 'lon1']], 
            )
            new_output_template = new_output_template.rename({'lat1': 'lat', 'lon1': 'lon'})
            new_output_template['lat'].attrs = output_template['lat'].attrs
            new_output_template['lon'].attrs = output_template['lat'].attrs
            new_output_template = new_output_template.assign_coords(
                lat=np.arange(new_output_template.dims['lat'], dtype=np.float64),
                lon=np.arange(new_output_template.dims['lon'], dtype=np.float64),
            )
            output_template = new_output_template
        else:
            raise ValueError("GC-Classic restart resizing not implemented. Please provide a restart file template with the proper resolution.")
    else:
        output_shape = output_template_shape

    output_dims = [('lat', 'lon'), output_shape]
    logging.info("Regridding the input restart file")
    transform = sparselt.esmf.load_weights(weights, input_dims, output_dims)
    ds = sparselt.xr.apply(transform, ds, output_template)
    return ds


def update_encoding(ds):
    logging.info(f"Updating encoding")
    for name in ds.data_vars:
        ds[name].encoding.update({'dtype': 'float32'})
        if 'missing_value' in ds[name].encoding and '_FillValue' in ds[name].encoding:
            del ds[name].encoding['missing_value']
    return ds


def check_for_nans(ds):
    nan_vars = []
    for name in ds.data_vars:
        if ds[name].isnull().any().item():
            nan_vars.append(name)
    if len(nan_vars) > 0:
        logging.warning(f"Dataset has {len(nan_vars)}/{len(ds.data_vars)} variables with NaN values")
        logging.debug(f"Variables with NaN values: {nan_vars}")


def regrid_restart_file(input_restart, regrid_weights, output_restart_template):
    logging.info(f"Input restart file: {input_restart}")
    logging.info(f"Regridding weights: {regrid_weights}")
    logging.info(f"Output template restart file: {output_restart_template}")
        
    ds = open_dataset(input_restart)
    check_for_nans(ds)
    output_template = open_dataset(output_restart_template)

    input_is_gchp_restart = is_gchp_restart_file(ds)
    output_is_gchp_restart = is_gchp_restart_file(output_template)
    logging.info(f"Input restart file type is '{'GCHP' if input_is_gchp_restart else 'GC-Classic'}'")
    logging.info(f"Output restart file type is '{'GCHP' if output_is_gchp_restart else 'GC-Classic'}'")
    is_conversion = (input_is_gchp_restart != output_is_gchp_restart)
    if is_conversion:
        to_gchp = output_is_gchp_restart
        ds = rename_variables(ds, to_gchp)
        ds = reverse_lev(ds)

    ds, output_template = drop_variables(ds, output_template)
    ds = regrid(ds, output_template, weights_file=regrid_weights)
    ds = update_encoding(ds)
    check_for_nans(ds)
    ds.to_netcdf('new_restart_file.nc')
    logging.info(f"Wrote 'new_restart_file.nc' with {len(ds.data_vars)} variables")
    cleanup_tempfile()


if __name__ == '__main__':
    logging.basicConfig(level=os.environ.get('LOGLEVEL', 'INFO').upper())
    if len(sys.argv) != 4:
        logging.error("This program has 3 required arguments:  input_restart regrid_weights output_restart_template")
        exit(1)
    input_restart = sys.argv[1]
    regrid_weights = sys.argv[2]
    output_restart_template = sys.argv[3]
    regrid_restart_file(input_restart, regrid_weights, output_restart_template)
