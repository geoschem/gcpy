import argparse
import os
import numpy as np
import xarray as xr
try:
    import xesmf as xe
    from distutils.version import LooseVersion
    if LooseVersion(xe.__version__) < LooseVersion("0.2.1"):
        raise ImportError(
            "file_regrid.py requires xESMF version 0.2.1 or higher.")
except ImportError as e:
    print('file_regrid.py requires xESMF version 0.2.1 or higher!\n\nSee the installation ' + \
          'instructions here: https://xesmf.readthedocs.io/en/latest/installation.html\n')
import pandas as pd

from gcpy.grid import get_input_res, get_vert_grid, get_grid_extents
from gcpy.regrid import make_regridder_S2S, reformat_dims, make_regridder_L2S, \
    make_regridder_C2L, make_regridder_L2L
from gcpy.util import reshape_MAPL_CS

def file_regrid(
        fin, fout, dim_format_in, dim_format_out, cs_res_out=0,
        ll_res_out='0x0', sg_params_in=[1.0, 170.0, -90.0],
        sg_params_out=[1.0, 170.0, -90.0], vert_params_out=[[], []]):
    """
    Regrids an input file to a new horizontal grid specification and saves it
    as a new file.

    Args:
        fin: str
            The input filename
        fout: str
            The output filename (file will be overwritten if it already exists)
        dim_format_in: str
            Format of the input file's dimensions (choose from: classic, 
            checkpoint, diagnostic), where classic denotes lat/lon and 
            checkpoint / diagnostic are cubed-sphere formats
        dim_format_out: str
            Format of the output file's dimensions (choose from: classic,
            checkpoint, diagnostic), where classic denotes lat/lon
            and checkpoint / diagnostic are cubed-sphere formats

    Keyword Args (optional):
        cs_res_out: int
            The cubed-sphere resolution of the output dataset. 
            Not used if dim_format_out is classic
            Default value: 0
        ll_res_out: str
            The lat/lon resolution of the output dataset. 
            Not used if dim_format_out is not classic
            Default value: '0x0'
        sg_params_in: list[float, float, float]
            Input grid stretching parameters 
            [stretch-factor, target longitude, target latitude].
            Not used if dim_format_in is classic
            Default value: [1.0, 170.0, -90.0] (No stretching)
        sg_params_out: list[float, float, float]
            Output grid stretching parameters 
            [stretch-factor, target longitude, target latitude].
            Not used if dim_format_out is classic
            Default value: [1.0, 170.0, -90.0] (No stretching)
        vert_params_out: list(list, list) of list-like types
            Hybrid grid parameter A in hPa and B (unitless) in [AP, BP] format.
            Needed for lat/lon output if not using full 72-level or 47-level grid
            Default value: [[], []]

    """

    # Load dataset
    ds_in = xr.open_dataset(fin, decode_cf=False)
    ds_in = ds_in.load()
    time = ds_in.time
    cs_res_in = 0
    if dim_format_in != 'classic':
        # Reformat dimensions to T, Z, F, Y, X
        ds_in = reformat_dims(ds_in, format=dim_format_in, towards_common=True)

        # Drop variables that don't look like fields
        non_fields = [
            v for v in ds_in.variables.keys()
            if len(set(ds_in[v].dims) - {'T', 'Z', 'F', 'Y', 'X'}) > 0 or
            len(ds_in[v].dims) == 0]
        ds_in = ds_in.drop(non_fields)

        # Transpose to T, Z, F, Y, X
        ds_in = ds_in.transpose('T', 'Z', 'F', 'Y', 'X')

        assert ds_in.dims['X'] == ds_in.dims['Y']
        cs_res_in = ds_in.dims['X']

    elif dim_format_in == 'classic' and dim_format_out != 'classic':
        ds_in = drop_and_rename_classic_vars(ds_in)

    # save type of data for later restoration
    original_dtype = np.dtype(ds_in[list(ds_in.data_vars)[0]])

    oface_files=[]
    if cs_res_in == cs_res_out and all(
            [v1 == v2 for v1, v2 in zip(sg_params_in, sg_params_out)]):
        print('Skipping regridding since grid parameters are identical')
        ds_out = ds_in

    elif dim_format_in != 'classic' and dim_format_out != 'classic':
        # CS/SG to CS/SG
        # Make regridders
        regridders = make_regridder_S2S(
            cs_res_in,
            cs_res_out,
            sf_in=sg_params_in[0],
            tlon_in=sg_params_in[1],
            tlat_in=sg_params_in[2],
            sf_out=sg_params_out[0],
            tlon_out=sg_params_out[1],
            tlat_out=sg_params_out[2])
        # Save temporary output face files to minimize RAM usage
        oface_files = [os.path.join('.',fout+str(x)) for x in range(6)]
        # For each output face, sum regridded input faces
        oface_datasets = []
        for oface in range(6):
            oface_regridded = []
            for iface, regridder in regridders[oface].items():
                ds_iface = ds_in.isel(F=iface)
                if 'F' in ds_iface.coords:
                    ds_iface = ds_iface.drop('F')
                oface_regridded.append(regridder(ds_iface, keep_attrs=True))
            oface_regridded = xr.concat(
                oface_regridded,
                dim='intersecting_ifaces').sum(
                'intersecting_ifaces',
                keep_attrs=True).expand_dims({'F':[oface]})
            oface_regridded.to_netcdf(
                oface_files[oface],
                format='NETCDF4_CLASSIC'
            )
        ds_out=xr.open_mfdataset(oface_files, combine='by_coords', concat_dim='F',engine='netcdf4')
        # Put regridded dataset back into a familiar format
        ds_out = ds_out.rename({
            'y': 'Y',
            'x': 'X',
        })
        # lat, lon are from xESMF which we don't want
        ds_out = ds_out.drop(['lat', 'lon'])

    elif dim_format_in == 'classic' and dim_format_out != 'classic':
        # LL to SG/CS
        llres_in = get_input_res(ds_in)[0]
        # make regridders
        regridders = make_regridder_L2S(
            llres_in, cs_res_out, sg_params=sg_params_out)
        ds_out = xr.concat([regridders[face](ds_in, keep_attrs=True)
                            for face in range(6)], dim='nf')
        # flip vertical
        ds_out = ds_out.reindex(lev=ds_out.lev[::-1])
        ds_out = ds_out.rename({
            'y': 'Ydim',
            'x': 'Xdim',
        })
        # lat, lon are from xESMF which we don't want
        ds_out = ds_out.drop(['lat', 'lon'])

        if dim_format_out == 'checkpoint':
            # convert to checkpoint format
            ds_out = reshape_MAPL_CS(ds_out)
            mi = pd.MultiIndex.from_product([
                np.linspace(1, 6, 6),
                np.linspace(1, cs_res_out, cs_res_out)
            ])
            ds_out = ds_out.assign_coords({'lat': mi})
            ds_out = ds_out.unstack('lat')

            ds_out = ds_out.stack(lat=['lat_level_0', 'lat_level_1'])
            ds_out = ds_out.assign_coords({
                'lat': np.linspace(1, 6 * cs_res_out, 6 * cs_res_out),
                'lon': np.linspace(1, ds_out.lon.size, ds_out.lon.size),
                'lev': np.linspace(ds_out.lev.size, 1, ds_out.lev.size)
            })
            ds_out = ds_out.transpose('time', 'lev', 'lat', 'lon')
        else:
            # convert to diagnostic format
            ds_out = ds_out.transpose('time', 'lev', 'nf', 'Ydim', 'Xdim')
            ds_out = ds_out.assign_coords({
                'nf': np.linspace(1, 6, 6),
                'lev': np.linspace(1, 72, 72)})
            print(
                'WARNING: xarray coordinates are not fully implemented for diagnostic format')

    elif dim_format_in != 'classic' and dim_format_out == 'classic':
        # SG/CS to LL
        regridders = make_regridder_C2L(
            cs_res_in, ll_res_out, sg_params=sg_params_in)
        ds_out = xr.concat(
            [regridders[face](ds_in.isel(F=face),
                              keep_attrs=True) for face in range(6)],
            dim='F').sum(
            'F', keep_attrs=True)
        ds_out = ds_out.rename({
            'T': 'time',
            'Z': 'lev'})
        ds_out = drop_and_rename_classic_vars(ds_out, towards_gchp=False)
        ds_out = ds_out.reindex(lev=ds_out.lev[::-1])
        _, lev_coords, _ = get_vert_grid(ds_out, *vert_params_out)
        ds_out = ds_out.assign_coords({'lev': lev_coords})
        ds_out['lat'].attrs = {'long_name': 'Latitude',
                               'units': 'degrees_north',
                               'axis': 'Y'}
        ds_out['lon'].attrs = {'long_name': 'Longitude',
                               'units': 'degrees_east',
                               'axis': 'X'}
    elif dim_format_in == 'classic' and dim_format_out == 'classic':
        # ll to ll
        in_extent = get_grid_extents(ds_in)
        out_extent = in_extent
        ll_res_in = get_input_res(ds_in)[0]
        [lat_in, lon_in] = list(map(float, ll_res_in.split('x')))
        [lat_out, lon_out] = list(map(float, ll_res_out.split('x')))

        if lat_in == lat_out and lon_in == lon_out:
            print('Skipping regridding since grid parameters are identical')
            ds_out = ds_in
        else:
            lon_attrs = ds_in.lon.attrs
            lat_attrs = ds_in.lat.attrs
            # drop non-regriddable variables
            non_fields = [v for v in ds_in.variables.keys(
            ) if 'lat' not in ds_in[v].dims and 'lon' not in ds_in[v].dims]
            non_fields_ds = ds_in[non_fields]
            ds_in = ds_in.drop(non_fields)

            regridder = make_regridder_L2L(
                ll_res_in,
                ll_res_out,
                reuse_weights=True,
                in_extent=in_extent,
                out_extent=out_extent)
            ds_out = regridder(ds_in, keep_attrs=True)
            ds_out = ds_out.merge(non_fields_ds)
            ds_out['lon'].attrs = lon_attrs
            ds_out['lat'].attrs = lat_attrs
            ds_out = ds_out.transpose('time', 'lev', 'ilev', 'lat', 'lon')

    if dim_format_in != 'classic' and dim_format_out != 'classic':
        # Reformat dimensions to desired output format
        ds_out = reformat_dims(
            ds_out,
            format=dim_format_out,
            towards_common=False)
    if dim_format_out != 'classic':
        # Store stretched-grid parameters as metadata
        ds_out.attrs['stretch_factor'] = sg_params_out[0]
        ds_out.attrs['target_longitude'] = sg_params_out[1]
        ds_out.attrs['target_latitude'] = sg_params_out[2]
        ds_out.attrs['cs_res'] = cs_res_out

    ds_out = ds_out.assign_coords({'time': time})
    # correct precision changes (accidental 32-bit to 64-bit)
    # save attributes (no longer needed in xarray >=0.16.1)
    attrs = ds_out.attrs
    data_attrs = {var : ds_out[str(var)].attrs for var in list(ds_out.variables)}
    ds_out = ds_out.astype(original_dtype)
    for var in list(ds_out.variables):
        ds_out[str(var)].attrs = data_attrs[var]
    ds_out.attrs = attrs
    # Write dataset
    ds_out.to_netcdf(
        fout,
        format='NETCDF4_CLASSIC'
    )
    # Print the resulting dataset
    print(ds_out)
    # Remove any temporary files
    for f in oface_files: os.remove(f)


def rename_restart_variables(ds, towards_gchp=True):
    """
    Renames restart variables according to GEOS-Chem Classic and GCHP conventions.

    Args:
        ds: xarray.Dataset
            The input dataset

    Keyword Args (optional):
        towards_gchp: bool
            Whether renaming to (True) or from (False) GCHP format
            Default value: True

    Returns:
        xarray.Dataset
            Input dataset with variables renamed
    """

    if towards_gchp:
        old_str = 'SpeciesRst'
        new_str = 'SPC'
    else:
        old_str = 'SPC'
        new_str = 'SpeciesRst'
    return ds.rename({name: name.replace(old_str, new_str, 1)
                      for name in list(ds.data_vars)
                      if name.startswith(old_str)})


def drop_and_rename_classic_vars(ds, towards_gchp=True):
    """
    Renames and drops certain restart variables according to GEOS-Chem Classic
    and GCHP conventions.

    Args:
        ds: xarray.Dataset
            The input dataset

    Keyword Args (optional):
        towards_gchp: bool
            Whether going to (True) or from (False) GCHP format
            Default value: True

    Returns:
        xarray.Dataset
            Input dataset with variables renamed and dropped
    """

    if towards_gchp:
        ds = ds.rename(
            {name: name.replace('Met_', '', 1).replace('Chem_', '', 1)
             for name in list(ds.data_vars)
             if name.startswith('Met_') or name.startswith('Chem_')})
        if 'DELPDRY' in list(ds.data_vars): ds = ds.rename({'DELPDRY': 'DELP_DRY'})
        ds = ds.drop_vars(['P0',
                           'hyam',
                           'hybm',
                           'hyai',
                           'hybi',
                           'AREA',
                           'ilev',
                           'PS1DRY',
                           'PS1WET',
                           'TMPU1',
                           'SPHU1',
                           'StatePSC'],
                          errors='ignore')
    else:
        renames = {'DELP_DRY': 'Met_DELPDRY',
                   'BXHEIGHT': 'Met_BXHEIGHT',
                   'TropLev': 'Met_TropLev',
                   'DryDepNitrogen': 'Chem_DryDepNitrogen',
                   'WetDepNitrogen': 'Chem_WetDepNitrogen',
                   'H2O2AfterChem': 'Chem_H2O2AfterChem',
                   'SO2AfterChem': 'Chem_SO2AfterChem',
                   'KPPHvalue': 'Chem_KPPHvalue'}
        data_vars = list(ds.data_vars)
        new_renames = renames.copy()
        for key in renames.keys():
            if key not in data_vars:
                del(new_renames[key])
        ds = ds.rename(new_renames)

    return rename_restart_variables(ds, towards_gchp=towards_gchp)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='General cubed-sphere to cubed-sphere regridder.')
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
    parser.add_argument(
        '--sg_params_in', metavar='P', type=float, nargs=3,
        default=[1.0, 170.0, -90.0],
        help='input grid stretching parameters (stretch-factor, target longitude, target latitude)')
    parser.add_argument(
        '--sg_params_out', metavar='P', type=float, nargs=3,
        default=[1.0, 170.0, -90.0],
        help='output grid stretching parameters (stretch-factor, target longitude, target latitude)')
    parser.add_argument('--cs_res_out',
                        metavar='RES',
                        type=int,
                        required=False,
                        help='output grid\'s cubed-sphere resolution')
    parser.add_argument(
        '--ll_res_out',
        metavar='RES',
        type=str,
        required=False,
        help='output grid\'s lat/lon resolution in \'latxlon\' format')
    parser.add_argument(
        '--dim_format_in',
        metavar='WHICH',
        type=str,
        choices=[
            'checkpoint',
            'diagnostic',
            'classic'],
        required=True,
        help='format of the input file\'s dimensions (choose from: checkpoint, diagnostic)')
    parser.add_argument(
        '--dim_format_out',
        metavar='WHICH',
        type=str,
        choices=[
            'checkpoint',
            'diagnostic',
            'classic'],
        required=True,
        help='format of the output file\'s dimensions (choose from: checkpoint, diagnostic)')
    parser.add_argument(
        '--vert_params_out',
        metavar='VERT',
        type=list,
        required=False,
        help='Hybrid grid parameter A in hPa and B (unitless) in [AP, BP] format')

    args = parser.parse_args()
    file_regrid(
        args.filein,
        args.fileout,
        args.dim_format_in,
        args.dim_format_out,
        args.cs_res_out,
        args.ll_res_out,
        args.sg_params_in,
        args.sg_params_out,
        args.vert_params_out)
