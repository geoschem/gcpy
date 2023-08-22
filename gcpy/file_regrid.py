"""
GCPy regridding utility: Regrids data between lat/lon, cubed-sphere,
and cubed-sphere stretched grids.
"""
import argparse
import os
import warnings
import numpy as np
import xarray as xr
try:
    import xesmf as xe
    from distutils.version import LooseVersion
    if LooseVersion(xe.__version__) < LooseVersion("0.2.1"):
        raise ImportError(
            "file_regrid.py requires xESMF version 0.2.1 or higher.")
except ImportError as exc:
    print('file_regrid.py requires xESMF version 0.2.1 or higher!\n\nSee the installation ' + \
          'instructions here: https://xesmf.readthedocs.io/en/latest/installation.html\n')

from gcpy.grid import get_input_res, get_vert_grid, get_grid_extents
from gcpy.regrid import make_regridder_S2S, reformat_dims, \
    make_regridder_L2S, make_regridder_C2L, make_regridder_L2L
from gcpy.util import verify_variable_type

# Ignore any FutureWarnings
warnings.simplefilter(action='ignore', category=FutureWarning)


def file_regrid(
        fin,
        fout,
        dim_format_in,
        dim_format_out,
        cs_res_out=0,
        ll_res_out='0x0',
        sg_params_in=None,
        sg_params_out=None,
        vert_params_out=None
):
    """
    Regrids an input file to a new horizontal grid specification
    and saves it as a new file.

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
    verify_variable_type(fin, str)
    verify_variable_type(fout, str)
    verify_variable_type(dim_format_in, str)
    verify_variable_type(dim_format_out, str)

    # Assign default values for optional keywords
    if sg_params_in is None:
        sg_params_in = [1.0, 170.0, -90.0]
    if sg_params_out is None:
        sg_params_out = [1.0, 170.0, -90.0]
    if vert_params_out is None:
        vert_params_out = [[], []]

    # Load dataset
    ds_in = xr.open_dataset(
        fin,
        decode_cf=False,
        engine='netcdf4'
    ).load()
    cs_res_in = 0

    # Make sure all xarray.Dataset global & variable attributes are kept
    with xr.set_options(keep_attrs=True):

        # ==============================================================
        # Regrid data
        # ==============================================================

        # save type of data for later restoration
        original_dtype = np.dtype(ds_in[list(ds_in.data_vars)[0]])

        oface_files=[]
        if cs_res_in == cs_res_out and all(
                [v1 == v2 for v1, v2 in zip(sg_params_in, sg_params_out)]):

            # ----------------------------------------------------------
            # Input CS/SG grid == Output CS/SG grid
            # ----------------------------------------------------------
            print('Skipping regridding since grid parameters are identical')
            ds_out = ds_in

        elif dim_format_in != 'classic' and dim_format_out != 'classic':

            # ----------------------------------------------------------
            # Input grid is CS/SG; Output grid is CS/SG
            # ----------------------------------------------------------
            ds_out = regrid_cssg_to_cssg(
                fout,
                ds_in,
                dim_format_in,
                sg_params_in,
                cs_res_out,
                dim_format_out,
                sg_params_out
            )

        elif dim_format_in == 'classic' and dim_format_out != 'classic':

            # ----------------------------------------------------------
            # Input grid is LL; Output grid is CS/SG
            # ----------------------------------------------------------
            ds_out = regrid_ll_to_cssg(
                ds_in,
                cs_res_out,
                dim_format_out,
                sg_params_out
            )

        elif dim_format_in != 'classic' and dim_format_out == 'classic':

            # ----------------------------------------------------------
            # Input grid is CS/SG; Output grid is LL
            # ----------------------------------------------------------
            ds_out = regrid_cssg_to_ll(
                ds_in,
                cs_res_in,
                dim_format_in,
                sg_params_in,
                ll_res_out,
                vert_params_out=vert_params_out
            )

        elif dim_format_in == 'classic' and dim_format_out == 'classic':

            # ----------------------------------------------------------
            # Input grid is LL; Output grid is LL
            # ----------------------------------------------------------
            ds_out = regrid_ll_to_ll(
                ds_in,
                ll_res_out
            )

        # ==============================================================
        # Post-regridding stuff
        # ==============================================================

        # Correct precision changes (accidental 32-bit to 64-bit)
        ds_out = ds_out.astype(original_dtype)

        # Write dataset to file
        ds_out.to_netcdf(
            fout,
            format='NETCDF4'
        )

        # Print the resulting dataset
        print(ds_out)


def prepare_cssg_input_grid(
        ds_in,
        dim_format_in
):
    """
    Reformats cubed-sphere/stretched grid data to the universal
    format and drops non-regriddable fields.

    Args:
    -----
    ds_in : xr.Dataset
        Input grid (cubed-sphere or stretched grid)
    dim_format_in : str
        Either "checkpoint" (for restart files)
        or "diagnostic" (for History diagnostic files)

    Returns:
    --------
    ds_out : xr.Dataset
        Data with reformatted dimensions and dropped fields
    cs_res_in : int
        Cubed-sphere/stretched grid resolution
    """
    # Reformat dimensions to "common dimensions (T, Z, F, Y, X)
    ds_in = reformat_dims(
        ds_in,
        dim_format_in,
        towards_common=True
    )

    # Drop variables that don't look like fields
    non_fields = [
        v for v in ds_in.variables.keys()
        if len(set(ds_in[v].dims) - {'T', 'Z', 'F', 'Y', 'X'}) > 0
        or len(ds_in[v].dims) == 0]
    ds_in = ds_in.drop(non_fields)

    # Transpose to T, Z, F, Y, X
    ds_in = ds_in.transpose('T', 'Z', 'F', 'Y', 'X')

    assert ds_in.dims['X'] == ds_in.dims['Y']
    cs_res_in = ds_in.dims['X']

    return ds_in, cs_res_in


def regrid_cssg_to_cssg(
        fout,
        ds_in,
        dim_format_in,
        sg_params_in,
        cs_res_out,
        dim_format_out,
        sg_params_out,
):
    """
    Regrids from the cubed-sphere/stretched grid to a different
    cubed-sphere/stretched grid resolution.

    Args:
    -----
    fout : str
        File name template
    ds_in : xarray.Dataset
        Data on a cubed-sphere/stretched grid
    dim_format_in, dim_format_out : str
        Input & output grid format ("checkpoint", "diagnostic")
    cs_res_out : int
        Cubed-sphere grid resolution
    sg_params_in, sg_params_out: list[float, float, float]
        Input & output grid stretching parameters
        [stretch-factor, target longitude, target latitude].

    Returns:
    --------
    ds_out : xarray.Dataset
        Data regridded to the output lat-lon grid
    """
    with xr.set_options(keep_attrs=True):

        # Change CS/SG dimensions to universal format
        # and drop non-regriddable variables
        ds_in, cs_res_in = prepare_cssg_input_grid(
            ds_in,
            dim_format_in
        )

        # Make regridders
        regridders = make_regridder_S2S(
            cs_res_in,
            cs_res_out,
            sf_in=sg_params_in[0],
            tlon_in=sg_params_in[1],
            tlat_in=sg_params_in[2],
            sf_out=sg_params_out[0],
            tlon_out=sg_params_out[1],
            tlat_out=sg_params_out[2]
        )

        # Save temporary output face files to minimize RAM usage
        oface_files = [os.path.join('.',fout+str(x)) for x in range(6)]

        # For each output face, sum regridded input faces
        for oface in range(6):
            oface_regridded = []
            for (iface, regridder) in regridders[oface].items():
                ds_iface = ds_in.isel(F=iface)
                if 'F' in ds_iface.coords:
                    ds_iface = ds_iface.drop('F')
                oface_regridded.append(
                    regridder(
                        ds_iface,
                        keep_attrs=True
                    )
                )
            oface_regridded = xr.concat(
                oface_regridded,
                dim='intersecting_ifaces'
            ).sum(
                'intersecting_ifaces',
                keep_attrs=True).expand_dims({'F':[oface]})
            oface_regridded.to_netcdf(
                oface_files[oface],
                format='NETCDF4'
            )

        # Combine face files
        ds_out=xr.open_mfdataset(
            oface_files,
            combine='nested',
            concat_dim='F',
            engine='netcdf4'
        )

        # lat, lon are from xESMF which we don't want
        ds_out = drop_lon_and_lat(ds_out)

        # Put regridded dataset back into a familiar format
        ds_out = ds_out.rename({
            'y': 'Y',
            'x': 'X',
        })

        # Reformat dimensions from "common dimension format"
        # to CS/GG "checkpoint" or "diagnostics" format
        ds_out = reformat_dims(
            ds_out,
            dim_format_out,
            towards_common=False
        )

        # Save stretched-grid metadata
        ds_out = save_cssg_metadata(
            ds_out,
            cs_res_out,
            sg_params_out
        )

        # Remove any temporary files
        for oface in oface_files:
            os.remove(oface)

    return ds_out


def regrid_cssg_to_ll(
        ds_in,
        cs_res_in,
        dim_format_in,
        sg_params_in,
        ll_res_out,
        vert_params_out=None,
):
    """
    Regrids from the cubed-sphere/stretched grid to the lat-lon grid.

    Args:
    -----
    ds_in : xarray.Dataset
        Data on a cubed-sphere/stretched grid
    cs_res_in : int
        Cubed-sphere grid resolution
    sg_params_in: list[float, float, float]
        Input grid stretching parameters
        [stretch-factor, target longitude, target latitude].
    ll_res_out : str
        Output grid lat/lon resolution (e.g. "4x5")

    Returns:
    --------
    ds_out : xarray.Dataset
        Data regridded to the output lat-lon grid
    """
    if vert_params_out is None:
        vert_params_out = [[], []]

    with xr.set_options(keep_attrs=True):

        # Change CS/SG dimensions to universal format
        # and drop non-regriddable variables
        ds_in, cs_res_in = prepare_cssg_input_grid(
            ds_in,
            dim_format_in
        )

        # Regrid data
        regridders = make_regridder_C2L(
            cs_res_in,
            ll_res_out,
            sg_params=sg_params_in
        )
        ds_out = xr.concat(
            [regridders[face](ds_in.isel(F=face), keep_attrs=True)
             for face in range(6)],
            dim='F'
        ).sum('F', keep_attrs=True)

        # Update dimensions and attributes on the lat-lon grid
        ds_out = ds_out.rename({'T': 'time', 'Z': 'lev'})
        ds_out = drop_and_rename_classic_vars(ds_out, towards_gchp=False)
        ds_out = ds_out.reindex(lev=ds_out.lev[::-1])
        _, lev_coords, _ = get_vert_grid(ds_out, vert_params_out)
        ds_out = ds_out.assign_coords({'lev': lev_coords})
        ds_out['lat'].attrs = {
            'long_name': 'Latitude',
            'units': 'degrees_north',
            'axis': 'Y'
        }
        ds_out['lon'].attrs = {
            'long_name': 'Longitude',
            'units': 'degrees_east',
            'axis': 'X'
        }

    return ds_out


def regrid_ll_to_cssg(
        ds_in,
        cs_res_out,
        dim_format_out,
        sg_params_out,
):
    """
    Regrids from the lat-lon grid to the cubed-sphere/stretched grid.

    Args:
    -----
    ds_in : xarray.Dataset
        Data on a lat/lon grid
    cs_res_in : int
        Cubed-sphere grid resolution
    dim_format_out : str
        Either "checkpoint" (for restart files) or
        "diagnostic" (for History diagnostic files).
    sg_params_out: list[float, float, float]
        Output grid stretching parameters
        [stretch-factor, target longitude, target latitude].

    Returns:
    --------
    ds_out : xarray.Dataset
        Data regridded to the output cubed-sphere/stretched-grid
    """
    with xr.set_options(keep_attrs=True):

        # Drop non-regriddable variables when going from ll -> cs
        ds_in = drop_and_rename_classic_vars(ds_in)

        # Input lat/lon grid resolution
        llres_in = get_input_res(ds_in)[0]

        # Regrid data to CS/SG
        regridders = make_regridder_L2S(
            llres_in,
            cs_res_out,
            sg_params=sg_params_out
        )
        ds_out = xr.concat(
            [regridders[face](ds_in, keep_attrs=True) for face in range(6)],
            dim='nf'
        )

        # Flip vertical levels
        ds_out = ds_out.reindex(lev=ds_out.lev[::-1])

        # Drop lon & lat, which are from xESMF
        ds_out = drop_lon_and_lat(ds_out)

        # Rename dimensions to the "common dimension format"
        ds_out = ds_out.rename({
            'time': 'T',
            'lev': 'Z',
            'nf': 'F',
            'y': 'Y',
            'x': 'X',
        })

        # Reformat dimensions from "common dimension format"
        # to CS/GG "checkpoint" or "diagnostics" format
        ds_out = reformat_dims(
            ds_out,
            dim_format_out,
            towards_common=False
        )

        # Save stretched-grid metadata
        ds_out = save_cssg_metadata(
            ds_out,
            cs_res_out,
            sg_params_out
        )

    return ds_out


def regrid_ll_to_ll(
        ds_in,
        ll_res_out
):
    """
    Regrid from the lat/lon grid to the cubed-sphere/stretched grid.

    Args:
    -----
    ds_in : xarray.Dataset
        Data on a lat/lon grid
    ll_res_out : str
        Output grid lat-lon grid resolution (e.g. "4x5")

    Returns:
    --------
    ds_out : xarray.Dataset
        Data regridded to the output lat-lon grid.
    """
    # Keep all xarray global & variable attributes
    with xr.set_options(keep_attrs=True):

        # Get the input & output extents
        in_extent = get_grid_extents(ds_in)
        out_extent = in_extent
        ll_res_in = get_input_res(ds_in)[0]
        [lat_in, lon_in] = list(map(float, ll_res_in.split('x')))
        [lat_out, lon_out] = list(map(float, ll_res_out.split('x')))

        # Return if the output & input grids are the same
        if lat_in == lat_out and lon_in == lon_out:
            ds_out = ds_in
            return ds_out

        # Drop non-regriddable variables
        non_fields = [
            var for var in ds_in.variables.keys() \
            if 'lat' not in ds_in[var].dims       \
                and 'lon' not in ds_in[var].dims
            ]
        ds_in = ds_in.drop(["lat_bnds", "lon_bnds"])
        non_fields_ds = ds_in[non_fields]
        ds_in = ds_in.drop(non_fields)

        # Create the regridder
        regridder = make_regridder_L2L(
        ll_res_in,
            ll_res_out,
            reuse_weights=True,
            in_extent=in_extent,
            out_extent=out_extent
        )
        ds_out = regridder(ds_in, keep_attrs=True, verbose=False)

        # Add the non-regriddable fields back
        ds_out = ds_out.merge(non_fields_ds)

        # Change order of dimensions
        ds_out = ds_out.transpose(
            'time', 'lev', 'ilev', 'lat', 'lon', ...
        )

    return ds_out


def save_cssg_metadata(
        dset,
        cs_res_out,
        sg_params_out
):
    """
    Saves the stretched-grid metadata to an xarray.Dataset object
    containing cubed-sphere/stretched grid data

    Args:
    -----
    dset : xarray.Dataset
        Data on the stretched grid.
    cs_res_out : int
        Cubed-sphere grid resolution.
    sg_params_out: list[float, float, float]
        Output grid stretching parameters
        [stretch-factor, target longitude, target latitude]

    Returns:
    --------
    dset_out : xarray.Dataset
        The original data, plus stretched grid metadata.
    """
    with xr.set_options(keep_attrs=True):
        dset.attrs['stretch_factor'] = sg_params_out[0]
        dset.attrs['target_longitude'] = sg_params_out[1]
        dset.attrs['target_latitude'] = sg_params_out[2]
        dset.attrs['cs_res'] = cs_res_out

    return dset


def rename_restart_variables(dset, towards_gchp=True):
    """
    Renames restart variables according to GEOS-Chem Classic
    and GCHP conventions.

    Args:
        dset: xarray.Dataset
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
    return dset.rename({name: name.replace(old_str, new_str, 1)
                      for name in list(dset.data_vars)
                      if name.startswith(old_str)})


def drop_lon_and_lat(dset):
    """
    Drops lon and lat variables, which are added by xESMF.
    These are not needed for GCHP data files.

    Args:
        dset: xarray.Dataset
            The input data.

    Returns:
        dset: xarray.Dataset
            The input data minus "lat" and "lon" variables.
    """
    verify_variable_type(dset, xr.Dataset)

    with xr.set_options(keep_attrs=True):
        if "lat" in dset.variables.keys():
            dset = dset.drop("lat")
        if "lon" in dset.variables.keys():
            dset = dset.drop("lon")

    return dset


def drop_and_rename_classic_vars(dset, towards_gchp=True):
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
    with xr.set_options(keep_attrs=True):
        if towards_gchp:
            dset = dset.rename(
                {name: name.replace('Met_', '', 1).replace('Chem_', '', 1)
                 for name in list(dset.data_vars)
                 if name.startswith('Met_') or name.startswith('Chem_')})
            if 'DELPDRY' in list(dset.data_vars):
                dset = dset.rename({'DELPDRY': 'DELP_DRY'})
            dset = dset.drop_vars(
                ['P0',
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
                 'StatePSC',
                 'lon_bnds',
                 'lat_bnds'
                 ],
                errors='ignore'
            )
        else:
            renames = {
                'DELP_DRY': 'Met_DELPDRY',
                'BXHEIGHT': 'Met_BXHEIGHT',
                'TropLev': 'Met_TropLev',
                'DryDepNitrogen': 'Chem_DryDepNitrogen',
                'WetDepNitrogen': 'Chem_WetDepNitrogen',
                'H2O2AfterChem': 'Chem_H2O2AfterChem',
                'SO2AfterChem': 'Chem_SO2AfterChem',
                'KPPHvalue': 'Chem_KPPHvalue'
            }
            data_vars = list(dset.data_vars)
            new_renames = renames.copy()
            for items in renames.items():
                if items[0] not in data_vars:
                    del(new_renames[items[0]])
            dset = dset.rename(new_renames)

        return rename_restart_variables(dset, towards_gchp=towards_gchp)


def main():
    """
    Main program for file_regrid.  Parses command-line arguments and
    calls the file_regrid routine.

    Command-line arguments:
    -----------------------
    -i, --filein
        Input file, contains original data.

    -o --fileout
        Output file, contains regridded data.

    --sg-params-in
        Input grid stretching parameters (GCHP only).

    --sg-params-out
        Output grid stretching parameters (GCHP only).

    --dim-format-in
        Format of the input file's dimensions:
        ("checkpoint", "diagnostics". "classic")

    --dim-format-out
        Format of the output file's dimensions:
        ("checkpoint", "diagnostics", "classic")

    --cs_res_out
        Cubed-sphere resolution for the output file (e.g 24, 48, 360)

    --ll_res_out
        Resolution for the output file in 'latxlon` format

    --vert_params_out
        Hybrid grid parameter A in hPa and B (unitless) in [AP, BP] format
    """

    # Parse arguments from the command line
    parser = argparse.ArgumentParser(
        description='General cubed-sphere to cubed-sphere regridder.'
    )
    parser.add_argument(
        '-i', '--filein',
        metavar='FIN',
        type=str,
        required=True,
        help='input NetCDF file'
    )
    parser.add_argument(
        '-o', '--fileout',
        metavar='FOUT',
        type=str,
        required=True,
        help='name of output file'
    )
    parser.add_argument(
        '--sg_params_in',
        metavar='P',
        type=float,
        nargs=3,
        default=[1.0, 170.0, -90.0],
        help='input grid stretching parameters (stretch-factor, target longitude, target latitude)'
    )
    parser.add_argument(
        '--sg_params_out',
        metavar='P',
        type=float,
        nargs=3,
        default=[1.0, 170.0, -90.0],
        help='output grid stretching parameters (stretch-factor, target longitude, target latitude)'
    )
    parser.add_argument(
        '--cs_res_out',
        metavar='RES',
        type=int,
        required=False,
        help='output grid\'s cubed-sphere resolution'
    )
    parser.add_argument(
        '--ll_res_out',
        metavar='RES',
        type=str,
        required=False,
        help='output grid\'s lat/lon resolution in \'latxlon\' format'
    )
    parser.add_argument(
        '--dim_format_in',
        metavar='WHICH',
        type=str,
        choices=[
            'checkpoint',
            'diagnostic',
            'classic'],
        required=True,
        help='format of the input file\'s dimensions (choose from: checkpoint, diagnostic)'
    )
    parser.add_argument(
        '--dim_format_out',
        metavar='WHICH',
        type=str,
        choices=[
            'checkpoint',
            'diagnostic',
            'classic'],
        required=True,
        help='format of the output file\'s dimensions (choose from: checkpoint, diagnostic)'
    )
    parser.add_argument(
        '--vert_params_out',
        metavar='VERT',
        type=list,
        required=False,
        help='Hybrid grid parameter A in hPa and B (unitless) in [AP, BP] format'
    )
    args = parser.parse_args()

    # Regrid the file
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


# Only call when run as standalone
if __name__ == '__main__':
    main()
