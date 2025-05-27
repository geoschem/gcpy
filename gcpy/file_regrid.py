"""
Regrids data horizontally between lat/lon and/or cubed-sphere grids
(including stretched grids).
"""
import argparse
import os
import warnings
import numpy as np
import xarray as xr
from gcpy.grid import get_input_res, get_grid_extents, \
    get_ilev_coord, get_lev_coord
from gcpy.regrid import make_regridder_S2S, reformat_dims, \
    make_regridder_L2S, make_regridder_C2L, make_regridder_L2L
from gcpy.util import verify_variable_type
from gcpy.cstools import get_cubed_sphere_res, is_gchp_lev_positive_down

# Ignore any FutureWarnings
warnings.simplefilter(action="ignore", category=FutureWarning)


def file_regrid(
        filein,
        fileout,
        dim_format_in,
        dim_format_out,
        cs_res_out=0,
        ll_res_out="0x0",
        sg_params_in=None,
        sg_params_out=None,
        verbose=False,
        weightsdir="."
):
    """
    Regrids an input file to a new horizontal grid specification
    and saves it as a new file.

    Args:
    -----
    filein: str
        The input filename
    fileout: str
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
    ------------------------
    cs_res_out: int
        The cubed-sphere resolution of the output dataset.
        Not used if dim_format_out is classic
        Default value: 0
    ll_res_out: str
        The lat/lon resolution of the output dataset.
        Not used if dim_format_out is not classic
        Default value: "0x0"
    sg_params_in: list[float, float, float]
        Input grid stretching parameters
        [stretch-factor, target longitude, target latitude].
        Not used if dim_format_in is classic.
        Default value: [1.0, 170.0, -90.0] (No stretching)
    sg_params_out: list[float, float, float]
        Output grid stretching parameters
        [stretch-factor, target longitude, target latitude].
        Not used if dim_format_out is classic
        Default value: [1.0, 170.0, -90.0] (No stretching)
    verbose : bool
        Toggles verbose output on (True) or off (False).
    weightsdir : str
        Path to the directory containing regridding weights (or
        where weights will be created).  Default value: "."
    """
    verify_variable_type(filein, str)
    verify_variable_type(fileout, str)
    verify_variable_type(dim_format_in, str)
    verify_variable_type(dim_format_out, str)

    # TODO: Consider renaming checkpoint, classic, diagnostic,
    # which may be confusing to users.

    # Error check arguments
    valid_formats = ["checkpoint", "classic", "diagnostic"]
    if dim_format_in not in valid_formats:
        msg = f"Argument 'dim_format_in' must be one of: {valid_formats}!"
        raise ValueError(msg)
    if dim_format_out not in valid_formats:
        msg = f"Argument 'dim_format_out' must be one of: {valid_formats}!"
        raise ValueError(msg)

    # Assign default values for optional keywords
    if sg_params_in is None:
        sg_params_in = [1.0, 170.0, -90.0]
    if sg_params_out is None:
        sg_params_out = [1.0, 170.0, -90.0]

    # ------------------------------------------------------------------
    # There still seem to be a few issues with regridding to cubed-
    # sphere stretched grids.  For time time being, stop with error
    # if sg_params_in or sg_params_out do not equal the defaults.
    #  -- Bob Yantosca & Lizzie Lundgren (24 Oct 2023)
    if not np.array_equal(sg_params_in, [1.0, 170.0, -90.0]) or \
       not np.array_equal(sg_params_out, [1.0, 170.0, -90.0]):
        msg = "Regridding to or from cubed-sphere stretched grids is\n" + \
            "currently not supported.  Please use the offline regridding\n" + \
            "method described in the Regridding section of gcpy.readthedocs.io."
        raise RuntimeError(msg)
    # ------------------------------------------------------------------

    # Load dataset
    dset = xr.open_dataset(
        filein,
        decode_cf=False,
        engine="netcdf4"
    ).load()
    cs_res_in = get_cubed_sphere_res(dset)

    # Verbose printout of inputs
    if verbose:
        print("Inputs to file_regrid.py")
        print(f"filein         :  {filein}")
        print(f"dim_format_in  :  {dim_format_in}")
        if "classic" not in dim_format_in:
            print(f"sg_params_in   :  {sg_params_in}")
        print(f"fileout        :  {fileout}")
        print(f"dim_format_out :  {dim_format_out}")
        if "classic" in dim_format_out:
            print(f"ll_res_out     :  {ll_res_out}")
        else:
            print(f"cs_res_out     :  {cs_res_out}")
            print(f"sg_params_out  :  {sg_params_out}")
        print(f"verbose        :  {verbose}")
        print(f"weightsdir     :  {weightsdir}")

    # Make sure all xarray.Dataset global & variable attributes are kept
    with xr.set_options(keep_attrs=True):

        # ==============================================================
        # Regrid data
        # ==============================================================

        # Save type of data for later restoration
        # Avoid using the dtype of GCHP cubed-sphere grid variables
        dset_tmp = dset
        dtype_orig = np.dtype(dset[list(dset_tmp.data_vars.keys())[-1]])
        dset_tmp = xr.Dataset()

        if dim_format_in != "classic" and dim_format_out != "classic":

            # ----------------------------------------------------------
            # Input grid is CS/SG; Output grid is CS/SG
            # ----------------------------------------------------------
            dset = regrid_cssg_to_cssg(
                fileout,
                dset,
                dim_format_in,
                sg_params_in,
                cs_res_out,
                dim_format_out,
                sg_params_out,
                verbose=verbose,
                weightsdir=weightsdir
            )

        elif dim_format_in == "classic" and dim_format_out != "classic":

            # ----------------------------------------------------------
            # Input grid is LL; Output grid is CS/SG
            # ----------------------------------------------------------
            dset = regrid_ll_to_cssg(
                dset,
                cs_res_out,
                dim_format_out,
                sg_params_out,
                verbose=verbose,
                weightsdir=weightsdir
            )

        elif dim_format_in != "classic" and dim_format_out == "classic":

            # ----------------------------------------------------------
            # Input grid is CS/SG; Output grid is LL
            # ----------------------------------------------------------
            dset = regrid_cssg_to_ll(
                dset,
                cs_res_in,
                dim_format_in,
                sg_params_in,
                ll_res_out,
                verbose=verbose,
                weightsdir=weightsdir
            )

        elif dim_format_in == "classic" and dim_format_out == "classic":

            # ----------------------------------------------------------
            # Input grid is LL; Output grid is LL
            # ----------------------------------------------------------
            dset = regrid_ll_to_ll(
                dset,
                ll_res_out,
                verbose=verbose,
                weightsdir=weightsdir
            )

        # ==============================================================
        # Post-regridding stuff
        # ==============================================================

        # Correct precision changes (accidental 32-bit to 64-bit)
        # NOTE: Add a workaround to prevent the xr.DataArray.astype
        # function from overwriting the "lev" dimension.
        dset_tmp = dset.astype(
            dtype=dtype_orig,
            casting="same_kind",
            copy=False
        )
        dset = dset_tmp.assign_coords(lev=dset.lev)

        # Write dataset to file
        dset.to_netcdf(
            fileout,
            mode="w",
            format="NETCDF4",
            engine="netcdf4",
            unlimited_dims=["time"],
        )

        # Free memory of the temporary dataset
        dset_tmp = xr.Dataset()

        # Print the resulting dataset
        if verbose:
            print(dset)


def prepare_cssg_input_grid(
        dset,
        dim_format_in
):
    """
    Reformats cubed-sphere/stretched grid data to the universal
    format and drops non-regriddable fields.

    Args:
    -----
    dset : xr.Dataset
        Input grid (cubed-sphere or stretched grid)
    dim_format_in : str
        Either "checkpoint" (for restart files)
        or "diagnostic" (for History diagnostic files)

    Returns:
    --------
    dset : xr.Dataset
        Data with reformatted dimensions and dropped fields
    cs_res_in : int
        Cubed-sphere/stretched grid resolution
    """

    # Reformat dimensions to "common dimensions (T, Z, F, Y, X)
    dset = reformat_dims(
        dset,
        dim_format_in,
        towards_common=True
    )

    # Drop variables that don't look like fields
    # NOTE: Don't drop "lons" and "lats" if present.
    non_fields = [
        v for v in dset.variables.keys()
        if len(set(dset[v].dims) - {"T", "Z", "F", "Y", "X"}) > 0
        or len(dset[v].dims) == 0]
    dset_in = dset.drop(non_fields)

    # Transpose to T, Z, F, Y, X
    dset = dset_in.transpose("T", "Z", "F", "Y", "X")

    assert dset.dims["X"] == dset.dims["Y"]
    cs_res_in = dset.dims["X"]

    return dset, cs_res_in


def regrid_cssg_to_cssg(
        fileout,
        dset,
        dim_format_in,
        sg_params_in,
        cs_res_out,
        dim_format_out,
        sg_params_out,
        verbose=False,
        weightsdir="."
):
    """
    Regrids from the cubed-sphere/stretched grid to a different
    cubed-sphere/stretched grid resolution.

    Args:
    -----
    fileout : str
        File name template
    dset : xarray.Dataset
        Data on a cubed-sphere/stretched grid
    dim_format_in, dim_format_out : str
        Input & output grid format ("checkpoint", "diagnostic")
    cs_res_out : int
        Cubed-sphere grid resolution
    sg_params_in, sg_params_out: list[float, float, float]
        Input & output grid stretching parameters
        [stretch-factor, target longitude, target latitude].

    Keyword Args (optional):
    ------------------------
    verbose : bool
        Toggles verbose output on (True) or off (False).
    weightsdir : str
        Path to the directory containing regridding weights (or
        where weights will be created).  Default value: "."

    Returns:
    --------
    dset : xarray.Dataset
        Data regridded to the output lat-lon grid
    """
    if verbose:
        print("file_regrid.py: Regridding from CS/SG to CS/SG")

    # Keep all xarray attributes
    with xr.set_options(keep_attrs=True):

        # Flip vertical levels (if necessary) and
        # set the lev:positive attribute accordingly
        dset = flip_lev_coord_if_necessary(
            dset,
            dim_format_in=dim_format_in,
            dim_format_out=dim_format_out
        )

        # Change CS/SG dimensions to universal format
        # and drop non-regriddable variables
        dset, cs_res_in = prepare_cssg_input_grid(
            dset,
            dim_format_in
        )

        # ==============================================================
        # Only regrid if the cubed-sphere grids are similar
        # (i.e. same resolution & stretched-grid parameters)
        # ==============================================================
        if cs_res_in == cs_res_out and \
           np.array_equal(sg_params_in, sg_params_out) and \
           dim_format_in == dim_format_out:
            print("Skipping regridding since grid parameters are identical")

            # Put regridded dataset back into a familiar format
            dset = dset.rename({
                "y": "Y",
                "x": "X",
            })

            return dset

        # Make regridders
        regridders = make_regridder_S2S(
            cs_res_in,
            cs_res_out,
            sf_in=sg_params_in[0],
            tlon_in=sg_params_in[1],
            tlat_in=sg_params_in[2],
            sf_out=sg_params_out[0],
            tlon_out=sg_params_out[1],
            tlat_out=sg_params_out[2],
            weightsdir=weightsdir
        )

        # Save temporary output face files to minimize RAM usage
        oface_files = [os.path.join(".",fileout+str(x)) for x in range(6)]

        # For each output face, sum regridded input faces
        for oface in range(6):
            oface_regridded = []
            for (iface, regridder) in regridders[oface].items():
                dset_iface = dset.isel(F=iface)
                if "F" in dset_iface.coords:
                    dset_iface = dset_iface.drop("F")
                oface_regridded.append(
                    regridder(
                        dset_iface,
                        keep_attrs=True
                    )
                )
            oface_regridded = xr.concat(
                oface_regridded,
                dim="intersecting_ifaces"
            ).sum(
                "intersecting_ifaces",
                keep_attrs=True).expand_dims({"F":[oface]})
            oface_regridded.to_netcdf(
                oface_files[oface],
                format="NETCDF4",
                engine="netcdf4",
                mode="w"
            )

        # Combine face files
        dset = xr.open_mfdataset(
            oface_files,
            combine="nested",
            concat_dim="F",
            engine="netcdf4"
        )

        # Remove any temporary files
        for oface in oface_files:
            os.remove(oface)

        # ==============================================================
        # Reshape the data if necessary
        # ==============================================================

        # Put regridded dataset back into a familiar format
        dset = dset.rename({
            "y": "Y",
            "x": "X",
        })

        # Reformat dimensions from "common dimension format"
        # to CS/GG "checkpoint" or "diagnostics" format
        dset = reformat_dims(
            dset,
            format=dim_format_out,
            towards_common=False
        )

        # Rename variables if we are going between different grid types
        if "checkpoint" in dim_format_in and "diagnostic" in dim_format_out:
            dset = rename_restart_variables(
                dset,
                towards_gchp=False
            )
        if "diagnostic" in dim_format_in and "checkpoint" in dim_format_out:
            dset = rename_restart_variables(
                dset,
                towards_gchp=True
            )

        # Fix names and attributes of of coordinate variables depending
        # on the format of the ouptut grid (checkpoint or diagnostic).
        dset = adjust_cssg_grid_and_coords(
            dset,
            dim_format_in,
            dim_format_out
        )

        # Save stretched-grid metadata as global attrs
        dset = save_cssg_metadata(
            dset,
            cs_res_out,
            dim_format_out,
            sg_params_out,
            verbose=verbose
        )

    return dset


def regrid_cssg_to_ll(
        dset,
        cs_res_in,
        dim_format_in,
        sg_params_in,
        ll_res_out,
        verbose=False,
        weightsdir="."
):
    """
    Regrids from the cubed-sphere/stretched grid to the lat-lon grid.

    Args:
    -----
    dset : xarray.Dataset
        Data on a cubed-sphere/stretched grid
    cs_res_in : int
        Cubed-sphere grid resolution
    dim_format_in : str
        Input grid format ("checkpoint", "diagnostic")
    sg_params_in: list[float, float, float]
        Input grid stretching parameters
        [stretch-factor, target longitude, target latitude].
    ll_res_out : str
        Output grid lat/lon resolution (e.g. "4x5")

    Keyword Args (optional):
    ------------------------
    verbose: bool
        Toggles verbose printout on (True) or off (False)
    weightsdir : str
        Path to the directory containing regridding weights (or
        where weights will be created).  Default value: "."

    Returns:
    --------
    dset : xarray.Dataset
        Data regridded to the output lat-lon grid
    """
    if verbose:
        print("file_regrid.py: Regridding from CS/SG to LL")

    with xr.set_options(keep_attrs=True):

        # Flip vertical levels (if necessary) and
        # set the lev:positive attribute accordingly
        dset = flip_lev_coord_if_necessary(
            dset,
            dim_format_in=dim_format_in,
            dim_format_out="classic"
        )

        # Drop non-regriddable variables (if any)
        dset = drop_classic_vars(
            dset,
            towards_gchp=False
        )

        # Change CS/SG dimensions to universal format
        # and drop non-regriddable variables
        dset, cs_res_in = prepare_cssg_input_grid(
            dset,
            dim_format_in
        )

        # Regrid data
        regridders = make_regridder_C2L(
            cs_res_in,
            ll_res_out,
            sg_params=sg_params_in,
            weightsdir=weightsdir
        )
        dset = xr.concat(
            [regridders[face](dset.isel(F=face), keep_attrs=True)
             for face in range(6)],
            dim="F"
        ).sum("F", keep_attrs=True)

        # Update dimensions and attributes on the lat-lon grid
        dset = dset.rename({
            "T": "time",
            "Z": "lev"
        })

        # If regridding from a GCHP checkpoint/restart file, then
        # rename variables to adhere GCClassic name conventions.
        if "checkpoint" in dim_format_in:
            dset = rename_restart_variables(
                dset,
                towards_gchp=False
            )

        # Save lat/lon coordinate metadata
        dset = save_ll_metadata(
            dset,
            verbose=verbose
        )

        # Drop cubed-sphere variables
        if "lons" in dset.data_vars:
            dset = dset.drop_vars(["lons"])
        if "lats" in dset.data_vars:
            dset = dset.drop_vars(["lats"])

    return dset


def regrid_ll_to_cssg(
        dset,
        cs_res_out,
        dim_format_out,
        sg_params_out,
        verbose=False,
        weightsdir="."
):
    """
    Regrids from the lat-lon grid to the cubed-sphere/stretched grid.

    Args:
    -----
    dset : xarray.Dataset
        Data on a lat/lon grid
    cs_res_in : int
        Cubed-sphere grid resolution
    dim_format_out : str
        Either "checkpoint" (for restart files) or
        "diagnostic" (for History diagnostic files).
    sg_params_out: list[float, float, float]
        Output grid stretching parameters
        [stretch-factor, target longitude, target latitude].

    Keyword Args (optional):
    ------------------------
    verbose : bool
        Toggles verbose output on (True) or off (False).
    weightsdir : str
        Path to the directory containing regridding weights (or
        where weights will be created).  Default value: "."

    Returns:
    --------
    dset : xarray.Dataset
        Data regridded to the output cubed-sphere/stretched-grid
    """
    if verbose:
        print("file_regrid.py: Regridding from LL to CS/SG")

    with xr.set_options(keep_attrs=True):

        # Flip vertical levels (if necessary) and set lev:positive
        dset = flip_lev_coord_if_necessary(
            dset,
            dim_format_in="classic",
            dim_format_out=dim_format_out
        )

        # Drop non- regriddable variables when going from ll -> cs
        dset = drop_classic_vars(dset)

        # If regridding to a GCHP checkpoint/restart file, then
        # rename variables to adhere to GCHP naming conventions.
        if "checkpoint" in dim_format_out:
            dset = rename_restart_variables(
                dset,
                towards_gchp=True
            )

        # Input lat/lon grid resolution
        llres_in = get_input_res(dset)[0]

        # Regrid data to CS/SG
        regridders = make_regridder_L2S(
            llres_in,
            cs_res_out,
            sg_params=sg_params_out,
            weightsdir=weightsdir
        )
        dset = xr.concat(
            [regridders[face](dset, keep_attrs=True) for face in range(6)],
            dim="nf"
        )

        # Rename dimensions to the "common dimension format"
        dset = dset.rename({
            "time": "T",
            "lev": "Z",
            "nf": "F",
            "y": "Y",
            "x": "X",
            "lat": "Y",
            "lon": "X"
        })

        # Reformat dims from "common dimension format" to "diagnostic"
        # (we will convert to "checkpoint" later)
        dset = reformat_dims(
            dset,
            format="diagnostic",
            towards_common=False
        )

        # Fix names and attributes of of coordinate variables depending
        # on the format of the ouptut grid (checkpoint or diagnostic).
        # Also convert the "diagnostic" grid to the "checkpoint" grid
        # if "checkpoint" output are requested.
        dset = adjust_cssg_grid_and_coords(
            dset,
            dim_format_in="diagnostic",
            dim_format_out=dim_format_out
        )

        # Save stretched-grid metadata as global attrs
        dset = save_cssg_metadata(
            dset,
            cs_res_out,
            dim_format_out,
            sg_params_out,
            verbose=verbose
        )

    return dset


def regrid_ll_to_ll(
        dset,
        ll_res_out,
        verbose=False,
        weightsdir="."
):
    """
    Regrid from the lat/lon grid to the cubed-sphere/stretched grid.

    Args:
    -----
    dset : xarray.Dataset
        Data on a lat/lon grid
    ll_res_out : str
        Output grid lat-lon grid resolution (e.g. "4x5")

    Keyword Args (optional):
    ------------------------
    verbose : bool
        Toggles verbose output on (True) or off (False).
    weightsdir : str
        Path to the directory containing regridding weights (or
        where weights will be created).  Default value: "."

    Returns:
    --------
    dset : xarray.Dataset
        Data regridded to the output lat-lon grid.
    """
    if verbose:
        print("file_regrid.py: Regridding from LL to LL")

    with xr.set_options(keep_attrs=True):

        # Get the input & output extents
        in_extent = get_grid_extents(dset)
        out_extent = in_extent
        ll_res_in = get_input_res(dset)[0]
        [lat_in, lon_in] = list(map(float, ll_res_in.split("x")))
        [lat_out, lon_out] = list(map(float, ll_res_out.split("x")))

        # Return if the output & input grids are the same
        if lat_in == lat_out and lon_in == lon_out:
            print("Skipping regridding since grid parameters are identical")
            return dset

        # Drop non-regriddable variables
        non_fields = [
            var for var in dset.variables.keys() \
            if "lat" not in dset[var].dims       \
                and "lon" not in dset[var].dims
            ]
        if "lat_bnds" in dset.data_vars:
            dset = dset.drop(["lat_bnds"])
        if "lon_bnds" in dset.data_vars:
            dset = dset.drop(["lon_bnds"])
        non_fields = dset[non_fields]
        dset = dset.drop(non_fields)

        # Set the lev:positive attribute accordingly
        dset = flip_lev_coord_if_necessary(
            dset,
            dim_format_in="classic",
            dim_format_out="classic"
        )

        # Decide if we are regridding a data file or a mask
        # by testing for the variable name "MASK"
        method = "conservative"
        if "MASK" in dset.data_vars:
            method = "nearest_s2d"

        # Create the regridder and regrid the data
        regridder = make_regridder_L2L(
        ll_res_in,
            ll_res_out,
            reuse_weights=True,
            in_extent=in_extent,
            out_extent=out_extent,
            weightsdir=weightsdir,
            method=method,
        )
        dset = regridder(
            dset,
            keep_attrs=True
        )

        # Add the non-regriddable fields back
        dset = dset.merge(non_fields)

        # Change order of dimensions
        dset = dset.transpose(
            "time", "lev", "ilev", "lat", "lon", ...
        )

        # Save lat/lon coordinate metadata
        dset = save_ll_metadata(
            dset,
            verbose=verbose
        )

    return dset


def flip_lev_coord_if_necessary(
        dset,
        dim_format_in,
        dim_format_out
):
    """
    Flips the "lev" and "ilev" coords of an xarray.Dataset in the
    vertical depending on the values of dim_format_in and
    dim_format_out.  Also sets the attributes "lev:positive" and
    "ilev:positive" accordingly.

    Args:
    -----
    dset : xarray.Dataset
        The input dataset.
    dim_format_in : str
        Input grid format ("classic", "checkpoint", "diagnostic").
    dim_format_out : str
        Output grid format ("classic", "checkpoint", "diagnostic").

    Args:
    -----
    dset : xarray.Dataset
        The modified dataset.

    Remarks:
    --------
    (1) classic    : lev is in ascending order  (lev:positive="up"  )
    (2) diagnostic : lev is in ascending order* (lev:positive="up"  )
    (3) checkpoint : lev is in descending order (lev:positive="down")

    *Except for the Emissions collection, which has lev arranged
     in descending order.

    TODO: Make ths function more robust for all cases, since GCHP
    diagnostics may or may not have lev:positive="up".
    """
    verify_variable_type(dset, xr.Dataset)
    verify_variable_type(dim_format_in, str)
    verify_variable_type(dim_format_out, str)

    # ==================================================================
    # Case 1: checkpoint/diagnostic to classic
    # lev, ilev need to be in ascending order
    # ==================================================================
    if dim_format_in != "classic" and dim_format_out == "classic":

        # Flip lev and set to eta values at midpoints (if necessary)
        if "ilev" in dset.coords:
            if is_gchp_lev_positive_down(dset):
                dset = dset.reindex(ilev=dset.ilev[::-1])
            coord = get_ilev_coord(
                n_lev=dset.dims["ilev"],
                top_down=False
            )
            dset = dset.assign_coords({"ilev": coord})
            dset.ilev.attrs["positive"] = "up"

        # Flip lev and set to eta values at midpoints (if necessary)
        if "lev" in dset.coords:
            if is_gchp_lev_positive_down(dset):
                dset = dset.reindex(lev=dset.lev[::-1])
            coord = get_lev_coord(
                n_lev=dset.dims["lev"],
                top_down=False
            )
            dset = dset.assign_coords({"lev": coord})
            dset.lev.attrs["positive"] = "up"

        return dset

    # ==================================================================
    # Case 2: classic/diagnostic to checkpoint
    # lev needs to be in descending order (with ascending indices)
    #
    # TODO: Check for Emissions diagnostic (not a common use case)
    # ==================================================================
    if dim_format_in != "checkpoint" and dim_format_out == "checkpoint":

        if "lev" in dset.coords:
            if not is_gchp_lev_positive_down(dset):
                dset = dset.reindex(lev=dset.lev[::-1])
            coord = get_lev_coord(
                n_lev=dset.dims["lev"],
                gchp_indices=True
            )
            dset = dset.assign_coords({"lev": coord})
            dset.lev.attrs["positive"] = "down"

        return dset

    # ==================================================================
    # Case 3: classic/checkpoint to diagnostic:
    # lev, ilev need to be in ascending order
    #
    # TODO: Check for Emissions diagnostic (not a common use case)
    # ==================================================================
    if dim_format_in != "diagnostic" and dim_format_out == "diagnostic":

        if "lev" in dset.coords:
            if is_gchp_lev_positive_down(dset):
                dset = dset.reindex(lev=dset.lev[::-1])
            coord = get_lev_coord(
                n_lev=dset.dims["lev"],
                gchp_indices=True
            )
            dset = dset.assign_coords({"lev": coord})
            dset.lev.attrs["positive"] = "up"

        return dset

    # ==================================================================
    # Case 4: checkpoint to checkpoint
    # No flipping needed, but add lev:positive="down"
    # ==================================================================
    if dim_format_in == "checkpoint" and dim_format_out == "checkpoint":

        if "lev" in dset.coords:
            dset.lev.attrs["positive"] = "down"
        return dset

    return dset


def save_ll_metadata(
        dset,
        verbose=False,
):
    """
    Updates the lat-lon coordinate metadata in an xarray.Dataset object.

    Args:
    -----
    dset : xarray.Dataset
        The input data (on lat-lon grid).

    Keyword Arguments:
    ------------------
    verbose : bool
        Toggles verbose printout on (True) or off (False)

    Returns:
    --------
    dset : xarray.Dataset
        Original data plus updated coordinate metadata.
    """
    with xr.set_options(keep_attrs=True):

        dset.time.attrs = {
            "axis": "T"
        }

        dset.lat.attrs = {
            "long_name": "Latitude",
            "units": "degrees_north",
            "axis": "Y"
        }

        dset.lon.attrs = {
            "long_name": "Longitude",
            "units": "degrees_east",
            "axis": "X"
        }

        if "ilev" in dset.coords:
            dset.ilev.attrs["long_name"] = \
                "hybrid level at interfaces ((A/P0)+B)"
            dset.ilev.attrs["units"] = "level"
            dset.ilev.attrs["axis"] = "Z"

        if "lev" in dset.coords:
            dset.lev.attrs["long_name"] = \
                "hybrid level at midpoints ((A/P0)+B)"
            dset.lev.attrs["units"] = "level"
            dset.lev.attrs["axis"] = "Z"

    if verbose:
        print("file_regrid.py: In routine save_ll_metadata:")
        print(dset.coords)

    return dset


def save_cssg_metadata(
        dset,
        cs_res_out,
        dim_format_out,
        sg_params_out,
        verbose=False
):
    """
    Saves the stretched-grid metadata to an xarray.Dataset object
    containing cubed-sphere/stretched grid data.

    Args:
    -----
    dset : xarray.Dataset
        Data on the stretched grid.
    cs_res_out : int
        Cubed-sphere grid resolution.
    dim_format_out : str
        Either "checkpoint" (for restart files) or
        "diagnostic" (for History diagnostic files).
    sg_params_out: list[float, float, float]
        Output grid stretching parameters
        [stretch-factor, target longitude, target latitude].
    verbose : bool
        Toggles verbose printout on (True) or off (False).

    Returns:
    --------
    dset : xarray.Dataset
        The original data, plus stretched grid metadata.
    """
    if verbose:
        print("file_regrid.py: Saving CS/SG coordinate metadata")

    with xr.set_options(keep_attrs=True):

        # Stretched-grid global attrs
        dset.attrs["stretch_factor"] = sg_params_out[0]
        dset.attrs["target_longitude"] = sg_params_out[1]
        dset.attrs["target_latitude"] = sg_params_out[2]
        dset.attrs["cs_res"] = cs_res_out

        # Special handling for "checkpoint" format
        if "checkpoint" in dim_format_out:
            if "lon" in dset.dims:
                dset.lon.attrs = {
                    "standard_name": "longitude",
                    "long_name": "Longitude",
                    "units": "degrees_east",
                    "axis": "X"
                }
            if "lat" in dset.dims:
                dset.lat.attrs = {
                    "standard_name": "latitude",
                    "long_name": "Latitude",
                    "units": "degrees_north",
                    "axis": "Y"
                }

        # Special handling for "checkpoint" format
        if "diagnostic" in dim_format_out:
            if "lons" in dset.dims:
                dset.lons.attrs = {
                    "standard_name": "longitude",
                    "long_name": "Longitude",
                    "units": "degrees_east",
                    "axis": "X"
                }
            if "lats" in dset.dims:
                dset.lats.attrs = {
                    "standard_name": "la7titude",
                    "long_name": "Latitude",
                    "units": "degrees_north",
                    "axis": "Y"
                }

        # ilev:positive is set by flip_lev_coord_if_necessary
        if "ilev" in dset.coords:
            dset.ilev.attrs["long_name"] = \
                "hybrid level at interfaces ((A/P0)+B)"
            dset.ilev.attrs["units"] = "level"
            dset.ilev.attrs["axis"] = "Z"

        # lev:positive is set by flip_lev_coord_if_necessary
        if "lev" in dset.coords:
            dset.lev.attrs["long_name"] = \
                "hybrid level at midpoints ((A/P0)+B)"
            dset.lev.attrs["units"] = "level"
            dset.lev.attrs["axis"] = "Z"

    return dset


def rename_restart_variables(
        dset,
        towards_gchp=True
):
    """
    Renames restart variables according to GEOS-Chem Classic
    and GCHP conventions.

    Args:
    -----
    dset : xarray.Dataset
        The input dataset.

    Keyword Args (optional):
    ------------------------
    towards_gchp: bool
        Whether renaming to (True) or from (False) GCHP format
        Default value: True

    Returns:
    --------
    dset : xarray.Dataset
       The modified dataset.
    """
    verify_variable_type(dset, xr.Dataset)

    # Keep all xarray attribute settings
    with xr.set_options(keep_attrs=True):

        # Dictionary for name replacements
        old_to_new = {}

        # ==============================================================
        # classic/diagnostic -> checkpoint
        # ==============================================================
        if towards_gchp:
            for var in dset.data_vars.keys():
                if var.startswith("Met_"):
                    if "DELPDRY" in var:
                        old_to_new[var] = "DELP_DRY"
                    else:
                        old_to_new[var] = var.replace("Met_", "")
                if var.startswith("Chem_"):
                    old_to_new[var] = var.replace("Chem_", "")
                if var.startswith("SpeciesRst_"):
                    old_to_new[var] = var.replace("SpeciesRst_", "SPC_")
                if var.startswith("SpeciesConcVV_"):
                    old_to_new[var] = var.replace("SpeciesConcVV_", "SPC_")

            return dset.rename(old_to_new)

        # ==============================================================
        # checkpoint -> classic/diagnostic
        # ==============================================================
        for var in dset.data_vars.keys():
            if var in ("DELP_DRY", "DELPDRY"):
                old_to_new[var] = "Met_DELPDRY"
            if var == "BXHEIGHT":
                old_to_new[var] = "Met_BXHEIGHT"
            if var == "StatePSC":
                old_to_new[var] = "Chem_StatePSC"
            if var == "KPPHvalue":
                old_to_new[var] = "Chem_KPPHvalue"
            if var == "DryDepNitrogen":
                old_to_new[var] = "ChemDryDepNitrogen"
            if var == "WetDepNitrogen":
                old_to_new[var] = "Chem_WetDepNitrogen"
            if var == "SO2AfterChem":
                old_to_new[var] = "Chem_SO2AfterChem"
            if var == "JNO2":
                old_to_new[var] = "Chem_JNO2"
            if var == "JOH":
                old_to_new[var] = "Chem_JOH"
            if var == "H2O2AfterChem":
                old_to_new[var] = "Chem_H2O2AfterChem"
            if var == "ORVCsesq":
                old_to_new[var] = "Chem_ORVCsesq"
            if var == "AeroH2OSNA":
                old_to_new[var] = "Chem_AeroH2OSNA"
            if var.startswith("SPC_"):
                old_to_new[var] = var.replace("SPC_", "SpeciesRst_")

        return dset.rename(old_to_new)


def adjust_cssg_grid_and_coords(
        dset,
        dim_format_in,
        dim_format_out,
):
    """
    Adjusts cubed-sphere/stretched-grid coordinate names and attributes.

    Args:
    -----
    dset : xarray.Dataset
        The input data
    dim_format_in, dim_format_out: str
        Either "checkpoint" (for checkpoint/restart files) or
        "diagnostic" (for History diagnostic files).

    Returns:
    --------
    dset : xarray.Dataset
       The input data with updated coordinate names & attributes.

    Remarks:
    --------
    "diagnostic" dimension format: (time, lev, nf, Ydim, Xdim)
    "checkpoint" dimension format: (time, lev, lat, lon); lat = 6*lon
    """
    # Keep all xarray attributes intact
    with xr.set_options(keep_attrs=True):

        # ==============================================================
        # Rename coordinates returned by the xESMF regridding to
        # the "lons" and "lats" coordinates as saved out by MAPL.
        # ==============================================================
        if "diagnostic" in dim_format_in:
            if "Xdim" in dset.variables:
                dset = dset.rename_vars({"Xdim": "lons"})
            if "Ydim" in dset.variables:
                dset = dset.rename_vars({"Ydim": "lats"})

        if "checkpoint" in dim_format_in:
            if "lon" in dset.variables:
                dset = dset.rename_vars({"lon": "lons"})
            if "lat" in dset.variables:
                dset = dset.rename_vars({"lat": "lats"})

        if "lons" in dset.variables:
            dset.lons.attrs = {
                "standard_name": "longitude",
                "long_name": "Longitude",
                "units": "degrees_east"
            }

        if "lats" in dset.variables:
            dset.lats.attrs = {
                "standard_name": "latitude",
                "long_name": "latitude",
                "units": "degrees_north"
            }

        # ==================================================================
        # For "diagnostic" dimension format only
        # ==================================================================
        if "diagnostic" in dim_format_out:

            # Add "fake" Xdim and Ydim coordinates as done by MAPL,
            # which is needed for the GMAO GrADS visualization software.
            # NOTE: Use .values to convert to numpy.ndarray type in
            # order to avoid xarray from trying to redefileine dim "nf".
            if "lons" in dset.coords and "lats" in dset.coords:
                dset = dset.assign_coords({
                    "Xdim": dset.lons.isel(nf=0, Ydim=0).values,
                    "Ydim": dset.lats.isel(nf=0, Xdim=0).values
                })
            elif "lon" in dset.variables and "lat" in dset.variables:
                dset = dset.assign_coords({
                    "Xdim": dset.lon.isel(nf=0, Ydim=0).values,
                    "Ydim": dset.lat.isel(nf=0, Xdim=0).values
                })
            dset.Xdim.attrs = {
                "long_name": "Fake Longitude for GrADS Compatibility",
                "units": "degrees_east"
            }
            dset.Ydim.attrs = {
                "long_name": "Fake Latitude for GrADS Compatibility",
                "units": "degrees_north"
            }

            # Drop dimensions that may be left over from regridding
            if "lon" in dset.variables:
                dset = dset.drop_vars("lon")
            if "lat" in dset.variables:
                dset = dset.drop_vars("lat")

        # ==================================================================
        # For "checkpoint" dimension format only
        # ==================================================================
        if "checkpoint" in dim_format_out:

            # Reshape the grid from (time, lev, nf, Xdim, Ydim) dimensions
            # to (time, lev, lat, lon) dimensions (where lat/lon = 6)
            # Also drop any unnecessary variables
            dset = reshape_cssg_diag_to_chkpt(dset)
            if "lons" in dset.variables:
                dset = dset.drop_vars("lons")
            if "lats" in dset.variables:
                dset = dset.drop_vars("lats")

    return dset


def drop_classic_vars(
        dset,
        towards_gchp=True
):
    """
    Renames and drops certain restart variables according to
    GEOS-Chem Classic and GCHP conventions.

    Args:
    -----
    dset : xarray.Dataset
        The input dataset.

    Keyword Args (optional):
    ------------------------
    towards_gchp: bool
        Whether going to (True) or from (False) GCHP format.
        Default value: True

    Returns:
    --------
    dset : xarray.Dataset
        The modified dataset.
    """
    with xr.set_options(keep_attrs=True):
        if towards_gchp:
            dset = dset.drop_vars(
                ["P0",
                 "hyam",
                 "hybm",
                 "hyai",
                 "hybi",
                 "AREA",
                 "ilev",
                 "PS1DRY",
                 "PS1WET",
                 "TMPU1",
                 "SPHU1",
                 "StatePSC",
                 "lon_bnds",
                 "lat_bnds"],
                errors="ignore"
            )

        return dset


def order_dims_time_lev_lat_lon(dset):
    """
    Transposes dims of an Dataset to be in (time, lev, lat, lon) order.
    This corresponds to Fortran column-major ordering.


    Args:
    -----
    dset : xarray.Dataset
        The input dataset.

    Returns:
    --------
    dset : xarray.Dataset
        The modified dataset.
    """
    verify_variable_type(dset, xr.Dataset)

    if "lev" in dset.dims and "time" in dset.dims:
        dset = dset.transpose("time", "lev", "lat", "lon")
    elif "lev" in dset.dims:
        dset = dset.transpose("lev", "lat", "lon")
    elif "time" in dset.dims:
        dset = dset.transpose("time", "lat", "lon")
    else:
        dset = dset.transpose("lat", "lon")

    return dset


def reshape_cssg_diag_to_chkpt(
        dset,
        verbose=False
):
    """
    Reshapes a dataset from diagnostic to checkpoint dimension format.

    Args:
    -----
    dset : xarray.Dataset
        Dataset with dimensions (time, lev, nf, Xdim, Ydim).

    Keyword Args (optional)
    -----------------------
    verbose : bool
        Toggles verbose output on (True) or off (False).

    Returns:
    --------
    dset : xarray.Dataset
        Dataset wtih dimensions (time, lev, lat, lon), where lat/lon=6.
    """
    verify_variable_type(dset, xr.Dataset)

    if verbose:
        print("file_regrid.py: reshyaping diagnostic to checkpoint")

    # Keep xarray attributes unchanged
    with xr.set_options(keep_attrs=True):

        # ==============================================================
        # Get the size of the Xdim/YDim or lons/lats coords
        # ==============================================================
        if "Xdim" in dset.dims and "Ydim" in dset.dims:
            xdim = dset.dims["Xdim"]
            ydim = dset.dims["Ydim"]
        elif "lon" in dset.dims and "lat" in dset.dims:
            xdim = dset.dims["lon"]
            ydim = dset.dims["lat"]
        else:
            msg = "Dimensions (Xdim, Ydim) or (lon,lat)' not found!"
            raise ValueError(msg)

        # ==============================================================
        # Create the "lon" coord as a 1-D vector of values
        # ==============================================================
        if "Xdim" in dset.dims:
            dset = dset.rename_dims({"Xdim": "lon"})
        elif "lon" in dset.coords:
            dset = dset.drop_vars("lon")
        dset = dset.assign_coords({
            "lon": np.linspace(1, xdim, xdim, dtype=np.float64)
        })

        # ==============================================================
        # The dset,stack operation combines the nf and Ydim
        # dimensions into a MultiIndex (i.e. a list of tuples,
        # where each tuple is (face number, cell number).
        # We then have to unpack that into a linear list that
        # ranges from 1..nf*ydim.
        # ==============================================================
        if "nf" in dset.dims and "Ydim" in dset.dims:
            dset = dset.stack(lat=("nf", "Ydim"))
            multi_index_list = dset.lat.values
            lats = np.zeros(6 * ydim)             # 6 cubed-sphere faces
            for i, tpl in enumerate(multi_index_list):
                lats[i] = (tpl[1] + (tpl[0] * ydim)) + 1
            dset = dset.assign_coords({"lat": lats})

        # ==============================================================
        # Transpose dimensions
        # ==============================================================
        dset = order_dims_time_lev_lat_lon(dset)

        # ==============================================================
        # Drop coordinates not needed in checkpoint format files
        # ==============================================================
        if "lons" in dset.variables:
            dset = dset.drop_vars("lons")
        if "lats" in dset.variables:
            dset = dset.drop_vars("lats")

        return dset


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

    --verbose
        Toggles verbose printout on (True) or off (False).

    -w --weightsdir
        Directory where regridding weights are stored (or will be created)
    """

    # Tell parser which arguments to expect
    parser = argparse.ArgumentParser(
        description="General cubed-sphere to cubed-sphere regridder."
    )
    parser.add_argument(
        "-i", "--filein",
        metavar="FILEIN",
        type=str,
        required=True,
        help="input NetCDF file"
    )
    parser.add_argument(
        "-o", "--fileout",
        metavar="FILEOUT",
        type=str,
        required=True,
        help="name of output file"
    )
    parser.add_argument(
        "--sg_params_in",
        metavar="P",
        type=float,
        nargs=3,
        default=[1.0, 170.0, -90.0],
        help="input grid stretching parameters (stretch-factor, target longitude, target latitude)"
    )
    parser.add_argument(
        "--sg_params_out",
        metavar="P",
        type=float,
        nargs=3,
        default=[1.0, 170.0, -90.0],
        help="output grid stretching parameters (stretch-factor, target longitude, target latitude)"
    )
    parser.add_argument(
        "--cs_res_out",
        metavar="RES",
        type=int,
        required=False,
        help="output grid\"s cubed-sphere resolution"
    )
    parser.add_argument(
        "--ll_res_out",
        metavar="RES",
        type=str,
        required=False,
        help="output grid\"s lat/lon resolution in \"latxlon\" format"
    )
    parser.add_argument(
        "--dim_format_in",
        metavar="WHICH",
        type=str,
        choices=[
            "checkpoint",
            "diagnostic",
            "classic"],
        required=True,
        help="format of the input file's dimensions (choose from: checkpoint, diagnostic)"
    )
    parser.add_argument(
        "--dim_format_out",
        metavar="WHICH",
        type=str,
        choices=[
            "checkpoint",
            "diagnostic",
            "classic"],
        required=True,
        help="format of the output file's dimensions (choose from: checkpoint, diagnostic)"
    )
    parser.add_argument(
        "--verbose",
        metavar="VERB",
        type=bool,
        default=False,
        help="Toggles verbose output on (True) or off (False)"
    )
    parser.add_argument(
        "-w", "--weightsdir",
        metavar="WGT",
        type=str,
        default=".",
        help="Directory where regridding weights are found (or will be created)"
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
        args.verbose,
        args.weightsdir
    )


# Only call when run as standalone
if __name__ == "__main__":
    main()
