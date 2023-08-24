"""
Contains functions to make sure that data files to be read by
HEMCO adhere to the COARDS netCDF conventions.
"""
from os.path import join
from copy import deepcopy as dc
import xarray as xr
import numpy as np
import pandas as pd
from gcpy.util import verify_variable_type


def format_hemco_dimensions(
        dset,
        start_time="2000-01-01 00:00:00",
        lev_long_name="level",
        lev_units="level",
        lev_formula_terms=None,
        gchp=False
):
    """
    Formats time, lat, lon, and lev (optionally) attributes for coards
    compliance (HEMCO compatibility).

    Args:
        dset: xarray Dataset
            Dataset containing at least latitude and longitude
            variables, which must be named lat and lon, respectively.

    Keyword Args (optional):
        start_time: string of the format "YYYY-MM-DD HH:mm:ss"
            String containing the start time of the dataset for
            the purposes of encoding the time dimension. For GCHP
            compliance, the first time value must be 0 time units
            from the beginning of the unit. The default value is
            January 1, 2000.
        lev_long_name: string
            A detailed description of the level attribute. Examples
            include "level", "GEOS-Chem levels", "Eta centers", or
            "Sigma centers". Default is "level."
        lev_units: string
            The unit of the vertical levels, which should be "level",
            "eta_level", or "sigma_level". Setting both lev_units and
            lev_long_name to "level" allows HEMCO to regrid between
            vertical grids. Default is "level".
        lev_formula_terms: string or None
            If data is used that is not on the model vertical grid, the
            data must contain surface pressure values and the hybrid
            coefficients of the coordinate system together with the
            terms in the formula(e.g., ”ap: hyam b: hybm ps: PS”).
            Default is None.
        gchp: boolean
            Boolean identifying whether this file is for use in
            GCHP (True) or GEOS-Chem Classic (False). This is primarily
            used to set the lev attributes. The default value is
            False.

    Returns:
        dset: xarray Dataset
            An updated version of dset with encoding and attributes
            set to be coards/HEMCO compliant.
    """
    # Require that dset is an xarray Dataset object
    verify_variable_type(dset, xr.Dataset)

    # Check that latitude and longitude are found in the dataset
    ## First force all dimension names to be lowercase:
    dset = dset.rename_dims({k : k.lower() for k in dset.dims.keys()
                         if k != k.lower()})

    # Check and format each of the required dimensions
    dset = _format_lat(dset)
    dset = _format_lon(dset)
    dset = _format_time(dset, start_time)

    # If level is included in the dimensions, set its attributes
    if "lev" in dset.coordset:
        # Note: this is relatively untested (2023/08/21 HON)
        dset = _format_lev(dset, lev_long_name, lev_units,
                         lev_formula_terms, gchp)

    # Require data order to be time, lat, lon (optionally lev)
    dset = dset.transpose("time", "lat", "lon", ...)

    # Return the dataset
    return dset


def _update_variable_attributes(
        var_attrs,
        coards_attrs
):
    """
    Adds COARDS conforming variable attributes and/or replaces
    existing variable attributes with COARDS-conforming values.

    Args:
        var_attrs : dict
            Dictionary of variable attributes.
        coards_attrs : dict
            Dictionary of COARDS-conforming variable attributes.

    Returns
        var_attrs : dict
           Modified dictionary of variable attributes
    """
    verify_variable_type(var_attrs, dict)
    verify_variable_type(coards_attrs, dict)

    # Test if each COARDS-conforming attribute is
    # present in the list of variable attributes.
    found = {}
    for (name, _) in coards_attrs.items():
        found[name] = name in var_attrs.keys()

    # If the variable attribute has a COARDS-conforming name,
    # then replace it with a COARDS-conforming attribute value.
    #
    # If the variable attribute is missing, then add the
    # COARDS-conforming attribute to the list of variable attrs.
    #
    # This makes sure that we add/replace variable attrs
    # but do not clobber any other existing variable attrs.
    for (name, value) in coards_attrs.items():
        if found[name]:
            var_attrs.update({name: value})
        else:
            var_attrs[name] = value

    return var_attrs


def _format_lat(dset):
    '''
    Formats the latitude dimension for coards compliance.
    See define_HEMCO_dimensions for argument listings.
    '''
    # If there is a dimension is called latitude, rename it
    # (This function assumes ds has dimension names that are
    # all lower case)
    if "latitude" in dset.dims.keys():
        dset = dset.rename_dims({"latitude" : "lat"})

    # Require that lat is a monotonically increasing dimension
    _check_required_dim(dset, "lat")

    # Update attributes to be COARDS-conforming
    dset["lat"].attrs = _update_variable_attributes(
        dset["lat"].attrs,
        coards_attrs={
            "long_name": "latitude",
            "units": "degrees_north",
            "axis" : "Y"
        }
    )

    return dset


def _format_lon(
        dset
):
    '''
    Formats the longitude dimension for coards compliance.
    See define_HEMCO_dimensions for argument listings.
    '''
    # If there is a dimension is called longitude, rename it
    # (This function assumes dset has dimension names that are
    # all lower case)
    if "longitude" in dset.dims.keys():
        dset = dset.rename_dims({"longitude" : "lon"})

    # Require that lon is a monotonically increasing dimension
    _check_required_dim(dset, "lat")

    # Update attributes to be COARDS-conforming
    dset["lon"].attrs = _update_variable_attributes(
        dset["lon"].attrs,
        coards_attrs={
            "long_name": "longitude",
            "units": "degrees_east",
            "axis" : "X"
        }
    )

    return dset


def _format_time(
        dset,
        start_time
):
    '''
    Formats the time dimension for COARDS compliance.
    See define_HEMCO_dimensions for argument listings.
    '''
    if "time" not in dset.coordset:
        # If time isn't already in the coordset, create a dummy variable
        dset = dset.assign_coordset(time=pd.to_datetime(start_time))
        dset = dset.expand_dims("time")
    else:
        # Otherwise, update start_time to match the first time in the file,
        # consistent with GCHP requirements
        new_start_time = pd.to_datetime(dset["time"][0].values)
        new_start_time = new_start_time.strftime("%Y-%m-%d %H:%M:%S")
        print("Updating the reference start time from")
        print(f"{start_time} to {new_start_time}")
        print("so that time(0) = 0, consistent with GCHP requirements.")
        start_time = new_start_time

    # Now check that time is a monotonically increasing dimension
    _check_required_dim(dset, "time")

    # Set attributes and make sure they are COARDS conforming.
    dset["time"].encoding= {
        "units" : f"hours since {start_time}",
        "calendar" : "standard"
    }
    dset["time"].attrs = _update_variable_attributes(
        dset["time"].attrs,
        coards_attrs={
            "long_name": "Time",
            "axis" : "T"
        }
    )

    return dset


def _format_lev(
        dset,
        lev_long_name,
        lev_units,
        lev_formula_terms,
        gchp
):
    '''
    Formats the level dimension for COARDS compliance.
    See define_HEMCO_dimensions for argument listings.
    '''
    ## HON 2023/08/22: This is relatively untested

    # If there a dimension called level, rename it
    if "level" in dset.dims.keys():
        dset = dset.rename_dims({"level" : "lev"})

    # If formula is provided, check that the components of the
    # formula are included.
    if lev_formula_terms is not None:
        terms = lev_formula_terms.split(": ")
        terms = [term for i, term in enumerate(terms) if i % 2 == 1]
        for term in terms:
            if term not in dset.data_vars.keys():
                raise ValueError(
                    f"{term} is in lev_formula_terms and could \
                    not be found."
                )

    # If unit is level, require that the levels are integers
    if lev_units == "level" and \
       (dset["lev"] != dset["lev"].astype(int)).any():
        raise ValueError("lev has units of level but dimension values \
                            are not integers.")

    # Set attributes
    ## Set positive to match the GCHP/GEOS-Chem conventions
    positive = "up"
    if gchp:
        positive = "down"

    ## Set attributes and make sure they are COARDS-conforming.
    ## Setting both long_name and units to "level" allows HEMCO
    ## to regrid between vertical gridset (e.g., 47 -> 72 levels).
    dset["lev"].attrs = _update_variable_attributes(
        dset["lev"].attrs,
        coards_attrs={
            "long_name" : lev_long_name,
            "units" : lev_units,
            "positive" : positive,
            "axis" : "Z"
        }
    )
    if lev_formula_terms is not None:
        dset["lev"].attrs.update({
            "formula_terms" : lev_formula_terms
        })

    return dset


def _check_required_dim(
        dset,
        dim
):
    '''
    Checks required dimensions (time, latitude, and longitude)
    for COARDS compliance (that the dimension exists and is
    monotonically increasing).

    Args:
        dset: xarray Dataset
        dim: string ("time", "lat", or "lon")
            A string corresponding to the required dimension
    '''
    if dim not in ["time", "lat", "lon"]:
        raise ValueError(f"{dim} is not a required dimension.")

    # Check that the dim is included in
    if dim not in dset.dims.keys():
        raise ValueError(f"{dim} is not included in the dimensions.")

    # Require that the variable is monotonically increasing
    if np.any(np.diff(dset[dim]).astype("float") < 0):
        raise ValueError(f"{dim} is not monotonically increasing.")

    return dset


def format_hemco_variable(
        dset,
        var,
        long_name,
        units,
        **kwargs
):
    """
    Formats attributes for non-standard variables for COARDS compliance
    (HEMCO compatibility).

    Args:
        dset: xarray Dataset
            Dataset containing HEMCO input data.
        var: string
            The name of the non-standard variable to be formatted.
        long_name: string
            A required HEMCO attribute, a more descriptive name for
            var.
        units: string
            A required HEMCO attribute giving the units of var. See
            https://hemco.readthedocs.io/en/stable/hco-ref-guide/input-file-format.html
            for more information.
        **kwargs : dict
            Any other attributes wanted for the variable.

    Returns:
        dset: xarray Dataset
            An updated version of dset with variable attributes
            set to be COARDS/HEMCO compliant.
    """
    verify_variable_type(dset, xr.Dataset)
    verify_variable_type(var, str)
    verify_variable_type(long_name, str)
    verify_variable_type(units, str)

    # Add extra attributes if passed via **kwargs
    if len(kwargs) != 0:
        for (_, att_dict) in kwargs.items():
            dset[var].attrs.update(att_dict)

    # Update variable attributes to be COARDS-conforming
    # without clobbering any pre-existing attributes
    dset[var].attrs = _update_variable_attributes(
        dset[var].attrs,
        coards_attrs={
            "long_name" : long_name,
            "units" : units
        }
    )
    return dset


def save_hemco_netcdf(
        dset,
        save_dir,
        save_name,
        dtype="float",
        **kwargs
):
    """
    Saves COARDS compliant (HEMCO compatible) netcdf.

    Args:
        dset: xarray Dataset
            Dataset containing HEMCO input data.
        save_dir: string
            The directory where the data will be saved.
        save_name: string
            The name the file will be named under.

    Keyword Args (optional):
        dtype: data type
            The data type the data will be saved as. Default is
            float32 to minimize memory usage.
        kwargs: dictionary
            Any other attributes to be passed to the xarray
            to_netcdf function.
    """
    verify_variable_type(dset, xr.Dataset)
    verify_variable_type(save_dir, str)
    verify_variable_type(save_name, str)

    # Check that the save_name endset in .nc
    if save_name.split(".")[-1][:2] != "nc":
        save_name = f"{save_name}.nc"

    # Get time encoding before overwriting
    time_units = dset["time"].encoding["units"]
    calendar = dset["time"].encoding["calendar"]

    # Set default encoding and dtype for all variables and coordinates
    encoding = {"_FillValue" : None, "dtype" : dtype}
    var = {k : dc(encoding) for k in dset.keys()}
    coord = {k : dc(encoding) for k in dset.coordset}

    # Manually update the time encoding, which is often overwritten
    # by xarray defaults
    coord["time"]["units"] = time_units
    coord["time"]["calendar"] = calendar
    var.update(coord)

    # Save out
    dset.to_netcdf(
        join(save_dir, save_name),
        encoding=var,
        unlimited_dims=["time"],
        **kwargs
    )

    print("-"*70)
    print("Saved to", join(save_dir, save_name))
    print("-"*70)
