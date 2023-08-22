
import xarray as xr
import numpy as np
import pandas as pd
from copy import deepcopy as dc
from os.path import join

def format_HEMCO_dimensions(ds, 
                            start_time="2000-01-01 00:00:00",
                            lev_long_name="level", 
                            lev_units="level",
                            lev_formula_terms=None,
                            gchp=False):
    """
    Formats time, lat, lon, and lev (optionally) attributes for coards 
    compliance (HEMCO compatibility).
    
    Args:
        ds: xarray Dataset
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
            data must contain surface pressure values and the hybrid coefficients
            of the coordinate system together with the terms in the formula
            (e.g., ”ap: hyam b: hybm ps: PS”). Default is None.
        gchp: boolean
            Boolean identifying whether this file is for use in 
            GCHP (True) or GEOS-Chem Classic (False). This is primarily
            used to set the lev attributes. The default value is 
            False.
    
    Returns:
        ds: xarray Dataset
            An updated version of ds with encoding and attributes
            set to be coards/HEMCO compliant.
    """
    # Require that ds is an xarray Dataset object
    if not isinstance(ds, xr.Dataset):
        raise TypeError("The ds argument must be an xarray Dataset.")

    # Check that latitude and longitude are found in the dataset
    ## First force all dimension names to be lowercase:
    ds = ds.rename_dims({k : k.lower() for k in ds.dims.keys()
                         if k != k.lower()})

    # Check and format each of the required dimensions
    ds = _format_lat(ds)
    ds = _format_lon(ds)
    ds = _format_time(ds, start_time)

    # If level is included in the dimensions, set its attributes
    if "lev" in ds.coords:
        # Note: this is relatively untested (2023/08/21 HON)
        ds = _format_lev(ds, lev_long_name, lev_units,
                         lev_formula_terms, gchp)
    
    # Require data order to be time, lat, lon (optionally lev)
    ds = ds.transpose("time", "lat", "lon", ...)

    # Return the dataset
    return ds


def _format_lat(ds):
    '''
    Formats the latitude dimension for coards compliance.
    See define_HEMCO_dimensions for argument listings.
    '''
    # If there is a dimension is called latitude, rename it
    # (This function assumes ds has dimension names that are 
    # all lower case)
    if "latitude" in ds.dims.keys():
        ds = ds.rename_dims({"latitude" : "lat"})

    # Require that lat is a monotonically increasing dimension
    _check_required_dim(ds, "lat")

    # Set attributes
    ds["lat"].attrs = {"long_name": "latitude", 
                       "units": "degrees_north",
                       "axis" : "Y"}

    return ds


def _format_lon(ds):
    '''
    Formats the longitude dimension for coards compliance.
    See define_HEMCO_dimensions for argument listings.
    '''
    # If there is a dimension is called longitude, rename it
    # (This function assumes ds has dimension names that are 
    # all lower case)
    if "longitude" in ds.dims.keys():
        ds = ds.rename_dims({"longitude" : "lon"})

    # Require that lon is a monotonically increasing dimension
    _check_required_dim(ds, "lat")

    # Set attributes
    ds["lon"].attrs = {"long_name": "longitude", 
                       "units": "degrees_east",
                       "axis" : "X"}
    
    return ds


def _format_time(ds, start_time):
    '''
    Formats the time dimension for coards compliance.
    See define_HEMCO_dimensions for argument listings.
    '''
    if "time" not in ds.coords:
        # If time isn't already in the coords, create a dummy variable
        ds = ds.assign_coords(time=pd.to_datetime(start_time))
        ds = ds.expand_dims("time")
    else:
        # Otherwise, update start_time to match the first time in the file,
        # consistent with GCHP requirements
        new_start_time = pd.to_datetime(ds["time"][0].values)
        new_start_time = new_start_time.strftime("%Y-%m-%d %H:%M:%S")
        print(f"Updating the reference start time from")
        print(f"{start_time} to {new_start_time}")
        print(f"so that time(0) = 0, consistent with GCHP requirements.")
        start_time = new_start_time

    # Now check that time is a monotonically increasing dimension
    _check_required_dim(ds, "time")

    # Set attributes
    ds["time"].encoding= {"units" : f"hours since {start_time}",
                          "calendar" : "standard"}
    ds["time"].attrs = {"long_name" : "Time", "axis" : "T"}

    return ds


def _format_lev(ds, lev_long_name, lev_units, lev_formula_terms, gchp):
    '''
    Formats the level dimension for coards compliance.
    See define_HEMCO_dimensions for argument listings.
    '''
    ## HON 2023/08/22: This is relatively untested

    # If there a dimension called level, rename it
    if "level" in ds.dims.keys():
        ds = ds.rename_dims({"level" : "lev"})
    
    # If formula is provided, check that the components of the 
    # formula are included.
    if lev_formula_terms is not None:
        terms = lev_formula_terms.split(": ")
        terms = [t for i, t in enumerate(terms) if i % 2 == 1]
        for t in terms:
            if t not in ds.data_vars.keys():
                raise ValueError(f"{t} is in lev_formula_terms and could \
                                    not be found.")
    
    # If unit is level, require that the levels are integers
    if (lev_units == "level") and (ds["lev"] != ds["lev"].astype(int)).any():
        raise ValueError("lev has units of level but dimension values \
                            are not integers.")

    # Set attributes
    ## Set positive to match the GCHP/GEOS-Chem conventions
    if gchp:
        positive = "down"
    else:
        positive = "up"

    ## Setting both long_name and units to "level" allows HEMCO
    ## to regrid between vertical grids (e.g., 47 -> 72 levels).
    lev_attrs = {"long_name" : lev_long_name,
                 "units" : lev_units,
                 "positive" : positive,
                 "axis" : "Z"}
    if lev_formula_terms is not None:
        lev_attrs.update({"formula_terms" : lev_formula_terms})
    
    ## Set the attributes
    ds["lev"].attrs = lev_attrs

    return ds


def _check_required_dim(ds, dim):
    '''
    Checks required dimensions (time, latitude, and longitude)
    for COARDS compliance (that the dimension exists and is
    monotonically increasing).

    Args:
        ds: xarray Dataset
        dim: string ("time", "lat", or "lon")
            A string corresponding to the required dimension
    '''
    if dim not in ["time", "lat", "lon"]:
        raise ValueError(f"{dim} is not a required dimension.")

    # Check that the dim is included in 
    if dim not in ds.dims.keys():
        raise ValueError(f"{dim} is not included in the dimensions.")

    # Require that the variable is monotonically increasing
    if np.any(np.diff(ds[dim]).astype("float") < 0):
        raise ValueError(f"{dim} is not monotonically increasing.")
    
    return ds


def format_HEMCO_variable(ds, var, long_name, units, **kwargs):
    """
    Formats attributes for non-standard variables for coards compliance 
    (HEMCO compatibility).
    
    Args:
        ds: xarray Dataset
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
        kwargs: dictionary
            Any other attributes wanted for the variable.
    
    Returns:
        ds: xarray Dataset
            An updated version of ds with variable attributes
            set to be coards/HEMCO compliant.
    """
    ds[var].attrs = {"long_name" : long_name, "units" : units,
                       **kwargs}
    return ds


def save_HEMCO_netcdf(ds, save_dir, save_name, dtype="float", **kwargs):
    """
    Saves coards compliant (HEMCO compatible) netcdf.
    
    Args:
        ds: xarray Dataset
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
    # Check that the save_name ends in .nc
    if save_name.split(".")[-1][:2] != "nc":
        save_name = f"{save_name}.nc"
    
    # Get time encoding before overwriting
    time_units = ds["time"].encoding["units"]
    calendar = ds["time"].encoding["calendar"]

    # Set default encoding and dtype for all variables and coordinates
    encoding = {"_FillValue" : None, "dtype" : dtype}
    var = {k : dc(encoding) for k in ds.keys()}
    coord = {k : dc(encoding) for k in ds.coords}
    
    # Manually update the time encoding, which is often overwritten
    # by xarray defaults
    coord["time"]["units"] = time_units
    coord["time"]["calendar"] = calendar
    var.update(coord)

    # Save out
    ds.to_netcdf(join(save_dir, save_name), encoding=var,
                 unlimited_dims=["time"], **kwargs)
    
    print("-"*70)
    print("Saved to", join(save_dir, save_name))
    print("-"*70)