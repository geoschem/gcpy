"""
Contains functions to make sure that data files to be read by
HEMCO adhere to the COARDS netCDF conventions.
"""
from os.path import join
from copy import deepcopy as dc
import warnings
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
    Formats time, lat, lon, and lev attributes for COARDS compliance
    (HEMCO compatibility).

    Parameters
    ----------
    dset : xarray.Dataset
        Dataset containing at least latitude and longitude variables,
        which must be named ``lat`` and ``lon`` respectively.
    start_time : str, optional
        Start time of the dataset used to encode the time dimension,
        formatted as ``"YYYY-MM-DD HH:mm:ss"``. For GCHP compliance,
        the first time value must be 0 time units from the start of
        the unit. Default: ``"2000-01-01 00:00:00"``
    lev_long_name : str, optional
        Descriptive name for the level attribute. Examples include
        ``"level"``, ``"GEOS-Chem levels"``, ``"Eta centers"``, or
        ``"Sigma centers"``. Default: ``"level"``
    lev_units : str, optional
        Units of the vertical levels. Should be one of ``"level"``,
        ``"eta_level"``, or ``"sigma_level"``. Setting both
        ``lev_units`` and ``lev_long_name`` to ``"level"`` allows
        HEMCO to regrid between vertical grids.
        Default: ``"level"``
    lev_formula_terms : str or None, optional
        If data is not on the model vertical grid, this must contain
        the surface pressure values and hybrid coefficients of the
        coordinate system together with the formula terms (e.g.
        ``"ap: hyam b: hybm ps: PS"``). Default: ``None``
    gchp : bool, optional
        Whether this file is intended for use in GCHP (``True``) or
        GEOS-Chem Classic (``False``). Primarily used to set the
        ``lev`` attributes. Default: ``False``

    Returns
    -------
    dset : xarray.Dataset
        Updated dataset with encoding and attributes set for
        COARDS/HEMCO compliance.
    """
    # Require that dset is an xarray Dataset object
    verify_variable_type(dset, xr.Dataset)

    # Force all dimension names to be lowercase
    dset = dset.rename_dims(
        {k: k.lower() for k in dset.dims.keys() if k != k.lower()}
    )

    # Check and format each of the required dimensions
    dset = _format_lat(dset)
    dset = _format_lon(dset)
    dset = _format_time(dset, start_time)

    # If level is included in the dimensions, set its attributes
    if "lev" in dset.coords or "level" in dset.coords:
        # Note: this is relatively untested (2023/08/21 HON)
        dset = _format_lev(
            dset, lev_long_name, lev_units, lev_formula_terms, gchp
        )

    # Require data order to be time, lat, lon (optionally lev)
    if "lev" in dset.coords:
        dset = dset.transpose("time", "lev", "lat", "lon", ...)
    else:
        dset = dset.transpose("time", "lat", "lon", ...)

    return dset


def _update_variable_attributes(var_attrs, coards_attrs):
    """
    Adds COARDS-conforming variable attributes and/or replaces
    existing variable attributes with COARDS-conforming values.

    Existing attributes that are not listed in ``coards_attrs`` are
    preserved (i.e. not clobbered).

    Parameters
    ----------
    var_attrs : dict
        Dictionary of existing variable attributes.
    coards_attrs : dict
        Dictionary of COARDS-conforming variable attributes to add
        or update.

    Returns
    -------
    var_attrs : dict
        Modified dictionary of variable attributes with COARDS
        values added or updated.
    """
    verify_variable_type(var_attrs, dict)
    verify_variable_type(coards_attrs, dict)

    # Test if each COARDS-conforming attribute is present in the
    # list of variable attributes.
    found = {}
    for (name, _) in coards_attrs.items():
        found[name] = name in var_attrs.keys()

    # If the variable attribute has a COARDS-conforming name,
    # replace it with a COARDS-conforming value.
    # If the attribute is missing, add the COARDS-conforming
    # attribute without clobbering any other existing attributes.
    for (name, value) in coards_attrs.items():
        if found[name]:
            if var_attrs[name] != value:
                print(f"Updating attribute value for {name}:")
                print(f"    Original value : {var_attrs[name]}")
                print(f"    New value:     : {value}")
            var_attrs.update({name: value})
        else:
            var_attrs[name] = value

    return var_attrs


def _format_lat(dset):
    """
    Formats the latitude dimension for COARDS compliance.

    Renames ``"latitude"`` to ``"lat"`` if necessary, checks that
    ``"lat"`` is monotonically increasing, and sets COARDS-conforming
    attributes.

    Parameters
    ----------
    dset : xarray.Dataset
        Dataset with dimension names already converted to lowercase.

    Returns
    -------
    dset : xarray.Dataset
        Updated dataset with COARDS-conforming latitude attributes.
    """
    # If there is a dimension called latitude, rename it
    if "latitude" in dset.dims.keys():
        dset = dset.rename_dims({"latitude": "lat"})

    # Require that lat is a monotonically increasing dimension
    _check_required_dim(dset, "lat")

    # Update attributes to be COARDS-conforming
    dset["lat"].attrs = _update_variable_attributes(
        dset["lat"].attrs,
        coards_attrs={
            "long_name": "latitude",
            "units": "degrees_north",
            "axis": "Y",
        }
    )

    return dset


def _format_lon(dset):
    """
    Formats the longitude dimension for COARDS compliance.

    Renames ``"longitude"`` to ``"lon"`` if necessary, checks that
    ``"lon"`` is monotonically increasing, and sets COARDS-conforming
    attributes.

    Parameters
    ----------
    dset : xarray.Dataset
        Dataset with dimension names already converted to lowercase.

    Returns
    -------
    dset : xarray.Dataset
        Updated dataset with COARDS-conforming longitude attributes.
    """
    # If there is a dimension called longitude, rename it
    if "longitude" in dset.dims.keys():
        dset = dset.rename_dims({"longitude": "lon"})

    # Require that lon is a monotonically increasing dimension
    _check_required_dim(dset, "lat")

    # Update attributes to be COARDS-conforming
    dset["lon"].attrs = _update_variable_attributes(
        dset["lon"].attrs,
        coards_attrs={
            "long_name": "longitude",
            "units": "degrees_east",
            "axis": "X",
        }
    )

    return dset


def _format_time(dset, start_time):
    """
    Formats the time dimension for COARDS compliance.

    If no time coordinate is present, a dummy time coordinate is
    created from ``start_time``. Otherwise, ``start_time`` is updated
    to match the first time in the dataset, consistent with GCHP
    requirements (``time(0) = 0``).

    Parameters
    ----------
    dset : xarray.Dataset
        Dataset to be updated.
    start_time : str
        Reference start time formatted as
        ``"YYYY-MM-DD HH:mm:ss"``.

    Returns
    -------
    dset : xarray.Dataset
        Updated dataset with COARDS-conforming time encoding and
        attributes.
    """
    if "time" not in dset.coords:
        # If time isn't already in the coords, create a dummy variable
        print(
            f"Assigning time coordinate from input start_time {start_time}."
        )
        dset = dset.assign_coords(time=pd.to_datetime(start_time))
        dset = dset.expand_dims("time")
    else:
        # Otherwise, update start_time to match the first time in the
        # file, consistent with GCHP requirements
        new_start_time = pd.to_datetime(dset["time"][0].values)
        new_start_time = new_start_time.strftime("%Y-%m-%d %H:%M:%S")
        print("Updating the reference start time from")
        print(f"{start_time} to {new_start_time}")
        print("so that time(0) = 0, consistent with GCHP requirements.")
        start_time = new_start_time

    # Check that time is a monotonically increasing dimension
    _check_required_dim(dset, "time")

    # Set encoding and attributes to be COARDS-conforming
    dset["time"].encoding = {
        "units": f"hours since {start_time}",
        "calendar": "standard",
    }
    dset["time"].attrs = _update_variable_attributes(
        dset["time"].attrs,
        coards_attrs={
            "long_name": "Time",
            "axis": "T",
        }
    )

    return dset


def _format_lev(dset, lev_long_name, lev_units, lev_formula_terms, gchp):
    """
    Formats the level dimension for COARDS compliance.

    Renames ``"level"`` to ``"lev"`` if necessary, validates
    ``lev_formula_terms``, and sets COARDS-conforming level
    attributes.

    Parameters
    ----------
    dset : xarray.Dataset
        Dataset to be updated.
    lev_long_name : str
        Descriptive name for the level attribute (e.g. ``"level"``,
        ``"GEOS-Chem levels"``, ``"Eta centers"``).
    lev_units : str
        Units of the vertical levels. Must be one of ``"level"``,
        ``"eta_level"``, or ``"sigma_level"``.
    lev_formula_terms : str or None
        Formula terms string for hybrid-sigma coordinates (e.g.
        ``"ap: hyam b: hybm ps: PS"``), or ``None``.
    gchp : bool
        Whether the file is for GCHP (``True``) or GEOS-Chem Classic
        (``False``). Sets ``lev:positive`` to ``"down"`` for GCHP
        and ``"up"`` for GEOS-Chem Classic.

    Returns
    -------
    dset : xarray.Dataset
        Updated dataset with COARDS-conforming level attributes.

    Raises
    ------
    ValueError
        If ``lev_units`` is not one of ``"level"``, ``"eta_level"``,
        or ``"sigma_level"``.

    Warns
    -----
    UserWarning
        If both ``lev_formula_terms`` and ``lev['formula_terms']``
        are provided simultaneously.
    UserWarning
        If neither ``lev_formula_terms`` nor
        ``lev['formula_terms']`` are provided.
    UserWarning
        If any term listed in ``lev_formula_terms`` is not found
        in the dataset variables.

    Notes
    -----
    This function is relatively untested as of 2023/08/22 (HON).
    """
    # If there is a dimension called level, rename it
    if "level" in dset.dims.keys():
        dset = dset.rename_dims({"level": "lev"})

    # Check whether both lev_formula_terms and lev["formula_terms"]
    # are present -- if so, warn and use the argument value.
    if (lev_formula_terms is not None) and \
            ("formula_terms" in dset["lev"].attrs):
        warnings.warn(
            "Both lev_formula_terms and lev['formula_terms'] are provided."
            " The provided lev_formula_term is being used."
        )
    elif (lev_formula_terms is None) and \
            ("formula_terms" not in dset["lev"].attrs):
        warnings.warn(
            "Neither lev_formula_terms nor lev['formula_terms'] are provided."
            " Skipping lev_formula_terms formatting."
        )
    elif "formula_terms" in dset["lev"].attrs:
        lev_formula_terms = dset["lev"].attrs["formula_terms"]

    # If lev_formula_terms is now defined, check that all referenced
    # variables are present in the dataset
    if lev_formula_terms is not None:
        terms = lev_formula_terms.split(" ")
        terms = [term for i, term in enumerate(terms) if i % 2 == 1]
        failed_terms = [
            term for term in terms if term not in dset.data_vars.keys()
        ]
        if len(failed_terms) > 0:
            warnings.warn(
                f"The following values are in lev_formula_terms and could"
                f" not be found: {failed_terms}"
            )

    # Validate lev_units
    if lev_units not in ["level", "eta_level", "sigma_level"]:
        raise ValueError(
            f"lev has units of {lev_units}. Please set it to one "
            "of level, eta_level, or sigma_level."
        )

    # Set positive attribute following GCHP/GEOS-Chem conventions
    positive = "down" if gchp else "up"

    # Setting both long_name and units to "level" allows HEMCO to
    # regrid between vertical grids (e.g. 47 -> 72 levels).
    dset["lev"].attrs = _update_variable_attributes(
        dset["lev"].attrs,
        coards_attrs={
            "long_name": lev_long_name,
            "units": lev_units,
            "positive": positive,
            "axis": "Z",
        }
    )
    if lev_formula_terms is not None:
        dset["lev"].attrs.update({"formula_terms": lev_formula_terms})

    return dset


def _check_required_dim(dset, dim):
    """
    Checks that a required dimension exists and is monotonically
    increasing.

    Parameters
    ----------
    dset : xarray.Dataset
        Dataset to be validated.
    dim : str
        Name of the required dimension. Must be one of
        ``"time"``, ``"lat"``, or ``"lon"``.

    Returns
    -------
    dset : xarray.Dataset
        The input dataset, unchanged if all checks pass.

    Raises
    ------
    ValueError
        If ``dim`` is not one of ``"time"``, ``"lat"``, or ``"lon"``.
    ValueError
        If ``dim`` is not present in the dataset dimensions.
    ValueError
        If ``dim`` is not monotonically increasing.
    """
    if dim not in ["time", "lat", "lon"]:
        raise ValueError(f"{dim} is not a required dimension.")

    if dim not in dset.dims.keys():
        raise ValueError(f"{dim} is not included in the dimensions.")

    if np.any(np.diff(dset[dim]).astype("float") < 0):
        raise ValueError(f"{dim} is not monotonically increasing.")

    return dset


def check_hemco_variables(dset):
    """
    Checks that all variables in a dataset have the attributes
    required for HEMCO compliance (``units`` and ``long_name``).

    Parameters
    ----------
    dset : xarray.Dataset
        Dataset whose variables are to be checked.

    Raises
    ------
    ValueError
        If any variable is missing the ``units`` or ``long_name``
        attribute.
    """
    verify_variable_type(dset, xr.Dataset)

    print("Checking dataset variables for HEMCO compliance.")
    required_attrs = ["units", "long_name"]
    missing = False
    for (name, _) in dset.items():
        attr_names = [name for (name, _) in dset[name].attrs.items()]
        missing_attrs = [
            name for name in required_attrs if name not in attr_names
        ]
        if len(missing_attrs) > 0:
            missing = True
            print(f"  {name} missing {missing_attrs}")

    if missing:
        raise ValueError(
            "Required units missing from dataset variables."
        )
    print("Dataset variables are HEMCO compliant.")


def format_hemco_variable(
        dset,
        var,
        long_name=None,
        units=None,
        **kwargs
):
    """
    Formats attributes for a non-standard variable for COARDS
    compliance (HEMCO compatibility).

    Parameters
    ----------
    dset : xarray.Dataset
        Dataset containing HEMCO input data.
    var : str
        Name of the variable to be formatted.
    long_name : str, optional
        Descriptive name for ``var``. Required by HEMCO unless
        already present as a variable attribute.
        Default: ``None``
    units : str, optional
        Units of ``var``. Required by HEMCO unless already present
        as a variable attribute. See the `HEMCO input file format
        documentation
        <https://hemco.readthedocs.io/en/stable/hco-ref-guide/input-file-format.html>`_
        for more information. Default: ``None``
    **kwargs : dict
        Additional attributes to set on the variable.

    Returns
    -------
    dset : xarray.Dataset
        Updated dataset with COARDS/HEMCO-conforming variable
        attributes.

    Raises
    ------
    ValueError
        If ``long_name`` or ``units`` is not provided and cannot be
        found in the existing variable attributes.
    """
    verify_variable_type(dset, xr.Dataset)
    verify_variable_type(var, str)

    # Check required attributes
    coards_attrs = {"long_name": long_name, "units": units}
    for name, value in coards_attrs.items():
        if value is not None:
            verify_variable_type(value, str)
        elif name in dset[var].attrs:
            coards_attrs[name] = dset[var].attrs[name]
        else:
            raise ValueError(f"{name} is not defined for {var}")

    # Update variable attributes to be COARDS-conforming without
    # clobbering any pre-existing attributes
    dset[var].attrs = _update_variable_attributes(
        dset[var].attrs,
        coards_attrs=coards_attrs
    )

    # Add extra attributes if passed via **kwargs
    if len(kwargs) != 0:
        for (_, att_dict) in kwargs.items():
            dset[var].attrs.update(att_dict)

    return dset


def save_hemco_netcdf(
        dset,
        save_dir,
        save_name,
        dtype="float",
        **kwargs
):
    """
    Saves a COARDS-compliant (HEMCO-compatible) netCDF file.

    Parameters
    ----------
    dset : xarray.Dataset
        Dataset containing HEMCO input data.
    save_dir : str
        Directory in which the file will be saved.
    save_name : str
        Filename for the output file. A ``.nc`` extension will be
        appended if not already present.
    dtype : str or numpy.dtype, optional
        Data type used when writing data to disk. Defaults to
        ``"float"`` (float32) to minimise file size.
    **kwargs : dict
        Additional keyword arguments passed to
        :func:`xarray.Dataset.to_netcdf`.

    Notes
    -----
    The time encoding (``units`` and ``calendar``) is preserved
    explicitly to prevent xarray from overwriting it with its own
    defaults during the save step.
    """
    verify_variable_type(dset, xr.Dataset)
    verify_variable_type(save_dir, str)
    verify_variable_type(save_name, str)

    # Ensure the filename ends with .nc
    if save_name.split(".")[-1][:2] != "nc":
        save_name = f"{save_name}.nc"

    # Preserve time encoding before it can be overwritten
    time_units = dset["time"].encoding["units"]
    calendar   = dset["time"].encoding["calendar"]

    # Set default encoding and dtype for all variables and coordinates
    encoding = {"_FillValue": None, "dtype": dtype}
    var   = {k: dc(encoding) for k in dset.keys()}
    coord = {k: dc(encoding) for k in dset.coords}

    # Manually restore the time encoding, which xarray often overwrites
    coord["time"]["units"]    = time_units
    coord["time"]["calendar"] = calendar
    var.update(coord)

    # Save the dataset
    dset.to_netcdf(
        join(save_dir, save_name),
        encoding=var,
        unlimited_dims=["time"],
        **kwargs
    )

    print("-" * 70)
    print("Saved to", join(save_dir, save_name))
    print("-" * 70)
