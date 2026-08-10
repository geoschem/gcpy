"""
Contains methods for converting the units of data.
Mainly used for model benchmarking purposes.
"""
import numpy as np
import xarray as xr
import gcpy.constants as physconsts
from gcpy.cstools import is_cubed_sphere_diag_grid, is_cubed_sphere_rst_grid
from gcpy.util import reshape_MAPL_CS


def adjust_units(units):
    """
    Creates a consistent unit string used in unit conversion routines.

    Parameters
    ----------
    units : str
        Input unit string.

    Returns
    -------
    adjusted_units : str
        Output unit string, adjusted to a consistent value.

    Raises
    ------
    TypeError
        If units is not of type str.

    Notes
    -----
    Unit list is incomplete -- currently geared to units from common
    model diagnostics (e.g. kg/m2/s, kg, and variants).

    Examples
    --------
    >>> adjust_units("kg m-2 s-1")
    'kg/m2/s'
    >>> adjust_units("kgC/m2/s")
    'kgC/m2/s'
    """
    # Error check arguments
    if not isinstance(units, str):
        raise TypeError("Units must be of type str!")

    # Strip all spaces in the unit string
    units_squeezed = units.replace(" ", "")

    if units_squeezed in ["kg/m2/s", "kgm-2s-1", "kgm^-2s^-1"]:
        unit_desc = "kg/m2/s"

    elif units_squeezed in [
            "kgC/m2/s",
            "kgCm-2s-1",
            "kgCm^-2s^-1",
            "kgc/m2/s",
            "kgcm-2s-1",
            "kgcm^-2s^-1",
    ]:
        unit_desc = "kgC/m2/s"

    elif units_squeezed in ["molec/cm2/s", "moleccm-2s-1", "moleccm^-2s^-1"]:
        unit_desc = "molec/cm2/s"

    else:
        unit_desc = units_squeezed

    return unit_desc


def convert_kg_to_target_units(data_kg, target_units, kg_to_kgC):
    """
    Converts a data array from kg to one of several target unit types.

    Parameters
    ----------
    data_kg : numpy.ndarray
        Input data array in units of kg.
    target_units : str
        Name of the target units. Supported values include:
        ``'Tg'``, ``'Tg C'``, ``'Gg'``, ``'Gg C'``, ``'Mg'``,
        ``'Mg C'``, ``'kg'``, ``'kg C'``, ``'g'``, ``'g C'``.
    kg_to_kgC : float
        Conversion factor from kg species to kg carbon.

    Returns
    -------
    data : numpy.ndarray
        Output data array converted to the units specified by
        ``target_units``.

    Raises
    ------
    ValueError
        If ``target_units`` is not a supported unit string.

    Notes
    -----
    Only unit conversions corresponding to the GEOS-Chem benchmarks
    have been implemented. This is an internal routine intended to be
    called from :func:`convert_units`.

    Examples
    --------
    >>> import numpy as np
    >>> data_kg = np.array([1.0e9, 2.0e9])
    >>> convert_kg_to_target_units(data_kg, 'Tg', 1.0)
    array([1., 2.])
    """
    # Convert to target unit
    if target_units == "Tg":
        data = data_kg * 1e-9

    elif target_units == "Tg C":
        data = data_kg * kg_to_kgC * 1.0e-9

    elif target_units == "Gg":
        data = data_kg * 1e-6

    elif target_units == "Gg C":
        data = data_kg * kg_to_kgC * 1.0e-6

    elif target_units == "Mg":
        data = data_kg * 1e-3

    elif target_units == "Mg C":
        data = data_kg * kg_to_kgC * 1.0e-3

    elif target_units == "kg":
        data = data_kg

    elif target_units == "kg C":
        data = data_kg * kg_to_kgC

    elif target_units == "g":
        data = data_kg * 1e3

    elif target_units == "g C":
        data = data_kg * kg_to_kgC * 1.0e3

    else:
        msg = "Target units {} are not yet supported!".format(target_units)
        raise ValueError(msg)

    return data


def convert_units(
        dr,
        species_name,
        species_properties,
        target_units,
        interval=[2678400.0],
        area_m2=None,
        delta_p=None,
        box_height=None,
):
    """
    Converts data in an xarray DataArray from native units to target units.

    Parameters
    ----------
    dr : xarray.DataArray
        Data to be converted from native units to target units.
    species_name : str
        Name of the species corresponding to the data stored in ``dr``.
    species_properties : dict
        Dictionary containing species properties such as molecular
        weights and other metadata for the given species.
    target_units : str
        Units to which the data will be converted.
    interval : list of float, optional
        Length of the averaging period in seconds.
        Default: ``[2678400.0]``
    area_m2 : xarray.DataArray, optional
        Surface area in square metres.
        Default: ``None``
    delta_p : xarray.DataArray, optional
        Delta-pressure between top and bottom edges of grid box
        (dry air) in hPa.
        Default: ``None``
    box_height : xarray.DataArray, optional
        Grid box height in metres.
        Default: ``None``

    Returns
    -------
    dr_new : xarray.DataArray
        Data converted to ``target_units``. Identical to ``dr`` except
        for the data values and the ``units`` attribute.

    Raises
    ------
    ValueError
        If ``MW_g`` is not found in ``species_properties``.
    ValueError
        If conversion from ``molmol-1dry`` is requested but
        ``area_m2`` or ``delta_p`` are not provided.
    ValueError
        If a mass unit (containing ``'g'``) is requested but
        ``MW_g`` is ``None``.
    ValueError
        If ``interval`` has length > 1 but ``dr`` has no time dimension.
    ValueError
        If the native units of ``dr`` are not supported.

    Notes
    -----
    Only unit conversions corresponding to common GEOS-Chem benchmark
    output have been implemented. When ``molmol-1`` is present as the
    unit, dry air is assumed.

    Examples
    --------
    >>> dr_converted = convert_units(
    ...     dr, 'O3', species_props, 'Tg',
    ...     interval=[2678400.0], area_m2=area
    ... )
    """
    # Get species molecular weight information
    if "MW_g" in species_properties.keys():
        mw_g = species_properties.get("MW_g")
    else:
        msg = "Cannot find molecular weight MW_g for species {}".format(
            species_name)
        msg += "!\nPlease add the MW_g field for {}".format(species_name)
        msg += " to the species_database.yml file."
        raise ValueError(msg)

    # If the species metadata does not contain EmMW_g, use MW_g instead
    if "EmMW_g" in species_properties.keys():
        emitted_mw_g = species_properties.get("EmMW_g")
    else:
        emitted_mw_g = mw_g

    # If the species metadata does not contain MolecRatio, use 1.0 instead
    if "MolecRatio" in species_properties.keys():
        moles_C_per_mole_species = species_properties.get("MolecRatio")
    else:
        moles_C_per_mole_species = 1.0

    # EDGE CASE: For GCHP vs GCHP mass tables, we often need to multiply
    # data taken from a GCHP checkpoint file by the surface area taken
    # from the StateMet History file.  If this is the case, then reshape
    # the area array to the GCHP checkpoint grid before multiplication.
    if is_cubed_sphere_rst_grid(dr) and is_cubed_sphere_diag_grid(area_m2):
        area_m2 = reshape_MAPL_CS(area_m2)

    # ==============================
    # Compute conversion factors
    # ==============================

    # Physical constants
    Avo = physconsts.AVOGADRO
    mw_air = physconsts.MW_AIR_g
    g0 = physconsts.G

    # Get a consistent value for the units string
    units = adjust_units(dr.units)

    # Error checks
    if units == "molmol-1dry" and area_m2 is None:
        raise ValueError(
            "Conversion from {} to {} for {} requires area_m2 as input".format(
                units, target_units, species_name
            )
        )
    if units == "molmol-1dry" and delta_p is None:
        raise ValueError(
            "Conversion from {} to {} for {} requires delta_p as input".format(
                units, target_units, species_name
            )
        )
    if "g" in target_units and mw_g is None:
        raise ValueError(
            "Conversion from {} to {} for {} requires MW_g definition "
            "in species_database.yml".format(units, target_units, species_name)
        )

    # Conversion factor for kg species to kg C
    kg_to_kgC = (emitted_mw_g * moles_C_per_mole_species) / mw_g

    # Mass of dry air in kg (required when converting from v/v)
    if 'molmol-1' in units:
        air_mass = delta_p.values * 100.0 / g0 * area_m2.values

        if "g" in target_units:
            vv_to_kg = air_mass / mw_air * mw_g

        if "molec" in target_units:
            vv_to_MND = (
                air_mass / mw_air * Avo / (area_m2 * box_height) / 1e6
            )

    # ================================================
    # Get number of seconds per time in dataset
    # ================================================
    numsec = interval

    if len(interval) > 1:
        if 'time' in dr.dims:
            numnewdims = len(dr.dims) - (dr.dims.index('time') + 1)
            for _ in range(numnewdims):
                numsec = numsec[:, np.newaxis]
        else:
            raise ValueError(
                'Interval passed to convert_units has length greater than '
                'one but data array has no time dimension'
            )

    # ==============================
    # Compute target units
    # ==============================

    if units == "kg/m2/s":
        data_kg = dr * area_m2
        data_kg = data_kg.values * numsec
        data = convert_kg_to_target_units(data_kg, target_units, kg_to_kgC)

    elif units == "kgC/m2/s":
        data_kg = dr * area_m2 / kg_to_kgC
        data_kg = data_kg.values * numsec
        data = convert_kg_to_target_units(data_kg, target_units, kg_to_kgC)

    elif units == "kg":
        data_kg = dr.values
        data = convert_kg_to_target_units(data_kg, target_units, kg_to_kgC)

    elif units == "kgC":
        data_kg = dr.values / kg_to_kgC
        data = convert_kg_to_target_units(data_kg, target_units, kg_to_kgC)

    elif 'molmol-1' in units:
        if "g" in target_units:
            data_kg = dr.values * vv_to_kg
            data = convert_kg_to_target_units(data_kg, target_units, kg_to_kgC)
        elif "molec" in target_units:
            data = dr.values * vv_to_MND

    else:
        raise ValueError(
            "Units ({}) in variable {} are not supported".format(
                units, species_name
            )
        )

    # ==============================
    # Return result
    # ==============================
    dr_new = xr.DataArray(
        data, name=dr.name, coords=dr.coords, dims=dr.dims, attrs=dr.attrs
    )
    dr_new.attrs["units"] = target_units

    return dr_new


def check_units(ref_da, dev_da, enforce_units=True):
    """
    Ensures the units of two xarray DataArrays are the same.

    Parameters
    ----------
    ref_da : xarray.DataArray
        First data array containing a ``units`` attribute.
    dev_da : xarray.DataArray
        Second data array containing a ``units`` attribute.
    enforce_units : bool, optional
        Whether to raise an AssertionError if ref and dev units do
        not match.
        Default: ``True``

    Returns
    -------
    units_match : bool
        ``True`` if the units of ``ref_da`` and ``dev_da`` match,
        ``False`` otherwise.

    Raises
    ------
    AssertionError
        If ``enforce_units`` is ``True`` and the units do not match.

    Examples
    --------
    >>> match = check_units(ref_da, dev_da, enforce_units=False)
    >>> print(match)
    True
    """
    units_ref = ref_da.units.strip()
    units_dev = dev_da.units.strip()
    if units_ref != units_dev:
        units_match = False
        print("WARNING: ref and dev concentration units do not match!")
        print("Ref units: {}".format(units_ref))
        print("Dev units: {}".format(units_dev))
        if enforce_units:
            assert units_ref == units_dev, \
                "Units do not match: ref {} and dev {}!".format(
                    units_ref, units_dev)
    else:
        units_match = True
    return units_match


def data_unit_is_mol_per_mol(da):
    """
    Checks whether the units of an xarray DataArray are mol/mol.

    Parameters
    ----------
    da : xarray.DataArray
        Data array containing a ``units`` attribute.

    Returns
    -------
    is_molmol : bool
        ``True`` if the units attribute of ``da`` is a recognised
        mol/mol unit string, ``False`` otherwise.

    Notes
    -----
    Recognised mol/mol unit strings are:
    ``'mol mol-1 dry'``, ``'mol/mol'``, ``'mol mol-1'``.

    Examples
    --------
    >>> data_unit_is_mol_per_mol(da)
    True
    """
    conc_units = ["mol mol-1 dry", "mol/mol", "mol mol-1"]
    is_molmol = da.units.strip() in conc_units
    return is_molmol
