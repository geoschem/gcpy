"""
Contains methods for converting the units of data.
Mainly used for model benchmarking purposes.
"""


import numpy as np
import xarray as xr
import gcpy.constants as physconsts


def adjust_units(units):
    """
    Creates a consistent unit string that will be used in the unit
    conversion routines below.

    Args:
        units: str
            Input unit string.

    Returns:    
        adjusted_units: str
            Output unit string, adjusted to a consistent value.

    Remarks:    
        Unit list is incomplete -- currently is geared to units from
        common model diagnostics (e.g. kg/m2/s, kg, and variants).
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
    Converts a data array from kg to one of several types of target units.

    Args:    
        data_kg: numpy ndarray
            Input data array, in units of kg.

        target_units: str
            String containing the name of the units to which the "data_kg"
            argument will be converted.  Examples: 'Tg', 'Tg C', 'Mg',
            'Mg C', 'kg, 'kg C', etc.

        kg_to_kg_C: float
            Conversion factor from kg to kg carbon.

    Returns:     
        data: numpy ndarray
            Ouptut data array, converted to the units specified
            by the 'target_units' argument.

    Remarks:     
        At present, only those unit conversions corresponding to the
        GEOS-Chem benchmarks have been implemented.

        This is an internal routine, which is meant to be called
        directly from convert_units.
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

    # Return converted data
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
    Converts data stored in an xarray DataArray object from its native
    units to a target unit.

    Args:    
        dr: xarray DataArray
            Data to be converted from native units to target units.
        species_name: str
            Name of the species corresponding to the data stored in "dr".
        species_properties: dict
            Dictionary containing species properties (e.g. molecular
            weights and other metadata) for the given species.
        target_units: str
            Units to which the data will be converted.

    Keyword Args (optional):    
        interval: float
            The length of the averaging period in seconds.
            Default value: [2678400.0]
        area_m2: xarray DataArray
            Surface area in square meters
            Default value: None
        delta_p: xarray DataArray
            Delta-pressure between top and bottom edges of grid box (dry air)
            in hPa
            Default value: None
        box_height: xarray DataArray
            Grid box height in meters
            Default value: None

    Returns:    
        dr_new: xarray DataArray
            Data converted to target units.

    Remarks:    
        At present, only certain types of unit conversions have been
        implemented (corresponding to the most commonly used unit
        conversions for model benchmark output).

        When molmol-1 is present as unit, assumes dry air.
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

    # If the species metadata does not containe MolecRatio, use 1.0 instead
    if "MolecRatio" in species_properties.keys():
        moles_C_per_mole_species = species_properties.get("MolecRatio")
    else:
        moles_C_per_mole_species = 1.0

    # ==============================
    # Compute conversion factors
    # ==============================

    # Physical constants
    Avo = physconsts.AVOGADRO       # molecules/mole
    mw_air = physconsts.MW_AIR_g    # g/mole
    g0 = physconsts.G               # m/s2

    # Get a consistent value for the units string
    # (ignoring minor differences in formatting)
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
            "Conversion from {} to {} for {} requires MW_g definition in species_database.yml".format(
                units, target_units, species_name))

    # Conversion factor for kg species to kg C
    kg_to_kgC = (emitted_mw_g * moles_C_per_mole_species) / mw_g

    # Mass of dry air in kg (required when converting from v/v)
    if 'molmol-1' in units:
        air_mass = delta_p * 100.0 / g0 * area_m2

        # Conversion factor for v/v to kg
        # v/v * kg dry air / g/mol dry air * g/mol species = kg species
        if "g" in target_units:
            vv_to_kg = air_mass / mw_air * mw_g

        # Conversion factor for v/v to molec/cm3
        # v/v * kg dry air * mol/g dry air * molec/mol dry air /
        #  (area_m2 * box_height ) * 1m3/10^6cm3 = molec/cm3
        if "molec" in target_units:
            vv_to_MND = air_mass / mw_air * Avo / (area_m2 * box_height) / 1e6

    # ================================================
    # Get number of seconds per time in dataset
    # ================================================

    # Number of seconds is passed via the interval argument
    numsec = interval

    # Special handling is required if multiple times in interval (for
    # broadcast)
    if len(interval) > 1:
        if 'time' in dr.dims:
            # Need to right pad the interval array with new axes up to the
            # time dim of the dataset to enable broadcasting
            numnewdims = len(dr.dims) - (dr.dims.index('time') + 1)
            for _ in range(numnewdims):
                numsec = numsec[:, np.newaxis]
        else:
            # Raise an error if no time in dataset but interval has length > 1
            raise ValueError(
                'Interval passed to convert_units has length greater than one but data array has no time dimension')

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

    #    elif units == 'molec/cm2/s':
    #        # Implement later

    #    elif units == 'atomsC/cm2/s':
    #         implement later

    elif 'molmol-1' in units:

        if "g" in target_units:
            data_kg = dr.values * vv_to_kg
            data = convert_kg_to_target_units(data_kg, target_units, kg_to_kgC)

        elif "molec" in target_units:
            data = dr.values * vv_to_MND

    else:
        raise ValueError(
            "Units ({}) in variable {} are not supported".format(
                units, species_name))

    # ==============================
    # Return result
    # ==============================

    # Create a new DataArray.  This will be exactly the same as the old
    # DataArray, except that the data will have been converted to the
    # target_units, and the units string will have been adjusted accordingly.
    dr_new = xr.DataArray(
        data, name=dr.name, coords=dr.coords, dims=dr.dims, attrs=dr.attrs
    )
    dr_new.attrs["units"] = target_units

    # Return to calling routine
    return dr_new


def check_units(ref_da, dev_da, enforce_units=True):
    """
    Ensures the units of two xarray DataArrays are the same.

    Args:    
        ref_da: xarray DataArray
            First data array containing a units attribute.
        dev_da: xarray DataArray
            Second data array containing a units attribute.

    Keyword Args (optional):    
        enforce_units: bool
            Whether to stop program if ref and dev units do not match
            Default value: True
    Returns:    
        units_match: bool
    """
    units_ref = ref_da.units.strip()
    units_dev = dev_da.units.strip()
    if units_ref != units_dev:
        units_match = False
        print("WARNING: ref and dev concentration units do not match!")
        print("Ref units: {}".format(units_ref))
        print("Dev units: {}".format(units_dev))
        if enforce_units:
            # if enforcing units, stop the program if
            # units do not match
            assert units_ref == units_dev, \
                "Units do not match: ref {} and dev {}!".format(
                    units_ref, units_dev)
    else:
        units_match = True
    return units_match


def data_unit_is_mol_per_mol(da):
    """
    Check if the units of an xarray DataArray are mol/mol based on a set
    list of unit strings mol/mol may be.

    Args:    
        da: xarray DataArray
            Data array containing a units attribute

    Returns:    
        is_molmol: bool
            Whether input units are mol/mol
    """
    conc_units = ["mol mol-1 dry", "mol/mol", "mol mol-1"]
    is_molmol = False
    if da.units.strip() in conc_units:
        is_molmol = True
    return is_molmol
