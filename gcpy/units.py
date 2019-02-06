'''
Contains methods for converting the units of data.
Mainly used for model benchmarking purposes.
'''

import numpy as np
import xarray as xr


def adjust_units(units):
    '''
    Creates a consistent unit string that will be used in the unit
    conversion routines below.

    Args:
        units : str
            Input unit string.

    Returns:
        adjusted_units: str
            Output unit string, adjusted to a consistent value.

    Remarks:
        Unit list is incomplete -- currently is geared to units from
        common model diagnostics (e.g. kg/m2/s, kg, and variants).

    Examples:
        >>> import gcpy
        >>> print(gcpy.adjust_units('kg/m2/s'))
        kg/m2/s
        >>> print(gcpy.adjust_units('kg m-2 s-1'))
        kg/m2/s
        >>> print(gcpy.adjust_units('kg m^-2 s^-1))
        kg/m2/s
    '''

    # Strip all spaces in the unit string
    units_squeezed = units.replace(' ', '')

    if units_squeezed in ['kg/m2/s', 'kgm-2s-1', 'kgm^-2s^-1']:
        unit_desc = 'kg/m2/s'

    elif units_squeezed in ['kgC/m2/s', 'kgCm-2s-1', 'kgCm^-2s^-1',
                            'kgc/m2/s', 'kgcm-2s-1', 'kgcm^-2s^-1']:
        unit_desc = 'kgC/m2/s'

    elif units_squeezed in [ 'molec/cm2/s', 'moleccm-2s-1', 'moleccm^-2s^-1']:
        unit_desc = 'molec/cm2/s'

    else:
        unit_desc = units_squeezed

    return unit_desc


def convert_kg_to_target_units(data_kg, target_units, kg_to_kgC):
    '''
    Converts a data array from kg to one of several types of target units.

    Args:
        data_kg : numpy ndarray
            Input data array, in units of kg.

        target_units : str
            String containing the name of the units to which the "data_kg"
            argument will be converted.  Examples: 'Tg', 'Tg C', 'Mg', 
            'Mg C', 'kg, 'kg C', etc.

        kg_to_kg_C : float
            Conversion factor from kg to kg carbon.

     Returns:
        data : numpy ndarray
            Ouptut data array, converted to the units specified
            by the 'target_units' argument.

     Remarks:
        At present, only those unit conversions corresponding to the
        GEOS-Chem benchmarks have been implemented.

        This is an internal routine, which is meant to be called
        directly from convert_units.
    '''

    # Convert to target unit
    if target_units == 'Tg':
        data = data_kg * 1e-9

    elif target_units == 'Tg C':
        data = data_kg * kg_to_kgC * 1.0e-9

    elif target_units == 'Gg':
        data = data_kg * 1e-6

    elif target_units == 'Gg C':
        data = data_kg * kg_to_kgC * 1.0e-6

    elif target_units == 'Mg':
        data = data_kg * 1e-3

    elif target_units == 'Mg C':
        data = data_kg * kg_to_kgC * 1.0e-3

    elif target_unit == 'kg':
        data = data_kg

    elif target_unit == 'kg C':
        data = data_kg * kg_to_kgC

    elif target_unit == 'g':

        data = data_kg * 1e3
    elif target_unit == 'g C':

        data = data_kg * kg_to_kgC * 1.0e3

    else:
        msg = 'Target units {} are not yet supported! {}'.format(target_units)
        raise ValueError(msg)

    # Return converted data
    return data


def convert_units(dr, species_name, species_properties,
                  target_units, interval, area_m2):
    '''
    Converts data stored in an xarray DataArray object from its native
    units to a target unit.

    Args:
        dr : xarray DataArray
            Data to be converted from native units to target units.

        species_name : str
            Name of the species corresponding to the data stored in "dr".

        species_properties : dict
            Dictionary containing species properties (e.g. molecular
            weights and other metadata) for the given species.

        target_units : str
            Units to which the data will be converted.

        interval : float
            The length of the averaging period in seconds.

        area_m2 : xarray DataArray
            Surface area in square meters

    Returns:
        dr_new : xarray DataArray
            Data converted to target units.

    Remarks:
        At present, only certain types of unit conversions have been
        implemented (corresponding to the most commonly used unit
        conversions for model benchmark output).

    Example:
        >>> import.gcpy
        >>> import xarray as xr
        >>> import json
        >>> ds = xr.open_dataset("myfile.nc")
        >>> dr = ds["CO"]
        >>> properties = json.load(open(species_database.json))
        >>> dr_new = convert_units(dr, "CO", properties.get("CO"),
                     "Tg", interval=86400.0, ds["AREA"])

    '''

    # Get species molecular weight information
    mol_wt_g = species_properties.get('MW_g')
    emitted_mol_wt_g = species_properties.get('EmMW_g')
    moles_C_per_mole_species = species_properties.get('MolecRatio')

    # Ratio for converting kg species to kg C
    kg_to_kgC = (emitted_mol_wt_g * moles_C_per_mole_species) / mol_wt_g

    # Get a consistent value for the units string
    # (ignoring minor differences in formatting)
    units = adjust_units(dr.units)

    # ==============================
    # First compute units to kg
    # ==============================

    if units == 'kg/m2/s':
        data_kg = dr.values * area_m2.values * interval

    elif units == 'kgC/m2/s':
        data_kg = dr.values * area_m2.values * interval / kg_to_kgC

    elif units == 'kg':
        data_kg = dr.values

    elif units == 'kgC':
        data_kg = dr.values / kg_to_kgC

#    elif units == 'molec/cm2/s':
#        # Implement later

#    elif units == 'atomsC/cm2/s':
#         implement later
    else:
        raise ValueError('Units () are not supported'.format(units))

    # ===============================
    # Then compute to target units
    # ===============================

    # Convert
    data = convert_kg_to_target_units(data_kg, target_units, kg_to_kgC)

    # ==============================
    # Return result
    # ==============================

    # Create a new DataArray.  This will be exactly the same as the old
    # DataArray, except that the data will have been converted to the
    # target_units, and the units string will have been adjusted accordingly.
    dr_new = xr.DataArray(data, name=dr.name, coords=dr.coords,
                          dims=dr.dims, attrs=dr.attrs)
    dr_new.attrs['units'] = target_units

    # Return to calling routine
    return dr_new

