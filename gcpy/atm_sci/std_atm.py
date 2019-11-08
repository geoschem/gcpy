1;95;0c1""" Standard atmosphere functions """

import numpy as np
from numpy import asarray

import xarray as xr

from collections import namedtuple

from scipy.interpolate import interp1d

from .. constants import R_GAS_UNIV, AVOGADRO, G

# For estimation of temperature from altitude
__ISA_SETUP_DONE = False
__ISA_INTERPOLATOR_FROM_Z = None
__ISA_INTERPOLATOR_TO_Z = None
__ISA_MGR = None
__T_REF_VEC = None
__Z_REF_VEC = None

def isa_from_alt(z):
    """ Calculates properties of the International Standard Atmosphere.

    Parameters
    ----------
    z:  float
        Array of altitudes, m

    Returns
    -------
    #TODO

    Notes
    -----
    #TODO

    Examples
    --------
    #TODO

    See Also
    --------
    [NONE]

    """
    # "Convert" to a form we're comfortable with
    zarr = asarray(z)

    if not __ISA_SETUP_DONE:
        isa_setup()

    # Interpolate temperature based on z (now in km)
    temperature = __ISA_INTERPOLATOR_FROM_Z(zarr/1000.0)

    # Now get the rest
    # This is very ugly!
    zvec = zarr.reshape(zarr.size)
    tvec = temperature.reshape(temperature.size)
    pvec = tvec.copy()

    # Make local copy in meters
    zref_m = __Z_REF_VEC*1000.0
    lapse_m = __LAPSE_REF_VEC/1000.0

    # Force initial setup
    i_lo = None
    for i_point, z_curr in enumerate(zvec):
        t_curr = tvec[i_point]

        # Are we still increasing in altitude?
        if i_lo is None or z_curr <= zref_m[i_lo]:
            i_lo = 0
            alpha_temp = 0.0
            # Exponential offset
            p_base = 101325.0 * np.exp(-1.0*__ISA_MGR*zref_m[i_lo]/
                                        __T_REF_VEC[i_lo])

        while z_curr > zref_m[i_lo+1]:
            # Update the pressure to the new altitude
            if np.abs(alpha_temp) > 0:
                p_base = p_base*np.power(__T_REF_VEC[i_lo+1]/__T_REF_VEC[i_lo],
                                    __ISA_MGR/(-1.0*alpha_temp))
            else:
                p_base = p_base * np.exp(__ISA_MGR*(zref_m[i_lo]-zref_m[i_lo+1])
                                        /__T_REF_VEC[i_lo])

            i_lo += 1
            alpha_temp = lapse_m[i_lo]

        if np.abs(alpha_temp) > 0:
            p_temp = p_base*np.power(t_curr/__T_REF_VEC[i_lo],
                                        __ISA_MGR/(-1.0*alpha_temp))
        else:
            p_temp = p_base*np.exp(__ISA_MGR*(zref_m[i_lo]-z_curr)
                                            /__T_REF_VEC[i_lo])
        pvec[i_point] = p_temp

    pressure = pvec.reshape(temperature.shape)

    result = dict(temperature=temperature,pressure=pressure)
    """
    result = xr.Dataset({   'temperature': temperature,
                            'pressure': pressure,
                            'z': z})
                            """
    return result

def isa_temperature(z):
    ISA_data = isa_from_alt(z)
    return ISA_data["temperature"]

def isa_pressure(z):
    ISA_data = isa_from_alt(z)
    return ISA_data["pressure"]

def isa_from_pressure(p):
    if not __ISA_SETUP_DONE:
        isa_setup()
    z_m = __ISA_INTERPOLATOR_TO_Z(asarray(p))
    return z_m

def isa_setup():
    # Check if setup has already been set
    global __ISA_SETUP_DONE
    global __ISA_INTERPOLATOR_FROM_Z
    global __ISA_INTERPOLATOR_TO_Z
    global __ISA_MGR
    global __Z_REF_VEC
    global __T_REF_VEC
    global __LAPSE_REF_VEC

    if __ISA_SETUP_DONE:
        return

    # Reference altitudes, km
    __Z_REF_VEC = np.array([-0.1,0,11,20,32,47,51,71,84.8520,1e6])
    # Reference temperature lapse rates, K/km
    __LAPSE_REF_VEC = np.array([0.0,-6.5,0.0,1.0,2.8,0.0,-2.8,-2.0,0.0])
    z_ref_diff = np.diff(__Z_REF_VEC)

    # Change in temperature between each level
    t_ref_diff = z_ref_diff*__LAPSE_REF_VEC

    # Sum up from a reference temperature of 288.15 K (15 C)
    __T_REF_VEC = np.cumsum(np.insert(t_ref_diff,0,288.15))

    # Set up the interpolator to estimate temperature
    __ISA_INTERPOLATOR_FROM_Z = interp1d(__Z_REF_VEC, __T_REF_VEC)

    # For calculating gas composition
    MW_gas = [    28.0134,
                31.9988,
                39.948,
                44.00995,
                20.183,
                4.0026,
                83.80,
                131.30,
                16.04303,
                2.01594]
    MW_gas = np.array(MW_gas)

    gas_frac = [    0.78084,
                0.209476,
                0.00934,
                0.000314,
                0.00001818,
                0.00000524,
                0.00000114,
                0.000000087,
                0.000002,
                0.0000005]
    gas_frac = np.array(gas_frac)

    # Normalize, to be 100% safe
    gas_frac = gas_frac/np.sum(gas_frac)

    global __ISA_MGR
    __ISA_MGR = np.sum(gas_frac * MW_gas)*1.0e-3*G/R_GAS_UNIV

    # Flag setup as complete and return
    __ISA_SETUP_DONE = True

    # Set up the interpolator to estimate temperature
    z_interp = np.arange(0e3,82.0e3,1.0)
    p_interp = isa_pressure(z_interp)
    __ISA_INTERPOLATOR_TO_Z = interp1d(p_interp,z_interp)

    return
