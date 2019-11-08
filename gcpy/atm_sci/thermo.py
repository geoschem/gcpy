""" Thermodynamics and equation of state calculations """

import numpy as np
from numpy import asarray

from .. constants import G, R_EARTH


def airdens(pressure, temperature=None):
    """ Compute air mass density for a given temperature and pressure.

    If temperature is not provided, a reference temperature is estimated
    using the US Standard Atmosphere.

    Parameters
    ----------
    pressure : float
        Ambient pressure in hPa (mb)
    temperature : float
        Ambient temperature in K

    Returns
    -------
    Air mass density in molecules/cm^3.

    Examples
    --------

    >>> from gcpy.atm_sci.std_atm import airdens

    Print air density at STP/temperature

    >>> airdens(1013.25, 273.15)
    2.6900000e19

    Print air density using reference temperatures from the US Standard
    Atmosphere

    >>> pressures = [(i+1)*100 + 500 for i in range(5)]
    >>> airdens(pressures)
    np.array([  1.67413740e+10,   1.89029148e+10,   2.10103673e+10,
                2.30998231e+10,   2.52037568e+10])

    See Also
    --------
    ussa_alt
    ussa_temp

    """

    pressure = asarray(pressure)

    if temperature is None:
        # TODO: Implement US Std Atm lookup for temperatures based on pressure
        alt = ussa_alt(pressure)
        temperature = ussa_temp(alt)
        #temperature = 273.15
    temperature = asarray(temperature)

    air_density = 2.69e10 * (273.15 / temperature) * (pressure / 1013.25)
    air_density = asarray(air_density)

    # Mask out densities where temperature and pressure were invalid (< 0)
    mask = (temperature <= 0) | (pressure <= 0)
    air_density[mask] = np.nan

    return air_density


def e_h2o(temperature, ice_ref=False, minval=-1e-3):
    """ Calculate water vapor pressure for a given temperature.

    The parameterization used here has been taken from the NASA GTE project
    data description:

    .. math::

        e = 10^{(a - b/T)} T^{-c}

    where :math:`T` is the dew or frostpoint temperature in K and
    :math:`a, b, c` are empirically derived constants.

    Parameters
    ----------
    temperature : float
        Dew or frostpoint reading in K; if you supply the dry air temperature
        (or static air temperature), you will get a value for the water vapor
        saturation pressure.
    ice_ref : logical (default = False)
        Flag indicating that saturation vapor pressure should be calculated
        with respect to ice instead of water.
    minval : float (default = -0.001)
        Minimum valid data value.

    Returns
    -------
    Water vapor pressure in hPa or NaN if `temperature` was less than `min_val`

    Notes
    -----
    N/A

    Examples
    --------

    >>> from gcpy.atm_sci.std_atm import e_h2o

    Calculate water vapor pressure for a dewpoint reading of 266 K

    >>> ph2o = e_h2o(266.)
    >>> ph2o
    3.6174179338886723

    Compute relative humidity by dividing saturation pressure of dry
    temperature

    >>> rh = ph2o / e_h2o(283.)
    0.29460701099796127

    """

    temperature = asarray(temperature)

    if ice_ref:
        # Use constants for frostpoint
        a = 11.4816
        b = 2705.21
        c = 0.32286
    else:
        a = 23.5518
        b = 2937.4
        c = 4.9283

    wvp_h2o = np.power(10., a - (b / temperature)) * np.power(temperature, -c)
    wvp_h2o[temperature < minval] = np.nan

    return wvp_h2o


def hystat(alt, temp, psurf=1013.25):
    """ Compute atmospheric pressures in hydrostatic equilibrium. This
    function is adapted from the Harvard photochemical point model (CHEM1D)

    Parameters
    ----------
    alt : array of floats
        Altitude in km
    temp : array of floats
        Ambient temperature at altitude, in K, corresponding to the given
        altitudes
    psurf : float (default=1013.25)
        Surface presssure, mb

    Returns
    -------
    Pressures corresponding to hydrostatic equilibrium, in mb

    """

    pass


def model_to_airnumdens():
    # TODO
    pass


def numden():
    # TODO
    pass


def rh(dewpoint, temperature, ice=False):
    """ Calculate the relative humidity from temperature and dew/frost point.

    Parameters
    ----------
    dewpoint : float
        Dewpoint (or frostpoint) temperature in K
    temperature : float
        Dry (or static) air temperature in K
    ice : logical
        If true, calculate relative humidity with respect to liquid water
        (dewpoint): else, with respect to ice (frostpoint)

    Returns
    -------
    Relative humidity in %

    Notes
    -----
    N/A

    Example
    -------

    >>> from gcpy.atm_sci.thermo import rh

    Print the RH for dewpoint=366K and temp=278K

    >>> rh(266., 278.)
    41.4736

    """

    if (dewpoint > temperature):
        raise ValueError("ERROR in rh: dew-/frostpoint greater than static "
                         "air temperature")

    # Calculate water vapor pressure and saturation vapor pressure
    e = e_h2o(dewpoint, ice)
    esat = e_h2o(temperature, ice)

    return (e/esat) * 100.0


def ussa_alt(pressure):
    """ Compute altitude in km for a given pressure corresponding to the US
    Standard Atmosphere

    Parameters
    ----------
    pressure : float
        Ambient pressure in hPa (mb); must correspond to an altitude of less
        than 100 km.

    Returns
    -------
    Altitude in kilometers

    Notes
    -----

    Computes approximate altitudes (logp fit to US Standard Atmosphere); tested
    vs interpolated values, 5th degree polynomial gives good results (ca. 1%
    for 0-100 km, ca. 0.5% below 30 km)

    Examples
    --------

    >>> from gcpy.atm_sci.std_atm import ussa_alt
    >>> ussa_alt([1000., 800., 600., 400., 200.])
    np.array([0.106510, 1.95628, 4.20607, 7.16799, 11.8405])

    See Also
    --------
    ussa_press
    ussa_temp

    """

    pressure = asarray(pressure, dtype=float)

    # Mask pressures where P < 0.0003 mb - correspond to about 100 km)
    min_pressure = 3e-4

    # Construct a 5th-degree polynomial in descending order with given
    # coefficients
    coeffs = [ 48.0926, -17.5703, 0.278656, 0.485718, -0.0493698, -0.0283890 ]
    P = np.poly1d(coeffs[::-1], r=False)
    logp = np.log10(pressure)

    alts = P(logp)

    # Apply
    alts[pressure < min_pressure] = np.nan

    return alts


def ussa_press():
    # TODO
    pass


def ussa_temp(altitude):
    """ Compute reference temperature for a given altitude corresponding to the
    US Standard Atmosphere.

    Parameters
    ----------
    altitude : float
        Ambient altitude in kilometers; must lie in the range 0-50 km.

    Returns
    -------
    Temperature at given altitude in K.

    Notes
    -----

    This function evaluates a 6th-order polynomial which had been fitted to
    USSA data from 0-50 km. Accuracy is on the order of 2 K below 8 km, and 5 K
    from 8-50 km. Note that this is less than the actual variation in
    atmospheric temperatures.

    This routine was designed to assign a tempearture value to CTM grid boxes
    in order to allow conversion from mixing ratios to concentrations and vice
    versa.

    Examples
    --------

    >>> from gcpy.atm_sci.std_atm import ussa_temp
    >>> ussa_temp([0, 10, 20, 30])
    np.array([288.283, 226.094, 216.860, 229.344])

    See Also
    --------
    ussa_press
    ussa_alt

    """

    altitude = asarray(altitude, dtype=float)

    # Mask altitudes above 50 km
    max_altitude = 50.0

    # Construct a 6th-degree polynomial in descending order with given
    # coefficients
    coeffs = [  2.88283e2,   -5.20534, -6.75992e-1, 8.75339e-2,
              -3.62036e-3, 6.57752e-5, -4.43960e-7  ]
    P = np.poly1d(coeffs[::-1], r=False)

    pressure = P(altitude)
    pressure[altitude > max_altitude] = np.nan

    return pressure
