"""
Contains definitions of physical and chemical constants,
as well as other needed global variables.
"""

# ======================================================================
# Physical constant definitions (NIST, 2014)
# ======================================================================

#: Avogadro's number [mol-1]
AVOGADRO = 6.022140857e+23

# Boltzmann's constant [J/K]
BOLTZ = 1.38064852e-23

#: Acceleration due to gravity [m s-2]
G = 9.80665

#: "Equal area" radius of the Earth [km] and [m]
# "Gives the correct total surface area when modeled as a sphere
R_EARTH_m = 6371.0072e+3
R_EARTH_km = 6371.0072

# Typical molar mass of air [kg mol-1] and [g mol-1]
MW_AIR_g = 28.9644
MW_AIR_kg = 28.9644e-3

# Molar mass of water [kg mol-1]
MW_H2O_kg = 18.016e-3

# Gas constant in dry air [J/K/kg]
RD = 287.0

# Molar gas constant [J/K/mol]
RSTARG = 8.3144598

# Gas constant for water vapor [J/K/kg]
RV = 461.0

# ======================================================================
# netCDF variables that we should skip reading
# (these are from cubed-sphere data sets)
# ======================================================================
skip_these_vars = [
    "anchor",
    "ncontact",
    "orientation",
    "contacts",
    "cubed_sphere"
]
