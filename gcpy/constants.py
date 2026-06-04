"""
Contains definitions of physical and chemical constants,
as well as other needed global variables.

Notes
-----
Physical constant definitions are taken from NIST, 2014.
"""

#: Avogadro's number [mol-1]
AVOGADRO = 6.022140857e+23

#: Boltzmann's constant [J/K]
BOLTZ = 1.38064852e-23

#: Acceleration due to gravity [m s-2]
G = 9.80665

#: "Equal area" radius of the Earth [km] and [m]
#: "Gives the correct total surface area when modeled as a sphere
R_EARTH_m = 6371.0072e+3
R_EARTH_km = 6371.0072

#: Typical molar mass of air [g mol-1]
MW_AIR_g = 28.9644

#: Typical molar mass of air [kg mol-1]
MW_AIR_kg = 28.9644e-3

#: Molar mass of water [kg mol-1]
MW_H2O_kg = 18.016e-3

#: Gas constant in dry air [J/K/kg]
RD = 287.0

#: Molar gas constant [J/K/mol]
RSTARG = 8.3144598

#: Gas constant for water vapor [J/K/kg]
RV = 461.0

#: netCDF variables that we should skip reading
#: (these are from cubed-sphere data sets)
SKIP_THESE_VARS = [
    "anchor",
    "ncontact",
    "orientation",
    "contacts",
    "cubed_sphere"
]

#: Width for emissions & mass tables
TABLE_WIDTH = 105

#: Column width for emissions & mass tables
COL_WIDTH = 20

#: Default encoding for file I/O
ENCODING = "UTF-8"

#: Default global lat/lon extent
GLOBAL_LL_EXTENT = [-180, 180, -90, 90]

#: Default cubed-sphere tretched-grid stretch factor (i.e. no stretching)
NO_STRETCH_SG_STRETCH_FACTOR = 1

#: Default cubed-sphere stretched-grid target longitude (i.e. no stretching)
NO_STRETCH_SG_TARGET_LON = 170

#: Default cubed-sphere stretched-grid target latitude (i.e. no stretching)
NO_STRETCH_SG_TARGET_LAT = -90

#: Default stretched-grid parameters list (stretch factor, target lon, target lat)
NO_STRETCH_SG_PARAMS = [
    NO_STRETCH_SG_STRETCH_FACTOR,
    NO_STRETCH_SG_TARGET_LON,
    NO_STRETCH_SG_TARGET_LAT
]
