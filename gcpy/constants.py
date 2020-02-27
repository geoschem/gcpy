""" 
Contains definitions of physical and chemical constants,
as well as other needed global variables.
"""

#: Acceleration due to gravity [m s-2]
G = 9.80665

#: "Equal area" radius of the Earth [km]
# "Gives the correct total surface area when modeled as a sphere
R_EARTH = 6371.0072

#: Avogadro's number [mol-1]
AVOGADRO = 6.022140857e+23

#: Universal gas contant [J mol-1 K-1]
R_GAS_UNIV = 8.314

#: Boltzmann's constant [J K-1]
K_BOLTZMANN = 1.381e-23

# Typical molar mass of air [kg mol-1]
MW_AIR = 28.9644e-3

# Molar mass of water [kg mol-1]
MW_H2O = 18.016e-3

# netCDF variables that we should skip reading
skip_these_vars = ["anchor", 
                   "ncontact", 
                   "orientation", 
                   "contacts", 
                   "cubed_sphere"]

