#!/usr/bin/env python

"""
Computes the budget of Pb and Be7 from the TransportTracers benchmarks.

NOTE: This works for GC-Classic, but may need modifications for GCHP.
 -- Bob Yantosca (21 Jan 2020)
"""

# ======================================================================
# Imports etc.
# ======================================================================

from calendar import monthrange
import numpy as np
import os
from os.path import join
import gcpy.constants as constants
from gcpy.benchmark import get_troposphere_mask
import warnings
import xarray as xr

# Suppress harmless run-time warnings (mostly about underflow in division)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# ======================================================================
# GLOBAL VARIABLES: Configurables (MUST EDIT)
# ======================================================================

# Define a class for passing local data
class _GlobVars:
    """
    Private class _GlobVars contains global data that needs to be
    shared among the methods in this module.
    """
    def __init__(self, devstr, devdir, plotsdir, year, overwrite):
        """
        Initializes the _GlobVars class.
        """
        # ------------------------------
        # Arguments from outside
        # ------------------------------
        self.devstr = devstr
        self.devdir = devdir
        self.plotsdir = plotsdir
        self.overwrite = overwrite
        
        # ------------------------------
        # Benchmark year
        # ------------------------------
        self.y0 = year
        self.y1 = self.y0 + 1
        self.y0_str = "{}".format(self.y0)
        self.y1_str = "{}".format(self.y1)

        # ------------------------------
        # Collection file lists
        # ------------------------------
        datadir = join(devstr, "OutputDir")
        Aerosols = join(datadir, 
                        "GEOSChem.Aerosols.{}*.nc4".format(self.y0_str))
        StateMet = join(datadir, 
                        "GEOSChem.StateMet.{}*.nc4".format(self.y0_str))
        SpeciesConc = join(datadir, 
                           "GEOSChem.SpeciesConc.{}*.nc4".format(self.y0_str))

        # ------------------------------        
        # Read data collections
        # ------------------------------

        # Diagnostics
        self.ds_aer = xr.open_mfdataset(Aerosols)
        self.ds_cnc = xr.open_mfdataset(SpeciesConc)
        self.ds_met = xr.open_mfdataset(StateMet)

        # Area and troposphere mask
        self.area_m2 = self.ds_met["AREA"].isel(time=0)
        self.area_cm2 = self.area_m2 * 1.0e4
        self.tropmask = get_troposphere_mask(self.ds_met)

        # ------------------------------       
        # Months and days
        # ------------------------------
        self.N_MONTHS = 12
        self.N_MONTHS_FLOAT = self.N_MONTHS * 1.0
        
        # Days per month in the benchmark year
        self.d_per_mon = np.zeros(self.N_MONTHS)
        for t in range(self.N_MONTHS):
            self.d_per_mon[t] = monthrange(self.y0, t+1)[1] * 1.0
            
        # Days in the benchmark year
        self.d_per_yr = np.sum(self.d_per_mon)
            
        # Fraction of year occupied by each month
        self.frac_of_yr = np.zeros(self.N_MONTHS)
        for t in range(self.N_MONTHS):
            self.frac_of_yr[t] = self.d_per_mon[t] / self.d_per_yr

        # ------------------------------
        # Species info
        # ------------------------------

        # List of species (and subsets for the trop & strat)
        self.species_list = ["BCPI", "SO4", "DST1", "SALA", "SALC" ]

        # Molecular weights
        # NOTE: In future, get these from species database (bmy, 2/25/20)
        self.mw = { "BCPI": 12.0, 
                    "SO4": 96.0, 
                    "DST1": 29.0, 
                    "SALA": 31.4,
                    "SALC": 31.4,
                    "Air" : 28.9644}

        # ------------------------------
        # Conversion factors
        # ------------------------------
        self.kg_per_mol = {}
        self.vv_to_Tg = {}
        self.kg_s_to_g_d = {}

        for spc in self.species_list:

            # kg/mole for each species
            self.kg_per_mol[spc] = constants.AVOGADRO / (self.mw[spc]  * 1e-3)

            # v/v dry --> Tg
            self.vv_to_Tg[spc] = self.ds_met["Met_AD"].values  \
                               * (self.mw[spc] / self.mw["Air"]) * 1e-9

# ======================================================================
# Functions
# ======================================================================

def diff(globvars, dict0, dict1):
    """
    Function to take the difference of two dict objects.
    Assumes that both objects have the same keys.

    Args:
    -----
        globvars : obj of type _GlobVars
            Global variables needed for budget computations.
        dict0, dict1 : dict
            Dictionaries to be subtracted (dict1 - dict0)
    """
    result = {}
    for key, value in dict0.items():
        result[key] = dict1[key] - dict0[key]

    return result


def annual_average(globvars, ds, collection, conv_factor):
    """
    Computes the annual average of budgets or fluxes.

    Args:
    -----
        globvars : obj of type _GlobVars
            Global variables needed for budget computations.
        ds : xarray Dataset
            Data to be averaged
        collection : str
            Name of the diagnostic collection.
        conv_factor : str
            Conversion factor to be applied.
    """

    # Initialize
    q = {}
    q_sum_f = np.zeros(globvars.N_MONTHS)
    q_sum_t = np.zeros(globvars.N_MONTHS)
    q_sum_s = np.zeros(globvars.N_MONTHS)
    result = {}
    
    for spc in globvars.species_list:
        
        # Whole-atmosphere and trop-only quantities [g]
        # NOTE: DryDep is by nature trop-only
        varname = collection.strip() + "_" + spc
        q[spc + "_f"] = ds[varname].values * conv_factor[spc]
        q[spc + "_t"] = np.ma.masked_array(q[spc + "_f"], globvars.tropmask)

        # Compute monthly averages, weighted by # of days in month
        # Special handling for Drydep, which is a 3-D array
        for t in range(globvars.N_MONTHS):
            q_sum_f[t] = np.sum(q[spc + "_f"][t,:,:,:]) * globvars.d_per_mon[t]
            q_sum_t[t] = np.sum(q[spc + "_t"][t,:,:,:]) * globvars.d_per_mon[t] 
            q_sum_s[t] = q_sum_f[t] - q_sum_t[t]
            
        # Take annual averages
        result[spc + "_f"] = np.sum(q_sum_f) / globvars.d_per_yr
        result[spc + "_t"] = np.sum(q_sum_t) / globvars.d_per_yr
        result[spc + "_s"] = np.sum(q_sum_s) / globvars.d_per_yr

    return result


def print_aerosol_burdens(globvars, data):
    """
    Prints burdens for aerosols species from a
    1-year FullChemBenchmark simulation.
    
    Args:
    -----
        globvars: object of type _GlobVars
            Global variables needed for budget computations.
        data: dict
            Nested dictionary containing budget info.
    """

    # Directory in which budgets & burdens tables. will be created
    table_dir = "{}/Aerosols".format(globvars.plotsdir)

    # Create table_dir if it doesn't already exist (if overwrite=True)
    if os.path.isdir(table_dir) and not globvars.overwrite:
        err_str = "Pass overwrite=True to overwrite files in that directory"
        print("Directory {} exists. {}".format(table_dir, err_str))
        return
    elif not os.path.isdir(table_dir):
        os.mkdir(table_dir)
        
    # File name
    filename = "{}/Aerosol_Burdens_{}.txt".format(table_dir, globvars.devstr)

    # Open file and print budgets
    with open(filename, "w+") as f:

        # Print header
        print("%"*79, file=f)
        print(" Annual average global aerosol burdens for {} in {}".format(
            globvars.y0, globvars.devstr), file=f) 
        print("%"*79, file=f)
        print(file=f)

        for spc in globvars.species_list:
            print("{} burden [Tg]".format(spc), file=f)
            print( "    Stratosphere :   {:10.8f}".format(
                data["burden"][spc + "_s"]), file=f)
            print( "    Troposphere  :   {:10.8f}".format(
                data["burden"][spc + "_t"]), file=f)
            print( "    Strat + Trop :   {:10.8f}".format(
                data["burden"][spc + "_f"]), file=f)
            print(file=f)

        f.close()
        
        
def aerosol_budgets_and_burdens(devstr, devdir,
                                plotsdir, year, overwrite=True):
    """
    Compute FullChemBenchmark aerosol budgets & burdens

    Args:
    -----
        devdir : str
            Benchmark directory (containing links to data).
        devstr : str
            Denotes the "Dev" benchmark version.
        plotsdir : str 
            Directory where budget tables will be created.
        year : int
            The year of the benchmark simulation (e.g. 2016). 
        overwrite : bool
            Overwrite burden & budget tables? (default=True)
    """

    # ==================================================================
    # Initialize
    # ==================================================================

    # Dictionary containing output data
    data = {}
    
    # Internal class with Get global variables necessary for computations
    globvars = _GlobVars(devstr, devdir, plotsdir, year, overwrite)

    # ==================================================================
    # Aerosol burdens [Tg]
    # ==================================================================

    # Compute aerosol burdens
    data["burden"] = annual_average(globvars,
                                    ds=globvars.ds_cnc,
                                    collection='SpeciesConc',
                                    conv_factor=globvars.vv_to_Tg)


    # Print aerosol burdens
    print_aerosol_burdens(globvars, data)
    
    # ==================================================================
    # Annual average AOD's [Tg]
    # ==================================================================


if __name__ == "__main__":

    # Make sure we have enough arguments
    if  len(sys.argv) != 5:
        err_msg = "Usage: budgets_tt.py devstr, maindir, plotsdir, year"
        raise ValueError(err_msg)
    
    # Call the driver program
    aerosol_budgets_and_burdens(sys.argv[1:4])
