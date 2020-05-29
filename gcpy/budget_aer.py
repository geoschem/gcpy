#!/usr/bin/env python

"""
Computes the aerosol budgets and burdens from 1-year
FullChemBenchmark simulations.
"""

# ======================================================================
# Imports etc.
# ======================================================================

from calendar import monthrange
import numpy as np
import os
from os.path import join
import gcpy.constants as constants
from gcpy.grid import get_troposphere_mask
import warnings
import xarray as xr
from yaml import load as yaml_load_file

# Suppress harmless run-time warnings (mostly about underflow in division)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# ======================================================================
# Define a class for passing global variables to the methods below
# ======================================================================

class _GlobVars:
    """
    Private class _GlobVars contains global data that needs to be
    shared among the methods in this module.
    """
    def __init__(self, devstr, devdir, dst, year, overwrite):
        """
        Initializes the _GlobVars class.

        Args:
        -----
            devstr : str
                Label denoting the "Dev" version.
            devdir : str
                Directory where benchmark diagnostic files are found.
            dst : str
                Directory where plots & tables will be created.
            year : int
                Year of the benchmark simulation.
            overwrite : bool
                Denotes whether to ovewrite existing budget tables.
        """        """
        Initializes the _GlobVars class.
        """
        # ------------------------------
        # Arguments from outside
        # ------------------------------
        self.devstr = devstr
        self.devdir = devdir
        self.dst = dst
        self.overwrite = overwrite
        
        # ------------------------------
        # Benchmark year
        # ------------------------------
        self.y0 = year
        self.y1 = self.y0 + 1
        self.y0_str = "{}".format(self.y0)
        self.y1_str = "{}".format(self.y1)

        # ------------------------------
        # Species info
        # ------------------------------

        # Directory where diagnostic files are found
        datadir = devdir

        # List of species (and subsets for the trop & strat)
        self.species_list = ["BCPI", "OCPI", "SO4", "DST1", "SALA", "SALC" ]

        # Read the species database
        try:
            path = join(datadir, "species_database.yml")
            spcdb = yaml_load_file(open(path))
        except FileNotFoundError:
            path = join(os.path.dirname(__file__), "species_database.yml")
            spcdb = yaml_load_file(open(path))

        # Molecular weights [g mol-1], as taken from the species database
        self.mw = {}
        for v in self.species_list:
            self.mw[v] = spcdb[v]["MW_g"]
        self.mw["Air"] = constants.MW_AIR * 1.0e3

        # Get the list of relevant AOD diagnostics from a YAML file
        path = join(os.path.dirname(__file__), "aod_species.yml")
        aod = yaml_load_file(open(path))
        self.aod_list = [v for v in aod.keys() if "Dust" in v or "Hyg" in v]

        # Descriptive names
        self.spc2name = {"BCPI": "Black Carbon",
                         "DST1": "Dust",
                         "OCPI": "Organic Carbon",
                         "SO4" : "Sulfate",
                         "SALA": "Sea Salt (accum)",
                         "SALC": "Sea Salt (coarse)"}

        # ------------------------------
        # Collection file lists
        # ------------------------------
        Aerosols = join(datadir, 
                        "*.Aerosols.{}*.nc4".format(self.y0_str))
        StateMet = join(datadir, 
                        "*.StateMet.{}*.nc4".format(self.y0_str))
        SpeciesConc = join(datadir, 
                           "*.SpeciesConc.{}*.nc4".format(self.y0_str))
        
        # ------------------------------        
        # Read data collections
        # ------------------------------

        # Diagnostics
        skip_vars = constants.skip_these_vars
        self.ds_aer = xr.open_mfdataset(Aerosols, data_vars=self.aod_list)
        self.ds_cnc = xr.open_mfdataset(SpeciesConc, drop_variables=skip_vars)
        self.ds_met = xr.open_mfdataset(StateMet, drop_variables=skip_vars)

        # Troposphere mask
        self.tropmask = get_troposphere_mask(self.ds_met)

        # Number of vertical levels
        self.N_LEVS  = self.ds_cnc.dims["lev"]

        # Set a flag to denote if this data is from GCHP
        self.is_gchp = "nf" in self.ds_cnc.dims.keys()

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

        # --------------------------------
        # Surface area
        # (kludgey but it works)
        # --------------------------------
        if self.is_gchp:   
            area = self.ds_met["Met_AREAM2"].values
            a = area.shape
            self.area_m2 = np.zeros([a[0], N_LEVS, a[1], a[2], a[3]])
            for t in range(self.N_MONTHS):
                for k in range(self.N_LEVS):
                    self.area_m2[t,k,:,:,:] = area[t,:,:,:]
            self.total_area_m2 = np.sum(self.area_m2[0,0,:,:,:])
        else:
            area = self.ds_met["AREA"].values
            a = area.shape
            self.area_m2 = np.zeros([a[0], self.N_LEVS, a[1], a[2]])
            for t in range(self.N_MONTHS):
                for k in range(self.N_LEVS):
                    self.area_m2[t,k,:,:] = area[t,:,:]
            self.total_area_m2 = np.sum(self.area_m2[0,0,:,:])


        # ------------------------------
        # Conversion factors
        # ------------------------------

        # v/v dry --> Tg
        self.vv_to_Tg = {}
        for spc in self.species_list:
            self.vv_to_Tg[spc] = self.ds_met["Met_AD"].values  \
                               * (self.mw[spc] / self.mw["Air"]) * 1e-9
            
# ======================================================================
# Methods
# ======================================================================

def annual_average(globvars):
    """
    Computes the annual average budget of aerosol species.

    Args:
    -----
        globvars : obj of type _GlobVars
            Global variables needed for budget computations.

    Returns:
    --------
        result : dict
            Contains annual average budgets.
    """

    # Initialize
    q = {}
    q_sum_f = np.zeros(globvars.N_MONTHS)
    q_sum_t = np.zeros(globvars.N_MONTHS)
    q_sum_s = np.zeros(globvars.N_MONTHS)
    result = {}

    # Define the axes we need to sum over to make monthly sums
    if globvars.is_gchp:
        sum_axes = (1,2,3,4)
    else:
        sum_axes = (1,2,3)

    # Loop over species
    for spc in globvars.species_list:
        
        # Whole-atmosphere and trop-only quantities [g]
        # NOTE: DryDep is by nature trop-only
        varname = "SpeciesConc_" + spc
        q[spc + "_f"] = globvars.ds_cnc[varname].values \
                      * globvars.vv_to_Tg[spc]
        q[spc + "_t"] = np.ma.masked_array(q[spc + "_f"], globvars.tropmask)

        # Compute monthly sums, weighted by the number of days per month
        q_sum_f = np.sum(q[spc + "_f"], axis=sum_axes) * globvars.d_per_mon
        q_sum_t = np.sum(q[spc + "_t"], axis=sum_axes) * globvars.d_per_mon
        q_sum_s = q_sum_f - q_sum_t

        # Compute annual averages
        result[spc + "_f"] = np.sum(q_sum_f) / globvars.d_per_yr
        result[spc + "_t"] = np.sum(q_sum_t) / globvars.d_per_yr
        result[spc + "_s"] = np.sum(q_sum_s) / globvars.d_per_yr

    return result


def annual_average_aod(globvars):
    """
    Computes the annual average AODs (weighted by area and month).

    Args:
    -----
        globvars : obj of type _GlobVars
            Global variables needed for budget computations.

    Returns:
    --------
        result : dict
            Contains annual average AODs.
    """

    # Initialize
    q = {}
    q_sum_f = np.zeros(globvars.N_MONTHS)
    q_sum_t = np.zeros(globvars.N_MONTHS)
    q_sum_s = np.zeros(globvars.N_MONTHS)
    result = {}

    # Define axes to sum over, and total surface area
    if globvars.is_gchp:
        sum_axes = (1,2,3,4)
    else:
        sum_axes = (1,2,3)

    # Loop over AOD variables
    for varname in globvars.aod_list:
        
        # Get the corresponding species name 
        if "Dust" in varname:
            spc = "DST1"
        else:
            spc = varname.split("_")[1]
            
        # Whole-atmosphere AOD [1]
        q[spc + "_f"] = globvars.ds_aer[varname].values

        # Tropospheric-only AOD [1]
        q[spc + "_t"] = np.ma.masked_array(q[spc + "_f"], globvars.tropmask)

        # Create monthly sums, weighted by the number of days per month
        q_sum_f = np.sum(q[spc + "_f"] * globvars.area_m2, axis=sum_axes) \
                * globvars.d_per_mon
        q_sum_t = np.sum(q[spc + "_t"] * globvars.area_m2, axis=sum_axes) \
                * globvars.d_per_mon
        q_sum_s = q_sum_f - q_sum_t
        
        # Take annual averages
        result[spc + "_f"] = np.sum(q_sum_f)         \
                           / globvars.total_area_m2  \
                           / globvars.d_per_yr
        result[spc + "_t"] = np.sum(q_sum_t)         \
                           / globvars.total_area_m2  \
                           / globvars.d_per_yr
        result[spc + "_s"] = np.sum(q_sum_s)         \
                           / globvars.total_area_m2  \
                           / globvars.d_per_yr

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
    # Create the plot directory hierarchy if it doesn't already exist
    if os.path.isdir(globvars.dst) and not globvars.overwrite:
        err_str = "Pass overwrite=True to overwrite files in that directory"
        print("Directory {} exists. {}".format(globvars.dst, err_str))
        return
    elif not os.path.isdir(globvars.dst):
        os.makedirs(globvars.dst)
        
    # File name
    filename = "{}/Aerosol_Burdens_{}.txt".format(globvars.dst, globvars.devstr)

    # Open file and print budgets
    with open(filename, "w+") as f:

        # Print top header
        print("%"*79, file=f)
        print(" Annual average global aerosol burdens for {} in {}".format(
            globvars.y0, globvars.devstr), file=f) 
        print(" (weighted by the number of days per month)", file=f)
        print("%"*79, file=f)
        line = "\n                        Strat         Trop    Strat+Trop\n"
        line += "                    ----------   ----------   ----------"
        print(line, file=f)

        # Print burdens
        for spc in globvars.species_list:
            line = "{} burden [Tg] :  {:10.8f}   {:10.8f}   {:10.8f}\n".format(
                spc.ljust(4),
                data[spc + "_s"],
                data[spc + "_t"],
                data[spc + "_f"])
            print(line, file=f)

        # Close file
        f.close()


def print_annual_average_aod(globvars, data):
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
    # Create table_dir if it doesn't already exist (if overwrite=True)
    if os.path.isdir(globvars.dst) and not globvars.overwrite:
        err_str = "Pass overwrite=True to overwrite files in that directory"
        print("Directory {} exists. {}".format(globvars.dst, err_str))
        return
    elif not os.path.isdir(globvars.dst):
        os.makedirs(globvars.dst)
        
    # File name
    filename = "{}/Global_Mean_AOD_{}.txt".format(globvars.dst, globvars.devstr)

    # Open file and print budgets
    with open(filename, "w+") as f:

        # Print header
        print("%"*79, file=f)
        print(" Annual average global AODs for {} in {} (@ 550nm)".format(
            globvars.y0, globvars.devstr), file=f) 
        print(" (weighted by surface area and unumber of days per month",
              file=f)
        print("%"*79, file=f)
        line = "\n                                        Strat         Trop   Strat+Trop\n"
        line += "                                  -----------   ----------   ----------"
        print(line, file=f)

        # Print burdens
        for spc in globvars.species_list:
            line = "{} mean AOD [1] :  {:11.9f}   {:10.8f}   {:10.8f}\n".format(
                globvars.spc2name[spc].ljust(17),
                data[spc + "_s"],
                data[spc + "_t"], 
                data[spc + "_f"])
            print(line, file=f)

        # Close file
        f.close()


def aerosol_budgets_and_burdens(devstr, devdir, year, 
                                dst='./1yr_benchmark', overwrite=True):
    """
    Compute FullChemBenchmark aerosol budgets & burdens

    Args:
    -----
        devdir : str
            Benchmark directory (containing links to data).
        devstr : str
            Denotes the "Dev" benchmark version.
        year : int
            The year of the benchmark simulation (e.g. 2016). 

    Keyword Args (optional):
    ------------------------
        dst : str 
            Directory where budget tables will be created.
        overwrite : bool
            Overwrite burden & budget tables? (default=True)
    """
    # Initialize a private class with required global variables
    globvars = _GlobVars(devstr, devdir, dst, year, overwrite)

    # Aerosol burdens [Tg]
    burdens = annual_average(globvars)
    print_aerosol_burdens(globvars, burdens)
    
    # Annual average AOD's [Tg]
    aods = annual_average_aod(globvars)
    print_annual_average_aod(globvars, aods)


if __name__ == "__main__":

    # Make sure we have enough arguments
    if  len(sys.argv) != 5:
        err_msg = "Usage: budgets_tt.py devstr, maindir, dst, year"
        raise ValueError(err_msg)
    
    # Call the driver program
    aerosol_budgets_and_burdens(sys.argv[1:4])
