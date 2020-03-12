#!/usr/bin/env python

# ======================================================================
# Imports etc.
# ======================================================================

from calendar import monthrange
import numpy as np
import os
from os.path import join
import pandas as pd
import warnings
import xarray as xr

# Tell matplotlib not to look for an X-window
os.environ['QT_QPA_PLATFORM']='offscreen'

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
    def __init__(self, devstr, files, dst, bmk_type,
                 label, interval, species, overwrite):
        """
        Initializes the _GlobVars class.

        Args:
        -----
            devstr : str
                Label denoting the "Dev" version.
            files : list of str
                Files to be read from disk.
            dst : str
                Directory where plots & tables will be created.
            bmk_type : str
                "TransportTracersBenchmark" or "FullChemBenchmark".
            label : str
                Label for the plot (e.g. "2016" or "Apr2016", etc.)
            interval : float
                Number of seconds in the budget period.
            overwrite : bool
                Denotes whether to ovewrite existing budget tables.
        """
        # -----------------------------------------
        # Arguments
        # ------------------------------------------
        self.devstr = devstr
        self.files = files
        self.dst = dst
        self.bmk_type = bmk_type
        self.label = label.strip()
        self.interval = interval
        self.species = species
        self.overwrite = overwrite

        # -----------------------------------------
        # Other variables
        # ------------------------------------------

        # Is this an annual total?
        self.annual = self.interval > 3.0e7
                    
        # ------------------------------------------
        # Read data from disk
        # ------------------------------------------

        # Read the Budget collection
        self.ds_bdg = xr.open_mfdataset(files)

        # Variable list
        self.varlist = [v for v in self.ds_bdg.data_vars.keys() \
                        if "Strat" not in v and "PBL" not in v  \
                        and "_" in v]
        
        # Reduce the dataset to size
        self.ds_bdg = self.ds_bdg[self.varlist]

        # ------------------------------------------
        # Species info
        # ------------------------------------------

        # Default species if none are passed
        if self.species is None:
            self.species = [v.split("_")[1] for v in self.varlist \
                            if "Full" in v]
            
        # Also exclude the stratospheric species
        self.species = [s for s in self.species if "Strat" not in s]

        # Sort the species alphabetically
        self.species.sort()

        # -----------------------------------------
        # Conversion factors and units
        # -----------------------------------------
        self.conv_factor = {}
        self.units ={}
        for spc in self.species:
            if "TransportTracers" in bmk_type:
                if "Tracer" in spc:
                    if self.annual:
                        # Passive species, annual
                        self.conv_factor[spc] = self.interval * 1e-9  
                        self.units[spc] = '[Tg/yr]'
                    else:
                        # Passive species, specific month
                        self.conv_factor[spc] = self.interval * 1e-6
                        self.units[spc] = '[Gg]'
                else:
                    if self.annual:
                        # Radionuclide, annual
                        self.conv_factor[spc] = self.interval
                        self.units[spc] = '[kg/yr]'
                    else:
                        # Radionuclide, specific month
                        self.conv_factor[spc] = self.interval
                        self.units[spc] = '[kg]'
            else:
                if self.annual:
                    # FullChem species, annual
                    self.conv_factor[spc] = self.interval * 1e-9
                    self.units[spc] = '[Tg/yr]'
                else:
                    # FullChem species, specific month
                    self.conv_factor[spc] = self.interval * 1e-6
                    self.units[spc] = '[Gg]'
        
        # ------------------------------------------
        # Other parameters
        # ------------------------------------------

        # Number of characters to pad
        self.pad = 27

        # List of operations
        self.ops = ["Chemistry", "Convection", "EmisDryDep",
                    "Mixing", "Transport", "WetDep"]

        # Atmospheric regimes
        self.regimes = ["Strat", "Trop", "Strat+Trop"]

# ======================================================================
# Methods
# ======================================================================

def short_name(spc):
    """
    Shortens a species name for display in the budget table.

    Args:
    -----
        spc : str
            Species name
    
    Returns:
    --------
        spc_short : str
            Abbreviated species name
    """
    spc_short = spc.replace("Uniform", "Un")
    spc_short = spc_short.replace("Tracer", "Trac")
    spc_short = spc_short.replace("Anthro", "Anth")
    spc_short = spc_short.replace("Emis", "Em")
    spc_short = spc_short.replace("day", "d")
    return spc_short

def create_dataframe(globvars, spc):
    """
    Returns budget information for a given species
    as a DataFrame object.

    Args:
    -----
        globvars : _Globvars 
            Global variables

    Returns:
    --------
        df : pandas DataFrame
    """
    # Operations will be the DataFrame index
    index = []
    for op in globvars.ops:
        index.append("{} {}".format(spc, op).ljust(globvars.pad))
    index.append("{} TOTAL".format(spc).ljust(globvars.pad))
    n_index = len(index)
                     
    # Atmospheric regimes will be DataFrame series
    df_dict = {}
    for r in globvars.regimes:
        df_dict[r] = np.zeros(n_index)

    return pd.DataFrame(df_dict, index=index)


def compute_operations_budgets(globvars):
    """
    Computes operations budgets for species.

    Args:
    -----
        globvars : _Globvars 
            Global variables

    Returns:
    --------
        frames : dict of DataFrame
            Operations budgets (one DataFrame per species)
    """
    # Create a dictionary of DataFrames (one per species)
    # Passive tracer species names will be shortened
    frames = {}
    for spc in globvars.species:
        frames[spc] = create_dataframe(globvars, short_name(spc))
        
    # Loop over all variables in the data file
    for v in globvars.varlist:

        # Extract DataFrame coordinates (index, regime) from the var name
        name = v.replace("Budget", "")
        substrs = name.split("_")
        op = substrs[0]
        spc = substrs[1]
        if spc not in globvars.species:
            continue
        spc_short = short_name(spc)
        if "Trop" in op:
            op = op.replace("Trop", "")
            regime = "Trop"
        elif "Full" in op:
            op = op.replace("Full", "")
            regime = "Strat+Trop"
        index = "{} {}".format(spc_short, op).ljust(globvars.pad)

        # Compute the total change in mass from each operation [kg/s]
        # and convert to [kg/yr] or [Tg/yr] depending on the species
        total = np.sum(globvars.ds_bdg[v].values) * globvars.conv_factor[spc]

        # Add total into the proper DataFrame element
        frames[spc].loc[index, regime] = total

    # Compute stratospheric values
    for k in frames.keys():
        frames[k]["Strat"] = frames[k]["Strat+Trop"] - frames[k]["Trop"]

    # Compute the column sum over each atmospheric regime
    for k in frames.keys():
        index = "{} TOTAL".format(short_name(k)).ljust(globvars.pad)
        for regime in globvars.regimes:
            frames[k].loc[index,regime] = frames[k][regime].sum()

    return frames


def print_operations_budgets(globvars, dataframes):
    """
    Args:
    -----
        globvars : _Globvars 
            Global variables
        dataframes : dict of DataFrame
            Operations budgets (one DataFrame per species)
    """
    # Set numeric format to be 16 chars wide with 8 decimals
    pd.options.display.float_format = '{:16.8f}'.format

    # Create the plot directory hierarchy if it doesn't already exist
    if os.path.isdir(globvars.dst) and not globvars.overwrite:
        err_str = "Pass overwrite=True to overwrite files in that directory"
        print("Directory {} exists. {}".format(dst, err_str))
        return
    elif not os.path.isdir(globvars.dst):
        os.makedirs(globvars.dst)

    # Filename to contain budget info
    filename = "{}/{}_operations_budgets_{}.txt".format(
        globvars.dst, globvars.devstr, globvars.label)
    
    # Print budgets to the file
    with open(filename, "w+") as f:    
        for k in dataframes.keys():

            # Header
            print("{} budget diagnostics from {} for {}".format(
                globvars.devstr, k, globvars.label), file=f)
            print("Units : {}".format(globvars.units[k]), file=f)

            # Data
            print(dataframes[k], file=f)
            index = "{} TOTAL".format(short_name(k)).ljust(globvars.pad)
            total = 0.0
            for regime in globvars.regimes:
                total += dataframes[k].loc[index, regime]
            index = index.replace("TOTAL", "GLOBAL")
            print("-"*79, file=f)
            print("{}{:16.8f}".format(index, total), file=f)
            print("\n", file=f)

        f.close()


def make_operations_budget_table(devstr, files, bmk_type, label,
                                 dst=None, interval=None,
                                 species=None, overwrite=True):
    """
    Prints the "operations budget" (i.e. change in mass after
    each operation) from a GEOS-Chem benchmark simulation.

    Args:
    -----
        devstr : str
            Label denoting the "Dev" version.
        files : list of str
            List of files to read
        bmk_type : str
            "TransportTracersBenchmark" or "FullChemBenchmark".
        year : int
            Year of the benchmark simulation.

    Keyword Args (optional):
    ------------------------
        dst : str
            Directory where plots & tables will be created.
        interval : float
            Number of seconds in the diagnostic interval.
        species : list of str
            List of species for which budgets will be created.
        overwrite : bool
            Denotes whether to ovewrite existing budget tables.
    """
    # Initialize a class to hold global variables
    globvars = _GlobVars(devstr, files, dst, bmk_type,
                         label, interval, species, overwrite)

    # Compute operations budgets.  Return a dictionary of
    # DataFrame objects (one DataFrame per species).
    dataframes = compute_operations_budgets(globvars)

    # Print the operations budgets
    print_operations_budgets(globvars, dataframes)
