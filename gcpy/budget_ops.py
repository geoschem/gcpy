#!/usr/bin/env python

# ======================================================================
# Imports etc.
# ======================================================================

from calendar import monthrange
from gcpy.util import add_missing_variables, compare_varnames
import gcpy.constants as constants
from joblib import Parallel, delayed, cpu_count, parallel_backend
import numpy as np
import os
from os.path import join
import pandas as pd
import warnings
import xarray as xr
from joblib import Parallel, delayed, cpu_count, parallel_backend

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
    def __init__(self, refstr, reffiles, devstr, devfiles, dst,
                 bmk_type, label, interval, species, overwrite):
        """
        Initializes the _GlobVars class.

        Args:
        -----
            refstr, devsr : str
                Labels denoting the "Ref" and "Dev" versions.
            reffiles, devfiles : list of str
                Lists of data files for the "Ref" and "Dev" versions.
            dst : str
                Directory where plots & tables will be created.
            bmk_type : str
                "TransportTracersBenchmark" or "FullChemBenchmark".
            label : str
                Label for the plot (e.g. "2016" or "Apr2016", etc.)
            interval : float
                Number of seconds in the budget period.
            species : list of str
                Restrict computation of budgets to certain species.
                Default: print budgets for all common species in
                both Ref and Dev.
            overwrite : bool
                Denotes whether to ovewrite existing budget tables.
        """
        # -----------------------------------------
        # Arguments
        # ------------------------------------------
        self.refstr = refstr
        self.reffiles = reffiles
        self.devstr = devstr
        self.devfiles = devfiles
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
        skip_vars = constants.skip_these_vars
        self.ref_bdg = xr.open_mfdataset(reffiles, drop_variables=skip_vars)
        self.dev_bdg = xr.open_mfdataset(devfiles, drop_variables=skip_vars)
        # Add NaNs for missing variables in each dataset
        [self.ref_bdg, self.dev_bdg] = \
            add_missing_variables(self.ref_bdg, self.dev_bdg)
        # Get the list of common names
        vardict = compare_varnames(self.ref_bdg, self.dev_bdg, quiet=True)
        self.varlist = vardict["commonvars2D"]
        self.varlist = [v for v in self.varlist \
                        if "Strat" not in v and "_" in v]

        # Reduce the dataset to size
        self.ref_bdg = self.ref_bdg[self.varlist]
        self.dev_bdg = self.dev_bdg[self.varlist]

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
        self.conv_fac = {}
        self.units ={}
        for spc in self.species:
            if "TransportTracers" in bmk_type:
                if "Tracer" in spc:
                    if self.annual:
                        # Passive species, annual
                        self.conv_fac[spc] = self.interval * 1e-9
                        self.units[spc] = '[Tg/yr]'
                    else:
                        # Passive species, specific month
                        self.conv_fac[spc] = self.interval * 1e-6
                        self.units[spc] = '[Gg]'
                else:
                    if self.annual:
                        # Radionuclide, annual
                        self.conv_fac[spc] = self.interval
                        self.units[spc] = '[kg/yr]'
                    else:
                        # Radionuclide, specific month
                        self.conv_fac[spc] = self.interval
                        self.units[spc] = '[kg]'
            else:
                if self.annual:
                    # FullChem species, annual
                    self.conv_fac[spc] = self.interval * 1e-9
                    self.units[spc] = '[Tg/yr]'
                else:
                    # FullChem species, specific month
                    self.conv_fac[spc] = self.interval * 1e-6
                    self.units[spc] = '[Gg]'

        # ------------------------------------------
        # Other parameters
        # ------------------------------------------

        # Number of characters to pad
        self.pad = 20

        # Atmospheric regimes for which columns are computed
        self.FULL = "Full"
        self.TROP = "Trop"
        self.PBL  = "PBL"
        self.STRAT = "Strat"
        self.regimes = [self.FULL, self.TROP, self.PBL]

        # List of budget operations (rows of the dataframe)
        self.ACCUM = "ACCUMULATION"
        self.ops = ["Chemistry", "Convection", "EmisDryDep",
                     "Mixing", "Transport", "WetDep", self.ACCUM]

        # Columns of the dataframe
        self.REF = "Ref"
        self.DEV = "Dev"
        self.DIFF = "Dev - Ref"
        self.PCTDIFF = "% diff"
        self.columns = [self.REF, self.DEV, self.DIFF, self.PCTDIFF]

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


def create_dataframe(globvars):
    """
    Returns budget information for a given species
    as a DataFrame object.

    Args:
    -----
        globvars : _Globvars
            Global variables
        spc : str
            Species name

    Returns:
    --------
        df : pandas DataFrame
    """
    # Operations (e.g. chemistry, wetdep) will be rows
    rows = []
    for op in globvars.ops:
        rows.append("{}".format(op).ljust(globvars.pad))
    n_rows = len(rows)

    # Columns will be Ref, Dev, Dev-Ref, %diff
    df_dict = {}
    for column in globvars.columns:
        df_dict[column] = np.zeros(n_rows)

    return pd.DataFrame(df_dict, index=rows)


def compute_operations_budgets(globvars, varlist, regime):
    """
    Computes operations budgets for species.

    Args:
    -----
        globvars : _Globvars
            Global variables
        regime : str
            One of "Full", "Trop", or "PBL".

    Returns:
    --------
        frames : dict of DataFrame
            Operations budgets (one DataFrame per species)
    """
    # Create a dictionary of DataFrames (one per species)
    # Passive tracer species names will be shortened
    frames = {}
    for spc in globvars.species:
        frames[spc] = create_dataframe(globvars)
    # Loop over all variables in the data file
    for v in varlist:

        # Extract DataFrame coordinates (row, column) from the var name
        name = v.replace("Budget", "")
        substrs = name.split("_")
        op = substrs[0]
        spc = substrs[1]
        if spc not in globvars.species:
            continue
        spc_short = short_name(spc)

        # Create the row title (same as the variable name, but
        # with the words "Budget" and "Full/Trop/PBL" stripped out
        if regime in op:
            op = op.replace(regime, "")
        else:
            continue
        row = "{}".format(op).ljust(globvars.pad)

        # Compute the total change in mass from each operation [kg/s]
        # and convert to [kg/yr] or [Tg/yr] depending on the species
        totref = globvars.ref_bdg[v].sum().values * globvars.conv_fac[spc]
        totdev = globvars.dev_bdg[v].sum().values * globvars.conv_fac[spc]
        diff = totdev - totref
        pctdiff = ( diff / totref ) * 100.0
        
        # Add total into the proper DataFrame element
        frames[spc].loc[row, globvars.REF] = totref
        frames[spc].loc[row, globvars.DEV] = totdev
        frames[spc].loc[row, globvars.DIFF] = diff
        frames[spc].loc[row, globvars.PCTDIFF] = pctdiff
    # Compute the column sum over each atmospheric regime
    for k in frames.keys():
        row = globvars.ACCUM.ljust(globvars.pad)
        frames[k].loc[row, "Ref"] = frames[k]["Ref"].sum()
        frames[k].loc[row, "Dev"] = frames[k]["Dev"].sum()
        frames[k].loc[row, "Dev - Ref"] = frames[k]["Dev - Ref"].sum()
        frames[k].loc[row, "% diff"] = ( frames[k].loc[row, "Dev - Ref"] 
                                         / frames[k].loc[row, "Ref"] ) * 100.0
    return frames


def print_operations_budgets(globvars, dataframes, regime):
    """
    Args:
    -----
        globvars : _Globvars
            Global variables
        dataframes : dict of DataFrame
            Operations budgets (one DataFrame per species)
        regime : str
            One of "Full", "Trop", or "PBL".
    """

    # Create the plot directory hierarchy if it doesn't already exist
    if os.path.isdir(globvars.dst) and not globvars.overwrite:
        msg = "Directory {} exists. ".format(dst)
        msg += "Pass overwrite=True to overwrite files in that directory"
        raise ValueError(msg)
    elif not os.path.isdir(globvars.dst):
        os.makedirs(globvars.dst)

    # Define a dictionary for the regime descriptions
    desc = "{}Column".format(regime)
        
    # Filename to contain budget info
    filename = "{}/Budgets_After_Operations_{}_{}.txt".format(
        globvars.dst, regime, globvars.label)

    # Print budgets to the file
    with open(filename, "w+") as f:

        # Header
        print("#"*78, file=f)
        print(" {} budget diagnostics for {}".format(
            desc, globvars.label), file=f)
        print("\n", file=f)
        print("NOTES:", file=f)
        print(" - When using the non-local mixing scheme (default),", file=f)
        print("   'Mixing' includes emissions and dry deposition", file=f)
        print("   applied below the PBL. 'EmisDryDep' therefore", file=f)
        print("   only captures fluxes above the PBL.", file=f)
        print(" - When using full mixing, 'Mixing' and 'EmisDryDep'", file=f)
        print("   are fully separated.", file=f)
        print("#"*78, file=f)
        print(file=f)

        # Data
        for k in dataframes.keys():
            print("{} budgets (Ref={}; Dev={})".format(
                k, globvars.refstr, globvars.devstr), file=f)
            print("Units : {}".format(globvars.units[k]), file=f)
            print(dataframes[k], file=f)
            print("\n", file=f)

        f.close()


def make_operations_budget_table(
    refstr,
    reffiles,
    devstr,
    devfiles,
    bmk_type,
    label,
    dst=None,
    interval=None,
    species=None,
    overwrite=True,
    pd_float_format="{:13.6f}"
):
    """
    Prints the "operations budget" (i.e. change in mass after
    each operation) from a GEOS-Chem benchmark simulation.

    Args:
    -----
        refstr, devstr : str
            Labels denoting the "Ref" and "Dev" versions
        reffiles, devfiles : list of str
            Lists of files to read from the "Ref" and "Dev" version.
        bmk_type : str
            "TransportTracersBenchmark" or "FullChemBenchmark".
        label : str
            Contains the date or date range for each dataframe title.
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
        pd_float_format : str
            Floating-point format for Pandas Dataframe.
            Default value: "{:13.6f}"
    """

    # ==================================================================
    # Initialization
    # ==================================================================
    
    # Set numeric format for Pandas Dataframe
    pd.options.display.float_format = pd_float_format.format
    # Initialize a class to hold global variables
    globvars = _GlobVars(refstr, reffiles, devstr, devfiles, dst,
                         bmk_type, label, interval, species, overwrite)
    # Nested dictionary (dictonary of dictionary of dataframes)
    bdg = {}

    # ==================================================================
    # Compute operations budgets for Full, Trop, PBL regimes
    # TODO: Parallelize it
    # ==================================================================
    def call_compute(regime):
        # Restrict variables to the current atmospheric regime
        # (so we don't waste time looping over variables we won't use here)
        varlist = [v for v in globvars.varlist if regime in v]
        # Return a dictionary of dataframes with all species included
        # and store into the bdg dictonary.
        bdg[regime] = compute_operations_budgets(globvars, varlist, regime)
        # Print the budgets for the current atmospheric regime to a file
        print_operations_budgets(globvars, bdg[regime], regime)
        return {regime : bdg[regime]}
    results = Parallel(n_jobs=-1) (delayed(call_compute)(regime) for regime in globvars.regimes)
    bdg = {list(result.keys())[0] : result[list(result.keys())[0]] for result in results}
    # ==================================================================
    # Compute the Strat budgets by diffing Full and Trop
    # ==================================================================
    # First create empty dataframes for the Strat regime

    dfs = {}
    for spc in globvars.species:
        dfs[spc] = create_dataframe(globvars)
    bdg[globvars.STRAT] = dfs

    # Then compute Strat as the difference of Full - Trop
    for spc in globvars.species:
        bdg[globvars.STRAT][spc] = bdg[globvars.FULL][spc] \
                                 - bdg[globvars.TROP][spc]

    # Print the budgets for the Strat regime to a file
    print_operations_budgets(globvars, bdg[globvars.STRAT], globvars.STRAT)
