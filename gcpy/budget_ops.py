#!/usr/bin/env python

# ======================================================================
# Imports etc.
# ======================================================================

from calendar import monthrange
from gcpy.benchmark import add_missing_variables
import gcpy.constants as constants
from gcpy.core import compare_varnames
from joblib import Parallel, delayed, cpu_count, parallel_backend
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
        self.pad = 27

        # Atmospheric regimes for which columns are computed
        self.FULL = "Full"
        self.TROP = "Trop"
        self.PBL  = "PBL"
        self.STRAT = "Strat"
        self.regimes = [self.FULL, self.TROP, self.PBL]

        # List of budget operations (rows of the dataframe)
        self.ops = ["Chemistry", "Convection", "EmisDryDep",
                    "Mixing", "Transport", "WetDep"]

        # Columns of the dataframe
        self.REF = "Ref"
        self.DEV = "Dev"
        self.DIFF = "Dev - Ref"
        self.PCTDIFF = "% diff"
        self.columns = [self.REF, self.DEV, self.DIFF, self.PCTDIFF]

        # Title for the accumulation term
        self.ACCUM = "{} ACCUMULATION"

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
        spc : str
            Species name

    Returns:
    --------
        df : pandas DataFrame
    """
    # Operations (e.g. chemistry, wetdep) will be rows
    rows = []
    for op in globvars.ops:
        rows.append("{} {}".format(spc, op).ljust(globvars.pad))
    rows.append(globvars.ACCUM.format(spc).ljust(globvars.pad))
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
        frames[spc] = create_dataframe(globvars, short_name(spc))

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
        row = "{} {}".format(spc_short, op).ljust(globvars.pad)

        # Compute the total change in mass from each operation [kg/s]
        # and convert to [kg/yr] or [Tg/yr] depending on the species
        totref = np.sum(globvars.ref_bdg[v].values) * globvars.conv_fac[spc]
        totdev = np.sum(globvars.dev_bdg[v].values) * globvars.conv_fac[spc]
        diff = totdev - totref
        pctdiff = ( diff / totref ) * 100.0

        # Add total into the proper DataFrame element
        frames[spc].loc[row, globvars.REF] = totref
        frames[spc].loc[row, globvars.DEV] = totdev
        frames[spc].loc[row, globvars.DIFF] = diff
        frames[spc].loc[row, globvars.PCTDIFF] = pctdiff

    # Compute the column sum over each atmospheric regime
    for k in frames.keys():
        row = globvars.ACCUM.format(short_name(k)).ljust(globvars.pad)
        for column in globvars.columns:
            frames[k].loc[row, column] = frames[k][column].sum()

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

    # Filename to contain budget info
    filename = "{}/{}_operations_budgets_{}_{}.txt".format(
        globvars.dst, globvars.devstr, regime, globvars.label)

    # Define a dictionary for the regime descriptions
    desc = { "Full": "Full atm",  "Trop": "Trop only",
             "PBL": "PBL only", "Strat": "Strat only"}

    # Print budgets to the file
    with open(filename, "w+") as f:
        for k in dataframes.keys():
            print("{} budget diagnostics ({}) from {} for {}".format(
                globvars.devstr, desc[regime], k, globvars.label), file=f)
            print("Units : {}".format(globvars.units[k]), file=f)
            print(dataframes[k], file=f)
            print("\n", file=f)
        f.close()


def make_operations_budget_table(refstr, reffiles, devstr, devfiles,
                                 bmk_type, label,
                                 dst=None, interval=None,
                                 species=None, overwrite=True,
                                 pd_float_format="{:13.6f}"):
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

    # Set numeric format to be 16 chars wide with 8 decimals
    pd.options.display.float_format = pd_float_format.format

    # Initialize a class to hold global variables
    globvars = _GlobVars(refstr, reffiles, devstr, devfiles, dst,
                         bmk_type, label, interval, species, overwrite)

    # Nested dictionary (dictonary of dictionary of dataframes)
    bdg = {}

    # ==================================================================
    # Compute operations budgets for Full, Trop, PBL regimes
    # TODO: Parallelize
    # ==================================================================
    for regime in globvars.regimes:

        # Restrict variables to the current atmospheric regime
        # (so we don't waste time looping over variables we won't use here)
        varlist = [v for v in globvars.varlist if regime in v]

        # Return a dictionary of dataframes with all species included
        # and store into the bdg dictonary.
        bdg[regime] = compute_operations_budgets(globvars, varlist, regime)

        # Print the budgets for the current atmospheric regime to a file
        print_operations_budgets(globvars, bdg[regime], regime)

    # ==================================================================
    # Compute the Strat budgets by diffing Full and Trop
    # ==================================================================

    # First create empty dataframes for the Strat regime
    dfs = {}
    for spc in globvars.species:
        dfs[spc] = create_dataframe(globvars, short_name(spc))
    bdg[globvars.STRAT] = dfs

    # Then compute Strat as the difference of Full - Trop
    for spc in globvars.species:
        bdg[globvars.STRAT][spc] = bdg[globvars.FULL][spc] \
                                 - bdg[globvars.TROP][spc]

    # Print the budgets for the Strat regime to a file
    print_operations_budgets(globvars, bdg[globvars.STRAT], globvars.STRAT)
