#!/usr/bin/env python
"""
Computes strat-trop exchange (taken as species flux across 100 hPa)
for selected benchmark species.
"""

# ======================================================================
# Imports etc.
# ======================================================================

from calendar import monthrange, month_abbr
import datetime
import gcpy.constants as constants
import numpy as np
import os
from os.path import join
import pandas as pd
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
    def __init__(self, devstr, files, plotsdir, year,
                 bmk_type, species=["O3"], overwrite=True):
        """
        Initializes the _GlobVars class.

        Args:
        -----
            devstr : str
                Label denoting the "Dev" version.
            devdir : str
                Directory where benchmark diagnostic files are found.
            plotsdir : str
                Directory where plots & tables will be created.
            year : int
                Year of the benchmark simulation.
            species : list of str
                Species for which STE fluxes are desired.
            overwrite : bool
                Denotes whether to ovewrite existing budget tables.
        """
        # ------------------------------
        # Arguments from outside
        # ------------------------------
        self.devstr = devstr
        self.files = files
        self.plotsdir = plotsdir
        self.overwrite = overwrite
        self.species = species
        self.is_TransportTracers = "TransportTracers" in bmk_type

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

        # Read the species database
        #try:
        #    path = join(devdir, "species_database.yml")
        #    spcdb = yaml_load_file(open(path))
        #except FileNotFoundError:
        #    path = join(os.path.dirname(__file__), "species_database.yml")
        #    spcdb = yaml_load_file(open(path))

        # ------------------------------
        # Read data collections
        # ------------------------------

        # Variable names
        self.data_vars = {}
        for spc in self.species:
            self.data_vars[spc] = "AdvFluxVert_" + spc

        # Vertical flux diagnostics
        self.ds_flx = xr.open_mfdataset(files)

        # Set a flag to denote if this data is from GCHP
        self.is_gchp = "nf" in self.ds_flx.dims.keys()

        # 100 hPa level (assumes GEOS-FP/MERRA2 72-level grid)
        self.hPa_100_lev = 35

        # ------------------------------
        # Months and days
        # ------------------------------
        self.N_MONTHS = 12
        self.N_MONTHS_FLOAT = self.N_MONTHS * 1.0

        # Days per month in the benchmark year
        self.d_per_mon = np.zeros(self.N_MONTHS)
        for t in range(self.N_MONTHS):
            self.d_per_mon[t] = monthrange(self.y0, t + 1)[1] * 1.0

        # Month names
        self.mon_name = []
        for t in range(self.N_MONTHS):
            self.mon_name.append("{} {}".format(
                month_abbr[t + 1], self.y0_str))
        self.mon_name.append("Annual Mean")

        # Days in the benchmark year
        self.d_per_yr = np.sum(self.d_per_mon)

        # ------------------------------
        # Conversion factors
        # ------------------------------

        # kg/s --> Tg/yr
        kg_s_to_Tg_yr = 1.0e-9 * 86400.0 * self.d_per_yr

        # kg/s --> g/yr
        kg_s_to_g_yr = 1000.0 * 86400.0 * self.d_per_yr

        # Define the conversion factor dict
        self.conv_factor = {}
        for spc in species:
            if self.is_TransportTracers:
                self.conv_factor[spc] = kg_s_to_g_yr
            else:
                self.conv_factor[spc] = kg_s_to_Tg_yr


# ======================================================================
# Methods
# ======================================================================

def compute_ste(globvars):
    """
    Computes the strat-trop-exchange, taken as species flux
    across the 100hPa pressure level.

    Args:
    -----
        globvars : obj of type _GlobVars
            Global variables needed for budget computations.

    Returns:
    --------
        result: Pandas DataFrame
            Strat-trop fluxes
    """

    # 100 hPa level index (convert to Python array notation)
    lev = globvars.hPa_100_lev - 1

    # Create a dictionary to define a DataFrame
    df_dict = {}
    for spc in globvars.species:
        df_dict[spc] = np.zeros(globvars.N_MONTHS + 1)

    # Define the DataFrame: species * months
    df = pd.DataFrame(df_dict, index=globvars.mon_name)

    # Compute monthly strat-trop exchange fluxes,
    # taken as the vertical flux across the 100hPa level
    for t in range(globvars.N_MONTHS):
        month = globvars.mon_name[t]
        for spc in globvars.species:
            var = globvars.data_vars[spc]
            data = globvars.ds_flx[var].isel(time=t, lev=lev).values
            df.loc[month, spc] = np.sum(data) * globvars.conv_factor[spc]

            # Compute annual mean, weighted by the number of days per month
            data = df.loc[month, spc] * globvars.d_per_mon[t]
            df.loc["Annual Mean", spc] += np.sum(data)

    # Divide annual mean by days per year
    df.loc["Annual Mean"] /= globvars.d_per_yr
    return df


def print_ste(globvars, df):
    """
    Prints the strat-trop exchange table.

    Args:
    -----
        globvars : _GlobVars
            Global variables
        df : pandas DataFrame
            Strat-trop exchange table
    """

    # Create table_dir if it doesn't already exist (if overwrite=True)
    if os.path.isdir(globvars.plotsdir) and not globvars.overwrite:
        err_str = "Pass overwrite=True to overwrite files in that directory"
        print("Directory {} exists. {}".format(table_dir, err_str))
        return
    elif not os.path.isdir(globvars.plotsdir):
        os.mkdir(globvars.plotsdir)

    # Save the file in the Tables folder of plotsdir
    filename = "{}/{}.strat_trop_exchange_{}.txt".format(
        globvars.plotsdir, globvars.devstr, globvars.y0_str)

    # Set numeric format to be 11 chars wide with 4 decimals
    pd.options.display.float_format = '{:11.4f}'.format

    # Open filename for output
    with open(filename, "w+") as f:

        # Print header
        print("%"*79, file=f)
        if globvars.is_TransportTracers:
            print(" Table 4. Strat-trop exchange in {} for year {}".format(
                globvars.devstr, globvars.y0_str), file=f)
            print("          (i.e. species flux across 100 hPa)", file=f)
            print("\n Units: g/yr", file=f)
        else:
            print(" Strat-trop exchange in {} for year {}".format(
                globvars.devstr, globvars.y0_str), file=f)
            print(" (i.e. species flux across 100 hPa)", file=f)
            print("\n Units: Tg/yr", file=f)
        print(file=f)
        print(" Annual mean is weighted by the number of days per month",
              file=f)
        print("%"*79, file=f)
        print(file=f)

        # Print the DataFrame
        print(df, file=f)


def make_benchmark_ste_table(devstr, files, plotsdir, year,
                             bmk_type="FullChemBenchmark",
                             species=["O3"],
                             overwrite=True):
    """
    Driver program.  Computes and prints strat-trop exchange for
    the selected species and benchmark year.

    Args:
    -----
        devstr : str
            Label denoting the "Dev" version.
        files : str
            List of files containing vertical fluxes.
        plotsdir : str
            Directory where plots & tables will be created.
        year : int
            Year of the benchmark simulation.
        species : list of str
            Species for which STE fluxes are desired.
        overwrite : bool
            Denotes whether to ovewrite existing budget tables.
    """

    # Initialize a private class with required global variables
    globvars = _GlobVars(devstr, files, plotsdir, year,
                         bmk_type, species, overwrite)

    # Compute the STE fluxes
    df = compute_ste(globvars)

    # Print the STE fluxes
    print_ste(globvars, df)
