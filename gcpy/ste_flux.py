#!/usr/bin/env python
"""
Computes strat-trop exchange (taken as species flux across 100 hPa)
for selected benchmark species.
"""

# ======================================================================
# Imports etc.
# ======================================================================

import os
from calendar import monthrange, month_abbr
import warnings
import numpy as np
import pandas as pd
import xarray as xr
import gcpy.constants as physconsts

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

    def __init__(self, devstr, files, dst, year,
                 bmk_type, species, overwrite, month):
        """
        Initializes the _GlobVars class.

        Args:
            devstr: str
                Label denoting the "Dev" version.
            devdir: str
                Directory where benchmark diagnostic files are found.
            year: int
                Year of the benchmark simulation.

        Keyword Args (optional):
            dst: str
                Directory where plots & tables will be created.
            bmk_type: str
                FullChemBenchmark or TransportTracersBenchmark.
            species: list of str
                Species for which STE fluxes are desired.
            overwrite: bool
                Denotes whether to ovewrite existing budget tables.
        """
        # ------------------------------
        # Arguments from outside
        # ------------------------------
        self.devstr = devstr
        self.files = files
        self.dst = dst
        self.overwrite = overwrite
        self.species = species
        self.month = month
        self.is_TransportTracers = "TransportTracers" in bmk_type

        # ------------------------------
        # Benchmark year
        # ------------------------------
        self.y0 = year
        self.y1 = self.y0 + 1
        self.y0_str = "{}".format(self.y0)
        self.y1_str = "{}".format(self.y1)

        # ------------------------------
        # Read data collections
        # ------------------------------

        # Variable names
        self.data_vars = {}
        for spc in self.species:
            self.data_vars[spc] = "AdvFluxVert_" + spc

        # Vertical flux diagnostics
        # Return a single Dataset containing data from all files.
        # NOTE: Need to add combine="nested" and concat_dim="time"
        # for xarray 0.15 and higher!!!
        v = xr.__version__.split(".")
        if int(v[0]) == 0 and int(v[1]) >= 15:
            try:
                self.ds_flx = xr.open_mfdataset(
                    files,
                    drop_variables=physconsts.skip_these_vars,
                    combine="nested",
                    concat_dim="time"
                )
            except FileNotFoundError:
                msg = "Could not find one or more files in {}".format(files)
                raise FileNotFoundError(msg)
        else:
            try:
                self.ds_flx = xr.open_mfdataset(
                    files,
                    drop_variables=physconsts.skip_these_vars,
                )
            except FileNotFoundError:
                msg = "Could not find one or more files in {}".format(files)
                raise FileNotFoundError(msg)

        # Set a flag to denote if this data is from GCHP
        self.is_gchp = "nf" in self.ds_flx.dims.keys()

        # 100 hPa level (assumes GEOS-FP/MERRA2 72-level grid)
        self.hPa_100_lev = 35

        # Set a flag to denote if this is a 1-year benchmark
        self.is_1yr = self.month is None

        # Months and days
        if self.is_1yr:

            # -----------------------------------
            # Months and days: 1-year benchmark
            # -----------------------------------
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

        else:

            # -----------------------------------
            # Months and days: 1-month benchmark
            # -----------------------------------
            self.N_MONTHS = 1
            self.N_MONTHS_FLOAT = self.N_MONTHS * 1.0

            # Days per month in the benchmark year
            self.d_per_mon = [monthrange(self.y0, self.month)[1] * 1.0]

            # Month name
            self.mon_name = ["{} {}".format(month_abbr[self.month],
                                            self.y0_str)]

            # Days in benchmark year
            self.d_per_yr = 0.0
            for t in range(12):
                self.d_per_yr += monthrange(self.y0, t + 1)[1] * 1.0

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
        globvars: obj of type _GlobVars
            Global variables needed for budget computations.

    Returns:
        result: Pandas DataFrame
            Strat-trop fluxes
    """

    # 100 hPa level index (convert to Python array notation)
    lev = globvars.hPa_100_lev - 1

    # 1-year benchmarks also include the annual mean row
    if globvars.is_1yr:
        n_rows = globvars.N_MONTHS + 1
    else:
        n_rows = globvars.N_MONTHS

    # Create a dictionary to define a DataFrame
    df_dict = {}
    for spc in globvars.species:
        df_dict[spc] = np.zeros(n_rows)

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
            if globvars.is_1yr:
                df.loc["Annual Mean", spc] += np.sum(data)

    # Divide annual mean by days per year
    if globvars.is_1yr:
        df.loc["Annual Mean"] /= globvars.d_per_yr

    return df


def print_ste(globvars, df):
    """
    Prints the strat-trop exchange table.

    Args:
        globvars: _GlobVars
            Global variables
        df: pandas DataFrame
            Strat-trop exchange table
    """
    # Create plot directory hierarchy if necessary
    if os.path.isdir(globvars.dst) and not globvars.overwrite:
        err_str = "Pass overwrite=True to overwrite files in that directory"
        print("Directory {} exists. {}".format(globvars.dst, err_str))
        return
    elif not os.path.isdir(globvars.dst):
        os.makedirs(globvars.dst)

    # Save the file in the Tables folder of dst
    filename = "{}/Strat_trop_exchange.txt".format(globvars.dst)

    # Set numeric format to be 11 chars wide with 4 decimals
    pd.options.display.float_format = '{:11.4f}'.format

    # Open filename for output
    with open(filename, "w+") as f:

        # Print header
        print("%" * 79, file=f)
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
        if globvars.is_1yr:
            print(" Annual mean is weighted by the number of days per month",
                  file=f)
        print("%" * 79, file=f)
        print(file=f)

        # Print the DataFrame
        print(df, file=f)


def make_benchmark_ste_table(devstr, files, year,
                             dst='./1yr_benchmark',
                             bmk_type="FullChemBenchmark",
                             species=["O3"],
                             overwrite=True,
                             month=None):
    """
    Driver program.  Computes and prints strat-trop exchange for
    the selected species and benchmark year.

    Args:
        devstr: str
            Label denoting the "Dev" version.
        files: str
            List of files containing vertical fluxes.
        year: str
            Year of the benchmark simulation.

    Keyword Args (optional):
        dst: str
            Directory where plots & tables will be created.
        bmk_type: str
            FullChemBenchmark or TransportTracersBenchmark.
        species: list of str
            Species for which STE fluxes are desired.
        overwrite: bool
            Denotes whether to ovewrite existing budget tables.
        month: float
            If passed, specifies the month of a 1-month benchmark.
            Default: None (denotes a 1-year benchmark)
    """

    # Initialize a private class with required global variables
    globvars = _GlobVars(devstr, files, dst, int(year),
                         bmk_type, species, overwrite, month)

    # Compute the STE fluxes
    df = compute_ste(globvars)

    # Print the STE fluxes
    print_ste(globvars, df)
