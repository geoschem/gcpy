#!/usr/bin/env python

"""
Creates a table of mean OH taken from GEOS-Chem "Classic" 1-year
FullChemBenchmark log files.

NOTE: At present, GCHP cannot write mean OH to the log files, and
so this module can only be used with GEOS-Chem "Classic".
"""

# ======================================================================
# Imports etc.
# ======================================================================

from calendar import monthrange, month_abbr
import os
from os.path import join
import numpy as np
import pandas as pd

# ======================================================================
# Define a class for passing global variables to the methods below
# ======================================================================


class _GlobVars:
    """
    Private class _GlobVars contains global data that needs to be
    shared among the methods in this module.
    """

    def __init__(self, reflogdir, refstr, devlogdir, devstr,
                 year, dst, overwrite):
        """
        Initializes the _GlobVars class.

        Args:
            reflogdir, devlogdir: str
                Directory containing log files from Ref and Dev.
            refstr, devstr: str
                String label for the Ref and Def simulations.
            year: int
                Year of the Ref and Dev benchmark simulations.
            dst: str
                Folder in which the mean OH table will be printed.
            overwrite: bool
                If true, will overwrite the existing OH table in dst.
        """
        # ------------------------------
        # Arguments from outside
        # ------------------------------
        self.reflogdir = reflogdir
        self.refstr = refstr
        self.devlogdir = devlogdir
        self.devstr = devstr
        self.dst = dst
        self.year = year
        self.overwrite = overwrite

        # ------------------------------
        # Benchmark year
        # ------------------------------
        self.y0 = year
        self.y0_str = "{}".format(self.y0)

        # ------------------------------
        # Months and days
        # ------------------------------
        self.N_MONTHS = 12
        self.N_MONTHS_FLOAT = self.N_MONTHS * 1.0

        # Days per month in the benchmark year
        self.d_per_mon = np.zeros(self.N_MONTHS)
        for t in range(self.N_MONTHS):
            self.d_per_mon[t] = monthrange(self.y0, t + 1)[1] * 1.0

        # Days in the benchmark year
        self.d_per_yr = np.sum(self.d_per_mon)

        # Month names
        self.mon_name = []
        for t in range(self.N_MONTHS):
            self.mon_name.append("{} {}".format(
                month_abbr[t + 1], self.y0_str))
        self.mon_name.append("Annual Mean")


def find_mean_oh(filename):
    """
    Searches a GEOS-Chem "Classic" log file for the Mean OH value.

    Args:
        filename: str
            GEOS-Chem "Classic" log file.
    """
    # Read through each log file from the bottom up
    mean_oh = 0.0
    with open(filename) as f:
        lines = f.readlines()
        for line in reversed(lines):
            if "Mean OH =" in line:
                line = line.replace("Mean OH =", "")
                line = line.replace("[1e5 molec/cm3]", "")
                line = line.strip()
                mean_oh = float(line)
                break

        f.close()
        return mean_oh


def compute_mean_oh_from_logs(globvars):
    """
    Computes mean OH from GEOS-Chem FullChemBenchmark log files.

    Args:
        globvars: _GlobVars
            Global variables
    """

    # Some dataframe indices
    ref = "Ref"
    dev = "Dev"
    diff = "Dev - Ref"
    mean = globvars.mon_name[-1]

    # Create a dictionary to define a DataFrame
    df_dict = {}
    for col in [ref, dev, diff]:
        df_dict[col] = np.zeros(globvars.N_MONTHS + 1)

    # Define the DataFrame: species * months
    df = pd.DataFrame(df_dict, index=globvars.mon_name)

    # Save mean OH from each month to the DataFrame
    for t in range(globvars.N_MONTHS):

        # Indices for dataframe
        mon_yr_str = "{}".format(globvars.mon_name[t])

        # GEOS-Chem log file template
        log_name = "log.{}{}01".format(globvars.year, str(t + 1).zfill(2))

        # Mean OH in Ref log for this month
        ref_log = join(globvars.reflogdir, log_name)
        mean_oh = find_mean_oh(ref_log)
        df.loc[mon_yr_str, ref] = mean_oh
        df.loc[mean, ref] += (mean_oh * globvars.d_per_mon[t])

        # Mean OH in Dev log for this month
        dev_log = join(globvars.devlogdir, log_name)
        mean_oh = find_mean_oh(dev_log)
        df.loc[mon_yr_str, dev] = mean_oh
        df.loc[mean, dev] += (mean_oh * globvars.d_per_mon[t])

    # Annual mean is weighted by days per month
    df.loc[mean] /= globvars.d_per_yr

    # Take Dev - Ref
    df[diff] = df[dev] - df[ref]

    return df


def print_mean_oh_from_logs(globvars, df):
    """
    Prints the mean OH table from 1-year FullChemBenchmark log files.

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
    filename = "{}/{}.mean_OH_from_log_files_{}.txt".format(
        globvars.dst, globvars.devstr, globvars.y0_str)

    # Set numeric format to be 11 chars wide with 4 decimals
    pd.options.display.float_format = '{:13.4f}'.format

    # Open filename for output
    with open(filename, "w+") as f:

        # Print header
        print("%" * 79, file=f)
        print(" Mean OH concentration table for year {}".format(
            globvars.year), file=f)
        print(" \n Units: 1e5 molec cm-3", file=f)
        print(" \n Ref = {}, Dev = {}".format(
            globvars.refstr, globvars.devstr), file=f)
        print("\n Annual mean OH is weighted by the number of days per month",
              file=f)
        print("%" * 79, file=f)
        print(file=f)

        # Print the DataFrame
        print(df, file=f)


def make_benchmark_oh_from_logs(reflogdir, refstr, devlogdir, devstr,
                                year, dst='./1yr_benchmark', overwrite=True):
    """
    Creates the table of mean OH concentrations, as obtained from log files.

    Args:
        reflogdir, devlogdir: str
            Directory containing log files from  Ref and Dev simulations.
        refstr, devstr: str
            String label for the Ref and Def simulations.
        year: str
            Year of the Ref and Dev benchmark simulations.

    Keyword Args (optional):
        dst: str
            Folder in which the mean OH table will be printed.
        overwrite: bool
            If true, will overwrite the existing OH table in dst.
    """

    # Convert year to integer
    year = int(year)

    # Define global variables
    globvars = _GlobVars(reflogdir, refstr, devlogdir, devstr,
                         year, dst, overwrite)

    # Compute mean OH
    df = compute_mean_oh_from_logs(globvars)

    # Print mean OH
    print_mean_oh_from_logs(globvars, df)


if __name__ == "__main__":

    # For debug testing
    maindir = "/path/to/benchmark/main/dir"
    refstr = "gcc_ref_version_str"
    devstr = "gcc_dev_version_str"
    refdir = join(maindir, refstr)
    devdir = join(maindir, devstr)
    reflogdir = join(refdir, "logs")
    devlogdir = join(devdir, "logs")
    year = 2016
    dst = join(devdir, "Plots", "Tables")

    make_benchmark_oh_from_logs(reflogdir, refstr, devlogdir, devstr,
                                year, dst=dst, overwrite=True)
