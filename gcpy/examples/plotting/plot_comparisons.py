#!/usr/bin/env python
"""
Six Panel Comparison Plots
--------------------------------------
This example script demonstrates the comparitive plotting
capabilities of GCPy, including single level plots as well as
global zonal mean plots. These comparison plots are frequently
used to evaluate results from different runs / versions of
GEOS-Chem, but can also be used to compare results from different
points in one run that are stored in separate xarray datasets.

The example data described here is in lat/lon format, but the same
code works equally well for cubed-sphere (GCHP) data.
"""
import argparse
import xarray as xr
from matplotlib import use as mpl_use
import matplotlib.pyplot as plt
from gcpy import plot
from gcpy.constants import skip_these_vars
from gcpy.util import rename_and_flip_gchp_rst_vars

# X11 backend needed for plt.show()
mpl_use("tkagg")


def plot_comparisons(
        ref,
        dev,
        varname,
        level
):
    """
    Example function to create six-panel comparison plots.

    Args:
    -----
    ref     (str) : Path to the "Ref" data file.
    dev     (str) : Path to the "Dev" data file.
    varname (str) : Variable to plot
    level   (int) : Level to plot (for single-level comparisons only).
    """
    # xarray allows us to read in any NetCDF file, the format of
    # GEOS-Chem diagnostics, #as an xarray Dataset
    #
    # The skip_these_vars list avoids trying to read certain
    # GCHP variables that cause data read issues.
    ref_ds = xr.open_dataset(
        ref,
        drop_variables=skip_these_vars
    )
    dev_ds = xr.open_dataset(
        dev,
        drop_variables=skip_these_vars
    )

    # Special handling is needed for GCHP restart files
    ref_ds = rename_and_flip_gchp_rst_vars(ref_ds)
    dev_ds = rename_and_flip_gchp_rst_vars(dev_ds)

    # ==================
    # Single level plots
    # ==================

    # compare_single_level generates sets of six panel plots for
    # data at a specified level in your datasets.  By default, the
    # level at index 0 (likely the surface) is plotted.
    #
    # You likely want to look at the same variables across both of
    # your datasets. If a variable is in one dataset but not the other,
    # the plots will show NaN values for the latter.  You can pass
    # variable names in a list to these comparison plotting functions
    # (otherwise all variables will plot).
    #
    # NOTE: For simplicity, we will just restrict the comparisons
    # to a single variable.  But you can add as many variables as
    # you like to varlist.
    varlist = [varname]

    # compare_single_level has many arguments which can be optionally
    # specified. The first four arguments are required.  They specify
    # your first xarray Dataset, the name of your first dataset,
    # your second xarray Dataset, and the name of your second dataset.
    # Here we will also pass a specific level and the names of the
    # variables you want to plot.
    plot.compare_single_level(
        ref_ds,
        'Ref version',
        dev_ds,
        'Dev version',
        ilev=level,
        varlist=varlist
    )
    plt.show()

    # Using plt.show(), you can view the plots interactively.
    # You can also save out the plots to a PDF.
    plot.compare_single_level(
        ref_ds,
        'Ref version',
        dev_ds,
        'Dev version',
        ilev=level,
        varlist=varlist,
        pdfname='single_level.pdf'
    )

    # ==================
    # Zonal Mean Plots
    # ==================

    # compare_zonal_mean generates sets of six panel plots containing
    # zonal mean data across your dataset.  compare_zonal_mean shares
    # many of the same arguments as compare_single_level.  You can
    # specify pressure ranges in hPa for zonal mean plotting (by
    # default every vertical level is plotted)
    plot.compare_zonal_mean(
        ref_ds,
        'Ref version',
        dev_ds,
        'Dev version',
        pres_range=[0, 100],
        varlist=varlist,
        pdfname='zonal_mean.pdf'
    )


def main():
    """
    Parses command-line arguments and calls plot_comparisons
    """

    # Tell the parser which arguments to look for
    parser = argparse.ArgumentParser(
        description="Single-panel plotting example program"
    )
    parser.add_argument(
        "-r", "--ref",
        metavar="REF",
        type=str,
        required=True,
        help="path to NetCDF file for the Ref model"
    )
    parser.add_argument(
        "-d", "--dev",
        metavar="DEV",
        type=str,
        required=True,
        help="path to NetCDF file for the Dev model"
    )
    parser.add_argument(
        "-v", "--varname",
        metavar="VAR",
        type=str,
        required=True,
        help="Variable name to plot"
    )
    parser.add_argument(
        "-l", "--level",
        metavar="LEV",
        type=int,
        required=True,
        help="level to plot (single-level plots only), starting at 0"
    )

    # Parse command-line arguments
    args = parser.parse_args()

    # Call the plot_single_panel routine
    plot_comparisons(
        args.ref,
        args.dev,
        args.varname,
        args.level
    )


if __name__ == "__main__":
    main()
