#!/usr/bin/env python
"""
Global and Regional Single Panel Plots
--------------------------------------
This example script demonstrates the core single panel plotting
capabilities of GCPy, including global and regional single level plots
as well as global zonal mean plots.

The example data described here is in lat/lon format, but the same code
works equally well for cubed-sphere (GCHP) data.

For full documentation on the plotting capabilities of GCPy
(including full argument lists), please see the GCPy documentation
at https://gcpy.readthedocs.io.

NOTE: If you are using GCPy from a Mac, set the environment variable:

   export MPLBACKEND="MacOSX"

Otherwise set:

   export MPLBACKEND="tkagg"

This will set the proper X11 backend (which is needed to open a plot
window on the screen.  There is some incompatibility with the Tck/Tk
backend "tkagg" in MacOS X operating systems.
"""
import argparse
import xarray as xr
import matplotlib.pyplot as plt
from gcpy.plot.single_panel import single_panel
from gcpy.util import rename_and_flip_gchp_rst_vars


def plot_single_panel(infile, varname, level):
    """
    Example routine to create single panel plots.

    Args:
    -----
    infile  (str) : Name of netCDF file to read.
    varname (str) : Name of variable to plot
    level   (int) : Model level for single-panel plots
                    in Python notation (starting from 0)
    """

    # xarray allows us to read in any NetCDF file
    dset = xr.open_dataset(infile)

    # If the data is from a GCHP restart file, rename variables and
    # flip levels to match the GEOS-Chem Classic naming and level
    # conventions.  Otherwise no changes will be made.
    dset = rename_and_flip_gchp_rst_vars(dset)

    # You can easily view the variables available for plotting
    # using xarray.  Each of these variables has its own xarray
    # DataArray within the larger Dataset container.
    print(dset.data_vars)

    # Most variables have some sort of prefix; in this example all
    # variables are prefixed with 'SpeciesRst_'. We'll select the
    # DataArray for ozone.
    darr = dset[varname]

    # Printing a DataArray gives a summary of the dimensions and attributes
    # of the data.
    print(darr)

    # ==================
    # Single-level Plots
    # ==================

    # gcpy.single_panel is the core plotting function of GCPy, able to
    # create a one panel zonal mean or single level plot.  Here we will
    # create a single level plot. We must manually index into the level
    # (in Python notation, starting from 0).
    darr_single_level = darr.isel(lev=level)

    # single_panel has many arguments which can be optionally specified.
    # The only argument you must always pass to a call to single_panel is
    # the DataArray that you want to plot. By default, the created plot
    # includes a colorbar with units read from the DataArray, an
    # automatic title (the data variable name in the DataArray), and
    # an extent equivalent to the full lat/lon extent of the DataArray
    single_panel(
        darr_single_level,
        title=f"{varname} at level {level}"
    )
    plt.show()

    # You can specify a specific area of the globe you would like plotted
    # using the 'extent' argument, which uses the format [min_longitude,
    # max_longitude, min_latitude, max_latitude] with bounds
    # [-180, 180, -90, 90]
    single_panel(
        darr_single_level,
        extent=[50, -90, -10, 60],
        title=f"{varname} at level {level} over N. Pacific"
    )
    plt.show()

    # Other commonly used arguments include specifying a title and a
    # colormap (defaulting to a White-Green-Yellow-Red colormap)
    #You can find more colormaps at
    # https://matplotlib.org/tutorials/colors/colormaps.html
    single_panel(
        darr_single_level,
        title=f"{varname} at level {level} over N. Pacific, viridis colormap",
        comap=plt.get_cmap("viridis"),
        log_color_scale=True,
        extent=[80, -90, -10, 60]
    )
    plt.show()

    # ===================
    # Zonal Mean Plotting
    # ===================

    # Use the plot_type argument to specify zonal_mean plotting
    single_panel(
        darr,
        plot_type="zonal_mean",
        title=f"Zonal mean plot for {varname}, full atmosphere"
    )
    plt.show()

    # You can specify pressure ranges in hPa for zonal mean plot
    # (by default every vertical level is plotted)
    single_panel(
        darr,
        pres_range=[0, 100],
        log_yaxis=True,
        log_color_scale=True,
        plot_type="zonal_mean",
        title=f"Zonal mean plot for {varname}, stratosphere-only"
    )
    plt.show()


def main():
    """
    Parses command-line arguments and calls plot_single_panel.
    """

    # Tell the parser which arguments to look for
    parser = argparse.ArgumentParser(
        description="Single-panel plotting example program"
    )
    parser.add_argument(
        "-i", "--infile",
        metavar="INF",
        type=str,
        required=True,
        help="input NetCDF file"
    )
    parser.add_argument(
        "-v", "--varname",
        metavar="VARNAME",
        type=str,
        required=True,
        help="variable to plot"
    )
    parser.add_argument(
        "-l", "--level",
        metavar="LEV",
        type=int,
        required=True,
        help="level to plot (single-panel plots only), starting at 0"
    )

    # Parse command-line arguments
    args = parser.parse_args()

    # Call the plot_single_panel routine
    plot_single_panel(
        args.infile,
        args.varname,
        args.level
    )


if __name__ == "__main__":
    main()
