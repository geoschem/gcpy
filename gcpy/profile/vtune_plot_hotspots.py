#!/usr/bin/env python3
"""
Plots the top functions by exclusive time from an
Intel VTune profile of hotspots listed by function.
"""
from sys import argv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from gcpy.util import verify_variable_type
from gcpy.plot.core import text_to_data_units
from gcpy.profile.vtune_utils import vtune_read_hotspots_csv


def vtune_read_hotspots(filename):
    """
    Reads an Intel VTune hotspot data report and returns a
    pd.DataFrame object with the results.

    Args
    filename : str          : CSV file containing hotspot data

    Returns
    dframe   : pd.DataFrame : DataFrame w/ hotspot data
    """

    # Read the profiling data
    dframe = vtune_read_hotspots_csv(filename)

    # Drop entries for system calls that can be listed more
    # than once in the report, causing confusing output
    dframe = dframe[~dframe["Source File"].isin(
        ["wait.h", "mutex.c", "do_spin", "simple-bar.h"]
    )]

    # Hotspots listed by function
    if "Function" in dframe.columns:
        dframe = dframe.set_index("Function")
        return dframe

    # Hotspots listed by source line
    if "Source File" in dframe.columns and "Source Line" in dframe.columns:

        # Concatenate "Source File" and "Source Line" into "Function",
        # which will become the index
        dframe["Function"] = dframe["Source File"] + \
            "_L" + dframe["Source Line"]
        dframe = dframe.set_index("Function")

        return dframe

    return None


def vtune_plot_functions(dframe, filename, n_min, n_max):
    """
    Plots functions having the largest CPU time as listed in
    an Intel Vtune hotspots report.

    Args
    dframe   : pd.DataFrame : Profiling information
    filename : str          : Name of the file
    n_min    : int          : Display functions starting with this index
    n_max    : int          : Display functions ending with this index
    """
    verify_variable_type(dframe, pd.DataFrame)
    verify_variable_type(filename, str)
    verify_variable_type(n_min, int)
    verify_variable_type(n_max, int)

    # ------------------------------------------------------------------
    # Prepare the data
    # ------------------------------------------------------------------

    # Extract info from the data
    labels = list(dframe.index)
    times  = [var[0] for var in dframe[["CPU Time [s]"]].to_numpy()]
    times  = [np.float64(var) for var in times]

    # Take the top n_funcs values (excluding total) and reverse
    labels = labels[n_min:n_max][::-1]
    times = times[n_min:n_max][::-1]

    # Invert values for bars to extend to the left
    neg_values = [-var for var in times]

    # ------------------------------------------------------------------
    # Create horizontal bar plot
    # ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(11,8))
    hbars = ax.barh(labels, neg_values, color="skyblue")
    fig.canvas.draw()

    # ------------------------------------------------------------------
    # Create time and percent labels to be displayed
    # at the right of each histogram bar
    # ------------------------------------------------------------------

    # Initialize
    time_len = 0.0

    # Starting position of the time label at right of plot
    xmin = min(neg_values) * 1.02
    time_x = round(abs(xmin) * 0.02)

    # Create labels for CPU times [s] and also compute
    # the length of the longest label in data units
    for hbar, time in zip(hbars, times):
        text = ax.text(
            time_x,
            hbar.get_y() + hbar.get_height()/2,
            f"{time}",
            va="center",
            ha="left",
        )
        time_len = max(text_to_data_units(ax, text), 0.0)

    # ------------------------------------------------------------------
    # Set X axis parameters
    # ------------------------------------------------------------------
    xmax = time_x + 3.0*time_len
    ax.set_xlim(xmin, xmax)
    ax.set_xlabel(
        "CPU time [s]",
        fontsize=10
    )
    ax.set_xticks([])

    # ------------------------------------------------------------------
    # Set top title
    # ------------------------------------------------------------------
    ax.set_title(
        f"Functions {n_min} - {n_max} in {filename}",
        fontsize=15
    )

    # ------------------------------------------------------------------
    # Show the plot!
    # ------------------------------------------------------------------
    plt.tight_layout()
    plt.show()


def main():
    """
    Main program. Reads arguments and calls routines to
    read and plot hotspots from Intel VTune output.
    """

    # Raise an error if too few or too many arguments are passed
    if len(argv) != 4:
        msg = "Usage: python -m gcpy.examples.vtune_plot_hotspots "
        msg += " FILENAME N_MIN N_MAX"
        raise ValueError(msg)

    # Verify the variable types
    filename = argv[1]
    n_min = int(argv[2])
    n_max = int(argv[3])

    # Read the profiling data
    dframe = vtune_read_hotspots(filename)
    if dframe is None:
        raise ValueError("Could not read Intel VTune hotspots data!")

    # Plot the profiling data
    vtune_plot_functions(dframe, filename, n_min, n_max)


if __name__ == '__main__':
    main()
