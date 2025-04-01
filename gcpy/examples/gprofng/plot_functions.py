#!/usr/bin/env python3
"""
Plots the top functions by exclusive time from a gprofng profile.
"""
from sys import argv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from gcpy.constants import ENCODING


def gprofng_read_functions(filename):
    """
    Reads function profiling information from gprofng output.

    Args
    filename : str          : Name of file w/ profiling information

    Returns
    data     : pd.DataFrame : Profiling information read from filename
    """
    results = {}

    with open(filename, "r", encoding=ENCODING) as ifile:
        line_count = 0

        for line in ifile:

            # Skip the first 5 header lines
            line_count += 1
            if line_count < 6:
                continue

            # Split timing info into columns
            columns = [var for var in line.strip().split(" ") if len(var) > 0]
            if len(columns) != 5:
                continue

            # Store in a dict
            results[columns[4]] = {
                "Exclusive time": np.float64(columns[0]),
                "Exclusive %": np.float64(columns[1]),
                "Inclusive time": np.float64(columns[2]),
                "Inclusive %": np.float64(columns[3])
            }

    return pd.DataFrame(results).transpose()


def gprofng_text_to_data_units(ax, text):
    """
    Computes the width of a label in data units.

    Args
    ax     : mpl.Axes.Subplot : MatPlotLib Axes.Subplot object
    text   : ax.text          : Text that is being plotted

    Returns
    length : float            : Length in data units
    """

    # Get the extent of the text as a Bbox object
    bbox = text.get_window_extent()

    # Convert Bbox width from pixels to data units
    inv = ax.transData.inverted()
    pixel_to_data = inv.transform([[bbox.x0, bbox.y0], [bbox.x1, bbox.y0]])
    width_in_data_units = abs(pixel_to_data[1][0] - pixel_to_data[0][0])

    return width_in_data_units


def gprofng_plot_functions(dframe, filename, n_min, n_max):
    """
    Plots functions having the largest exclusive time from
    gprofng profiling output.

    Args
    dframe   : pd.DataFrame : Profiling information
    filename : str          : Name of the file
    n_min    : int          : Display functions starting with this index
    n_max    : int          : Display functions ending with this index
    """
    # ------------------------------------------------------------------
    # Prepare the data
    # ------------------------------------------------------------------

    # Extract info from the data
    labels = list(dframe.index)
    times  = [var[0] for var in dframe[["Exclusive time"]].to_numpy()]
    percents = [var[0] for var in dframe[["Exclusive %"]].to_numpy()]

    # Take the top n_funcs values (excluding total) and reverse
    labels = labels[n_min:n_max][::-1]
    times = times[n_min:n_max][::-1]
    percents = percents[n_min:n_max][::-1]

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
    percent_len = 0.0

    # Starting position of the time label at right of plot
    xmin = min(neg_values) * 1.02
    time_x = round(abs(xmin) * 0.02)

    # Create labels for exclusive times [s] and also compute
    # the length of the longest label in data units
    for hbar, time in zip(hbars, times):
        text = ax.text(
            time_x,
            hbar.get_y() + hbar.get_height()/2,
            f"{time}",
            va="center",
            ha="left",
        )
        time_len = max(gprofng_text_to_data_units(ax, text), 0.0)

    # Starting position for the percent labels at right of plpot
    percent_x = time_x + round(2.0*time_len)

    # Create labels for exclusive times [T] and also compute
    # the length of the longest label in data units
    for hbar, percent in zip(hbars, percents):
        text = ax.text(
            percent_x,
            hbar.get_y() + hbar.get_height()/2,
            f"{percent}%",
            va="center",
            ha="left",
        )
        percent_len = max(gprofng_text_to_data_units(ax, text), 0.0)

    # ------------------------------------------------------------------
    # Set X axis parameters
    # ------------------------------------------------------------------
    xmax = percent_x + 3.0*percent_len
    ax.set_xlim(xmin, xmax)
    ax.set_xlabel(
        "Exclusive time (seconds and percent)",
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
    read and plot function profiling data from gprofng output.
    """

    # Raise an error if too few or too many arguments are passed
    if len(argv) != 4:
        msg = "Usage: python -m gcpy.examples.gprofng.plot_functions "
        msg += " FILENAME N_MIN N_MAX"
        raise ValueError(msg)

    # Verify the variable types
    filename = argv[1]
    n_min = int(argv[2])
    n_max = int(argv[3])

    # Read the profiling data
    dframe = gprofng_read_functions(filename)

    # Plot the profiling data
    gprofng_plot_functions(dframe, filename, n_min, n_max)


if __name__ == '__main__':
    main()
