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


def gprofng_plot_functions(dframe, filename, n_funcs):
    """
    Plots functions having the largest exclusive time from
    gprofng profiling output.

    Args
    dframe   : pd.DataFrame : Profiling information
    filename : str          : Name of the file
    n_funcs  : int          : Number of functions to display

    """
    print(n_funcs)


    # Sample data
    labels = list(dframe.index)
    times  = [var[0] for var in dframe[["Exclusive time"]].to_numpy()]
    percents = [var[0] for var in dframe[["Exclusive %"]].to_numpy()]

    # Take the top n_funcs values (excluding total) and reverse
    labels = labels[1:n_funcs][::-1]
    times = times[1:n_funcs][::-1]
    percents = percents[1:n_funcs][::-1]

    # Invert values for bars to extend to the left
    neg_values = [-var for var in times]

    # Create horizontal bar plot
    _, ax = plt.subplots(figsize=(11,8))
    hbars = ax.barh(labels, neg_values, color="skyblue")

    # Add the labels at the right of each bar
    for hbar, time, percent in zip(hbars, times, percents):

        # Print exclusive time [s]
        ax.text(
            1,
            hbar.get_y() + hbar.get_height()/2,
            f"{time}",
            va="center",
            ha="left",
        )

        # Print exclusive time [%]
        ax.text(
            10,
            hbar.get_y() + hbar.get_height()/2,
            f"{percent}%",
            va="center",
            ha="left",
        )

    # Set X axis parameters.  Plot times as negative and add a cushion
    # of -5 beyond the largest value.  Disable x-axis tick values.
    ax.set_xlim(min(neg_values) - 5, 20)
    ax.set_xlabel(
        "Exclusive time (seconds and percent)",
        fontsize=10
    )
    ax.set_xticks([])

    # Set top title
    ax.set_title(
        f"Top {n_funcs} functions in {filename}",
        fontsize=15
    )

    # Show the plot!
    plt.tight_layout()
    plt.show()


def main():
    """
    Main program. Reads arguments and calls routines to
    read and plot function profiling data from gprofng output.
    """

    # Raise an error if too few or too many arguments are passed
    if len(argv) != 3:
        msg = "Usage: python -m gcpy.examples.gprofng.plot_functions "
        msg += " FILENAME N_FUNCS"
        raise ValueError(msg)

    # Verify the variable types
    filename = argv[1]
    n_funcs = int(argv[2])

    # Read the profiling data
    dframe = gprofng_read_functions(filename)

    # Plot the profiling data
    gprofng_plot_functions(dframe, filename, n_funcs)


if __name__ == '__main__':
    main()
