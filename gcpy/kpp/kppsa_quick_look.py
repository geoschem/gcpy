#!/usr/bin/env python3
"""
Creates a "quick-look" plot from KPP-Standalone box model output.
"""

# Imports
import argparse
import matplotlib.pyplot as plt
from gcpy.util import verify_variable_type
from gcpy.kpp.kppsa_utils import \
    kppsa_get_file_list, kppsa_get_unique_site_names, \
    kppsa_plot_single_site, kppsa_prepare_site_data, \
    kppsa_read_csv_files


def kppsa_make_quick_look_plot(file_list, label, species):
    """
    Creates a quick-look plot from KPP-Standalone box model output.

    Args
    file_list : list : List of KPP-Standalone log files
    site_name : str  : Name of the site that you wish to plot
    label     : str  : Descriptive label for the data
    species   : str  : Name of the species that you wish to plot
    """
    verify_variable_type(file_list, list)
    verify_variable_type(label, str)
    verify_variable_type(species, str)

    # Read data
    dframe = kppsa_read_csv_files(file_list)

    # Get the site name from the DataFrame
    site_name = kppsa_get_unique_site_names(dframe)[0]

    # Get the data for the given species and site
    site_data, site_title = kppsa_prepare_site_data(
        dframe,
        site_name,
        species,
    )

    # Figure setup
    plt.style.use("seaborn-v0_8-darkgrid")

    # Define a new matplotlib.figure.Figure object for this page
    # Landscape width: 11" x 8"
    fig = plt.figure(figsize=(11, 8))
    fig.tight_layout()

    # Figure setup
    plt.style.use("seaborn-v0_8-darkgrid")

    # Plot species vertical profile at a given site
    kppsa_plot_single_site(
        fig,
        rows_per_page=1,
        cols_per_page=1,
        subplot_index=0,
        subplot_title=site_title,
        ref_data=site_data,
        ref_label=label,
        dev_data=None,
        dev_label=None,
        species=species,
        font_scale=2.0,
    )

    # Add top-of-page legend
    plt.legend(
        ncol=3,
        bbox_to_anchor=(0.5, 0.98),
        bbox_transform=fig.transFigure,
        loc='upper center'
    )

    # Show the plot
    plt.show()

    # Reset the plot style (this prevents the seaborn style from
    # being applied to other model vs. obs plotting scripts)
    plt.style.use("default")


def main():
    """
    Parses arguments and calls function kppsa_make_quick_look_plot
    to generate a "quick-look" plot from KPP-Standalone box model
    output.

    Command-line arguments
    --dirname (or -d) : Folder containing KPP-Standalone output
    --label   (or -l) : Label for top-of-plot
    --pattern (or -p) : Look for files matching this pattern
    --species (or -s) : Name of the species to plot
    """

    # Tell the parser which arguments to look for
    parser = argparse.ArgumentParser(
        description="Single-panel plotting example program"
    )
    parser.add_argument(
        "-d", "--dirname",
        metavar="DIRNAME",
        type=str,
        required=True,
        help="Directory containing KPP-Standalone output files"
    )
    parser.add_argument(
        "-l", "--label",
        metavar="LABEL",
        type=str,
        required=False,
        help="Descriptive label",
        default="KPP-Standalone output"
    )
    parser.add_argument(
        "-p", "--pattern",
        metavar="PATTERN",
        type=str,
        required=False,
        help="Search for file names matching this pattern",
    )
    parser.add_argument(
        "-s", "--species",
        metavar="SPECIES",
        type=str,
        required=True,
        help="Species to plot"
    )

    # Parse command-line arguments
    args = parser.parse_args()

    # Get a list of KPP-Standalone files matching the criteria
    file_list = kppsa_get_file_list(
        args.dirname,
        args.pattern,
    )

    # Create the quick look plot from KPP-Standalone output
    kppsa_make_quick_look_plot(
        file_list,
        args.label,
        args.species
    )


if __name__ == '__main__':
    main()
