#!/usr/bin/env python3
"""
Utility functions for visualizing output from the
KPP-Standalone box model.
"""

# Imports
from os.path import expanduser, join
from glob import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from gcpy.constants import ENCODING
from gcpy.util import get_element_of_series, verify_variable_type


def kppsa_get_file_list(
        input_dir,
        pattern=""
):
    """
    Returns a list of KPP-Standalone log files matching
    a search criteria.

    Args
    input_dir : str  : Directory with KPP-Standalone log files
    pattern   : str  : Read files matching this pattern (Default = "")

    Returns
    file_list : list : List of files matching the criteria
    """
    return glob(join(expanduser(input_dir), f"*{pattern}*"))


def kppsa_read_one_csv_file(file_name):
    """
    Reads a single log file (in CSV format) from the KPP
    standalone box model into a pandas.DataFrame object.

    Args
    file_name : str          : File to be read

    Returns
    dframe    : pd.DataFrame : DataFrame with the results
    """
    verify_variable_type(file_name, str)

    # Initialize variables
    latitude = ""
    level = ""
    longitude = ""
    location = ""
    pressure = ""
    timestamp = ""

    # Read file header contents
    with open(file_name, "r", encoding=ENCODING) as ifile:

        # Find the number of rows to skip
        skiprows = int(ifile.readline().strip()) + 1

        # Read the rest of the header contents
        for line in ifile:
            line = line.strip()

            # Extract selected metadata fields
            if "Location" in line:
                location = line.split(":")[1].strip()
            if "Timestamp" in line:
                substrs = line.split(":")
                timestamp = substrs[1].strip() + ":" + substrs[2].strip()
                timestamp = timestamp.replace("/", "-").replace(" ", "T")
            if "Longitude" in line:
                longitude = line.split(":")[1].strip()
            if "Latitude" in line:
                latitude = line.split(":")[1].strip()
            if "GEOS-Chem Vertical Level" in line:
                level = line.split(":")[1].strip()
            if "Pressure (hPa)" in line:
                pressure = line.split(":")[1].strip()

            # Exit the loop
            if "Species Name" in line:
                break

    # Read the CSV into a DataFrame object
    dframe = pd.read_csv(
        file_name,
        skiprows=skiprows,
        delimiter=",",
        dtype={
            "Initial Concentration (molec/cm3)": np.float64,
            "Final Concentration (molec/cm3)": np.float64,
        },
        engine="c"
    )

    # Add series with metadata obtained from the header
    dframe = dframe.assign(Location=location)
    dframe = dframe.assign(DateTime=np.datetime64(timestamp))
    dframe = dframe.assign(Longitude=float(longitude))
    dframe = dframe.assign(Latitude=float(latitude))
    dframe = dframe.assign(Level=int(level))
    dframe = dframe.assign(Pressure=float(pressure))

    return dframe


def kppsa_read_csv_files(file_list):
    """
    Reads all KPP standalone log files for a given site
    in a given directory.

    Args
    input_dir  : str          : Directory to search
    site       : str          : KPP standalone site name

    Returns
    dframe_all : pd.DataFrame : Observations at all levels

    """
    dframe_all = pd.DataFrame()
    for file_name in file_list:
        dframe = kppsa_read_one_csv_file(file_name)
        dframe_all = pd.concat([dframe_all, dframe], ignore_index=True)
        del dframe

    return dframe_all


def kppsa_prepare_site_data(
        dframe,
        site_name,
        species,
):
    """
    Returns a pd.DataFrame object containing data for a given species,
    and observation site, as well as the corresponding top-of-plot
    title.  Species data is limited from the surface to 500 hPa.

    Args
    dframe     : pd.DataFrame : KPP-Standalone output data
    site_name  : str          : Name of site to plot
    species    : species      : Name of species to plot
    
    Returns
    site_data  : pd.DataFrame : Data for the given site & species
    site_title : str          : Corresponding plot title string
    """

    # Exit if the data frame is set to None
    if dframe is None:
        return None, None

    # Get data for a given species at all locations
    site_data = dframe.loc[dframe["Location"] == site_name]
    site_data = site_data.loc[site_data["Species Name"] == species]
    site_data = site_data.loc[site_data["Pressure"] >= 500]
    site_data = site_data.sort_values(by="Pressure", ascending=False)

    # Create the top title for the subplot for this observation site
    # (use integer lon & lat values and N/S lat and E/W lon notation)
    lat = int(round(get_element_of_series(site_data["Latitude"], 0)))
    lon = int(round(get_element_of_series(site_data["Longitude"], 0)))
    time = get_element_of_series(site_data["DateTime"], 0)
    ystr = "S"
    if lat >= 0:
        ystr = "N"
    xstr = "W"
    if lon >= 0:
        xstr = "E"
    lon = abs(lon)
    lat = abs(lat)
    site_title = \
        f"{site_name.strip()} ({lat}$^\\circ${ystr},{lon}$^\\circ${xstr})"
    site_title += f" at {time}"

    return site_data, site_title


def kppsa_plot_single_site(
        fig,
        rows_per_page,
        cols_per_page,
        subplot_index,
        subplot_title,
        ref_data,
        ref_label,
        dev_data,
        dev_label,
        species,
        font_scale,
):
    """
    Plots observation data vs. model data at a single station site.

    Args:
    fig            : mpl.figure.Figure : Figure object for the plot
    rows_per_page  : int               : # of rows to plot on a page
    cols_per_page  : int               : # of columns to plot on a page
    subplot_index  : int               : Index of each subplot
    subplot_title  : str               : Title for each subplot
    ref_data       : pd.DataFrame      : Observations at each station site
    ref_label      : str               : Label for the Ref model data
    dev_data       : pd.DataFrame      :
    dev_label      : str               : Label for the Dev model data
    site_name      : str               : Name of the station site
    species        : pd.Series         : Data from the Ref model version
    font_scale     : float             : Scale fac to increase font size
    """

    # Error checks
    if species not in list(ref_data["Species Name"]):
        raise ValueError("Species not found in Ref model output!")
    if dev_data is not None:
        if species not in list(dev_data["Species Name"]):
            raise ValueError("Species not found in Dev model output!")

    # Create matplotlib axes object for this subplot
    # axes_subplot is of type matplotlib.axes_.subplots.AxesSubplot
    axes_subplot = fig.add_subplot(
        rows_per_page,
        cols_per_page,
        subplot_index + 1,
    )

    # Title for each subplot
    axes_subplot.set_title(
        subplot_title,
        weight='bold',
        fontsize=8*font_scale,
    )

    # Initial concentration
    axes_subplot.plot(
        ref_data["Initial Concentration (molec/cm3)"].astype(np.float64),
        ref_data["Pressure"],
        color='k',
        marker='o',
        markersize=3*font_scale,
        lw=1,
        label=f"Initial {species}",
    )

    # Final concentration (Ref)
    axes_subplot.plot(
        ref_data["Final Concentration (molec/cm3)"].astype(np.float64),
        ref_data["Pressure"],
        color='r',
        marker='o',
        markersize=3*font_scale,
        lw=1,
        label=f"{ref_label}",
    )

    # Final concentration (Dev)
    if dev_data is not None:
        axes_subplot.plot(
            dev_data["Final Concentration (molec/cm3)"].astype(np.float64),
            dev_data["Pressure"],
            color='b',
            marker='o',
            markersize=3*font_scale,
            lw=1,
            label=f"{dev_label}",
        )


    # Set X and Y axis labels
    axes_subplot.set_xlabel(
        f"{species} (molec/cm3)",
        fontsize=8*font_scale,
    )

    # Apply y-axis label only if this is a leftmost plot panel
    if subplot_index == 0 or subplot_index % cols_per_page == 0:
        axes_subplot.set_ylabel(
            "Pressure (hPa)",
            fontsize=8*font_scale,
        )

    # Set Y-axis range
    axes_subplot.set_ylim(
        1020.0,
        500.0
    )
    axes_subplot.set_yticks(
        [1000, 900, 800, 700, 600, 500]
    )
    axes_subplot.tick_params(
        axis='both',
        which='major',
        labelsize=8*font_scale
    )


def kppsa_plot_one_page(
        pdf,
        site_names,
        ref_dframe,
        ref_label,
        dev_dframe,
        dev_label,
        species="O3",
        rows_per_page=3,
        cols_per_page=3,
        font_scale=1.0,
):
    """
    Plots a single page of models vs. observations.

    Args:
    pdf             : pdf          : PDF object
    ref_dframe      : pd.DataFrame : Observations at each station site.
    ref_label       : str          : Label for the observational data
    dev_dframe      : pd.DataFrame : Data from the Ref model version
    dev_label       : str          : Label for the Ref model data
    species         : str          : Name of the species to plot
    dev_dataarray   : xr.DataArray : Data from the Dev model version
    dev_label       : str          : Label for the Dev model data
    dev_cs_grid     : str|None     : Metadata for Dev cubed-sphere grid
    gc_levels       : pd.DataFrame : Metadata for model vertical levels
    rows_per_page   : int          : Number of rows to plot on a page
    cols_per_page   : int          : Number of cols to plot on a page
    font_scale      : float        : PDF output file name
    """

    # Define a new matplotlib.figure.Figure object for this page
    # Landscape width: 11" x 8"
    fig = plt.figure(figsize=(11, 8))
    fig.tight_layout()

    # Loop over all of the stations that fit on the page
    for subplot_index, site_name in enumerate(site_names):

        # Get the data for the given species and site
        ref_site_data, site_title = kppsa_prepare_site_data(
            ref_dframe,
            site_name,
            species,
        )
        dev_site_data, _ = kppsa_prepare_site_data(
            dev_dframe,
            site_name,
            species,
        )

        # Plot species vertical profile at a given site
        kppsa_plot_single_site(
            fig,
            rows_per_page,
            cols_per_page,
            subplot_index,
            site_title,
            ref_site_data,
            ref_label,
            dev_site_data,
            dev_label,
            species,
            font_scale,
        )

    # Add extra spacing around plots
    plt.subplots_adjust(
        hspace=0.4,
        top=0.9
    )

    # Add top-of-page legend
    plt.legend(
        ncol=3,
        bbox_to_anchor=(0.5, 0.98),
        bbox_transform=fig.transFigure,
        loc='upper center'
    )

    # Save this page to the PDF file
    pdf.savefig(fig)


def kppsa_get_unique_site_names(dframe):
    """
    Returns a list of unique sites where KPP-Standalone box model
    output has been archived.

    Args
    dframe     : pd.DataFrame : Object containing KPP-Standalone output

    Returns
    site_names : list of str  : List of unique site names
    """

    # Sort the sites from north to south
    site_names = dframe[["Location", "Latitude"]].sort_values(
        by="Latitude",
        ascending=False
    )

    # Get a list of unique site names, preserving the ordering
    unique_site_names = []
    seen = set()
    for site_name in site_names["Location"]:
        if site_name not in seen:
            seen.add(site_name)
            unique_site_names.append(site_name)

    return unique_site_names
