#!/usr/bin/env python3
"""
Script to visualize output from the KPP standalone box model.
"""

# Imports
from os.path import abspath, basename, join
from sys import argv
from glob import glob
from datetime import datetime, timedelta
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from gcpy.constants import ENCODING
from gcpy.util import verify_variable_type


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
    longitude = ""
    level    = ""
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


def kppsa_read_csv_files(input_dir, site="all", search_str=""):
    """
    Reads all KPP standalone log files for a given site
    in a given directory.

    Args
    input_dir  : str          : Directory to search
    site       : str          : KPP standalone site name

    Returns
    dframe_all : pd.DataFrame : Observations at all levels
    
    """
    verify_variable_type(input_dir, str)
    verify_variable_type(site, str)

    # Get a list of files that match the site name
    # Make sure to convert to absolute path
    if site == "all":
        file_list = glob(abspath(join(input_dir) + f"*{search_str}*.log"))
    else:
        file_list = glob(abspath(join(input_dir, site) + f"*{search_str}*.log"))
        
    # Read data from each level into a dataframe
    # and concatenate into a new dataframe
    dframe_all = pd.DataFrame()
    for file_name in file_list:
        dframe = kppsa_read_one_csv_file(file_name)
        dframe_all = pd.concat([dframe_all, dframe], ignore_index=True)
        del dframe
        
    return dframe_all


def get_element_of_series(series, element):
    """
    Returns the first element of a pd.Series object.

    Args
    serie   : pd.Series : A pd.Series object
    element : int       : Element of the pd.Series object to return

    Returns
    value   : various   : The returned element
    """
    verify_variable_type(series, pd.Series)
    verify_variable_type(element, int)

    return list(series)[element]


def prepare_site_data(
        dframe,
        site_name,
        species,
):
    """
    """
    
    # Get data for a given species at all locations
    site_data = dframe.loc[dframe["Location"] == site_name]
    site_data = site_data.loc[site_data["Species Name"] == species]
    site_data = site_data.loc[site_data["Pressure"] >= 500]
    site_data = site_data.sort_values(by="Pressure", ascending=False)
    
    # Create the top title for the subplot for this observation site
    # (use integer lon & lat values and N/S lat and E/W lon notation)
    lat = int(round(get_element_of_series(dframe["Latitude"], 0)))
    lon = int(round(get_element_of_series(dframe["Longitude"], 0)))
    #time = get_element_of_series(dframe["DateTime"], 0)
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

    return site_data, site_title


def plot_single_site(
        fig,
        rows_per_page,
        cols_per_page,
        subplot_index,
        subplot_title,
        dframe,
        location,
        species,
):
    """
    Plots observation data vs. model data at a single station site.

    Args:
    obs_dataframe  : mpl.figure.Figure : Figure object for the plot
    rows_per_page  : int               : # of rows to plot on a page
    cols_per_page  : int               : # of columns to plot on a page
    subplot_index  : int               : Index of each subplot
    subplot_title  : str               : Title for each subplot
    subplot_ylabel : str               : Y-axis title for each subplot
    obs_dataframe  : pd.DataFrame      : Observations at each station site
    obs_site_name  : str               : Name of the station site
    ref_series     : pd.Series         : Data from the Ref model version
    ref_label      : str               : Label for the Ref model data
    dev_dataarray  : pd.Series         : Data from the Dev model version
    dev_label      : str               : Label for the Dev model data
    """

    # Error check
    if species not in list(dframe["Species Name"]):
        raise ValueError("Species is not included in the mechaism!")

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
        fontsize=8
    )

    # Initial concentration
    axes_subplot.plot(
        dframe["Initial Concentration (molec/cm3)"].astype(np.float64),
        dframe["Pressure"],
        color='k',
        marker='^',
        markersize=3,
        lw=1,
        label=f"Initial {species}",
    )

    # Final concentration
    axes_subplot.plot(
        dframe["Final Concentration (molec/cm3)"].astype(np.float64),
        dframe["Pressure"],
        color='r',
        marker='o',
        markersize=3,
        lw=1,
        label=f"Final {species}",
    )

    # Set X and Y axis labels
    axes_subplot.set_xlabel(
        f"{species} (molec/cm3)",
        fontsize=8
    )

    # Apply y-axis label only if this is a leftmost plot panel
    if subplot_index == 0 or subplot_index % cols_per_page == 0:
        axes_subplot.set_ylabel(
            "Pressure (hPa)",
            fontsize=8
        )

    # Set Y-axis range
    axes_subplot.set_ylim(
        1000.0,
        500.0
    )
    axes_subplot.set_yticks(
        [1000, 900, 800, 700, 600, 500]
    )
    axes_subplot.tick_params(
        axis='both',
        which='major',
        labelsize=8
    )


def plot_one_page(
        pdf,
        site_names,
        dframe,
        species="O3",
        rows_per_page=3,
        cols_per_page=3,
        pdf_file=f"site.pdf",
):
    """
    Plots a single page of models vs. observations.

    Args:
    obs_dataframe   : pd.DataFrame : Observations at each station site.
    obs_label       : str          : Label for the observational data
    obs_site_coords : dict         : Coords (lon/lat/alt) at each site.
    obs_site_names  : list         : Names of station sites per page
    ref_dataarray   : xr.DataArray : Data from the Ref model version
    ref_label       : str          : Label for the Ref model data
    ref_cs_grid     : str|None     : Metadata for Ref cubed-sphere grid
    dev_dataarray   : xr.DataArray : Data from the Dev model version
    dev_label       : str          : Label for the Dev model data
    dev_cs_grid     : str|None     : Metadata for Dev cubed-sphere grid
    gc_levels       : pd.DataFrame : Metadata for model vertical levels
    rows_per_page   : int          : Number of rows to plot on a page
    varname         : str          : Variable name for model data
    """

    # Define a new matplotlib.figure.Figure object for this page
    # Landscape width: 11" x 8"
    fig = plt.figure(figsize=(11, 8))
    fig.tight_layout()
    
    # Loop over all of the stations that fit on the page
    for subplot_index, site_name in enumerate(site_names):

        # Get the data 
        site_data, site_title = prepare_site_data(
            dframe,
            site_name,
            species,
        )

        # Plot species vertical profile at a given site
        plot_single_site(
            fig,
            rows_per_page,
            cols_per_page,
            subplot_index,
            site_title,
            site_data,
            site_name,
            species,
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


def get_unique_site_names(dframe):
    """
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

    
def make_kpp_standalone_plots(search_str, species):
    """
    """

    # Read data
    dframe = kppsa_read_csv_files("./output4/", search_str=search_str)

    # Figure setup
    plt.style.use("seaborn-v0_8-darkgrid")
    rows_per_page = 3
    cols_per_page = 2
    plots_per_page = rows_per_page * cols_per_page
    
    # Open the plot as a PDF document
    pdf = PdfPages(f"{species}.pdf")

    # Get a list of site names sorted from N to S
    site_names = get_unique_site_names(dframe)
    
    # Loop over the number of obs sites that fit on a page
    for start in range(0, len(site_names), plots_per_page):
        end = start + plots_per_page - 1
        plot_one_page(
            pdf,
            site_names[start:end+1],
            dframe,
            species,
            rows_per_page=rows_per_page,
            cols_per_page=cols_per_page
        )

    pdf.close()

    # Reset the plot style (this prevents the seaborn style from
    # being applied to other model vs. obs plotting scripts)
    plt.style.use("default")


if __name__ == '__main__':
    make_kpp_standalone_plots(argv[1], argv[2])

