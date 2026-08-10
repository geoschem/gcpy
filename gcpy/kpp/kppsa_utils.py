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


def kppsa_get_file_list(input_dir, pattern=""):
    """
    Returns a list of KPP-Standalone log files matching a search
    pattern.

    Parameters
    ----------
    input_dir : str
        Directory containing KPP-Standalone log files.
    pattern : str, optional
        Glob pattern used to filter filenames. All files containing
        this string will be returned.
        Default: ``""`` (matches all files)

    Returns
    -------
    file_list : list of str
        List of file paths matching the search criteria.
    """
    return glob(join(expanduser(input_dir), f"*{pattern}*"))


def kppsa_read_one_csv_file(file_name):
    """
    Reads a single KPP-Standalone log file in CSV format into a
    pandas DataFrame.

    Metadata fields (location, timestamp, longitude, latitude,
    vertical level, and pressure) are extracted from the file
    header and appended as additional columns.

    Parameters
    ----------
    file_name : str
        Path to the CSV-format KPP-Standalone log file to be read.

    Returns
    -------
    dframe : pandas.DataFrame
        DataFrame containing species concentrations and associated
        metadata extracted from the file header. Added columns are:
        ``Location``, ``DateTime``, ``Longitude``, ``Latitude``,
        ``Level``, and ``Pressure``.
    """
    verify_variable_type(file_name, str)

    # Initialize variables
    latitude  = ""
    level     = ""
    longitude = ""
    location  = ""
    pressure  = ""
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

            # Exit the loop once the column header row is reached
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
    Reads multiple KPP-Standalone log files and concatenates them
    into a single pandas DataFrame.

    Parameters
    ----------
    file_list : list of str
        Paths to KPP-Standalone CSV log files to be read.

    Returns
    -------
    dframe_all : pandas.DataFrame
        Concatenated DataFrame containing data from all files in
        ``file_list``.
    """
    dframe_all = pd.DataFrame()
    for file_name in file_list:
        dframe = kppsa_read_one_csv_file(file_name)
        dframe_all = pd.concat([dframe_all, dframe], ignore_index=True)
        del dframe

    return dframe_all


def kppsa_prepare_site_data(dframe, site_name, species):
    """
    Returns a DataFrame and plot title for a given species and
    observation site.

    Data is filtered to include only observations at or below
    500 hPa (surface to lower stratosphere) and sorted by
    descending pressure.

    Parameters
    ----------
    dframe : pandas.DataFrame
        KPP-Standalone output data as returned by
        :func:`kppsa_read_csv_files`.
    site_name : str
        Name of the observation site to extract data for.
    species : str
        Name of the species to extract data for.

    Returns
    -------
    site_data : pandas.DataFrame or None
        DataFrame filtered to the given site and species, sorted by
        descending pressure, and limited to pressures >= 500 hPa.
        Returns ``None`` if ``dframe`` is ``None``.
    site_title : str or None
        Plot title string including the site name, cardinal
        coordinates, and timestamp (e.g.
        ``"Mace Head (53°N,10°W) at 2019-07-01T00:00"``).
        Returns ``None`` if ``dframe`` is ``None``.
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
    lat  = int(round(get_element_of_series(site_data["Latitude"],  0)))
    lon  = int(round(get_element_of_series(site_data["Longitude"], 0)))
    time = get_element_of_series(site_data["DateTime"], 0)
    ystr = "S" if lat < 0 else "N"
    xstr = "W" if lon < 0 else "E"
    lon  = abs(lon)
    lat  = abs(lat)
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
    Plots initial and final species concentrations as vertical
    profiles at a single observation site.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure object into which the subplot will be added.
    rows_per_page : int
        Number of subplot rows on the page.
    cols_per_page : int
        Number of subplot columns on the page.
    subplot_index : int
        Zero-based index of this subplot within the page grid.
    subplot_title : str
        Title string displayed above this subplot.
    ref_data : pandas.DataFrame
        KPP-Standalone output for the ``Ref`` model version at this
        site, as returned by :func:`kppsa_prepare_site_data`.
    ref_label : str
        Legend label for the ``Ref`` model final concentration.
    dev_data : pandas.DataFrame or None
        KPP-Standalone output for the ``Dev`` model version at this
        site. Pass ``None`` to omit the Dev profile.
    dev_label : str or None
        Legend label for the ``Dev`` model final concentration.
        Ignored if ``dev_data`` is ``None``.
    species : str
        Name of the species being plotted, used for axis labels
        and error checking.
    font_scale : float
        Multiplicative scale factor applied to all font sizes,
        allowing text to be enlarged or reduced uniformly.

    Raises
    ------
    ValueError
        If ``species`` is not found in ``ref_data``.
    ValueError
        If ``dev_data`` is not ``None`` and ``species`` is not
        found in ``dev_data``.

    Notes
    -----
    The Y-axis (pressure) is plotted from 1020 hPa at the bottom
    to 500 hPa at the top. The Y-axis label is only applied to
    leftmost panels (``subplot_index`` is 0 or a multiple of
    ``cols_per_page``).
    """
    # Error checks
    if species not in list(ref_data["Species Name"]):
        raise ValueError("Species not found in Ref model output!")
    if dev_data is not None:
        if species not in list(dev_data["Species Name"]):
            raise ValueError("Species not found in Dev model output!")

    # Create matplotlib axes object for this subplot
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

    # Set Y-axis range and ticks
    axes_subplot.set_ylim(1020.0, 500.0)
    axes_subplot.set_yticks([1000, 900, 800, 700, 600, 500])
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
    Plots a single page of Ref vs. Dev vertical profiles for a
    list of observation sites and saves it to a PDF file.

    Parameters
    ----------
    pdf : matplotlib.backends.backend_pdf.PdfPages
        Open PDF file object to which the page will be saved.
    site_names : list of str
        Names of the observation sites to plot on this page.
    ref_dframe : pandas.DataFrame
        KPP-Standalone output for the ``Ref`` model version, as
        returned by :func:`kppsa_read_csv_files`.
    ref_label : str
        Legend label for the ``Ref`` model data.
    dev_dframe : pandas.DataFrame
        KPP-Standalone output for the ``Dev`` model version, as
        returned by :func:`kppsa_read_csv_files`.
    dev_label : str
        Legend label for the ``Dev`` model data.
    species : str, optional
        Name of the species to plot.
        Default: ``"O3"``
    rows_per_page : int, optional
        Number of subplot rows on the page.
        Default: ``3``
    cols_per_page : int, optional
        Number of subplot columns on the page.
        Default: ``3``
    font_scale : float, optional
        Multiplicative scale factor applied to all font sizes.
        Default: ``1.0``

    Notes
    -----
    The page is landscape-oriented (11" × 8"). A shared legend is
    placed at the top of the page. Vertical spacing between subplots
    is adjusted automatically.
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
    plt.subplots_adjust(hspace=0.4, top=0.9)

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
    Returns a list of unique observation site names from
    KPP-Standalone output, ordered from north to south.

    Parameters
    ----------
    dframe : pandas.DataFrame
        KPP-Standalone output data as returned by
        :func:`kppsa_read_csv_files`.

    Returns
    -------
    unique_site_names : list of str
        Unique site names sorted by descending latitude
        (north to south).
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
