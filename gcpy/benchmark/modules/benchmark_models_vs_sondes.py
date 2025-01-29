#!/usr/bin/env python
# coding: utf-8

"""Plots ozone sonde data vs. GEOS-Chem data for 1-yr benchmarks"""

import os
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from gcpy.cstools import extract_grid
from gcpy.grid import get_nearest_model_data, get_vert_grid
from gcpy.util import make_directory, replace_whitespace, verify_variable_type
from gcpy.benchmark.modules.benchmark_utils import \
    read_ref_and_dev, rename_speciesconc_to_speciesconcvv

# Tell matplotlib not to look for an X-window
os.environ["QT_QPA_PLATFORM"] = "offscreen"


def get_ref_and_dev_model_data(
        ref_filepaths,
        dev_filepaths,
        varname="SpeciesConcVV_O3",
):
    """
    Returns GEOS-Chem model ozone data as a Pandas Dataframe object

    Args
    file_path : str          : Path to the data files

    Returns
    data      : xr.DataArray : Ozone data [ppbv]
    """
    verify_variable_type(ref_filepaths, (str,list))
    verify_variable_type(dev_filepaths, (str,list))

    # Get the model data
    ref_data, dev_data = read_ref_and_dev(
        ref_filepaths,
        dev_filepaths,
        multi_file=True
    )

    # Rename variables starting with "SpeciesConc_" to "SpeciesConcVV_",
    # for backwards compatibility with GC versions prior to 14.1.0.
    ref_data = rename_speciesconc_to_speciesconcvv(ref_data)
    dev_data = rename_speciesconc_to_speciesconcvv(dev_data)

    # Return the species of interest and convert to ppbv
    return ref_data[varname]*1.0e9, dev_data[varname]*1.0e9


def get_obs_ozone_data(file_path):
    """
    Returns the ozone sonde observations as a pandas DataFrame object.

    Args
    file_path : str          : Path to ozone sonde file

    Returns
    data      : pd.DataFrame : Sonde observations (ppbv)
    """
    verify_variable_type(file_path, str)

    data = pd.read_csv(file_path, delimiter=',', index_col=0)
    data = data.rename(columns={"Latitude": "lat", "Longitude": "lon"})

    return data


def get_site_coords(
        obs_data,
        obs_site_metadata,
        site_name,
):
    """
    Returns the lon, lat, and surface pressure at each sonde site.

    Args
    obs_data          : pd.DataFrame : Observational data
    obs_site_metadata : pd.DataFrame : Surface pressure at observation sites
    site_name         : str          : Observation site name

    Returns
    lat               : float        : Latitude at observation site
    lon               : float        : Longitude at observation site
    p_sfc             : float        : Surface pressure (hPa) at obs site
    """
    search = obs_data['Site'] == site_name
    lat = obs_data[search]['lat'].iloc[0]
    lon = obs_data[search]['lon'].iloc[0]

    search = obs_site_metadata['Site'] == site_name
    p_sfc = obs_site_metadata[search]["Surface Pressure (hPa)"].iloc[0]

    return lat, lon, p_sfc


def get_seasonal_means(
        obs_data,
        ref_data,
        dev_data,
        months,
        varname="SpeciesConcVV_O3",
):
    """
    Returns seasonally averaged data for the observations
    and models.

    Args
    obs_data      : pd.DataFrame : O3 observations (ppbv)
    ref_data      : pd.DataFrame : O3 from Ref model (ppbv)
    dev_data      : pd.DataFrame : O3 from Dev model (ppbv)
    months        : list         : Months in each season

    Keyword Args
    varname       : str          : GEOS-Chem diagnostic variable name

    Returns
    obs_seas_mean : pd.DataFrame : Seasonal mean O3 observations (ppbv)
    std_dev       : pd.DataFrame : Standard deviation in mean_obs (ppbv)
    ref_seas_mean : pd.DataFrame ; Seasonal mean O3 from Ref (ppbv)
    dev_seas_mean : pd.DataFrame : Seasonal mean O3 from Dev (ppbv)
    """
    verify_variable_type(obs_data, pd.DataFrame)
    verify_variable_type(ref_data, pd.DataFrame)
    verify_variable_type(dev_data, pd.DataFrame)
    verify_variable_type(months, list)

    # Filter data for season
    obs_this_season = obs_data[(obs_data['month'].isin(months))]
    ref_this_season = ref_data[(ref_data["month"].isin(months))]
    dev_this_season = dev_data[(dev_data["month"].isin(months))]

    # Take the seasonal means
    obs_seas_mean = obs_this_season.groupby('pressure')['o3_ppb'].mean()
    std_dev = obs_this_season.groupby('pressure')['o3_ppb_sd'].mean()
    ref_seas_mean = ref_this_season.groupby('pressure')[varname].mean()
    dev_seas_mean = dev_this_season.groupby('pressure')[varname].mean()

    return obs_seas_mean, std_dev, ref_seas_mean, dev_seas_mean


def get_nearest_model_data_to_obs(
        data,
        lon,
        lat,
        p_sfc,
        cs_grid=None,
):
    """
    Returns the nearest GEOS-Chem data to an observation.  Also
    inserts the GEOS-Chem pressure levels into the dataset.

    Args
    data    : xr.DataArray    : GEOS-Chem data
    lon     : float           : Longitude at obs site
    lat     : float           : Latitude at obs site
    p_sfc   : float           : Sfc pressure at obs site

    Keyword Args
    cs_grid : xr.Dataset|None : Cubed-sphere grid metadata
    """
    verify_variable_type(data, xr.DataArray)
    verify_variable_type(cs_grid, (xr.Dataset, type(None)))

    # Nearest data to the observation (pd.DataFrame)
    nearest = get_nearest_model_data(
        data,
        lon,
        lat,
        cs_grid
    ).reset_index()

    # Vertical pressure grid for GEOS-Chem.
    _, p_mid, _ = get_vert_grid(data, p_sfc=p_sfc)

    # Repeat surface pressure for 12 months of data
    p_mid = np.tile(p_mid, 12)

    # Add "Pressure" and "month" columns
    nearest.insert(4, "pressure", p_mid)
    nearest["month"] = nearest["time"].dt.month

    return nearest


def plot_one_site(
        axes_subplot,
        season,
        season_idx,
        site_idx,
        obs_seas_mean,
        std_dev,
        ref_label,
        ref_seas_mean,
        dev_label,
        dev_seas_mean,
):
    """
    Plots ozonesonde vs. model data for all seasons at a single site
    (i.e. one column of a page).

    Args
    axes_supblot  : mpl.Axes  : The current subplot
    season        : str       : Name of the season
    season_idx    : int       : Index of the current season
    site_idx      : int       : Index of the current site
    obs_seas_mean : pd.Series : Seasonal mean O3 obs (ppbv)
    std_dev       : pd.Series : Std dev of mean_obs (ppbv)
    ref_label     : str       : Label for the Ref version
    ref_seas_mean : pd.Series ; Seasonal mean O3 from Ref (ppbv)
    dev_label     : str       : Label for the Dev version
    dev_seas_mean : pd.Series ; Seasonal mean O3 from Dev (ppbv)
    """
    verify_variable_type(season, str)
    verify_variable_type(season_idx, int)
    verify_variable_type(site_idx, int)
    verify_variable_type(obs_seas_mean, pd.Series)
    verify_variable_type(std_dev, pd.Series)
    verify_variable_type(ref_label, str)
    verify_variable_type(ref_seas_mean, pd.Series)
    verify_variable_type(dev_label, str)
    verify_variable_type(dev_seas_mean, pd.Series)

    # Plotting observations with error bars
    axes_subplot.errorbar(
        obs_seas_mean,
        obs_seas_mean.index,
        xerr=std_dev,
        fmt='-o',
        color='black',
        markersize=3,
        label='Observations'
    )

    # Plotting model data
    axes_subplot.plot(
        ref_seas_mean,
        ref_seas_mean.index,
        color='r',
        marker='o',
        markersize=3,
        lw=1,
        label=ref_label
    )
    axes_subplot.plot(
        dev_seas_mean,
        dev_seas_mean.index,
        color='g',
        marker='o',
        markersize=3,
        lw=1,
        label=dev_label,
    )

    # Inverting y-axis and setting y-axis to start from 1000 hPa
    axes_subplot.invert_yaxis()
    axes_subplot.set_ylim(1000, 50)

    # Setting x-axis scale to 0-160 ppb
    axes_subplot.set_xlim(0, 160)

    # X/Y ticks
    axes_subplot.tick_params(axis='both', labelsize=12)

    # Get current x-axis tick values
    xticks = axes_subplot.get_xticks()

    # Exclude '0' from the x-axis tick labels
    xtick_labels = [f'{tick:.0f}' if tick != 0 else '' for tick in xticks]

    # Set modified tick labels
    axes_subplot.set_xticks(xticks)
    axes_subplot.set_xticklabels(xtick_labels)

    # Add a legend in the fourth column and first row
    if season_idx == 1 and site_idx == 4:
        axes_subplot.legend(
            loc='lower right',
            fontsize='large'
        )

    # Adding latitude and season as text inside the plot
    if site_idx == 1:
        axes_subplot.text(
            0.95,
            0.2,
            f'{season}',
            transform=axes_subplot.transAxes,
            fontsize=20,
            verticalalignment='top',
            horizontalalignment="right",
        )


def page_adjustments(fig):
    """
    Adjusts the page settings after all the subplots have been made.

    Args
    fig : mpl.Figure : Figure object
    """

    # Eliminating space between subplots completely
    plt.subplots_adjust(wspace=0, hspace=0)

    # Adjust subplot layout
    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)

    # Setting common labels
    fig.text(
        0.5,
        0.01,
        'Ozone Concentration (ppb)',
        ha='center',
        fontsize=20
    )
    fig.text(
        0.0,
        0.5,
        'Pressure (hPa)',
        va='center',
        rotation='vertical',
        fontsize=20
    )


def sort_sites_by_lat(
        obs_data
):
    """
    Returns a list of sonde sites sorted from N to S in latitude.

    Args
    obs_data   : pd.DataFrame : Observations from sondes (all sites)

    Returns
    site_names : list         : Sorted list of site names N to S
    """
    sites = obs_data[['Site', "lat"]].sort_values(by=["lat"], ascending=False)

    return sites["Site"].unique()


def plot_the_data(
        obs_data,
        obs_site_metadata,
        ref_label,
        ref_data,
        dev_label,
        dev_data,
        dst="./benchmark",
        varname="SpeciesConcVV_O3",
):
    """
    Creates plots of model data vs. ozonesonde data

    Args
    obs_data          : pd.DataFrame : Observational data
    obs_site_metadata : pd.DataFrame : Sfs pressure at observation sites
    ref_label         : str          : Label for Ref model data
    ref_data          : xr.DataArray : Ref model data
    dev_label         : str          : Label for Dev model data
    data_dev          : xr.DataArray : Dev model data
    pdf_path          : str          : Path to PDF output file

    Keyword Args
    dst               : str          : Destination folder
    varname           : str          : Name of variable to plot
    """
    # ==================================================================
    # Initialization
    # ==================================================================

    # Get cubed-sphere metadata (will be None for non-CS grids)
    # This is needed to find the nearest grid box to the observations
    ref_cs_grid = extract_grid(ref_data)
    dev_cs_grid = extract_grid(dev_data)

    # Define seasons
    seasons = {
        'DJF': [12, 1, 2],
        'JJA': [3, 4, 5],
        'MAM': [6, 7, 8],
        'SON': [9, 10, 11]
    }


    # List of site names from N to S latitude
    sorted_sites = sort_sites_by_lat(obs_data)

    # Creating the PDF with no space between subplots
    # Open the plot as a PDF document
    pdf_file = f"{dst}/models_vs_sondes.{varname.split('_')[1]}.pdf"
    pdf_pages = PdfPages(pdf_file)

    # ==================================================================
    # Loop over pages (4 sites per page)
    # ==================================================================
    for page_idx in range(0, len(sorted_sites), 4):
        fig, axes = plt.subplots(
            4, 4,
            figsize=(15, 15),
            sharex=True,
            sharey=True
        )

        # Selecting the next four sites for this page
        sites_for_page = sorted_sites[page_idx:page_idx+4]

        # ==============================================================
        # Loop over sites
        # ==============================================================
        for site_idx, site in enumerate(sites_for_page, 1):

            # Get coordinates of observation site
            lat, lon, p_sfc = get_site_coords(
                obs_data,
                obs_site_metadata,
                site,
            )

            # Get nearest model data to the observation as pd.DataFrame
            nearest_ref = get_nearest_model_data_to_obs(
                ref_data,
                lon,
                lat,
                p_sfc,
                cs_grid=ref_cs_grid,
            )
            nearest_dev = get_nearest_model_data_to_obs(
                dev_data,
                lon,
                lat,
                p_sfc,
                cs_grid=dev_cs_grid,
            )

            # Adding site names at the top of each column
            axes[0, site_idx-1].set_title(f'{site} ({lat}°)', size=15)

            # ==========================================================
            # Loop over seasons (DJF, JJA, MAM, SON)
            # ==========================================================
            for season_idx, (season, months) in enumerate(seasons.items(), 1):

                # Get the seasonal means of observations & model data
                obs_seas_mean, std_dev, \
                ref_seas_mean, dev_seas_mean = get_seasonal_means(
                    obs_data[obs_data["Site"] == site],
                    nearest_ref,
                    nearest_dev,
                    months
                )

                # Create a plot for a single site (all seasons)
                axes_subplot = axes[season_idx-1, site_idx-1]
                plot_one_site(
                    axes_subplot,
                    season,
                    season_idx,
                    site_idx,
                    obs_seas_mean,
                    std_dev,
                    ref_label,
                    ref_seas_mean,
                    dev_label,
                    dev_seas_mean,
                )

        # ================================================================
        # Global page adjustments
        # ================================================================
        page_adjustments(fig)
        pdf_pages.savefig()
        plt.close()

    # ================================================================
    # Save the PDF
    # ================================================================
    pdf_pages.close()


def make_benchmark_models_vs_sondes_plots(
        obs_data_file,
        obs_site_file,
        ref_filepaths,
        ref_label,
        dev_filepaths,
        dev_label,
        dst="./benchmark",
        overwrite=False,
        varname="SpeciesConcVV_O3",

    ):
    """
    Creates plots of sonde data vs. GEOS-Chem output.  For use in the
    1-year benchmark plotting workflow.

    Args
    obs_data_file : str      : File containing sonde data
    obs_site_file : str      : File containing sonde site metadata
    ref_filepaths : str|list : Files for the GEOS-Chem Ref version
    ref_label     : str      : GEOS-Chem Ref version label
    dev_filepaths : str|list : Files for the GEOS-Chem Dev version
    dev_label     : str      : GEOS-Chem Dev version label

    Keyword Args
    dst           : str      : Folder where PDF w/ plots will be created
    overwrite     : bool     : Overwrite contents of dst folder?
    varname       : str      : GEOS-Chem diagnostic variable name
    verbose       : bool     : Activate verbose printout?

    """
    verify_variable_type(obs_data_file, str)
    verify_variable_type(obs_site_file, str)
    verify_variable_type(ref_filepaths, (str, list))
    verify_variable_type(ref_label, str)
    verify_variable_type(dev_filepaths, (str, list))
    verify_variable_type(dev_label, str)

    # Replace whitespace in the ref and dev labels
    ref_label = replace_whitespace(ref_label)
    dev_label = replace_whitespace(dev_label)

    # Create the destination folder
    make_directory(
        dst,
        overwrite=overwrite
    )

    # Get the observational data and site metadata
    obs_data = get_obs_ozone_data(obs_data_file)
    obs_site_metadata = pd.read_csv(obs_site_file)

    # Get the model data from Ref and Dev versions
    ref_data, dev_data = get_ref_and_dev_model_data(
        ref_filepaths,
        dev_filepaths,
        varname=varname,
    )

    # Create plots
    plot_the_data(
        obs_data,
        obs_site_metadata,
        ref_label,
        ref_data,
        dev_label,
        dev_data,
        dst,
    )
