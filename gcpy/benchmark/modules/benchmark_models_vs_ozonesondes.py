#!/usr/bin/env python
# coding: utf-8

"""Plots ozone sonde data vs. GEOS-Chem data for 1-yr benchmarks"""

import os
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from gcpy.grid import _GEOS_72L_AP, _GEOS_72L_BP, \
    get_nearest_model_data, vert_grid
from gcpy.util import verify_variable_type
from gcpy.benchmark.modules.benchmark_utils import get_geoschem_level_metadata


def get_model_ozone_data(file_path):
    """
    Returns model ozone data as a Pandas Dataframe object

    Args:
    -----
    file_path (str) : Path to the data files

    Returns:
    --------
    data (xr.DataArray) : Ozone data [ppbv]
    """
    verify_variable_type(file_path, str)

    # Read data (if path has wildcard, will read multiple data files)
    data = xr.open_mfdataset(file_path)

    # Convert from v/v to ppb
    data *= 1.0e9

    # Only pull out the Ozone data
    if "SpeciesConcVV_O3" in data.data_vars:
        return data["SpeciesConcVV_O3"]
    if "SpeciesConc_O3" in data.data_vars:
        return data["SpeciesConc_O3"]


def get_ref_and_dev_model_data(
        ref_dir,
        dev_dir,
        collection="SpeciesConc",
        year="2019"
):
    """
    Returns the GEOS-Chem model data for Ref & Dev versions

    Args:
    -----
    ref_dir (str) : Path to the Ref model SpeciesConc files
    dev_dir (str) : Path to the Dev model SpeciesConc files

    Keyword Args (optional):
    ------------------------
    collection (str) : Diagnostic collection (default: SpeciesConc)
    year       (str) : Benchmark year (default: 2019)

    Returns:
    --------
    data_ref, data_dev (xr.DataArray) : Ref & Dev Ozone data
    """
    verify_variable_type(ref_dir, str)
    verify_variable_type(dev_dir, str)

    # Get the model data
    data_ref = get_model_ozone_data(
        f"{ref_dir}/GEOSChem.{collection}.{year}*.nc4"
    )
    data_dev = get_model_ozone_data(
        f"{dev_dir}/GEOSChem.{collection}.{year}*.nc4"
    )

    return data_ref, data_dev


def get_obs_ozone_data(file_path):
    """
    Returns the ozone sonde observations as a pandas DataFrame object.

    Args:
    -----
    file_path (str) : Path to ozone sonde file

    Returns:
    --------
    data (pd.DataFrame) : Ozonesonde observations at various sites (ppbv)
    """
    verify_variable_type(file_path, str)

    data = pd.read_csv(file_path, delimiter=',', index_col=0)
    data = data.rename(columns={"Latitude": "lat", "Longitude": "lon"})

    return data


def get_seasonal_means(obs, ref, dev, months):
    """
    Returns seasonally averaged data for the observations
    and Ref and Dev models.  Also returns the data
    """
    verify_variable_type(obs, pd.DataFrame)
    verify_variable_type(ref, pd.DataFrame)
    verify_variable_type(dev, pd.DataFrame)
    verify_variable_type(months, list)

    # Filter data for season
    seasonal_obs = obs[(obs['month'].isin(months))]
    seasonal_ref = ref[(ref["month"].isin(months))]
    seasonal_dev = dev[(dev["month"].isin(months))]

    # Take the seasonal means
    mean_obs = seasonal_obs.groupby('pressure')['o3_ppb'].mean()
    std_obs = seasonal_obs.groupby('pressure')['o3_ppb_sd'].mean()
    mean_ref = seasonal_ref.groupby('pressure')['o3_ppb'].mean()
    mean_dev = seasonal_dev.groupby('pressure')['o3_ppb'].mean()

    return mean_obs, std_obs, mean_ref, mean_dev


def get_nearest_model_data_to_obs(
        gc_data,
        lon,
        lat,
        gc_pressure,
        gc_cs_grid=None
):
    """
    Returns the nearest GEOS-Chem data to an observation.  Also
    inserts the GEOS-Chem pressure levels into the dataset.
    """
    verify_variable_type(gc_data, (xr.Dataset, xr.DataArray))
    verify_variable_type(gc_cs_grid, (xr.Dataset, type(None)))

    # Nearest data to the observation (pd.DataFrame)
    data = get_nearest_model_data(gc_data, lon, lat, gc_cs_grid).reset_index()

    # NOTE: This would more accurately represent the pressure coordinate.
    # Need to get approx surface pressure at each site.
    #pmid = vert_grid(_GEOS_72L_AP, _GEOS_72L_BP, p_sfc=X).p_mid()
    #pmid = np.tile(pmid, 12)

    # Add "Pressure" and "month" columns
    data.insert(4, "pressure", gc_pressure)
    data["month"] = data["time"].dt.month

    # Rename columns to be consistent w/ obs data
    data = data.rename(
        columns={"SpeciesConcVV_O3": "o3_ppb", "SpeciesConc_O3": "o3_ppb"}
    )

    return data


def plot_one_site(
        axes_subplot,
        season,
        season_idx,
        site_idx,
        mean_obs,
        std_obs,
        gc_ref_label,
        mean_ref,
        gc_dev_label,
        mean_dev
):
    """
    """

    # Plotting observations with error bars
    axes_subplot.errorbar(
        mean_obs,
        mean_obs.index,
        xerr=std_obs,
        fmt='-o',
        color='black',
        markersize=3,
        label='Observations'
    )

    # Plotting model data
    axes_subplot.plot(
        mean_ref,
        mean_ref.index,
        color='r',
        marker='o',
        markersize=3,
        lw=1,
        label=gc_ref_label
    )
    axes_subplot.plot(
        mean_dev,
        mean_dev.index,
        color='g',
        marker='o',
        markersize=3,
        lw=1,
        label=gc_dev_label,
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
            0.05,
            0.8,
            f'{season}',
            transform=axes_subplot.transAxes,
            fontsize=20,
            verticalalignment='top',
        )

def page_adjustments(fig):
    """
    Adjusts the page settings after all the subplots have been made.

    Args:
    -----
    fig (mpl.Figure) : Figure object
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

    #plt.tight_layout()


def plot_the_data(
        obs_data,
        gc_ref_label,
        gc_ref_data,
        gc_dev_label,
        gc_dev_data,
        gc_pressure,
        pdf_path,
        gc_cs_grid_ref=None,
        gc_cs_grid_dev=None,
):
    """
    Creates plots of model data vs. ozonesonde data

    Args:
    -----
    obs_data     (pd.DataFrame) : Observational data
    gc_ref_label (str         ) : Label for Ref model data
    gc_ref_data  (xr.DataArray) : Ref model data
    gc_dev_label (str         ) : Label for Dev model data
    data_dev     (xr.DataArray) : Dev model data
    gc_pressure  (np.ndarray  ) : GC pressures (hPa), tiled x 12 months
    pdf_path     (str         ) : Path to PDF file that will be created
    """

    # Define seasons
    seasons = {
        'DJF': [12, 1, 2],
        'JJA': [3, 4, 5],
        'MAM': [6, 7, 8],
        'SON': [9, 10, 11]
    }

    # Getting unique sites and sorting them
    sorted_sites = sorted(obs_data['Site'].unique())

    # Creating the PDF with no space between subplots
    pdf_pages = PdfPages(pdf_path)

    # ==================================================================
    # Loop over pages (4 sites per page)
    # ==================================================================
    for page_idx in range(1):
    #for page_idx in range(0, len(sorted_sites), 4):
        fig, axs = plt.subplots(
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

            # Get the lon & lat for each site
            lat = obs_data[obs_data['Site'] == site]['lat'].iloc[0]
            lon = obs_data[obs_data['Site'] == site]['lon'].iloc[0]

            # Get nearest model data to the observation as pd.DataFrame
            nearest_ref = get_nearest_model_data_to_obs(
                gc_ref_data,
                lon,
                lat,
                gc_pressure,
                gc_cs_grid=gc_cs_grid_ref
            )
            nearest_dev = get_nearest_model_data_to_obs(
                gc_dev_data,
                lon,
                lat,
                gc_pressure,
                gc_cs_grid=gc_cs_grid_dev
            )

            # Adding site names at the top of each column
            axs[0, site_idx-1].set_title(f'{site} ({lat}°)', size=15)

            # ==========================================================
            # Loop over seasons (DJF, JJA, MAM, SON)
            # ==========================================================
            for season_idx, (season, months) in enumerate(seasons.items(), 1):

                # Get the seasonal means of observations & model data
                mean_obs, std_obs, mean_ref, mean_dev = get_seasonal_means(
                    obs_data[obs_data["Site"] == site],
                    nearest_ref,
                    nearest_dev,
                    months
                )

                # Create a plot for a single site (all seasons)
                plot_one_site(
                    axes_subplot=axs[season_idx-1, site_idx-1],
                    season=season,
                    season_idx=season_idx,
                    site_idx=site_idx,
                    mean_obs=mean_obs,
                    std_obs=std_obs,
                    gc_ref_label=gc_ref_label,
                    mean_ref=mean_ref,
                    gc_dev_label=gc_dev_label,
                    mean_dev=mean_dev
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
    #print('The pdf has been closed')


def main():
    """
    Main program (for testing)
    """

    # Define the base directories for the two models (customized)
    gc_ref_label = '14.2.0-rc.2'
    gc_dev_label = '14.3.0-rc.0'
    model_root_dir = "/n/holyscratch01/jacob_lab/ryantosca/BM/1yr"
    obs_root_dir = "/n/jacob_lab/Lab/obs_data_for_bmk/atom_obs"

    # figure saveout path
    pdf_file = '/n/holyscratch01/jacob_lab/ryantosca/BM/1yr/ozonesondes.pdf'

    # Get the model data
    gc_ref_data, gc_dev_data = get_ref_and_dev_model_data(
        f"{model_root_dir}/{gc_ref_label}/GCClassic/FullChem/OutputDir",
        f"{model_root_dir}/{gc_dev_label}/GCClassic/FullChem/OutputDir"
    )

    # Get the observational data
    obs_data = get_obs_ozone_data(
        f"{obs_root_dir}/allozonesondes_2010-2019.csv"
    )

    # GEOS-Chem pressure levels, tiled for 12 months
    gc_levels = get_geoschem_level_metadata()
    gc_pressure = np.tile(gc_levels["Pressure (hPa)"].to_numpy(), 12)

    # Create plots
    plot_the_data(
        obs_data,
        gc_ref_label,
        gc_ref_data,
        gc_dev_label,
        gc_dev_data,
        gc_pressure,
        pdf_file,
    )

if __name__ == '__main__':
    main()
