#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gcpy/benchmark/modules/benchmark_model_vs_obs.py

Python functions to plot modeled data from 1-year fullchem benchmark
simulations against observations for the year 2019.  At present, only
O3 plots are supported, but this can be extended in the future.

Author: Matt Rowlinson <matthew.rowlinson@york.ac.uk>

Linted with PyLint and incorporated into GCPy
by Bob Yantosca <yantosca@seas.harvard.edu>
"""
import os
import glob
from datetime import datetime, timedelta
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import xarray as xr
from gcpy.constants import skip_these_vars
from gcpy.util import dataset_reader, make_directory, reshape_MAPL_CS

def read_nas(
        input_file,
        verbose=False,

):
    """
    Read NASA Ames data files from EBAS (https://ebas-data.nilu.no)
    Creates data frame of O3 values converted to ppb and dictionary
    with key site information (name, lat, lon, altitude)

    Args:
    -----
    input_file : str
        Path to data file with observational data (e.g. sonde data).

    Keyword Args:
    -------------
    verbose : bool
        Toggles verbose printout on (True) or off (False).
        Default value: False

    Returns:
    --------
    obs_dataframe : pandas DataFrame
        Dataframe containing observational data from input_file.

    obs_site_coords : dict
        Dictionary containing formatted site name: lon, lat and altitude.
    """
    if not isinstance(input_file, str):
        raise TypeError("Argument 'input_file' is not of type str!")

    if verbose:
        print(f"read_nas: Reading {input_file}")

    with open(input_file, encoding='UTF-8') as the_file:
        header = np.array(
            [next(the_file) for x in range(155) ]
        )
        n_hdr = int(
            header[0].split(' ')[0]
        )
        st_ymd = header[6].split(' ')
        st_ymd = list(
            filter(
                None,
                st_ymd
            )
        )
        start_date = datetime(
            int(st_ymd[0]),
            int(st_ymd[1]),
            int(st_ymd[2])
        )
        for line in header:
            if 'Station name' in line:
                site = line.split(':')[1:]
                site = '_'.join(site).replace('\n','').\
                    replace('  ',' ').replace('/','-')
                site = site.replace('Atmospheric Observatory','')
                site = site.replace(' Research Station','')
            elif 'Station longitude:' in line:
                lon = float(line.split(' ')[-1].replace('\n',''))
            elif 'Station latitude:' in line:
                lat = float(line.split(' ')[-1].replace('\n',''))
            elif 'Station altitude:' in line:
                alt = float(line.split(' ')[-2].replace('\n',''))

    file_hdr = np.loadtxt(
        input_file,
        skiprows=n_hdr
    )
    obs_dataframe = pd.DataFrame(
        file_hdr,
        index=file_hdr[:,0]
    )
    obs_dataframe, qcflag = find_times(
        obs_dataframe,
        start_date
    )
    obs_dataframe = pd.DataFrame(
        {
            'Value': obs_dataframe.values/1.99532748,
            'Flag': qcflag
        },
        index=obs_dataframe.index
    )
    obs_dataframe = obs_dataframe[obs_dataframe.Flag == 0.000]
    obs_dataframe = obs_dataframe.loc['2019']
    obs_dataframe = obs_dataframe.resample('H').mean()
    obs_dataframe = pd.DataFrame(
        {
            site: obs_dataframe.Value
        },
        index=obs_dataframe.index
    )
    obs_site_coords = { site:
          {
              'lon': lon,
              'lat': lat,
              'alt': alt
          }
    }

    return obs_dataframe, obs_site_coords


def read_observational_data(
        path,
        verbose
):
    """
    Reads the observational O3 data from EBAS
    (taken from https://ebas-data.nilu.no/ on 15/05/2023)

    Loops over all data files (in NASA/Ames format) within
    a folder and concatenates them into a single DataFrame.

    Args:
    -----
    path : str
        Path to the observational data directory

    verbose : bool
        Toggles verbose printout on (True) or off (False).
        Default value: False

    Returns:
    --------
    obs_dataframe : pandas DataFrame
        DataFrame object with the observational data (i.e. station
        names, data, metadata).

    obs_site_coords : dict
        Dictionary with coordinates of each observation site
    """
    if not isinstance(path, str):
        raise TypeError("The 'path' argument is not of type str!")

    first = True
    obs_site_coords = {}
    dataframe = None
    for infile in sorted(glob.glob(f"{path}/*nas")):
        obs_dataframe, xyz = read_nas(
            infile,
            verbose=verbose
        )
        if first:
            dataframe = obs_dataframe
            obs_site_coords.update(xyz)
            first = False
        else:
            dataframe = pd.concat(
                [dataframe, obs_dataframe],
                axis=1
            )
            obs_site_coords.update(xyz)

    # If dataframe0 is undefined, the loop didn't execute... so throw error
    if dataframe is None:
        raise ValueError(f"Could not find data in {path}!")

    obs_dataframe = dataframe.groupby(
        dataframe.columns,
        axis=1
    ).max()

    return obs_dataframe, obs_site_coords


def read_model_data(
        filepaths,
        varname,
        verbose=False,
):
    """
    Reads model data from a netCDF file.  Adds special handling to look
    for species concentrations variable names starting with either
    "SpeciesConcVV" or "SpeciesConc".  This is necessary for backwards
    compatitbility with GEOS-Chem output prior to version 14.1.0.

    Args:
    -----
    filepaths : list of str
        List of data files to read.

    varname : str or list of str
        Variable name(s) to read from data files.

    Keyword Args:
    -------------
    varbose : bool
        Toggles verbose output on (True) or off (False).
        Default value: False

    Returns:
    --------
    dataarray : xarray DataArray
        DataArray object containing data read from files
        specified by the filepaths argument.
    """

    # Read the Ref and Dev model data
    reader = dataset_reader(
        multi_files=True,
        verbose=verbose,
    )

    # Set temporary variable name for use below
    varname_tmp = varname

    # First try reading the data as-is
    try:
        dataset = reader(
            filepaths,
            drop_variables=skip_these_vars,
            data_vars=[varname_tmp]
        ).load()

    # If we encounter a ValueError, it may be because the data is
    # older # and may have e.g. SpeciesConc fields instead of
    # SpeciesConcVV fields.  Reset the varname_tmp and try again.
    except ValueError:
        varname_tmp = varname_tmp.replace("VV", "")
        dataset = reader(
            filepaths,
            drop_variables=skip_these_vars,
            data_vars=[varname_tmp]
        ).load()

        # Rename to the original name to avid confusion with data
        # from GEOS-Chem versions prior to 14.1.0
        with xr.set_options(keep_attrs=True):
            dataset = dataset.rename({varname_tmp: varname})

    # If we fail again, then throw an error!
    except [FileNotFoundError, OSError, IOError] as exc:
        msg = f"get_model_data: Could not read Ref data for {varname}!"
        raise exc(msg) from exc

    # Create a DataArray object and convert to ppbv (if necessary)
    # Reshape GCHP data so that it has "lon", "lat" dimensions,
    # which will facilitate data handling elsewhere in this module.
    with xr.set_options(keep_attrs=True):
        dataarray = dataset[varname]
        if "mol mol-1" in dataarray.attrs["units"]:
            dataarray.values *= 1.0e9
            dataarray.attrs["units"] = "ppbv"
        if "nf" in dataarray.dims:
            dataarray = reshape_MAPL_CS(
                dataarray,
                multi_index_lat=False
            )

    return dataarray


def find_times(
        obs_dataframe,
        start_time
):
    """
    Convert timestamps in nasa ames data files to python datetime
    objects  Set DataFrame index to the new datetime array

    Args:
    ----------
    obs_dataframe : pandas DataFrame
        DataFrame with O3 values from GAW site

    start_time : str
        Reference start time for timestamp taken from nasa ames file

    Returns
    ------
    obs_dataframe: pandas DataFrame
        O3 in ppbV with datetime index

    qcflag : pandas Dataframe
        QC flag with datetime index
    """
    end_time = obs_dataframe[obs_dataframe.columns[1]]
    time_x = []

    for index in range(len(end_time)):
        time_x.append(start_time + timedelta(days=end_time.values[index]))

    obs_dataframe.index = time_x
    qcflag =obs_dataframe[obs_dataframe.columns[-1]]
    obs_dataframe = obs_dataframe[obs_dataframe.columns[2]]

    return obs_dataframe, qcflag


def find_nearest_3d(
        gc_data,
        gc_level_alts_m,
        lon_value,
        lat_value,
        alt_value
):
    """
    Find GEOS-Chem gridbox closest to the observational dataset.
    Uses lat, lon and alt from obs to select most appropriate GC data.

    Args:
    -----
    gc_data : xarray DataSet
        GEOS-Chem output to be processed

    gc_level_alts_m: pandas Series
        Altitudes of GEOS-Chem levels in meters

    lon_value : float
        GAW site longitude

    lat_value : float
        GAW site latitude

    alt_value : float
        GAW site altitude

    Returns:
    --------
    x_idx, y_idx, z_idx : numpy.int64
        GEOS-Chem grid box indices for the single gridbox
        closest to GAW site specifications
    """

    x_idx=(
        np.abs(
            gc_data.lon.values - float(lon_value)
        )
    ).argmin()

    y_idx=(
        np.abs(
            gc_data.lat.values - float(lat_value)
        )
    ).argmin()

    z_idx=(
        np.abs(
            gc_level_alts_m.values - float(alt_value)
        )
    ).argmin()

    return x_idx, y_idx, z_idx


def get_geoschem_level_metadata(
        filename=None,
        search_key=None,
        verbose=False,
):
    """
    Reads a comma-separated variable (.csv) file with GEOS-Chem vertical
    level metadata and returns it in a pandas DataFrame object.

    Args:
    -----
    filename : str
        Name of the comma-separated variable to read.
        Default value: "__file__/GC_72_vertical_levels.csv"

    Keyword Args:
    -------------
    search_key : str
        If present, will return metadata that matches this value.
        Default: None

    verbose : bool
        Toggles verbose printout on (True) or off (False).
        Default value: True

    Returns:
    --------
    metadata : pandas DataFrame
        Metadata for each of the GEOS-Chem vertical levels.
    """
    if filename is None:
        filename = os.path.join(
            os.path.dirname(__file__),
            "GC_72_vertical_levels.csv"
        )

    try:
        if verbose:
            print(f"get_geoschem_level_metadata: Reading {filename}")
        metadata = pd.read_csv(filename)
    except (IOError, OSError, FileNotFoundError) as exc:
        msg = f"Could not read GEOS-Chem level metadata in {filename}!"
        raise exc(msg) from exc

    if search_key is None:
        return metadata
    return metadata[search_key]


def prepare_data_for_plot(
        obs_dataframe,
        obs_site_coords,
        obs_site_name,
        ref_dataarray,
        dev_dataarray,
        gc_level_alts_m,
        varname="SpeciesConcVV_O3"
):
    """
    Prepares data for passing to routine plot_single_frames as follows:

    (1) Computes the mean of observations at the given station site.
    (2) Returns the GEOS-Chem Ref and Dev data at the gridbox closest
         to the given station site.
    (3) Creates the top-of-plot title for the given station site.

    Args:
    -----
    obs_dataframe : pandas DataFrame
        Observations at each station site.

    obs_site_coords : dict
        Coordinates (lon, lat, alt) for each observation station site.

    obs_site_name : str
        Name of the observation station site.

    ref_dataarray, dev_dataarray : xarray DataArray
        Data from the Ref and Dev model versions.

    ref_label, dev_label: str
        Labels describing the Ref and Dev datasets (e.g. version numbers)

    gc_level_alts_m : pandas DataFrame
        Metadata pertaining to GEOS-Chem vertical levels

    Keyword Args (Optional)
    -----------------------
    varname : str
        GEOS-Chem diagnostic name for the Ref and Dev model data.
        Default value: "SpeciesConcVV_O3"

    Returns:
    --------
    obs_dataframe : pandas DataFrame
        Meanb observational data at the given station site.

    ref_series, dev_series : pandas Series
        Data from the Ref and Dev model versions at the
        closest grid box to the observation station site.

    subplot_title : str
        Plot title string for the given observation station site.

    subplot_ylabel : str
        Label for the Y-axis (e.g. species name).
    """

    # Get Ref model data nearest to the observation site
    x_idx, y_idx, z_idx = find_nearest_3d(
        ref_dataarray,
        gc_level_alts_m,
        lon_value=round(obs_site_coords[obs_site_name]['lon'], 2),
        lat_value=round(obs_site_coords[obs_site_name]['lat'], 2),
        alt_value=round(obs_site_coords[obs_site_name]['alt'], 1)
    )
    ref_dataframe = ref_dataarray.isel(
        lon=x_idx,
        lat=y_idx,
        lev=z_idx
    ).to_dataframe()

    # Get Dev model data nearest to the observation site
    x_idx, y_idx, z_idx = find_nearest_3d(
        dev_dataarray,
        gc_level_alts_m,
        lon_value=round(obs_site_coords[obs_site_name]['lon'], 2),
        lat_value=round(obs_site_coords[obs_site_name]['lat'], 2),
        alt_value=round(obs_site_coords[obs_site_name]['alt'], 1)
    )
    dev_dataframe = dev_dataarray.isel(
        lon=x_idx,
        lat=y_idx,
        lev=z_idx
    ).to_dataframe()

    # Take the monthly mean of observations for plotting
    # (since some observation sites have multiple months of data)
    obs_dataframe = obs_dataframe.resample('M').mean()

    # Create the top title for the subplot for this observation site
    # (use integer lon & lat values and N/S lat and E/W lon notation)
    lon = int(round(obs_site_coords[obs_site_name]['lon'], 0))
    lat = int(round(obs_site_coords[obs_site_name]['lat'], 0))
    ystr = "S"
    if lat >= 0:
        ystr = "N"
    xstr = "W"
    if lon >= 0:
        xstr = "E"
    lon = abs(lon)
    lat = abs(lat)
    subplot_title = \
        f"{obs_site_name.strip()} ({lat}$^\\circ${ystr},{lon}$^\\circ${xstr})"

    # Y-axis label (i.e. species name)
    subplot_ylabel = varname.split("_")[1] + " (ppbv)"

    return obs_dataframe, ref_dataframe[varname], dev_dataframe[varname], \
        subplot_title, subplot_ylabel


def plot_single_station(
        fig,
        rows_per_page,
        cols_per_page,
        subplot_index,
        subplot_title,
        subplot_ylabel,
        obs_dataframe,
        obs_site_name,
        ref_series,
        ref_label,
        dev_series,
        dev_label
):
    """
    Plots observation data vs. model data at a single station site.

    Args:
    -----
    fig : matplotlib.figure.Figure
        Matplotlib Figure object containing the plot.

    rows_per_page, cols_per_page : int
        Number of rows and columns on each page of the plot.

    subplot_index : int
        Index of the subplot on the page.  Runs from 0 to
        (cols_per_page * rows_per_page - 1).

    subplot_title, subplot_ylabel : str
        Top title and y-axis label for each subplot

    obs_dataframe : pandas DataFrame
        Observational data.

    obs_site_name: : str
        Name of the observation station site.

    ref_series, dev_series : pandas Series
        GEOS-Chem data at closest grid box to the observation
        station site for the Ref and Dev model versions.

    ref_label, dev_label : str
        Descriptive labels (e.g. version numbers) for the
        GEOS-Chem Ref and Dev model versions.
    """
    if not isinstance(fig, Figure):
        msg = "The 'fig' argument is not of type matplotlib.figure.Figure!"
    if not isinstance(cols_per_page, int):
        msg = "The 'cols_per_page' argument is not of type int!"
    if not isinstance(rows_per_page, int):
        msg = "The 'rows_per_page' argument is not of type int!"
    if not isinstance(obs_dataframe, pd.DataFrame):
        msg = "The 'obs_dataframe' argument is not of type pandas.DataFrame!"
        raise TypeError(msg)
    if not isinstance(ref_series, pd.Series):
        msg = "The 'ref_series' argument is not of type pandas.Series!"
        raise TypeError(msg)
    if not isinstance(dev_series, pd.Series):
        msg = "The 'ref_series' argument is not of type pandas.Series!"
        raise TypeError(msg)

    # Create matplotlib axes object for this subplot
    # axes_subplot is of type matplotlib.axes_.subplots.AxesSubplot
    axes_subplot = fig.add_subplot(
        rows_per_page,
        cols_per_page,
        subplot_index + 1,
    )

    # Set title for top of each frame
    axes_subplot.set_title(
        f"{subplot_title}",
        weight='bold',
        fontsize=7
        )

    ## Plot observational data
    axes_subplot.plot(
        obs_dataframe.index,
        obs_dataframe[obs_site_name],
        color='k',
        marker='^',
        markersize=4,
        lw=1,
        label='Observations'
    )

    # Plot model data
    axes_subplot.plot(
        obs_dataframe.index,
        ref_series,
        color='r',
        marker='o',
        markersize=3,
        lw=1,
        label=ref_label
    )
    axes_subplot.plot(
        obs_dataframe.index,
        dev_series,
        color='g',
        marker='s',
        markersize=3,
        lw=1,
        label=dev_label
    )

    # Apply y-axis label only if this is a leftmost plot panel
    if subplot_index == 0 or subplot_index % cols_per_page == 0:
        axes_subplot.set_ylabel(
            subplot_ylabel,
            fontsize=8
        )

    # Set X-axis and Y-axis ticks and labels
    axes_subplot.set_xticks(
        obs_dataframe.index
    )
    # NOTE: In newer versions of matplotlib you can pass the
    # xticklabels keyword to the set_xticks function.  But we need
    # to set the xticklabels separately for backwards compatibility
    # with older matplotlib versions. -- Bob Yantosca (06 Jul 2023)
    axes_subplot.set_xticklabels(
        ['J','F','M','A','M','J','J','A','S','O','N','D']
    )
    axes_subplot.set_ylim(
        0,
        100
    )
    axes_subplot.set_yticks(
        [0, 25, 50, 75, 100]
    )
    axes_subplot.tick_params(
        axis='both',
        which='major',
        labelsize=6
    )


def plot_one_page(
        pdf,
        obs_dataframe,
        obs_site_coords,
        obs_site_names,
        ref_dataarray,
        ref_label,
        dev_dataarray,
        dev_label,
        gc_level_alts_m,
        rows_per_page=3,
        cols_per_page=3,
        varname="SpeciesConcVV_O3",
):
    """
    Plots a single page of models vs. observations.

    Args:
    -----
    obs_dataframe : pandas DataFrame
        Observations at each station site.

    obs_site_coords : dict
        Coordinates (lon, lat, alt) for each observation station site.

    obs_site_names : list of str
        Names of observation station sites that fit onto a single page.

    ref_dataarray, dev_dataarray : xarray DataArray
        Data from the Ref and Dev model versions.

    ref_label, dev_label: str
        Labels describing the Ref and Dev datasets (e.g. version numbers)

    gc_level_alts_m : pandas DataFrame
        Metadata pertaining to GEOS-Chem vertical levels

    Keyword Args:
    -------------

    rows_per_page, cols_per_page : int
        Number of rows and columns to plot on a single page.
        Default values: 3 rows, 3 columns

    varname : str
        Variable name for GEOS-Chem diagnostic data.
        Default value: "SpeciesConcVV_O3"

    verbose : bool
        Toggles verbose printout on (True) or off (False).
        Default value: False
    """
    if not isinstance(obs_dataframe, pd.DataFrame):
        msg = "The 'obs_dataframe' argument is not of type pandas.DataFrame!"
        raise TypeError(msg)
    if not isinstance(obs_site_coords, dict):
        msg = "The 'obs_site_coords' argument is not of type dict!"
        raise TypeError(msg)
    if not isinstance(ref_dataarray, xr.DataArray):
        msg = "The 'ref_dataset' argument is not of type xarray.DataArray!"
        raise TypeError(msg)
    if not isinstance(ref_label, str):
        msg = "The 'ref_label' argument is not of type str!"
        raise TypeError(msg)
    if not isinstance(dev_dataarray, xr.DataArray):
        msg = "The 'ref_dataset' argument is not of type xarray.DataArray!"
        raise TypeError(msg)
    if not isinstance(dev_label, str):
        msg = "The 'dev_label' argument is not of type str!"
        raise TypeError(msg)
    if not isinstance(gc_level_alts_m, pd.Series):
        msg = "The 'gc_level_alts_m' argument is not of type pandas.Series!"
        raise TypeError(msg)

    # Define a new matplotlib.figure.Figure object for this page
    # Landscape width: 11" x 8"
    fig = plt.figure(figsize=(11, 8))
    fig.tight_layout()

    # Loop over all of the stations that fit on the page
    for subplot_index, obs_site_name in enumerate(obs_site_names):

        # Find the model Ref & Dev data closest to the observational
        # station site.  Also take monthly average of observations,
        obs_dataframe, \
        ref_series, dev_series, \
        subplot_title, subplot_ylabel \
        = prepare_data_for_plot(
            obs_dataframe,                # pandas.DataFrame
            obs_site_coords,              # dict
            obs_site_name,                # str
            ref_dataarray,                # xarray.DataArray
            dev_dataarray,                # xarray.DataArray
            gc_level_alts_m,              # pandas.Series
            varname=varname,              # str
        )

        # Plot models vs. observation for a single station site
        plot_single_station(
            fig,                          # matplotlib.figure.Figure
            rows_per_page,                # int
            cols_per_page,                # int
            subplot_index,                # int
            subplot_title,                # str
            subplot_ylabel,               # str
            obs_dataframe,                # pandas.Dataframe
            obs_site_name,                # str
            ref_series,                   # pandas.Series
            ref_label,                    # str
            dev_series,                   # pandas.Series
            dev_label                     # str
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


def plot_models_vs_obs(
        obs_dataframe,
        obs_site_coords,
        ref_dataarray,
        ref_label,
        dev_dataarray,
        dev_label,
        gc_level_alts_m,
        varname="SpeciesConcVV_O3",
        dst="./benchmark",
        verbose=False
):
    """
    Plots models vs. observations using a 3 rows x 3 column layout.

    Args:
    -----
    obs_dataframe : pandas DataFrame
        Observations at each station site.

    obs_site_coords : dict
        Coordinates (lon, lat, alt) for each observation station site.

    ref_dataarray, dev_dataarray : xarray DataArray
        Data from the Ref and Dev model versions.

    ref_label, dev_label: str
        Labels describing the Ref and Dev datasets (e.g. version numbers)

    gc_level_alts_m : pandas DataFrame
        Metadata pertaining to GEOS-Chem vertical levels

    Keyword Args:
    -------------
    varname : str
        Variable name for GEOS-Chem diagnostic data.
        Default value: "SpeciesConcVV_O3"

    dst : str
        Root folder where output will be created.
        Default value: "./benchmark"

    verbose : bool
        Toggles verbose printout on (True) or off (False).
        Default value: False
    """
    if not isinstance(obs_dataframe, pd.DataFrame):
        msg = "The 'obs_dataframe' argument is not of type pandas.DataFrame!"
        raise TypeError(msg)
    if not isinstance(ref_dataarray, xr.DataArray):
        msg = "The 'ref_dataset' argument is not of type xarray.DataArray!"
        raise TypeError(msg)
    if not isinstance(ref_label, str):
        msg = "The 'ref_label' argument is not of type str!"
        raise TypeError(msg)
    if not isinstance(dev_dataarray, xr.DataArray):
        msg = "The 'ref_dataset' argument is not of type xarray.DataArray!"
        raise TypeError(msg)
    if not isinstance(dev_label, str):
        msg = "The 'dev_label' argument is not of type str!"
        raise TypeError(msg)
    if not isinstance(gc_level_alts_m, pd.Series):
        msg = "The 'gc_level_alts_m' argument is not of type pandas.Series!"
        raise TypeError(msg)

    # Figure setup
    plt.style.use('seaborn-darkgrid')
    rows_per_page = 3
    cols_per_page = 3
    plots_per_page = rows_per_page * cols_per_page

    # Open the plot as a PDF document
    pdf_file = f"{dst}/models_vs_obs.surface.{varname.split('_')[1]}.pdf"
    pdf = PdfPages(pdf_file)

    # Sort station sites N to S latitude order according to:
    # https://www.geeksforgeeks.org/python-sort-nested-dictionary-by-key/
    # NOTE: obs_site_names will be a MultiIndex list (a list of tuples)
    obs_site_names = sorted(
        obs_site_coords.items(),
        key = lambda x: x[1]['lat'],
        reverse=True
    )

    # Convert obs_site_names from a MultiIndex list to a regular list
    obs_site_names = [list(tpl)[0] for tpl in obs_site_names]

    # Loop over the number of obs sites that fit on a page
    for start in range(0, len(obs_site_names), plots_per_page):
        end = start + plots_per_page - 1

        # Plot obs sites that fit on a single page
        plot_one_page(
            pdf,                          # PdfPages
            obs_dataframe,                # pandas.DataFrame
            obs_site_coords,              # dict
            obs_site_names[start:end+1],  # list of str
            ref_dataarray,                # xarray.DataArray
            ref_label,                    # str
            dev_dataarray,                # xarray.DataArray
            dev_label,                    # str
            gc_level_alts_m,              # pandas.Series
            rows_per_page=rows_per_page,  # int
            cols_per_page=cols_per_page,  # int
            varname=varname               # str
        )

    # Close the PDF file after all pages are plotted.
    pdf.close()


def make_benchmark_models_vs_obs_plots(
        obs_filepaths,
        ref_filepaths,
        ref_label,
        dev_filepaths,
        dev_label,
        varname="SpeciesConcVV_O3",
        dst="./benchmark",
        verbose=False,
        overwrite=False
):
    """
    Driver routine to create plots
    """
    if not isinstance(obs_filepaths, list) and \
       not isinstance(obs_filepaths, str):
        msg = "The 'obs_filepaths' argument is not of type 'list' or 'str'!"
        raise TypeError(msg)
    if not isinstance(ref_filepaths, list) and \
       not isinstance(ref_filepaths, str):
        msg = "The 'ref_filepaths' argument is not of type 'list' or 'str'!"
        raise TypeError(msg)
    if not isinstance(ref_label, str):
        msg = "The 'ref_label' argument is not of type 'str'!"
        raise TypeError(msg)
    if not isinstance(dev_filepaths, list) and \
       not isinstance(dev_filepaths, str):
        msg = "The 'dev_filepaths' argument is not of type 'list' or 'str'!"
        raise TypeError(msg)
    if not isinstance(ref_label, str):
        msg = "The 'dev_label' argument is not of type 'str'!"
        raise TypeError(msg)

    # Create the destination folder
    make_directory(
        dst,
        overwrite=overwrite
    )

    # Get altitude [m] of GEOS-Chem level edges
    gc_level_alts_m = \
        get_geoschem_level_metadata(
            search_key="Altitude (km)"
        ) * 1.0e3

    # Read the observational data
    obs_dataframe, obs_site_coords = read_observational_data(
        obs_filepaths,
        verbose=verbose
    )

    # Read the model data
    ref_dataarray = read_model_data(
        ref_filepaths,
        varname=varname
    )
    dev_dataarray = read_model_data(
        dev_filepaths,
        varname=varname
    )

    # Plot data vs observations
    plot_models_vs_obs(
        obs_dataframe,                    # pandas.DataFrame
        obs_site_coords,                  # dict
        ref_dataarray,                    # xarray.DataArray
        ref_label,                        # str
        dev_dataarray,                    # xarray.DataArray
        dev_label,                        # str
        gc_level_alts_m,                  # pandas.Series
        varname=varname,                  # str
        dst=dst,                          # str
        verbose=verbose                   # bool
    )
