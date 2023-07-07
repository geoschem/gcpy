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
from math import ceil
import glob
from datetime import datetime, timedelta
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import xarray as xr
from gcpy.constants import skip_these_vars
from gcpy.util import dataset_reader, make_directory

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
    dataframe : pandas DataFrame
        Dataframe containing observational data from input_file.

    coords_dict : dict
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
    dataframe = pd.DataFrame(
        file_hdr,
        index=file_hdr[:,0]
    )
    dataframe, flag = find_times(
        dataframe,
        start_date
    )
    dataframe = pd.DataFrame(
        {
            'Value': dataframe.values/1.99532748,
            'Flag': flag
        },
        index=dataframe.index
    )
    dataframe = dataframe[dataframe.Flag == 0.000]
    dataframe = dataframe.loc['2019']
    dataframe = dataframe.resample('H').mean()
    dataframe = pd.DataFrame(
        {
            site: dataframe.Value
        },
        index=dataframe.index
    )
    coords_dict = { site:
          {
              'lon': lon,
              'lat': lat,
              'alt': alt
          }
    }
    return dataframe, coords_dict


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
    dataframe : pandas DataFrame
        DataFrame object with the observational data (i.e. station
        names, data, metadata).

    coords_dict : dict
        Dictionary with coordinates of each station name.
    """
    if not isinstance(path, str):
        raise TypeError("The 'path' argument is not of type str!")

    first = True
    coords_dict = {}
    dataframe0 = None
    for infile in sorted(glob.glob(f"{path}/*nas")):
        dataframe, xyz = read_nas(
            infile,
            verbose=verbose
        )
        if first:
            dataframe0 = dataframe
            coords_dict.update(xyz)
            first = False
        else:
            dataframe0 = pd.concat(
                [dataframe0, dataframe],
                axis=1
            )
            coords_dict.update(xyz)

    # If dataframe0 is undefined, the loop didn't execute... so throw error
    if dataframe0 is None:
        raise ValueError(f"Could not find data in {path}!")

    dataframe = dataframe0.groupby(
        dataframe0.columns,
        axis=1
    ).max()

    return dataframe, coords_dict


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

    # Convert data to ppbv if necessary
    with xr.set_options(keep_attrs=True):
        dataarray = dataset[varname]
        units = dataarray.attrs["units"]
        if "mol mol-1" in units or "v/v" in units:
            dataarray.values *= 1.0e9
            dataarray.attrs["units"] = "ppbv"

    return dataarray


def find_times(
        dataframe,
        start_time
):
    """
    Convert timestamps in nasa ames data files to python datetime
    objects  Set DataFrame index to the new datetime array

    Args:
    ----------
    dataframe : pandas DataFrame
        DataFrame with O3 values from GAW site

    start_time : str
        Reference start time for timestamp taken from nasa ames file

    Returns
    ------
    dataframe: pandas DataFrame
        O3 in ppbV with datetime index

    qcflag : pandas Dataframe
        QC flag with datetime index
    """
    end_time = dataframe[dataframe.columns[1]]
    time_x = []

    for index in range(len(end_time)):
        time_x.append(start_time + timedelta(days=end_time.values[index]))

    dataframe.index = time_x
    qcflag = dataframe[dataframe.columns[-1]]
    dataframe = dataframe[dataframe.columns[2]]

    return dataframe, qcflag


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

    gc_level_alts_m: pandas DataFrame
        Altitudes of GEOS-Chem levels in meters

    lon_value : float
        GAW site longitude

    lat_value : float
        GAW site latitude

    alt_value : float
        GAW site altitude

    Returns:
    --------
    x_idx, y_idx, z_idx
        GEOS-Chem grid box indices for the single gridbox
        closest to GAW site specifications
    """
    x_idx=(
        np.abs(
            gc_data.lon - float(lon_value)
        )
    ).argmin()

    y_idx=(
        np.abs(
            gc_data.lat - float(lat_value)
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
        station_name,
        coords_dict,
        gc_level_alts_m,
        obs_dataframe,
        ref_dataarray,
        dev_dataarray,
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
    station_name : str
        Name of the observation station site.

    coords_dict : dict
        Coordinates (lon, lat, alt) for each observation station site.

    gc_level_alts_m : pandas DataFrame
        Metadata pertaining to GEOS-Chem vertical levels

    obs_dataframe : pandas DataFrame
        Observations at each station site.

    ref_dataarray, dev_dataarray : xarray DataArray
        Data from the Ref and Dev model versions.

    ref_label, dev_label: str
        Labels describing the Ref and Dev datasets (e.g. version numbers)

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

    plot_title : str
        Plot title string for the given observation station site.

    yaxis_label : str
        Label for the Y-axis (e.g. species name).
    """

    # Round the station lon, lat, alt  values
    sta_lon = round(coords_dict[station_name]['lon'], 2)
    sta_lat = round(coords_dict[station_name]['lat'], 2)
    sta_alt = round(coords_dict[station_name]['alt'], 1)

    # Find nearest model box and get model data at that box
    x_idx, y_idx, z_idx = find_nearest_3d(
        ref_dataarray,
        gc_level_alts_m,
        lon_value=sta_lon,
        lat_value=sta_lat,
        alt_value=sta_alt
    )
    ref_dataframe = ref_dataarray.isel(
        lon=x_idx,
        lat=y_idx,
        lev=z_idx
    ).to_dataframe()
    dev_dataframe = dev_dataarray.isel(
        lon=x_idx,
        lat=y_idx,
        lev=z_idx
    ).to_dataframe()

    # Take the mean of observations for plotting
    obs_dataframe = obs_dataframe.resample('M').mean()

    # Top-of-plot title for the observation station site
    plot_title = "{station_name} ({sta_lat}$^\\circ$N,{sta_lon}$^\\circ$E)"

    # Y-axis label (i.e. species name)
    yaxis_label = varname.split("_")[1]

    return obs_dataframe, ref_dataframe[varname], dev_dataframe[varname], \
        plot_title, yaxis_label


def plot_single_station(
        fig,
        n_rows,
        n_cols,
        n_station,
        station_name,
        plot_title,
        yaxis_label,
        obs_dataframe,
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
        Matplotlib Figure object containing the plot

    n_rows, n_cols: int
        Number of rows and columns in the plot.

    n_station : int
        Number of the observation station.

    station_name : str
        Name of the observation station site.

    plot_title, yaxis-label : str
        Top of plot title and y-axis label.

    obs_dataframe : pandas DataFrame
        Observational data.

    ref_series, dev_series : pandas Series
        GEOS-Chem data at closest grid box to the observation
        station site for the Ref and Dev model versions.

    ref_label, dev_label : str
        Descriptive labels (e.g. version numbers) for the
        GEOS-Chem Ref and Dev model versions.
    """
    if not isinstance(fig, Figure):
        msg = "The 'fig' argument is not of type matplotlib.figure.Figure!"
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
    axes_subplot = fig.add_subplot(n_rows, n_cols, n_station+1)

    # Set title for top of each frame
    axes_subplot.set_title(
        f"{plot_title}",
        weight='bold',
        fontsize=7
        )

    ## Plot observational data
    axes_subplot.plot(
        obs_dataframe.index,
        obs_dataframe[station_name],
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
    if n_station == 0 or n_station % n_cols == 0:
        axes_subplot.set_ylabel(
            yaxis_label,
            fontsize=8
        )

    # Set X-axis and Y-axis ticks
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
        75
    )
    axes_subplot.set_yticks(
        [0, 25, 50, 75]
    )
    axes_subplot.tick_params(
        axis='both',
        which='major',
        labelsize=6
    )


def plot_models_vs_obs(
        ref_dataarray,
        ref_label,
        dev_dataarray,
        dev_label,
        gc_level_alts_m,
        obs_dataframe,
        coords_dict,
        varname="SpeciesConcVV_O3",
        dst="./benchmark",
        verbose=False
):
    """
    Plots models vs. observations using a 3 rows x 3 column layout.

    Args:
    -----
    ref_dataarray, dev_dataarray : xarray DataArray
        Data from the Ref and Dev model versions.

    ref_label, dev_label: str
        Labels describing the Ref and Dev datasets (e.g. version numbers)

    gc_level_alts_m : pandas DataFrame
        Metadata pertaining to GEOS-Chem vertical levels

    obs_dataframe : pandas DataFrame
        Observations at each station site.

    coords_dict : dict
        Coordinates (lon, lat, alt) for each observation station site.

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
    # Error checks
    if not isinstance(ref_dataarray, xr.DataArray):
        msg = "The 'ref_dataset' argument is not of type xarray DataArray!"
        raise TypeError(msg)
    if not isinstance(ref_label, str):
        msg = "The 'ref_label' argument is not of type str!"
        raise TypeError(msg)
    if not isinstance(dev_dataarray, xr.DataArray):
        msg = "The 'ref_dataset' argument is not of type xarray DataArray!"
        raise TypeError(msg)
    if not isinstance(dev_label, str):
        msg = "The 'dev_label' argument is not of type str!"
        raise TypeError(msg)
    if not isinstance(obs_dataframe, pd.DataFrame):
        msg = "The 'obs_dataframe' argument is not of type pandas DataFrame!"
        raise TypeError(msg)
    if not isinstance(gc_level_alts_m, pd.Series):
        msg = "The 'gc_level_alts_m' argument is not of type pandas Series!"
        raise TypeError(msg)

    # Figure setup
    plt.style.use('seaborn-darkgrid')
    n_cols = 3
    n_rows = ceil(len(obs_dataframe.columns) / 3)

    # Compute nuimber of pages that are necessary
    #n_pages = ncols * nrows
    #if n_pages % (nrows * ncols) > 0:
    #    npages += 1
    #
    #
    #for page in range(npages):
    #    plot_single_page(

    fig = plt.figure(figsize=(8.27, 29))

    # Loop over the number of columns in the obs_data dataframe
    for n_station, station_name, in enumerate(obs_dataframe.columns):
        if verbose:
            print(f"Plotting data for station: {station_name}")


        # Find the model Ref & Dev data closest to the observational
        # station site.  Also take average of observations for plotting.
        obs_dataframe, ref_series, dev_series, plot_title, yaxis_label = \
        prepare_data_for_plot(
            station_name,
            coords_dict,
            gc_level_alts_m,         # pandas.DataFrame
            obs_dataframe,           # pandas.DataFrame
            ref_dataarray,           # xarray.DataArray
            dev_dataarray,           # xarray.DataArray
            varname=varname
        )

        # Plot models vs. observation for a single station site
        plot_single_station(
            fig,                     # matplotlib.figure.Figure
            n_rows,
            n_cols,
            n_station,
            station_name,
            plot_title,
            yaxis_label,
            obs_dataframe,           # pandas.Dataframe
            ref_series,              # pandas.Series
            ref_label,
            dev_series,              # pandas.Series
            dev_label
        )

    ## Adjust spacing and save to pdf file
    plt.tight_layout()
    plt.subplots_adjust(
        hspace=0.4,
        top=0.9
    )
    plt.legend(
        ncol=3,
        bbox_to_anchor=(1.27, 26.8)
    )
    fig.savefig(
        f"{dst}/models_vs_obs.surface.{varname}.pdf",
        bbox_inches='tight'
    )
    plt.close()


def make_benchmark_models_vs_obs_plots(
        ref_filepaths,
        ref_label,
        dev_filepaths,
        dev_label,
        obs_filepaths,
        varname="SpeciesConcVV_O3",
        dst="./benchmark",
        verbose=False,
        overwrite=False
):
    """
    Driver rtou
    """

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
    obs_dataframe, coords_dict = read_observational_data(
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
        ref_dataarray,
        ref_label,
        dev_dataarray,
        dev_label,
        gc_level_alts_m,
        obs_dataframe,
        coords_dict,
        varname=varname,
        dst=dst,
        verbose=verbose
    )
