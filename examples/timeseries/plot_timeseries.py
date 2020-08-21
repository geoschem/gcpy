#!/usr/bin/env python
'''
Example of plotting timeseries data from GEOS-Chem and saving
the output to a PDF file.  You can modify this for your particular
diagnostic output.  This also contains a good overview of

This example script creates a PDF file with 2 pages.

   Page 1:
   -------
       O3 from the first model layer (from the "SpeciesConc"
       diagnostic collection is) plotted in blue.

       O3 at 10 meter height (from the "SpeciesConc_10m"
       diagnostic collection) is plotted in red.

   Page 2:
   -------
       HNO3 from the first model layer (from the SpeciesConc
       diagnostic collection is) plotted in blue.

       HNO3 at 10 meter height (from the SpeciesConc_10m
       diagnostic collection) is plotted in red.

You can of course modify this for your own particular applications.

Author:
-------
Bob Yantosca
yantosca@seas.harvard.edu
23 Aug 2019
'''

# Imports
import gcpy.constants as gcon
import os
import numpy as np
import matplotlib.dates as mdates
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import xarray as xr
import warnings

# Tell matplotlib not to look for an X-window, as we are plotting to
# a file and not to the screen.  This will avoid some warning messages.
os.environ['QT_QPA_PLATFORM'] = 'offscreen'

# Suppress harmless run-time warnings (mostly about underflow in division)
warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', category=UserWarning)


def find_files_in_dir(path, substrs):
    '''
    Returns a list of all files in a directory that match one or more
    substrings.

    Args:
    -----
        path : str
            Path to the directory in which to search for files.

        substrs : list of str
            List of substrings used in the search for files.

    Returns:
    --------
        file_list : list of str
            List of files in the directory (specified by path)
            that match all substrings (specified in substrs).
    '''
    
    # Initialize
    file_list = []

    # Walk through the given data directory.  Then for each file found,
    # add it to file_list if it matches text in search_list.
    for root, directory, files in os.walk(path):
        for f in files:
            for s in substrs:
                if s in f:
                    file_list.append(os.path.join(root, f))

    # Return an alphabetically sorted list of files
    file_list.sort()
    return file_list


def find_value_index(seq, val):
    '''
    Finds the index of a numpy array that is close to a value.

    Args:
    -----
        seq : numpy ndarray
            An array of numeric values.

        val : number
            The value to search for in seq.

    Returns:
    --------
        result : integer
            The index of seq that has a value closest to val.

    Remarks:
    --------
    This algorithm was found on this page:
    https://stackoverflow.com/questions/48900977/find-all-indexes-of-a-numpy-array-closest-to-a-value
    '''
    r = np.where(np.diff(np.sign(seq - val)) != 0)
    idx = r + (val - seq[r]) / (seq[r + np.ones_like(r)] - seq[r])
    idx = np.append(idx, np.where(seq == val))
    idx = np.sort(idx)
    result = np.round(idx)

    # NOTE: xarray needs integer values, so convert here!
    return int(result[0])


def read_geoschem_data(path, collections):
    '''
    Returns an xarray Dataset containing timeseries data.

    Args:
    -----
        path : str
            Directory path where GEOS-Chem diagnostic output
            files may be found.

        collections: list of str
            List of GEOS-Chem collections.  Files for these 
            collections will be read into the xarray Dataset.
            
    Returns:
    --------
        ds : xarray Dataset
            A Dataset object containing the GEOS-Chem diagnostic
            output corresponding to the collections that were
            specified.
    '''

    # Get a list of variables that GCPy should not read.
    # These are mostly variables introduced into GCHP with the MAPL v1.0.0
    # update.  These variables contain either repeated or non-standard
    # dimensions that can cause problems in xarray when combining datasets.
    skip_vars = gcon.skip_these_vars
    
    # Find all files in the given 
    file_list = find_files_in_dir(path, collections) 

    # Return a single xarray Dataset containing data from all files
    # NOTE: Need to add combine="nested" for xarray 0.15 and higher
    v = xr.__version__.split(".")
    if int(v[0]) == 0 and int(v[1]) >= 15: 
        return xr.open_mfdataset(file_list,
                                 drop_variables=skip_vars,
                                 combine="nested",
                                 concat_dim=None)
    else:
        return xr.open_mfdataset(file_list,
                                 drop_variables=skip_vars)


def plot_timeseries_data(ds, site_coords):
    '''
    Plots a timseries of data at a given (lat,lon) location.

    Args:
    -----
        ds : xarray Dataset
            Dataset containing GEOS-Chem timeseries data.

        site_coords : tuple
            Contains the coordinate (lat, lon) of a site location
            at which the timeseries data will be plotted.
    '''

    # ----------------------------------------------------------------------
    # Get the GEOS-Chem data for O3 and HNO3 corresponding to the
    # location of the observational station.  We will save these into
    # xarray DataArray objects, which we'll need for plotting.
    #
    # YOU CAN EDIT THIS FOR YOUR OWN PARTICULAR APPLICATION!
    # ----------------------------------------------------------------------

    # Find the indices corresponding to the site lon and lat
    lat_idx = find_value_index(ds.lat.values, site_coords[0])
    lon_idx = find_value_index(ds.lon.values, site_coords[1])

    # Save O3 from the first level (~60m height) (ppb) into a DataArray
    O3_L1 = ds['SpeciesConc_O3'].isel(lon=lon_idx, lat=lat_idx, lev=0)
    O3_L1 *= 1.0e9
    O3_L1.attrs['units'] = 'ppbv'

    # Save O3 @ 10m height into a DataArray
    O3_10m = ds['SpeciesConc10m_O3'].isel(lon=lon_idx, lat=lat_idx)
    O3_10m *= 1.0e9
    O3_10m.attrs['units'] = 'ppbv'

    # Save HNO3 from the first level (~60m height) into a DataArray
    HNO3_L1 = ds['SpeciesConc_HNO3'].isel(lon=lon_idx, lat=lat_idx, lev=0)
    HNO3_L1 *= 1.0e9
    HNO3_L1.attrs['units'] = 'ppbv'

    # Save HNO3 @ 10m height into a DataArray
    HNO3_10m = ds['SpeciesConc10m_HNO3'].isel(lon=lon_idx, lat=lat_idx)
    HNO3_10m *= 1.0e9
    HNO3_10m.attrs['units'] = 'ppbv'

    # ----------------------------------------------------------------------
    # Create a PDF file of the plots
    # ----------------------------------------------------------------------

    # Get min & max days of the plot span (for setting the X-axis range).
    # To better center the plot, add a cushion of 12 hours on either end.
    time = ds['time'].values
    datemin = np.datetime64(time[0]) - np.timedelta64(12, 'h')
    datemax = np.datetime64(time[-1]) + np.timedelta64(12, 'h')

    # Define a PDF object so that we can save the plots to PDF
    pdf = PdfPages('O3_and_HNO3.pdf')

    # Loop over number of desired pages (in this case, 2)
    for i in range(0, 2):

        # Create a new figure: 1 plot per page, 2x as wide as high
        figs, ax0 = plt.subplots(1, 1, figsize=[12, 6])

        # -----------------------------
        # Plot O3 on the first page
        # -----------------------------
        if i == 0:

            # 1st model level
            O3_L1.plot.line(ax=ax0, x='time', color='blue',
                            marker='o', label='O3 from 1st model level',
                            linestyle='-')

            # 10 mheight
            O3_10m.plot.line(ax=ax0, x='time', color='red',
                             marker='x', label='O3 at 10m height',
                             linestyle='-')

            # Set title (has to be after the line plots are drawn)
            ax0.set_title('O3 from the 1st model level and at 10m height')

            # Set Y-axis minor tick marks at every 2 ppb (5 intervals)
            ax0.yaxis.set_minor_locator(mticker.AutoMinorLocator(5))

            # Set y-axis title
            ax0.set_ylabel('O3 (ppbv)')

        # -----------------------------
        # Plot HNO3 on the second page
        # -----------------------------
        if i == 1:

            # 1st model level
            HNO3_L1.plot.line(ax=ax0, x='time', color='blue',
                              marker='o', label='HNO3 from 1st model level',
                              linestyle='-')

            # 10m height
            HNO3_10m.plot.line(ax=ax0, x='time', color='red',
                               marker='x', label='HNO3 at 10m height',
                               linestyle='-')

            # Set title (has to be after the line plots are drawn
            ax0.set_title('HNO3 from the 1st model level and at 10m height')

            # Set Y-axis minor tick marks at every 0.05 ppb (4 intervals)
            ax0.yaxis.set_minor_locator(mticker.AutoMinorLocator(4))

            # Set y-axis title
            ax0.set_ylabel('HNO3 (ppbv)')

        # -----------------------------
        # Set general plot parameters
        # -----------------------------

        # Add the plot legend
        ax0.legend()

        # Set the X-axis range
        ax0.set_xlim(datemin, datemax)

        # Set the X-axis major tickmarks
        locator = mdates.DayLocator()
        formatter = mdates.DateFormatter('%d')
        ax0.xaxis.set_major_locator(locator)
        ax0.xaxis.set_major_formatter(formatter)

        # Set X-axis minor tick marks at noon of each day
        # (i.e. split up the major interval into 2 bins)
        ax0.xaxis.set_minor_locator(mticker.AutoMinorLocator(2))

        # Don't rotate the X-axis jtick labels
        ax0.xaxis.set_tick_params(rotation=0)

        # Center the X-axis tick labels
        for tick in ax0.xaxis.get_major_ticks():
            tick.label1.set_horizontalalignment('center')

        # Set X-axis and Y-axis labels
        ax0.set_xlabel('Day of July (and August) 2016')

        # -----------------------------
        # Save this page to PDF
        # -----------------------------
        pdf.savefig(figs)
        plt.close(figs)

    # ----------------------------------------------------------------------
    # Save the PDF file to disk
    # ----------------------------------------------------------------------
    pdf.close()


def main():
    '''
    Main program.
    '''
    # Path where the data files live
    # (YOU MUST EDIT THIS FOR YUR OWN PARTICULAR APPLICATION!)
    path_to_data = '/path/to/GEOS-Chem/diagnostic/data/files'

    # Get a list of files in the ConcAboveSfc and SpeciesConc collections
    # (YOU CAN EDIT THIS FOR YOUR OWN PARTICULAR APPLICATION!)
    collections = ['ConcAboveSfc', 'SpeciesConc']
    
    # Read GEOS-Chem data into an xarray Dataset
    ds = read_geoschem_data(path_to_data, collections)

    # Plot timeseries data at Centerville, AL (32.94N, 87.18W)
    # (YOU CAN EDIT THIS FOR YOUR OWN PARTICULAR APPLICATION!)
    site_coords = (32.94, -87.18)
    plot_timeseries_data(ds, site_coords)


if __name__ == "__main__":
    main()
