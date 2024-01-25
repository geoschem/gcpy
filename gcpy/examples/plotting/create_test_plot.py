#!/usr/bin/env python

"""
Creates a dummy xarray DataArray object for testing the single_panel
plotting routine.  This script can be a useful check to see whether the
GCPy mamba/conda environment has been successfully installed.
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import gcpy


def main():
    """
    Main program: Creates a dummy plot
    """
    
    # Get X and Y coordinate arrays
    n_lat = 181
    n_lon = 361
    lat_vals = np.linspace(-90, 90, n_lat)
    lon_vals = np.linspace(-180, 180, n_lon)
    
    # Create a dummy numpy array for plotting
    # Populate it with the distance from the center of the plot
    data = np.zeros([n_lat, n_lon])
    for lat in range(n_lat):
        for lon in range(n_lon):
            data[lat, lon] = np.sqrt(lon_vals[lon]**2 + lat_vals[lat]**2)
            
    # Create a Dataarray from the numpy array
    darr = xr.DataArray(
        data=data,
        dims=["lat", "lon"],
        coords=dict(
            lat=("lat", lat_vals),
            lon=("lon", lon_vals),
        ),
    )
    
    # Create a plot
    gcpy.plot.single_panel(
        darr,
        plot_type="single_level",
        title="Test plot to check GCPy import",
        extent=[-180, 180, -90, 90],
        comap=plt.get_cmap("RdBu_r")
    )
    plt.show()


if __name__ == '__main__':
    main()
