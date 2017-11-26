""" Specific utilities re-factored from the benchmarking utilities. """

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
from cartopy.mpl.geoaxes import GeoAxes  # for assertion

from .plot import WhGrYlRd, add_latlon_ticks

# change default fontsize (globally)
# http://matplotlib.org/users/customizing.html
mpl.rcParams['font.size'] = 12
mpl.rcParams['axes.titlesize'] = 20

cmap_abs = WhGrYlRd  # for absolute magnitude
cmap_diff = 'RdBu_r'  # for difference plot


def plot_layer(dr, ax, title='', unit='', diff=False):
    '''Plot 2D DataArray as a lat-lon layer

    Parameters
    ----------
    dr : xarray.DataArray
        Dimension should be [lat, lon]

    ax : Cartopy GeoAxes object with PlateCarree projection
        Axis on which to plot this figure

    title : string, optional
        Title of the figure

    unit : string, optional
        Unit shown near the colorbar

    diff : Boolean, optional
        Switch to the color scale for difference plot
    '''

    # imshow() is 6x faster than pcolormesh(),
    # but with some limitations:
    # - only works with PlateCarree projection
    # - the left map boundary can't be smaller than -180,
    #   so the leftmost box (-182.5 for 4x5 grid) is slightly out of the map

    assert isinstance(ax, GeoAxes), (
           "Input axis must be cartopy GeoAxes! "
           "Can be created by: \n"
           "plt.axes(projection=ccrs.PlateCarree()) \n or \n"
           "plt.subplots(n, m, subplot_kw={'projection': ccrs.PlateCarree()})"
           )
    assert ax.projection == ccrs.PlateCarree(), (
           'must use PlateCarree projection'
           )

    fig = ax.figure  # get parent figure

    if diff:
        vmax = np.max(np.abs(dr.values))
        vmin = -vmax
        cmap = cmap_diff
    else:
        vmax = np.max(dr.values)
        vmin = 0
        cmap = cmap_abs

    im = dr.plot.imshow(ax=ax, vmin=vmin, vmax=vmax, cmap=cmap,
                        transform=ccrs.PlateCarree(),
                        add_colorbar=False)

    # can also pass cbar_kwargs to dr.plot() to add colorbar
    # but it is easier to tweak colorbar afterwards
    cb = fig.colorbar(im, ax=ax, shrink=0.6, orientation='horizontal', pad=0.1)
    cb.set_label(unit)

    # xarray automatically sets a title which might contain dimension data
    # surpress it or use user-provided title
    ax.set_title(title)

    ax.coastlines()
    add_latlon_ticks(ax)  # add ticks and gridlines


def plot_zonal(dr, ax, title='', unit='', diff=False):
    '''Plot 2D DataArray as a zonal profile

    Parameters
    ----------
    dr : xarray.DataArray
        dimension should be [lev, lat]

    ax : matplotlib axes object
        Axis on which to plot this figure

    title : string, optional
        Title of the figure

    unit : string, optional
        Unit shown near the colorbar

    diff : Boolean, optional
        Switch to the color scale for difference plot
    '''

    xtick_positions = np.array([-90, -60, -30, 0, 30, 60, 90])
    xticklabels = ['90$\degree$S',
                   '60$\degree$S',
                   '30$\degree$S',
                   '0$\degree$',
                   '30$\degree$N',
                   '60$\degree$N',
                   '90$\degree$N'
                   ]

    fig = ax.figure  # get parent figure

    if diff:
        vmax = np.max(np.abs(dr.values))
        vmin = -vmax
        cmap = cmap_diff
    else:
        vmax = np.max(dr.values)
        vmin = 0
        cmap = cmap_abs

    im = dr.plot.imshow(ax=ax, vmin=vmin, vmax=vmax, cmap=cmap,
                        add_colorbar=False)

    ax.set_aspect(1.5)  # the ratio of x-unit/y-unit in screen-space

    ax.set_xticks(xtick_positions)
    ax.set_xticklabels(xticklabels)
    ax.set_xlabel('')
    ax.set_ylabel('Level')

    # can also pass cbar_kwargs to dr.plot() to add colorbar
    # but it is easier to tweak colorbar afterwards
    cb = fig.colorbar(im, ax=ax, shrink=0.6, orientation='horizontal', pad=0.1)
    cb.set_label(unit)

    ax.set_title(title)


def make_pdf(ds1, ds2, filename, on_map=True, diff=False,
             title1='DataSet 1', title2='DataSet 2', unit=''):
    '''Plot all variables in two 2D DataSets, and create a pdf.

    ds1 : xarray.DataSet
        shown on the left column

    ds2 : xarray.DataSet
        shown on the right column

    filename : string
        Name of the pdf file

    on_map : Boolean, optional
        If True (default), use plot_layer() to plot
        If False, use plot_zonal() to plot

    diff : Boolean, optional
        Switch to the color scale for difference plot

    title1, title2 : string, optional
        Title for each DataSet

    unit : string, optional
        Unit shown near the colorbar
    '''

    if on_map:
        plot_func = plot_layer
        subplot_kw = {'projection': ccrs.PlateCarree()}
    else:
        plot_func = plot_zonal
        subplot_kw = None

    # plot all the variables in ds1
    # assume ds2 also has those variables

    varname_list = list(ds1.data_vars.keys())

    n_var = len(varname_list)
    print('Benchmarking {} variables'.format(n_var))

    n_row = 3  # how many rows per page
    n_page = (n_var-1) // n_row + 1  # how many pages

    print('generating a {}-page pdf'.format(n_page))
    print('Page: ', end='')

    pdf = PdfPages(filename)

    for ipage in range(n_page):
        print(ipage, end=' ')
        fig, axes = plt.subplots(n_row, 2, figsize=[16, 16],
                                 subplot_kw=subplot_kw)

        # a list of variables names with only 3 elements
        sub_varname_list = varname_list[n_row*ipage:n_row*(ipage+1)]

        for i, varname in enumerate(sub_varname_list):
            for j, ds in enumerate([ds1, ds2]):
                plot_func(ds[varname], axes[i][j], unit=unit, diff=diff)

            axes[i][0].set_title(varname+'; '+title1)
            axes[i][1].set_title(varname+'; '+title2)

        pdf.savefig(fig)
        plt.close(fig)  # don't show in notebook!
    pdf.close()  # close it to save the pdf
    print('done!')
