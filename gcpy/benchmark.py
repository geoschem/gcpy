""" Specific utilities re-factored from the benchmarking utilities. """

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages

from .plot import WhGrYlRd, add_latlon_ticks

# change default fontsize (globally)
# http://matplotlib.org/users/customizing.html
mpl.rcParams['font.size'] = 15

cmap_abs = WhGrYlRd  # for absolute magnitude
cmap_diff = 'RdBu_r'  # for difference plot


def plot_layer(dr, ax, fig, title=''):
    '''plot 2D DataArray as a lat-lon layer
    '''

    # imshow() is 6x faster than pcolormesh(),
    # but with some limitations:
    # - only works with PlateCarree projection
    # - the left map boundary can't be smaller than -180,
    #   so the leftmost box (-182.5 for 4x5 grid) is slightly out of the map

    im = dr.plot.imshow(ax=ax, cmap=cmap_abs, transform=ccrs.PlateCarree(),
                        add_colorbar=False)

    try:
        unit = dr.attrs['units']
    except:
        unit = ''

    # can also pass cbar_kwargs to dr.plot() to add colorbar
    # but it is easier to tweak colorbar afterwards
    cb = fig.colorbar(im, ax=ax, shrink=0.6, orientation='horizontal', pad=0.1)
    cb.set_label(unit)

    # xarray automatically sets a title which might contain dimension data
    # surpress it or use user-provided title
    ax.set_title(title)

    ax.coastlines()
    add_latlon_ticks(ax)  # add ticks and gridlines


def pdf_two_layers(ds1, ds2, filename):
    '''plot all variables in a 2D DataSet (lat-lon layer)
    '''

    # plot all the variables in ds1
    # assume ds2 also has those variables

    varname_list = list(ds1.data_vars.keys())

    n_var = len(varname_list)
    print('Benchmarking {} variables'.format(n_var))

    n_row = 3  # how many rows per page
    n_page = (n_var-1) // n_row + 1  # how many pages

    print('generating a {}-page pdf'.format(n_page))

    pdf = PdfPages(filename)

    for ipage in range(n_page):
        print(ipage, end=' ')
        fig, axes = plt.subplots(n_row, 2, figsize=[16, 16],
                                 subplot_kw={'projection': ccrs.PlateCarree()})

        # a list of variables names with only 3 elements
        sub_varname_list = varname_list[n_row*ipage:n_row*(ipage+1)]

        for i, varname in enumerate(sub_varname_list):
            for j, ds in enumerate([ds1, ds2]):
                plot_layer(ds[varname], axes[i][j], fig)

            axes[i][0].set_title(varname+'; surface', fontsize=20)
            axes[i][1].set_title(varname+'; 500 hpa', fontsize=20)

        pdf.savefig(fig)
        plt.close(fig)  # don't show in notebook!
    pdf.close()  # close it to save the pdf
    print('done!')
