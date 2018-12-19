""" Specific utilities re-factored from the benchmarking utilities. """

import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

import matplotlib as mpl
from matplotlib import ticker
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap

from cartopy import crs
from cartopy.mpl.geoaxes import GeoAxes  # for assertion

from PyPDF2 import PdfFileWriter, PdfFileReader

from .plot import WhGrYlRd, add_latlon_ticks
from .grid.horiz import make_grid_LL, make_grid_CS
from .grid.regrid import make_regridder_C2L

# change default fontsize (globally)
# http://matplotlib.org/users/customizing.html
#mpl.rcParams['font.size'] = 12
#mpl.rcParams['axes.titlesize'] = 20

cmap_abs = WhGrYlRd  # for plotting absolute magnitude
cmap_diff = 'RdBu_r'  # for plotting difference


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

    # imshow() is 6x faster than pcolormesh(), but with some limitations:
    # - only works with PlateCarree projection
    # - the left map boundary can't be smaller than -180,
    #   so the leftmost box (-182.5 for 4x5 grid) is slightly out of the map
    im = dr.plot.imshow(ax=ax, vmin=vmin, vmax=vmax, cmap=cmap,
                        transform=ccrs.PlateCarree(),
                        add_colorbar=False)

    # can also pass cbar_kwargs to dr.plot() to add colorbar,
    # but it is easier to tweak colorbar afterwards
    cb = fig.colorbar(im, ax=ax, shrink=0.6, orientation='horizontal', pad=0.1)
    cb.set_label(unit)

    # xarray automatically sets a title which might contain dimension info.
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

    # assume global field from 90S to 90N
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

    # this code block largely duplicates plot_layer()
    # TODO: remove duplication
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

    # the ratio of x-unit/y-unit in screen-space
    # 'auto' fills current figure with data without changing the figrue size
    ax.set_aspect('auto')

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

    # get a list of all variable names in ds1
    # assume ds2 also has those variables
    varname_list = list(ds1.data_vars.keys())

    n_var = len(varname_list)
    print('Benchmarking {} variables'.format(n_var))

    n_row = 3  # how many rows per page. TODO: should this be input argument?
    n_page = (n_var-1) // n_row + 1  # how many pages

    print('generating a {}-page pdf'.format(n_page))
    print('Page: ', end='')

    pdf = PdfPages(filename)

    for ipage in range(n_page):
        print(ipage, end=' ')
        fig, axes = plt.subplots(n_row, 2, figsize=[16, 16],
                                 subplot_kw=subplot_kw)

        # a list of 3 (n_row) variables names
        sub_varname_list = varname_list[n_row*ipage:n_row*(ipage+1)]

        for i, varname in enumerate(sub_varname_list):

            # Get min/max for both datasets to have same colorbar (ewl)
            #vmin = min([ds1[varname].data.min(), ds2[varname].data.min()])
            #vmax = max([ds1[varname].data.max(), ds2[varname].data.max()])

            for j, ds in enumerate([ds1, ds2]):
                if on_map:
                    plot_func(ds[varname], axes[i][j], unit=ds[varname].units, 
                              diff=diff)
                else:
                    # For now, assume zonal mean if plotting zonal (ewl)
                    plot_func(ds[varname].mean(axis=2), axes[i][j], 
                              unit=ds[varname].units, diff=diff)


            # TODO: tweak varname, e.g. Trim "TRC_O3" to "O3"
            axes[i][0].set_title(varname+'; '+title1)
            axes[i][1].set_title(varname+'; '+title2)

            # TODO: unit conversion according to data range of each variable,
            # e.g. mol/mol -> ppmv, ppbv, etc...

        pdf.savefig(fig)
        plt.close(fig)  # don't show in notebook!
    pdf.close()  # close it to save the pdf
    print('done!')

# Add docstrings later. Use this function for benchmarking or general comparisons.
'''Notes: The two GC files must be at the same grid resolution. You can use this function to plot interactively or to generate a multi-page pdf of plots. It take#s a list of variable names to plot for a single collection only. You can plot for any level and any time slice in the file. By default the colorbars for the plots will have the same range, but you can turn this feature off. Also by default the colorbar of the fractional difference between the model outputs will be limited to +/-2, but you can change this as well via the passed parameters.'''
def compare_single_level(refdata, refstr, devdata, devstr, llres, varlist=None, ilev=0, itime=0, savepdf=False, 
                         pdfname='map.pdf', match_cbar=True, normalize_by_area=False, check_units=True):
    
    # If no varlist is passed, plot all (surface only for 3D)
    if varlist == None:
        [varlist, commonvars2D, commonvars3D] = compare_varnames(refdata, devdata)
        print('Plotting all common variables (surface only if 3D)')
    n_var = len(varlist)

    # Get lat-lon grid
    llgrid = make_grid_LL(llres)

    # Get lat/lon extents
    [minlon, maxlon] = [min(llgrid['lon_b']), max(llgrid['lon_b'])]
    [minlat, maxlat] = [min(llgrid['lat_b']), max(llgrid['lat_b'])]
    
    # Create pdf (if saving)
    if savepdf:
        print('\nCreating {} for {} variables'.format(pdfname,n_var))
        pdf = PdfPages(pdfname)

    # Loop over variables
    for ivar in range(n_var):
        if savepdf: print('{} '.format(ivar), end='')
        varname = varlist[ivar]
        
        # Do some checks: dimensions and units
        varndim_ref = refdata[varname].ndim
        varndim_dev = devdata[varname].ndim      
        
        units_ref = refdata[varname].units.strip()
        units_dev = devdata[varname].units.strip()
        if check_units: 
            assert units_ref == units_dev, 'Units do not match for {}!'.format(varname)
            
        # if normalizing by area, adjust units to be per m2, and adjust title string
        units = units_ref
        subtitle_extra = ''
                    
        # Slice the data for ref, allowing for possibility of no time dimension
        vdims = refdata[varname].dims
        if 'time' in vdims and 'lev' in vdims: 
            ds_ref = refdata[varname].isel(time=itime,lev=ilev)
        elif 'time' in vdims and 'lev' not in vdims: 
            ds_ref = refdata[varname].isel(time=itime)
        elif 'time' not in vdims and 'lev' in vdims:
            ds_ref = refdata[varname].isel(lev=ilev)   
        else:
            ds_ref = refdata[varname]
 
        # Slice the data for dev, allowing for possibility of no time dimension
        vdims = devdata[varname].dims
        if 'time' in vdims and 'lev' in vdims: 
            ds_dev = devdata[varname].isel(time=itime,lev=ilev)
        elif 'time' in vdims and 'lev' not in vdims:  
            ds_dev = devdata[varname].isel(time=itime)
        elif 'time' not in vdims and 'lev' in vdims:
            ds_dev = devdata[varname].isel(lev=ilev) 
        else:
            ds_dev = devdata[varname]
            
        # if normalizing by area, transform on the native grid and adjust units and subtitle string
        exclude_list = ['WetLossConvFrac','Prod_','Loss_']
        if normalize_by_area and not any(s in varname for s in exclude_list):
            ds_ref.values = ds_ref.values / refdata['AREAM2'].values
            ds_dev.values = ds_dev.values / devdata['AREAM2'].values
            units = '{} m-2'.format(units)
            subtitle_extra = ', Normalized by Area'
            
        # Regrid the slices (skip for gchp vs gchp for now)
        vmin_ref = ds_ref.min()
        vmin_dev = ds_dev.min()
        vmin_cmn = np.min([vmin_ref, vmin_dev])
        vmax_ref = ds_ref.max()
        vmax_dev = ds_dev.max()
        vmax_cmn = np.max([vmax_ref, vmax_dev])
        if match_cbar: [vmin, vmax] = [vmin_cmn, vmax_cmn]
        
        # Create 3x2 figure
        figs, ((ax0, ax1), (ax2, ax3), (ax4, ax5)) = plt.subplots(3, 2, figsize=[12,14], 
                                                      subplot_kw={'projection': crs.PlateCarree()})
        # Give the figure a title
        offset = 0.96
        fontsize=25
        
        if 'lev' in refdata[varname].dims and 'lev' in devdata[varname].dims:
            if ilev == 0: levstr = 'Surface'
            elif ilev == 22: levstr = '500 hPa'
            else: levstr = 'Level ' +  str(ilev-1)
            figs.suptitle('{}, {}'.format(varname,levstr), fontsize=fontsize, y=offset)
        elif 'lat' in refdata[varname].dims and 'lat' in devdata[varname] and 'lon' in refdata[varname].dims and 'lon' in devdata[varname]: 
            figs.suptitle('{}'.format(varname), fontsize=fontsize, y=offset)
        else:
            print('Incorrect dimensions for {}!'.format(varname))   
        
        # Subplot (0,0): Ref
        ax0.coastlines()
        if not match_cbar: [vmin, vmax] = [vmin_ref, vmax_ref]
        plot0 = ax0.imshow(ds_ref, extent=(minlon, maxlon, minlat, maxlat), 
                           cmap=WhGrYlRd,vmin=vmin, vmax=vmax)
        ax0.set_title('{} (Ref){}\n{}'.format(refstr,subtitle_extra,llres)) 
        cb = plt.colorbar(plot0, ax=ax0, orientation='horizontal', pad=0.10)
        if (vmax-vmin) < 0.1 or (vmax-vmin) > 100:
            cb.locator = ticker.MaxNLocator(nbins=4)
            cb.update_ticks()
        cb.set_label(units)
        
        # Subplot (0,1): Dev
        ax1.coastlines()
        if not match_cbar: [vmin, vmax] = [vmin_dev, vmax_dev]
        plot1 = ax1.imshow(ds_dev, extent=(minlon, maxlon, minlat, maxlat), 
                           cmap=WhGrYlRd,vmin=vmin, vmax=vmax)
        ax1.set_title('{} (Dev){}\n{}'.format(devstr,subtitle_extra,llres)) 
        cb = plt.colorbar(plot1, ax=ax1, orientation='horizontal', pad=0.10)
        if (vmax-vmin) < 0.1 or (vmax-vmin) > 100:
            cb.locator = ticker.MaxNLocator(nbins=4)
            cb.update_ticks()
        cb.set_label(units)

        # Calculate difference, get dynamic range, configure colorbar, use gray for NaNs
        absdiff = np.array(ds_dev) - np.array(ds_ref)
        diffabsmax = max([np.abs(np.nanmin(absdiff)), np.abs(np.nanmax(absdiff))])        
        cmap = mpl.cm.RdBu_r
        cmap.set_bad(color='gray')
            
        # Subplot (1,0): Difference, dynamic range
        [vmin, vmax] = [-diffabsmax, diffabsmax]
        ax2.coastlines()
        plot2 = ax2.imshow(absdiff, extent=(minlon, maxlon, minlat, maxlat), 
                           cmap=cmap,vmin=vmin, vmax=vmax)
        ax2.set_title('Difference\nDev - Ref, Dynamic Range') 
        cb = plt.colorbar(plot2, ax=ax2, orientation='horizontal', pad=0.10)
        if (vmax-vmin) < 0.1 or (vmax-vmin) > 100:
            cb.locator = ticker.MaxNLocator(nbins=4)
            cb.update_ticks()
        cb.set_label(units)
        if np.all(absdiff==0): 
            cb.ax.set_xticklabels(['0.0', '0.0', '0.0', '0.0', '0.0']) 
        
        # Subplot (1,1): Difference, restricted range
        [vmin, vmax] = [-1e-9, 1e-9] # Placeholder; need a dictional for this (species: limit)
        ax3.coastlines()
        plot3 = ax3.imshow(absdiff, extent=(minlon, maxlon, minlat, maxlat), 
                           cmap=cmap,vmin=vmin, vmax=vmax)
        ax3.set_title('Difference\nDev - Ref, Fixed Range') 
        cb = plt.colorbar(plot3, ax=ax3, orientation='horizontal', pad=0.10)
        if (vmax-vmin) < 0.1 or (vmax-vmin) > 100:
            cb.locator = ticker.MaxNLocator(nbins=4)
            cb.update_ticks()
        cb.set_label(units)
        if np.all(absdiff==0): 
            cb.ax.set_xticklabels(['0.0', '0.0', '0.0', '0.0', '0.0'])

        # Calculate fractional difference, get dynamic range, set 0/0 to Nan
        fracdiff = (np.array(ds_dev) - np.array(ds_ref)) / np.array(ds_ref)
        fracdiff[(ds_dev==0) & (ds_ref==0)] = np.nan
        fracdiffabsmax = max([np.abs(np.nanmin(fracdiff)), np.abs(np.nanmax(fracdiff))])
        
        # Subplot (2,0): Fractional Difference, full dynamic range
        [vmin, vmax] = [-fracdiffabsmax, fracdiffabsmax]
        ax4.coastlines()
        plot4 = ax4.imshow(fracdiff, extent=(minlon, maxlon, minlat, maxlat), vmin=vmin, vmax=vmax, cmap=cmap)
        ax4.set_title('Fractional Difference\n(Dev-Ref)/Ref, Dynamic Range') 
        cb = plt.colorbar(plot4, ax=ax4, orientation='horizontal', pad=0.10)
        if (vmax-vmin) < 0.1 or (vmax-vmin) > 100:
            cb.locator = ticker.MaxNLocator(nbins=4)
            cb.update_ticks()
            if np.all(absdiff==0): 
                cb.ax.set_xticklabels(['0.0', '0.0', '0.0', '0.0', '0.0'])
        cb.set_label('unitless')  
        
        # Subplot (2,1): Fractional Difference, restricted
        [vmin, vmax] = [-2, 2]
        #[vmin, vmax] = [-0.5, 2] # doesn't work with this colorbar. Need to customize one. Already in gamap?
                                  # Switch to this if change to ratios (dev/ref)
        ax5.coastlines()
        plot5 = ax5.imshow(fracdiff, extent=(minlon, maxlon, minlat, maxlat),
                           cmap=cmap,vmin=vmin, vmax=vmax)
        ax5.set_title('Fractional Difference\n(Dev-Ref)/Ref, Fixed Range') 
        cb = plt.colorbar(plot5, ax=ax5, orientation='horizontal', pad=0.10)
        if np.all(absdiff==0): 
            cb.locator = ticker.MaxNLocator(nbins=4)
            cb.update_ticks()
            cb.ax.set_xticklabels(['0.0', '0.0', '0.0', '0.0', '0.0'])
        cb.set_label('unitless') 
            
        if savepdf:    
            pdf.savefig(figs)
            plt.close(figs)
            
    if savepdf: pdf.close()

def compare_gchp_single_level(refdata, refstr, devdata, devstr, varlist=None, weightsdir='.', ilev=0, 
                         itime=0, savepdf=False, pdfname='gchp_vs_gchp_map.pdf', match_cbar=True, 
                         full_ratio_range=False, normalize_by_area=False, area_ref=None, 
                         area_dev=None, check_units=True, flip_vert=False):
    
    # If no varlist is passed, plot all (surface only for 3D)
    if varlist == None:
        [varlist, commonvars2D, commonvars3D] = compare_varnames(refdata, devdata)
        print('Plotting all common variables (surface only if 3D)')
    n_var = len(varlist)

    # Get cubed sphere grids and regridder 
    # for now, do not regrid for gchp vs gchp. Assume same grid.
    csres_ref = refdata['lon'].size
    [csgrid_ref, csgrid_list_ref] = make_grid_CS(csres_ref)
    csres_dev = devdata['lon'].size
    [csgrid_dev, csgrid_list_dev] = make_grid_CS(csres_dev)

    # Create pdf (if saving)
    if savepdf:
        print('\nCreating {} for {} variables'.format(pdfname,n_var))
        pdf = PdfPages(pdfname)

    # Loop over variables
    for ivar in range(n_var):
        if savepdf: print('{} '.format(ivar), end='')
        varname = varlist[ivar]
        
        # Do some checks: dimensions and units
        varndim_ref = refdata[varname].ndim
        varndim_dev = devdata[varname].ndim
        if check_units: 
            assert varndim_ref == varndim_dev, 'Dimensions do not agree for {}!'.format(varname)
        units_ref = refdata[varname].units
        units_dev = devdata[varname].units
        if check_units: 
            assert units_ref == units_dev, 'Units do not match for {}!'.format(varname)
            
        # if normalizing by area, adjust units to be per m2, and adjust title string
        units = units_ref
        varndim = varndim_ref
        subtitle_extra = ''
                    
        # Slice the data
        if varndim == 4: 
            if flip_vert: 
                ds_ref = refdata[varname].isel(time=itime,lev=71-ilev)
                ds_dev = devdata[varname].isel(time=itime,lev=71-ilev)
            else: 
                ds_ref = refdata[varname].isel(time=itime,lev=ilev)
                ds_dev = devdata[varname].isel(time=itime,lev=ilev)
        elif varndim == 3: 
            refds = refdata[varname].isel(time=itime)
            devds = devdata[varname].isel(time=itime)
            
        # if normalizing by area, transform on the native grid and adjust units and subtitle string
        exclude_list = ['WetLossConvFrac','Prod_','Loss_']
        if normalize_by_area and not any(s in varname for s in exclude_list):
            ds_ref.values = ds_ref.values / area_ref
            ds_dev.values = ds_dev.values / area_dev
            units = '{} m-2'.format(units)
            subtitle_extra = ', Normalized by Area'
            
        # Regrid the slices (skip for gchp vs gchp for now)
        csdata_ref = ds_ref.data.reshape(6,csres_ref,csres_ref)
        csdata_dev = ds_dev.data.reshape(6,csres_dev,csres_dev)
        vmin_ref = csdata_ref.min()
        vmin_dev = csdata_dev.min()
        vmin_cmn = np.min([vmin_ref, vmin_dev])
        vmax_ref = csdata_ref.max()
        vmax_dev = csdata_dev.max()
        vmax_cmn = np.max([vmax_ref, vmax_dev])
        if match_cbar: [vmin, vmax] = [vmin_cmn, vmax_cmn]
        
        # Create 2x2 figure
        figs, ((ax0, ax1), (ax2, ax3)) = plt.subplots(2, 2, figsize=[12,9], 
                                                      subplot_kw={'projection': crs.PlateCarree()})
        # Give the figure a title
        offset = 0.96
        fontsize=25
        if varndim == 4:
            if ilev == 0: levstr = 'Surface'
            elif ilev == 22: levstr = '500 hPa'
            else: levstr = 'Level ' +  str(ilev-1)
            figs.suptitle('{}, {}'.format(varname,levstr), fontsize=fontsize, y=offset)
        elif varndim == 3: 
            figs.suptitle('{}'.format(varname), fontsize=fontsize, y=offset)
        else:
            print('varndim is 2 for {}! Must be 3 or 4.'.format(varname))
            
        # Subplot (0,0): Ref
        ax0.coastlines()
        if not match_cbar: [vmin, vmax] = [vmin_ref, vmax_ref]        
        masked_csdata = np.ma.masked_where(np.abs(csgrid_ref['lon'] - 180) < 2, csdata_ref)
        for i in range(6):
            plot0 = ax0.pcolormesh(csgrid_ref['lon_b'][i,:,:], csgrid_ref['lat_b'][i,:,:], masked_csdata[i,:,:], 
                                   cmap=WhGrYlRd,vmin=vmin, vmax=vmax)
        ax0.set_title('{} (Ref){}\nC{}'.format(refstr,subtitle_extra,str(csres_ref)))
        cb = plt.colorbar(plot0, ax=ax0, orientation='horizontal', pad=0.10)
        if (vmax-vmin) < 0.1 or (vmax-vmin) > 100:
            cb.locator = ticker.MaxNLocator(nbins=4)
            cb.update_ticks()
        cb.set_label(units)     
        
        # Subplot (0,1): Dev
        ax1.coastlines()
        if not match_cbar: [vmin, vmax] = [vmin_dev, vmax_dev]        
        masked_csdata = np.ma.masked_where(np.abs(csgrid_dev['lon'] - 180) < 2, csdata_dev)
        for i in range(6):
            plot1 = ax1.pcolormesh(csgrid_dev['lon_b'][i,:,:], csgrid_dev['lat_b'][i,:,:], 
                                   masked_csdata[i,:,:], cmap=WhGrYlRd,vmin=vmin, vmax=vmax)
        ax1.set_title('{} (Dev){}\nC{}'.format(devstr,subtitle_extra,str(csres_dev)))
        cb = plt.colorbar(plot1, ax=ax1, orientation='horizontal', pad=0.10)
        if (vmax-vmin) < 0.1 or (vmax-vmin) > 100:
            cb.locator = ticker.MaxNLocator(nbins=4)
            cb.update_ticks()
        cb.set_label(units)  
            
        # Subplot (1,0): Difference
        gc_absdiff = csdata_dev - csdata_ref
        diffabsmax = max([np.abs(gc_absdiff.min()), np.abs(gc_absdiff.max())])
        [vmin, vmax] = [-diffabsmax, diffabsmax]
        ax2.coastlines()
        # assume the same grid for now in gchp vs gchp
        masked_csdata = np.ma.masked_where(np.abs(csgrid_dev['lon'] - 180) < 2, gc_absdiff)
        for i in range(6):
            plot2 = ax2.pcolormesh(csgrid_dev['lon_b'][i,:,:], csgrid_dev['lat_b'][i,:,:], 
                                   masked_csdata[i,:,:], cmap='RdBu_r',vmin=vmin, vmax=vmax)
        ax2.set_title('Difference\n(Dev - Ref)')   
        cb = plt.colorbar(plot2, ax=ax2, orientation='horizontal', pad=0.10)
        if (vmax-vmin) < 0.1 or (vmax-vmin) > 100:
            cb.locator = ticker.MaxNLocator(nbins=4)
            cb.update_ticks()
        cb.set_label(units)  
    
        # Subplot (1,1): Fractional Difference (restrict to +/-2)
        gc_fracdiff = (csdata_dev - csdata_ref) / csdata_ref
        if full_ratio_range: [vmin, vmax] = [None, None]
        else: [vmin, vmax] = [-2, 2]
        ax3.coastlines()
        # assume the same grid for now in gchp vs gchp
        masked_csdata = np.ma.masked_where(np.abs(csgrid_dev['lon'] - 180) < 2, gc_fracdiff)
        for i in range(6):
            plot3 = ax3.pcolormesh(csgrid_dev['lon_b'][i,:,:], csgrid_dev['lat_b'][i,:,:], 
                                   masked_csdata[i,:,:], cmap='RdBu_r',vmin=vmin, vmax=vmax)
        ax3.set_title('Fractional Difference\n(Dev - Ref)/Ref')   
        cb = plt.colorbar(plot3, ax=ax3, orientation='horizontal', pad=0.10)
        cb.set_clim(vmin=vmin, vmax=vmax)
        cb.set_label('unitless')     
            
        if savepdf:    
            pdf.savefig(figs)
            plt.close(figs)
            
    if savepdf: pdf.close()

# Add docstrings later. Use this function for benchmarking or general comparisons.
def compare_zonal_mean(refdata, refstr, devdata, devstr, llres, varlist=None, itime=0, savepdf=False, 
                       pdfname='zonalmean.pdf', match_cbar=True, normalize_by_area=False, check_units=True ):

    # If no varlist is passed, plot all 3D variables in the dataset
    if varlist == None:
        [commonvars, commonvars2D, varlist] = compare_varnames(refdata, devdata)
        print('Plotting all 3D variables')
    n_var = len(varlist)
    
    # Get lat-lon grid
    llgrid = make_grid_LL(llres)
    
    # Universal plot setup
    xtick_positions = np.arange(-90,91,30)
    xticklabels = ['{}$\degree$'.format(x) for x in xtick_positions]
    ytick_positions = np.arange(0,61,20)
    yticklabels = [str(y) for y in ytick_positions]
    
    # Create pdf (if saving)
    if savepdf:
        print('\nCreating {} for {} variables'.format(pdfname, n_var))
        pdf = PdfPages(pdfname)

    # Loop over variables
    for ivar in range(n_var):
        if savepdf: print('{} '.format(ivar), end='')
        varname = varlist[ivar]
        
        # Do some checks: dimensions and units
        varndim_ref = refdata[varname].ndim
        varndim_dev = devdata[varname].ndim
        nlev = 72
        #assert varndim_ref == varndim_dev, 'Dimensions do not agree for {}!'.format(varname)
        units_ref = refdata[varname].units.strip()
        units_dev = devdata[varname].units.strip()
        assert units_ref == units_dev, 'Units do not match for {}!'.format(varname)
        
        # Set plot extent
        extent=(-90,90,0,nlev)
        
        # if normalizing by area, adjust units to be per m2, and adjust title string
        units = units_ref
        varndim = varndim_ref
        subtitle_extra = ''            
        
        # Slice the data.  Need to handle different incoming number of dimensions and the bpch case
        # where no time dimension is includes.
        if 'time' in refdata[varname].dims:
            ds_ref = refdata[varname].isel(time=itime)
        else:
            ds_ref = refdata[varname]
        if 'time' in devdata[varname].dims:
            ds_dev = devdata[varname].isel(time=itime)
        else:
            ds_dev = devdata[varname]
            
        # if normalizing by area, transform on the native grid and adjust units and subtitle string
        exclude_list = ['WetLossConvFrac','Prod_','Loss_']
        if normalize_by_area and not any(s in varname for s in exclude_list):
            ds_ref.values = ds_ref.values / refdata['AREAM2'].values[np.newaxis,:,:]
            ds_dev.values = ds_dev.values / devdata['AREAM2'].values[np.newaxis,:,:]
            units = '{} m-2'.format(units)
            subtitle_extra = ', Normalized by Area'
            
        # Calculate zonal mean of the regridded data
        zm_ref = ds_ref.mean(dim='lon')
        zm_dev = ds_dev.mean(dim='lon')
        
        # Get min and max for colorbar limits
        [vmin_ref, vmax_ref] = [zm_ref.min(), zm_ref.max()]
        [vmin_dev, vmax_dev] = [zm_dev.min(), zm_dev.max()]
        vmin_cmn = np.min([vmin_ref, vmin_dev])
        vmax_cmn = np.max([vmax_ref, vmax_dev])
        if match_cbar: [vmin, vmax] = [vmin_cmn, vmax_cmn]
        
        # Create 2x2 figure
        figs, ((ax0, ax1), (ax2, ax3), (ax4, ax5)) = plt.subplots(3, 2, figsize=[12,15.3], 
                                                      subplot_kw={'projection': crs.PlateCarree()})
        # Give the page a title
        offset = 0.96
        fontsize=25
        figs.suptitle('{}, Zonal Mean'.format(varname), fontsize=fontsize, y=offset)

        # Subplot 0: Ref
        if not match_cbar: [vmin, vmax] = [vmin_ref, vmax_ref]
        plot0 = ax0.imshow(zm_ref, cmap=WhGrYlRd, extent=extent, vmin=vmin, vmax=vmax)
        ax0.set_title('{} (Ref){}\n{}'.format(refstr, subtitle_extra, llres ))
        ax0.set_aspect('auto')
        ax0.set_xticks(xtick_positions)
        ax0.set_xticklabels(xticklabels)
        ax0.set_yticks(ytick_positions)
        ax0.set_yticklabels(yticklabels)
        cb = plt.colorbar(plot0, ax=ax0, orientation='horizontal', pad=0.10)
        if (vmax-vmin) < 0.001 or (vmax-vmin) > 1000:
            cb.locator = ticker.MaxNLocator(nbins=4)
            cb.update_ticks()
        cb.set_label(units)
        
        # Subplot 1: Dev
        if not match_cbar: [vmin, vmax] = [vmin_dev, vmax_dev]
        plot1 = ax1.imshow(zm_dev, cmap=WhGrYlRd, extent=extent, vmin=vmin, vmax=vmax)
        ax1.set_title('{} (Dev){}\n{}'.format(devstr, subtitle_extra, llres))
        ax1.set_aspect('auto')
        ax1.set_xticks(xtick_positions)
        ax1.set_xticklabels(xticklabels)
        ax1.set_yticks(ytick_positions)
        ax1.set_yticklabels(yticklabels)
        cb = plt.colorbar(plot1, ax=ax1, orientation='horizontal', pad=0.10)
        if (vmax-vmin) < 0.001 or (vmax-vmin) > 1000:
            cb.locator = ticker.MaxNLocator(nbins=4)
            cb.update_ticks()
        cb.set_label(units)
        
        # Calculate zonal mean difference
        zm_diff = np.array(zm_dev) - np.array(zm_ref)
                
        # Subplot 2: Difference, dynamic range
        diffabsmax = max([np.abs(zm_diff.min()), np.abs(zm_diff.max())])
        [vmin, vmax] = [-diffabsmax, diffabsmax]
        plot2 = ax2.imshow(zm_diff, cmap='RdBu_r', extent=extent, vmin=vmin, vmax=vmax)
        ax2.set_title('Difference\nDev - Ref, Dynamic Range')
        ax2.set_aspect('auto')
        ax2.set_xticks(xtick_positions)
        ax2.set_xticklabels(xticklabels)
        ax2.set_yticks(ytick_positions)
        ax2.set_yticklabels(yticklabels)
        cb = plt.colorbar(plot2, ax=ax2, orientation='horizontal', pad=0.10)
        if (vmax-vmin) < 0.001 or (vmax-vmin) > 1000:
            cb.locator = ticker.MaxNLocator(nbins=4)
            cb.update_ticks()
        cb.set_label(units)
        
        # Subplot 3: Difference, restricted range
        [vmin, vmax] = [-1e-9, 1e-9] # need a dictional for this (species: limit)
        plot3 = ax3.imshow(zm_diff, cmap='RdBu_r', extent=extent, vmin=vmin, vmax=vmax)
        ax3.set_title('Difference\nDev - Ref, Fixed Range')
        ax3.set_aspect('auto')
        ax3.set_xticks(xtick_positions)
        ax3.set_xticklabels(xticklabels)
        ax3.set_yticks(ytick_positions)
        ax3.set_yticklabels(yticklabels)
        cb = plt.colorbar(plot3, ax=ax3, orientation='horizontal', pad=0.10)
        if (vmax-vmin) < 0.001 or (vmax-vmin) > 1000:
            cb.locator = ticker.MaxNLocator(nbins=4)
            cb.update_ticks()
        cb.set_label(units)
        
        # Zonal mean fractional difference
        zm_fracdiff = (np.array(zm_dev) - np.array(zm_ref)) / np.array(zm_ref)
            
        # Subplot 4: Fractional Difference, dynamic range
        fracdiffabsmax = max([np.abs(zm_fracdiff.min()), np.abs(zm_fracdiff.max())])
        [vmin, vmax] = [-fracdiffabsmax, fracdiffabsmax]
        plot4 = ax4.imshow(zm_fracdiff, vmin=vmin, vmax=vmax, cmap='RdBu_r', extent=extent)
        ax4.set_title('Fractional Difference\n(Dev-Ref)/Ref, Dynamic Range')
        ax4.set_aspect('auto')
        ax4.set_xticks(xtick_positions)
        ax4.set_xticklabels(xticklabels)
        ax4.set_yticks(ytick_positions)
        ax4.set_yticklabels(yticklabels)
        cb = plt.colorbar(plot4, ax=ax4, orientation='horizontal', pad=0.10)
        if (vmax-vmin) < 0.1 or (vmax-vmin) > 100:
            cb.locator = ticker.MaxNLocator(nbins=4)
            cb.update_ticks()
            if np.all(zm_fracdiff==0): 
                cb.ax.set_xticklabels(['0.0', '0.0', '0.0', '0.0', '0.0'])
        cb.set_clim(vmin=vmin, vmax=vmax)
        cb.set_label('unitless')   
        
        # Subplot 5: Fractional Difference, restricted range
        [vmin, vmax] = [-2, 2]
        plot5 = ax5.imshow(zm_fracdiff, vmin=vmin, vmax=vmax, cmap='RdBu_r', extent=extent)
        ax5.set_title('Fractional Difference\n(Dev-Ref)/Ref, Fixed Range')
        ax5.set_aspect('auto')
        ax5.set_xticks(xtick_positions)
        ax5.set_xticklabels(xticklabels)
        ax5.set_yticks(ytick_positions)
        ax5.set_yticklabels(yticklabels)
        cb = plt.colorbar(plot5, ax=ax5, orientation='horizontal', pad=0.10)
        if (vmax-vmin) < 0.1 or (vmax-vmin) > 100:
            cb.locator = ticker.MaxNLocator(nbins=4)
            cb.update_ticks()
            if np.all(zm_fracdiff==0): 
                cb.ax.set_xticklabels(['0.0', '0.0', '0.0', '0.0', '0.0'])
        cb.set_clim(vmin=vmin, vmax=vmax)
        cb.set_label('unitless') 
            
        if savepdf:    
            pdf.savefig(figs)
            plt.close(figs)
            
    if savepdf: pdf.close()

def compare_gchp_zonal_mean(refdata, refstr, devdata, devstr, varlist=None,
                            weightsdir='.', itime=0, llres_cmp='1x1.25', savepdf=False,
                            pdfname='gchp_vs_gchp_map.pdf', match_cbar=True, full_ratio_range=False,
                            normalize_by_area=False, area_ref=None, area_dev=None, flip_vert=False):

    # If no varlist is passed, plot all 3D variables in the dataset
    if varlist == None:
        [commonvars, commonvars2D, varlist] = compare_varnames(refdata, devdata)
        print('Plotting all 3D variables')
    n_var = len(varlist)
    
    # Get lat-lon grid
    llgrid_cmp = make_grid_LL(llres_cmp)

    # Get cubed sphere grid and regridder for first data set
    csres_ref = refdata['lon'].size
    [csgrid_ref, csgrid_list_ref] = make_grid_CS(csres_ref)
    cs_regridder_list_ref = make_regridder_C2L(csres_ref, llres_cmp, weightsdir=weightsdir, 
                                               reuse_weights=True)

    # Get cubed sphere grid and regridder for first data set
    csres_dev = devdata['lon'].size
    [csgrid_dev, csgrid_list_dev] = make_grid_CS(csres_dev)
    cs_regridder_list_dev = make_regridder_C2L(csres_dev, llres_cmp, weightsdir=weightsdir, 
                                               reuse_weights=True)
    
    # Universal plot setup
    xtick_positions = np.arange(-90,91,30)
    xticklabels = ['{}$\degree$'.format(x) for x in xtick_positions]
    ytick_positions = np.arange(0,61,20)
    yticklabels = [str(y) for y in ytick_positions]
    
    # Create pdf (if saving)
    if savepdf:
        print('\nCreating {} for {} variables'.format(pdfname, n_var))
        pdf = PdfPages(pdfname)

    # Loop over variables
    for ivar in range(n_var):
        if savepdf: print('{} '.format(ivar), end='')
        varname = varlist[ivar]
        
        # Do some checks: dimensions and units
        varndim_ref = refdata[varname].ndim
        varndim_dev = devdata[varname].ndim
        nlev = 72
        assert varndim_ref == varndim_dev, 'GCHP dimensions do not agree for {}!'.format(varname)
        units_ref = refdata[varname].units
        units_dev = devdata[varname].units
        assert units_ref == units_dev, 'GCHP units do not match for {}!'.format(varname)
        
        # Set plot extent
        extent=(-90,90,0,nlev)
        
        # if normalizing by area, adjust units to be per m2, and adjust title string
        units = units_ref
        varndim = varndim_ref
        subtitle_extra = ''
        
        # Slice the data
        ds_ref = refdata[varname].isel(time=itime)
        ds_dev = devdata[varname].isel(time=itime)

        # if normalizing by area, transform on the native grid and adjust units and subtitle string
        exclude_list = ['WetLossConvFrac','Prod_','Loss_']
        if normalize_by_area and not any(s in varname for s in exclude_list):
            ds_ref.values = ds_ref.values / area_ref.values[np.newaxis,:,:]
            ds_dev.values = ds_dev.values / area_dev.values[np.newaxis,:,:]
            units = '{} m-2'.format(units)
            subtitle_extra = ', Normalized by Area'
            
        # Regrid the slices
        if flip_vert: 
            ds_ref.data = ds_ref.data[::-1,:,:]
            ds_dev.data = ds_dev.data[::-1,:,:]
        csdata_ref = ds_ref.data.reshape(nlev,6,csres_ref,csres_ref).swapaxes(0,1)
        csdata_dev = ds_dev.data.reshape(nlev,6,csres_dev,csres_dev).swapaxes(0,1)
        lldata_ref = np.zeros([nlev, llgrid_cmp['lat'].size, llgrid_cmp['lon'].size])
        lldata_dev = np.zeros([nlev, llgrid_cmp['lat'].size, llgrid_cmp['lon'].size])
        for i in range(6):
            regridder_ref = cs_regridder_list_ref[i]
            lldata_ref += regridder_ref(csdata_ref[i])
            regridder_dev = cs_regridder_list_dev[i]
            lldata_dev += regridder_dev(csdata_dev[i])
        
        # Calculate zonal mean of the regridded data
        zm_ref = lldata_ref.mean(axis=2)
        zm_dev = lldata_dev.mean(axis=2)
            
        # Get min and max for colorbar limits
        [vmin_ref, vmax_ref] = [zm_ref.min(), zm_ref.max()]
        [vmin_dev, vmax_dev] = [zm_dev.min(), zm_dev.max()]
        vmin_cmn = np.min([vmin_ref, vmin_dev])
        vmax_cmn = np.max([vmax_ref, vmax_dev])
        if match_cbar: [vmin, vmax] = [vmin_cmn, vmax_cmn]
        
        # Create 2x2 figure
        figs, ((ax0, ax1), (ax2, ax3)) = plt.subplots(2, 2, figsize=[12,12], 
                                                      subplot_kw={'projection': crs.PlateCarree()})
        # Give the page a title
        offset = 0.96
        fontsize=25
        figs.suptitle('{}, Zonal Mean'.format(varname), fontsize=fontsize, y=offset)

        # Subplot 0: Ref
        if not match_cbar: [vmin, vmax] = [vmin_ref, vmax_ref]
        plot0 = ax0.imshow(zm_ref, cmap=WhGrYlRd, extent=extent, vmin=vmin, vmax=vmax)
        ax0.set_title('{} (Ref){}\n{} regridded from {}'.format(refstr, subtitle_extra, 
                                                                llres_cmp, csres_ref))
        ax0.set_aspect('auto')
        ax0.set_xticks(xtick_positions)
        ax0.set_xticklabels(xticklabels)
        ax0.set_yticks(ytick_positions)
        ax0.set_yticklabels(yticklabels)
        cb = plt.colorbar(plot0, ax=ax0, orientation='horizontal', pad=0.10)
        if (vmax-vmin) < 0.001 or (vmax-vmin) > 1000:
            cb.locator = ticker.MaxNLocator(nbins=4)
            cb.update_ticks()
        cb.set_label(units)
        
        # Subplot 1: Dev
        if not match_cbar: [vmin, vmax] = [vmin_dev, vmax_dev]
        plot1 = ax1.imshow(zm_dev, cmap=WhGrYlRd, extent=extent, vmin=vmin, vmax=vmax)
        ax1.set_title('{} (Dev){}\n{} regridded from {}'.format(devstr, subtitle_extra, 
                                                                llres_cmp, csres_dev))
        ax1.set_aspect('auto')
        ax1.set_xticks(xtick_positions)
        ax1.set_xticklabels(xticklabels)
        ax1.set_yticks(ytick_positions)
        ax1.set_yticklabels(yticklabels)
        cb = plt.colorbar(plot1, ax=ax1, orientation='horizontal', pad=0.10)
        if (vmax-vmin) < 0.001 or (vmax-vmin) > 1000:
            cb.locator = ticker.MaxNLocator(nbins=4)
            cb.update_ticks()
        cb.set_label(units)
            
        # Subplot 2: Difference
        zm_absdiff = zm_dev - zm_ref
        diffabsmax = max([np.abs(zm_absdiff.min()), np.abs(zm_absdiff.max())])
        [vmin, vmax] = [-diffabsmax, diffabsmax]
        plot2 = ax2.imshow(zm_absdiff, cmap='RdBu_r', extent=extent, vmin=vmin, vmax=vmax)
        ax2.set_title('Difference\n(Dev - Ref)')
        ax2.set_aspect('auto')
        ax2.set_xticks(xtick_positions)
        ax2.set_xticklabels(xticklabels)
        ax2.set_yticks(ytick_positions)
        ax2.set_yticklabels(yticklabels)
        cb = plt.colorbar(plot2, ax=ax2, orientation='horizontal', pad=0.10)
        if (vmax-vmin) < 0.001 or (vmax-vmin) > 1000:
            cb.locator = ticker.MaxNLocator(nbins=4)
            cb.update_ticks()
        cb.set_label(units)
        
        # Subplot 3: Fractional Difference (restrict to +/-2)
        zm_fracdiff = (zm_dev - zm_ref) / zm_ref
        if full_ratio_range: [vmin, vmax] = [None, None]
        else: [vmin, vmax] = [-2, 2]
        plot3 = ax3.imshow(zm_fracdiff, vmin=vmin, vmax=vmax, cmap='RdBu_r', extent=extent)
        ax3.set_title('Fractional Difference\n(Dev - Ref)/Ref')
        ax3.set_aspect('auto')
        ax3.set_xticks(xtick_positions)
        ax3.set_xticklabels(xticklabels)
        ax3.set_yticks(ytick_positions)
        ax3.set_yticklabels(yticklabels)
        cb = plt.colorbar(plot3, ax=ax3, orientation='horizontal', pad=0.10)
        cb.set_clim(vmin=vmin, vmax=vmax)
        cb.set_label('unitless')      
            
        if savepdf:    
            pdf.savefig(figs)
            plt.close(figs)
            
    if savepdf: pdf.close()
    
def add_bookmarks_to_pdf( pdfname, varlist, remove_prefix='' ):
    # Existing pdf
    pdfobj = open(pdfname,"rb")
    input = PdfFileReader(pdfobj)
    numpages = input.getNumPages()
    
    # Pdf write and new filename
    pdfname_tmp = pdfname+'_with_bookmarks.pdf'
    output = PdfFileWriter()
    
    # Loop over variables (pages) in the file, removing the diagnostic prefix
    varnamelist = [k.replace(remove_prefix,'') for k in varlist]
    for i, varname in enumerate(varnamelist):
        output.addPage(input.getPage(i))
        output.addBookmark(varname,i)
        output.setPageMode('/UseOutlines')
        
    # Write to new file
    outputstream = open(pdfname_tmp,'wb')
    output.write(outputstream) 
    outputstream.close()
    
    # Replace the old file with the new
    os.rename(pdfname_tmp, pdfname)

def add_hierarchical_bookmarks_to_pdf( pdfname, catdict, remove_prefix='' ):
    # Existing pdf
    pdfobj = open(pdfname,"rb")
    input = PdfFileReader(pdfobj)
    numpages = input.getNumPages()
    
    # Pdf write and new filename
    pdfname_tmp = pdfname+'_with_bookmarks.pdf'
    output = PdfFileWriter()
    
    # Loop over variables (pages) in the file, removing the diagnostic prefix
    varnamelist = [k.replace(remove_prefix,'') for k in varlist]
    for i, varname in enumerate(varnamelist):
        output.addPage(input.getPage(i))
        output.addBookmark(varname,i)
        output.setPageMode('/UseOutlines')
        
    # Write to new file
    outputstream = open(pdfname_tmp,'wb')
    output.write(outputstream) 
    outputstream.close()
    
    # Replace the old file with the new
    os.rename(pdfname_tmp, pdfname)


