""" Specific utilities re-factored from the benchmarking utilities. """

import os
import shutil
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import json
import copy
from json import load as json_load_file
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap
from cartopy import crs
from cartopy.mpl.geoaxes import GeoAxes  # for assertion
from PyPDF2 import PdfFileWriter, PdfFileReader
from .plot import WhGrYlRd, add_latlon_ticks
from .grid.horiz import make_grid_LL, make_grid_CS
from .grid.regrid import make_regridder_C2L, make_regridder_L2L
from .grid.gc_vertical import GEOS_72L_grid
from . import core
from .units import convert_units

cmap_abs = WhGrYlRd  # for plotting absolute magnitude
cmap_diff = 'RdBu_r'  # for plotting difference

aod_spc = 'aod_species.json'
spc_categories = 'benchmark_categories.json'
emission_spc = 'emission_species.json' 
emission_inv = 'emission_inventories.json' 

def compare_single_level(refdata, refstr, devdata, devstr, varlist=None,
                         ilev=0, itime=0,  weightsdir=None, pdfname='',
                         cmpres=None, match_cbar=True, normalize_by_area=False,
                         enforce_units=True, flip_ref=False, flip_dev=False,
                         use_cmap_RdBu=False, verbose=False, 
                         log_color_scale=False, sigdiff_list=[]):
    '''
    Create single-level 3x2 comparison map plots for variables common in two xarray
    datasets. Optionally save to PDF. 

    Args:
        refdata : xarray dataset
            Dataset used as reference in comparison

        refstr  : str
            String description for reference data to be used in plots
     
        devdata : xarray dataset
            Dataset used as development in comparison

        devstr  : str
            String description for development data to be used in plots
 
    Keyword Args (optional):
        varlist : list of strings
            List of xarray dataset variable names to make plots for
            Default value: None (will compare all common variables)

        ilev : integer
            Dataset level dimension index using 0-based system
            Default value: 0   

        itime : integer
            Dataset time dimension index using 0-based system
            Default value: 0

        weightsdir : str
            Directory path for storing regridding weights
            Default value: None (will create/store weights in current directory)

        pdfname : str
            File path to save plots as PDF
            Default value: Empty string (will not create PDF)

        cmpres : str
            String description of grid resolution at which to compare datasets
            Default value: None (will compare at highest resolution of ref and dev)

        match_cbar : boolean
            Logical indicating whether to use same colorbar bounds across plots
            Default value: True

        normalize_by_area : boolean
            Logical indicating whether to normalize raw data by grid area
            Default value: False

        enforce_units : boolean
            Logical to force an error if reference and development variables units differ
            Default value: True

        flip_ref : boolean
            Logical to flip the vertical dimension of reference dataset 3D variables
            Default value: False

        flip_dev : boolean
            Logical to flip the vertical dimension of development dataset 3D variables
            Default value: False

        use_cmap_RdBu : boolean
            Logical to used a blue-white-red colormap for plotting raw reference and
            development datasets.
            Default value: False

        verbose : boolean
            Logical to enable informative prints
            Default value: False

        log_color_scale: boolean         
            Logical to enable plotting data (not diffs) on a log color scale.
            Default value: False

        sigdiff_list: list of str
            Returns a list of all quantities having significant differences.
            The criteria is: |max(fractional difference)| > 0.1
            Default value: []

    Returns:
        Nothing

    Example:
        >>> import matplotlib.pyplot as plt
        >>> import xarray as xr
        >>> from gcpy import benchmark
        >>> refds = xr.open_dataset('path/to/ref.nc4')
        >>> devds = xr.open_dataset('path/to/dev.nc4')
        >>> varlist = ['SpeciesConc_O3', 'SpeciesConc_NO']
        >>> benchmark.compare_single_level( refds, '12.3.2', devds, 'bug fix', varlist=varlist )
        >>> plt.show()
    '''

    # TODO: refactor this function and zonal mean plot function. There is a lot of overlap and
    # repeated code that could be abstracted.

    # Error check arguments
    if not isinstance(refdata, xr.Dataset):
        raise ValueError('The refdata argument must be an xarray Dataset!')

    if not isinstance(devdata, xr.Dataset):
        raise ValueError('The devdata argument must be an xarray Dataset!')

    # If no varlist is passed, plot all (surface only for 3D)
    if varlist == None:
        [commonvars, commonvarsOther, commonvars2D, commonvars3D] = core.compare_varnames(refdata, devdata)
        varlist = commonvars3D+commonvars2D
        print('Plotting all common variables (surface only if 3D)')
    n_var = len(varlist)

    # If no weightsdir is passed, set to current directory in case it is needed
    if weightsdir == None:
        weightsdir = '.'

    # If no pdf name passed, then do not save to PDF
    savepdf = True
    if pdfname == '':
        savepdf = False

    ####################################################################
    # Determine input grid resolutions and types
    ####################################################################

    # GCC output and GCHP output using pre-v1.0.0 MAPL have lat and lon dims

    # ref
    vdims = refdata.dims
    if 'lat' in vdims and 'lon' in vdims:
        refnlat = refdata.sizes['lat']
        refnlon = refdata.sizes['lon']
        if refnlat == 46 and refnlon == 72:
            refres = '4x5'
            refgridtype = 'll'
        elif refnlat == 91 and refnlon == 144:
            refres = '2x2.5'
            refgridtype = 'll'
        elif refnlat/6 == refnlon:
            refres = refnlon
            refgridtype = 'cs'
        else:
            print('ERROR: ref {}x{} grid not defined in gcpy!'.format(refnlat,refnlon))
            return
    else:
        # GCHP data using MAPL v1.0.0+ has dims time, lev, nf, Ydim, and Xdim
        refres = refdata.dims['Xdim']
        refgridtype = 'cs'

    # dev
    vdims = devdata.dims
    if 'lat' in vdims and 'lon' in vdims:
        devnlat = devdata.sizes['lat']
        devnlon = devdata.sizes['lon']
        if devnlat == 46 and devnlon == 72:
            devres = '4x5'
            devgridtype = 'll'
        elif devnlat == 91 and devnlon == 144:
            devres = '2x2.5'
            devgridtype = 'll'
        elif devnlat/6 == devnlon:
            devres = devnlon
            devgridtype = 'cs'
        else:
            print('ERROR: dev {}x{} grid not defined in gcpy!'.format(refnlat,refnlon))
            return
    else:
        devres = devdata.dims['Xdim']
        devgridtype = 'cs'

    ####################################################################
    # Determine comparison grid resolution and type (if not passed)
    ####################################################################

    # If no cmpres is passed then choose highest resolution between ref and dev.
    # If both datasets are cubed sphere then default to 1x1.25 for comparison.
    if cmpres == None:
        if refres == devres and refgridtype == 'll':
            cmpres = refres
            cmpgridtype = 'll'
        elif refgridtype == 'll' and devgridtype == 'll':
            cmpres = min([refres, devres])
            cmpgridtype = 'll'
        elif refgridtype == 'cs' and devgridtype == 'cs':
            cmpres = max([refres, devres])
            cmpgridtype = 'cs'
        else:
            cmpres = '1x1.25'
            cmpgridtype = 'll'
    elif 'x' in cmpres:
        cmpgridtype = 'll'
    else:
        cmpgridtype = 'cs'
        
    # Determine what, if any, need regridding.
    regridref = refres != cmpres
    regriddev = devres != cmpres
    regridany = regridref or regriddev
    
    ####################################################################
    # Make grids (ref, dev, and comparison)
    ####################################################################

    # Ref
    if refgridtype == 'll':
        refgrid = make_grid_LL(refres)
    else:
        [refgrid, regrid_list] = make_grid_CS(refres)

    # Dev
    if devgridtype == 'll':
        devgrid = make_grid_LL(devres)
    else:
        [devgrid, devgrid_list] = make_grid_CS(devres)

    # Comparison    
    if cmpgridtype == 'll':
        cmpgrid = make_grid_LL(cmpres)
    else:
        [cmpgrid, cmpgrid_list] = make_grid_CS(cmpres)
        
    ####################################################################
    # Make regridders, if applicable
    ####################################################################

    if regridref:
        if refgridtype == 'll':
            refregridder = make_regridder_L2L(refres, cmpres, weightsdir=weightsdir, reuse_weights=True)
        else:
            refregridder_list = make_regridder_C2L(refres, cmpres, weightsdir=weightsdir, reuse_weights=True)
    if regriddev:
        if devgridtype == 'll':
            devregridder = make_regridder_L2L(devres, cmpres, weightsdir=weightsdir, reuse_weights=True)
        else:
            devregridder_list = make_regridder_C2L(devres, cmpres, weightsdir=weightsdir, reuse_weights=True)

    ####################################################################
    # Get lat/lon extents, if applicable
    ####################################################################
    
    if refgridtype == 'll':
        [refminlon, refmaxlon] = [min(refgrid['lon_b']), max(refgrid['lon_b'])]
        [refminlat, refmaxlat] = [min(refgrid['lat_b']), max(refgrid['lat_b'])]
    if devgridtype == 'll':
        [devminlon, devmaxlon] = [min(devgrid['lon_b']), max(devgrid['lon_b'])]
        [devminlat, devmaxlat] = [min(devgrid['lat_b']), max(devgrid['lat_b'])]
    if cmpgridtype == 'll':
        [cmpminlon, cmpmaxlon] = [min(cmpgrid['lon_b']), max(cmpgrid['lon_b'])]
        [cmpminlat, cmpmaxlat] = [min(cmpgrid['lat_b']), max(cmpgrid['lat_b'])]

    ####################################################################
    # Create pdf if saving to file
    ####################################################################
    
    if savepdf:
        print('\nCreating {} for {} variables'.format(pdfname,n_var))
        pdf = PdfPages(pdfname)
        
    ####################################################################
    # Loop over variables
    ####################################################################

    print_units_warning = True
    for ivar in range(n_var):
        if savepdf: print('{} '.format(ivar), end='')
        varname = varlist[ivar]
        varndim_ref = refdata[varname].ndim
        varndim_dev = devdata[varname].ndim      

        # If units are mol/mol then convert to ppb
        conc_units = ['mol mol-1 dry','mol/mol','mol mol-1']
        if refdata[varname].units.strip() in conc_units:
            refdata[varname].attrs['units'] = 'ppbv'
            refdata[varname].values = refdata[varname].values * 1e9
        if devdata[varname].units.strip() in conc_units:
            devdata[varname].attrs['units'] = 'ppbv'
            devdata[varname].values = devdata[varname].values * 1e9

        # Binary diagnostic concentrations have units ppbv. Change to ppb.
        if refdata[varname].units.strip() == 'ppbv':
            refdata[varname].attrs['units'] = 'ppb'
        if devdata[varname].units.strip() == 'ppbv':
            devdata[varname].attrs['units'] = 'ppb'

        # Check that units match
        units_ref = refdata[varname].units.strip()
        units_dev = devdata[varname].units.strip()
        if units_ref != units_dev:
            print_units_warning=True
            if print_units_warning:
                print('WARNING: ref and dev concentration units do not match!')
                print('Ref units: {}'.format(units_ref))
                print('Dev units: {}'.format(units_dev))
            if enforce_units:
            # if enforcing units, stop the program if units do not match
               assert units_ref == units_dev, 'Units do not match for {}!'.format(varname)
            else:
               # if not enforcing units, just keep going after only printing warning once 
               print_units_warning = False

               
        ################################################################
        # Slice the data, allowing for the
        # possibility of no time dimension (bpch)
        ################################################################

        # Ref
        vdims = refdata[varname].dims
        if 'time' in vdims and 'lev' in vdims: 
            if flip_ref:
                ds_ref = refdata[varname].isel(time=itime,lev=71-ilev)
            else:
                ds_ref = refdata[varname].isel(time=itime,lev=ilev)
        elif 'lev' in vdims:
            if flip_ref:
                ds_ref = refdata[varname].isel(lev=71-ilev)
            else:
                ds_ref = refdata[varname].isel(lev=ilev)
        elif 'time' in vdims: 
            ds_ref = refdata[varname].isel(time=itime)
        else:
            ds_ref = refdata[varname]

        # Dev
        vdims = devdata[varname].dims
        if 'time' in vdims and 'lev' in vdims: 
            if flip_dev:
                ds_dev = devdata[varname].isel(time=itime,lev=71-ilev)
            else:
                ds_dev = devdata[varname].isel(time=itime,lev=ilev)
        elif 'lev' in vdims:
            if flip_dev:
                ds_dev = devdata[varname].isel(lev=71-ilev)
            else:
                ds_dev = devdata[varname].isel(lev=ilev)
        elif 'time' in vdims: 
            ds_dev = devdata[varname].isel(time=itime)
        else:
            ds_dev = devdata[varname]

        ################################################################
        # Reshape cubed sphere data if using MAPL v1.0.0+
        # TODO: update function to expect data in this format
        ################################################################

        # ref
        vdims = refdata[varname].dims
        if 'nf' in vdims and 'Xdim' in vdims and 'Ydim' in vdims:
            ds_ref = ds_ref.stack(lat=('nf', 'Ydim'))
            ds_ref = ds_ref.rename({'Xdim':'lon'})
            ds_ref = ds_ref.transpose('lat', 'lon')

        # dev
        vdims = devdata[varname].dims
        if 'nf' in vdims and 'Xdim' in vdims and 'Ydim' in vdims:
            ds_dev = ds_dev.stack(lat=('nf', 'Ydim'))
            ds_dev = ds_dev.rename({'Xdim':'lon'})
            ds_dev = ds_dev.transpose('lat', 'lon')
            
        ################################################################
        # Area normalization, if any
        ################################################################

        # if normalizing by area, adjust units to be per m2,
        # and adjust title string
        units = units_ref
        subtitle_extra = ''
        varndim = varndim_ref

        # if regridding then normalization by area may be necessary depending on units. Either pass
        # normalize_by_area=True to normalize all, or include units that should always be normalized
        # by area below. GEOS-Chem Classic output files include area and so do not need to be passed.
        exclude_list = ['WetLossConvFrac','Prod_','Loss_']
        if regridany and ( ( units == 'kg' or units == 'kgC' ) or normalize_by_area ):
            if not any(s in varname for s in exclude_list):
                if 'AREAM2' in refdata.data_vars.keys() and 'AREAM2' in devdata.data_vars.keys():
                    ds_ref.values = ds_ref.values / refdata['AREAM2'].values
                    ds_dev.values = ds_dev.values / devdata['AREAM2'].values
                else:
                    print('ERROR: Variables AREAM2 needed for area normalization missing from dataset')
                    return
                units = '{}/m2'.format(units)
                units_ref = units
                units_dev = units
                subtitle_extra = ', Normalized by Area'

        ###############################################################
        # Get comparison data sets, regridding input slices if needed
        ###############################################################

        # Reshape ref/dev cubed sphere data, if any
        if refgridtype == 'cs':
            ds_ref_reshaped = ds_ref.data.reshape(6,refres,refres)
        if devgridtype == 'cs':
            ds_dev_reshaped = ds_dev.data.reshape(6,devres,devres)

        # Ref
        if regridref:
            if refgridtype == 'll':
                # regrid ll to ll
                ds_ref_cmp = refregridder(ds_ref)
            else:
                # regrid cs to ll
                ds_ref_cmp = np.zeros([cmpgrid['lat'].size, cmpgrid['lon'].size])
                for i in range(6):
                    regridder = refregridder_list[i]
                    ds_ref_cmp += regridder(ds_ref_reshaped[i])
        else:
            ds_ref_cmp = ds_ref

        # Dev
        if regriddev:
            if devgridtype == 'll':
                # regrid ll to ll
                ds_dev_cmp = devregridder(ds_dev)
            else:
                # regrid cs to ll
                ds_dev_cmp = np.zeros([cmpgrid['lat'].size, cmpgrid['lon'].size])
                for i in range(6):
                    regridder = devregridder_list[i]
                    ds_dev_cmp += regridder(ds_dev_reshaped[i])
        else:
            ds_dev_cmp = ds_dev

        # Reshape comparison cubed sphere data, if any
        if cmpgridtype == 'cs':
            ds_ref_cmp_reshaped = ds_ref_cmp.data.reshape(6,cmpres,cmpres)
            ds_dev_cmp_reshaped = ds_dev_cmp.data.reshape(6,cmpres,cmpres)

        ################################################################
        # Get min and max values for use in the colorbars
        ################################################################

        # Ref
        vmin_ref = ds_ref.data.min()
        vmax_ref = ds_ref.data.max()

        # Dev
        vmin_dev = ds_dev.data.min()
        vmax_dev = ds_dev.data.max()

        # Comparison
        if cmpgridtype == 'cs':
            vmin_ref_cmp = ds_ref_cmp.data.min()
            vmax_ref_cmp = ds_ref_cmp.data.max()
            vmin_dev_cmp = ds_dev_cmp.data.min()
            vmax_dev_cmp = ds_dev_cmp.data.max()
            vmin_cmp = np.min([vmin_ref_cmp, vmin_dev_cmp])
            vmax_cmp = np.max([vmax_ref_cmp, vmax_dev_cmp]) 
        else:
            vmin_cmp = np.min([ds_ref_cmp.min(), ds_dev_cmp.min()])
            vmax_cmp = np.max([ds_ref_cmp.max(), ds_dev_cmp.max()])

        # Take min/max across all
        vmin_abs = np.min([vmin_ref, vmin_dev, vmin_cmp])
        vmax_abs = np.max([vmax_ref, vmax_dev, vmax_cmp])
        if match_cbar:
            [vmin, vmax] = [vmin_abs, vmax_abs]

        if verbose:
            print('vmin_ref: {}'.format(vmin_ref))
            print('vmin_dev: {}'.format(vmin_dev))
            print('vmin_ref_cmp: {}'.format(vmin_ref_cmp))
            print('vmin_dev_cmp: {}'.format(vmin_dev_cmp))
            print('vmin_cmp: {}'.format(vmin_cmp))
            print('vmin_abs: {}'.format(vmin_abs))
            print('vmax_ref: {}'.format(vmax_ref))
            print('vmax_dev: {}'.format(vmax_dev))
            print('vmax_ref_cmp: {}'.format(vmax_ref_cmp))
            print('vmax_dev_cmp: {}'.format(vmax_dev_cmp))
            print('vmax_cmp: {}'.format(vmax_cmp))
            print('vmax_abs: {}'.format(vmax_abs))

        ################################################################
        # Create 3x2 figure
        ################################################################
        
        figs, ((ax0, ax1), (ax2, ax3), (ax4, ax5)) = plt.subplots(3, 2, figsize=[12,14], 
                                                      subplot_kw={'projection': crs.PlateCarree()})
        # Give the figure a title
        offset = 0.96
        fontsize=25
        
        if 'lev' in ds_ref.dims and 'lev' in ds_dev.dims:
            if ilev == 0: levstr = 'Surface'
            elif ilev == 22: levstr = '500 hPa'
            else: levstr = 'Level ' +  str(ilev-1)
            figs.suptitle('{}, {}'.format(varname,levstr), fontsize=fontsize, y=offset)
        elif 'lat' in ds_ref.dims and 'lat' in ds_dev.dims and 'lon' in ds_ref.dims and 'lon' in ds_dev.dims:
            figs.suptitle('{}'.format(varname), fontsize=fontsize, y=offset)
        else:
            print('Incorrect dimensions for {}!'.format(varname))   

        ################################################################
        # Set colormap for raw data plots (first row)
        ################################################################

        if use_cmap_RdBu:
            cmap1 = copy.copy(mpl.cm.RdBu_r) # Copy to avoid application of cmap.set_bad used later on
        else:
            cmap1 = WhGrYlRd
            
        ################################################################
        # Subplot (0,0): Ref, plotted on ref input grid
        ################################################################

        # Set colorbar min/max
        if use_cmap_RdBu:
            if vmin_ref == 0 and vmax_ref == 0:
                [vmin0, vmax0] = [vmin_ref, vmax_ref]
            else:
                if not match_cbar:
                    absmax_ref = max([np.abs(vmin_ref), np.abs(vmax_ref)])
                    [vmin0, vmax0] = [-absmax_ref, absmax_ref]
                else:
                    absmax = max([np.abs(vmin_ref), np.abs(vmax_ref),
                                  np.abs(vmin_dev), np.abs(vmax_dev)])
                    [vmin0, vmax0] = [-absmax, absmax]
        else:
            if not match_cbar:
                [vmin0, vmax0] = [vmin_ref, vmax_ref]
            else:
                [vmin0, vmax0] = [vmin_abs, vmax_abs]

        if verbose: print('Subplot (0,0) vmin0, vmax0: {}, {}'.format(
                vmin0, vmax0))

        # Plot data for either lat-lon or cubed-sphere grids.
        ax0.coastlines()
        if refgridtype == 'll':
            # Set top title for lat-lon plot
            ax0.set_title('{} (Ref){}\n{}'.format(
                refstr, subtitle_extra, refres))

            # Check if the reference data is all zeroes
            ref_is_all_zeroes = np.all(ds_ref == 0)

            # Normalize colors for the plot
            norm0 = normalize_colors(vmin0, vmax0,
                                     is_all_zeroes=ref_is_all_zeroes,
                                     log_color_scale=log_color_scale)

            # Plot the data
            plot0 = ax0.imshow(ds_ref, extent=(refminlon, refmaxlon,
                                               refminlat, refmaxlat),
                               cmap=cmap1, norm=norm0)
        else:
            # Set top title for cubed-sphere plot
            ax0.set_title('{} (Ref){}\nc{}'.format(
                refstr, subtitle_extra, refres))

            # Mask the reference data
            masked_refdata = np.ma.masked_where(
                np.abs(refgrid['lon'] - 180) < 2, ds_ref_reshaped)

            # Check if the reference data is all zeroes
            ref_is_all_zeroes = np.all(masked_refdata == 0)

            # Normalize colors (put in range [0..1]) for the plot
            norm0 = normalize_colors(vmin0, vmax0,
                                     is_all_zeroes=ref_is_all_zeroes,
                                     log_color_scale=log_color_scale)

            # Plot each face of the cubed-sphere
            for i in range(6):
                plot0 = ax0.pcolormesh(refgrid['lon_b'][i,:,:],
                                       refgrid['lat_b'][i,:,:],
                                       masked_refdata[i,:,:],
                                       cmap=cmap1, norm=norm0)

        # Define the colorbar for log or linear color scales
        # If all values of Ref = 0, then manually set a tickmark at 0.
        cb = plt.colorbar(plot0, ax=ax0, orientation='horizontal', pad=0.10)
        if ref_is_all_zeroes:
            cb.locator = mpl.ticker.FixedLocator(locs=[0.0])
        else:
            if log_color_scale:
                cb.formatter = mpl.ticker.LogFormatter(base=10)
            else:
                if (vmax0-vmin0) < 0.1 or (vmax0-vmin0) > 100:
                    cb.locator = mpl.ticker.MaxNLocator(nbins=4)
        cb.update_ticks()
        cb.set_label(units_ref)

        ################################################################
        # Subplot (0,1): Dev, plotted on dev input grid
        ################################################################

        # Set colorbar min/max
        if use_cmap_RdBu:
            if vmin_dev == 0 and vmax_dev == 0:
                [vmin1, vmax1] = [vmin_dev, vmax_dev]
            else:
                if not match_cbar:
                    absmax_dev = max([np.abs(vmin_dev), np.abs(vmax_dev)])
                    [vmin1, vmax1] = [-absmax_dev, absmax_dev]
                else:
                    absmax = max([np.abs(vmin_ref), np.abs(vmax_ref),
                                  np.abs(vmin_dev), np.abs(vmax_dev)])
                    [vmin1, vmax1] = [-absmax, absmax]
        else:
            if not match_cbar:
                [vmin1, vmax1] = [vmin_dev, vmax_dev]
            else:
                [vmin1, vmax1] = [vmin_abs, vmax_abs]

        if verbose: print('Subplot (0,1) vmin1, vmax1: {}, {}'.format(vmin1,vmax1))

        # Plot for either lat-lon or cubed-sphere
        ax1.coastlines()
        if devgridtype == 'll':
            # Set top title for lat-lon plot
            ax1.set_title('{} (Dev){}\n{}'.format(
                devstr,subtitle_extra,devres))

            # Check if Dev is all zeroes
            dev_is_all_zeroes = np.all(ds_dev == 0)

            # Normalize colors (put in range [0..1]) for the plot
            norm1 = normalize_colors(vmin1, vmax1,
                                     is_all_zeroes=dev_is_all_zeroes,
                                     log_color_scale=log_color_scale)

            # Plot the data!
            plot1 = ax1.imshow(ds_dev, extent=(devminlon, devmaxlon,
                                               devminlat, devmaxlat),
                               cmap=cmap1, norm=norm1)
        else:
            # Set top title for cubed-sphere plot
            ax1.set_title('{} (Dev){}\nc{}'.format(
                devstr,subtitle_extra,devres))

            # Mask the data
            masked_devdata = np.ma.masked_where(
                np.abs(devgrid['lon'] - 180) < 2, ds_dev_reshaped)

            # Check if Dev is all zereos
            dev_is_all_zeroes = np.all(masked_devdata == 0)

            # Normalize colors (put in range [0..1]) for the plot
            norm1 = normalize_colors(vmin1, vmax1,
                                     is_all_zeroes=dev_is_all_zeroes,
                                     log_color_scale=log_color_scale)

            # Plot each face of the cubed-sphere
            for i in range(6):
                plot1 = ax1.pcolormesh(devgrid['lon_b'][i,:,:],
                                       devgrid['lat_b'][i,:,:],
                                       masked_devdata[i,:,:],
                                       cmap=cmap1, norm=norm1)

        # Define the colorbar for log or linear color scales
        # If all values of Dev = 0, then manually set a tickmark at 0.
        cb = plt.colorbar(plot1, ax=ax1, orientation='horizontal', pad=0.10)
        if dev_is_all_zeroes:
            cb.locator = mpl.ticker.FixedLocator(locs=[0.0])
        else:
            if log_color_scale:
                cb.formatter = mpl.ticker.LogFormatter(base=10)
            else:
                if (vmax1-vmin1) < 0.1 or (vmax1-vmin1) > 100:
                    cb.locator = mpl.ticker.MaxNLocator(nbins=4)
        cb.update_ticks()
        cb.set_label(units_dev)

        ################################################################
        # Calculate difference
        ################################################################

        if cmpgridtype == 'll':
            absdiff = np.array(ds_dev_cmp) - np.array(ds_ref_cmp)
        else:
            absdiff_raw = ds_dev_cmp_reshaped - ds_ref_cmp_reshaped
            absdiff = np.ma.masked_where(np.abs(cmpgrid['lon'] - 180) < 2, absdiff_raw)
        diffabsmax = max([np.abs(np.nanmin(absdiff)), np.abs(np.nanmax(absdiff))])

        ################################################################
        # Use gray for NaNs if plotting on lat-lon grid 
        # (has strange effect for cubed-sphere)
        ################################################################

        cmap_nongray = copy.copy(mpl.cm.RdBu_r)
        cmap_gray = copy.copy(mpl.cm.RdBu_r)
        cmap_gray.set_bad(color='gray')
            
        ################################################################
        # Subplot (1,0): Difference, dynamic range
        ################################################################

        [vmin, vmax] = [-diffabsmax, diffabsmax]
        if verbose: print('Subplot (1,0) vmin, vmax: {}, {}'.format(vmin, vmax))

        ax2.coastlines()
        if cmpgridtype == 'll':
            plot2 = ax2.imshow(absdiff, extent=(cmpminlon, cmpmaxlon,
                                                cmpminlat, cmpmaxlat),
                               cmap=cmap_gray, vmin=vmin, vmax=vmax)
        else:
            for i in range(6):
                plot2 = ax2.pcolormesh(cmpgrid['lon_b'][i,:,:],
                                       cmpgrid['lat_b'][i,:,:],
                                       absdiff[i,:,:],
                                       cmap=cmap_nongray,
                                       vmin=vmin, vmax=vmax)
        if regridany:
            ax2.set_title('Difference ({})\nDev - Ref, Dynamic Range'.format(
                cmpres))
        else:
            ax2.set_title('Difference\nDev - Ref, Dynamic Range')

        # Define the colorbar for the plot
        # If all values of absdiff = 0, manually set a tickmark at 0.
        cb = plt.colorbar(plot2, ax=ax2, orientation='horizontal', pad=0.10)
        if np.all(absdiff==0): 
            cb.locator = mpl.ticker.FixedLocator(locs=[0.0])
        else:
            if (vmax-vmin) < 0.1 or (vmax-vmin) > 100:
                cb.locator = mpl.ticker.MaxNLocator(nbins=4)
        cb.update_ticks()
        cb.set_label(units)

        ################################################################
        # Subplot (1,1): Difference, restricted range
        ################################################################

        # placeholder: use 5 and 95 percentiles as bounds
        [pct5, pct95] = [np.percentile(absdiff,5), np.percentile(absdiff, 95)] 

        abspctmax = np.max([np.abs(pct5),np.abs(pct95)])
        [vmin,vmax] = [-abspctmax, abspctmax]
        if verbose: print('Subplot (1,1) vmin, vmax: {}, {}'.format(vmin, vmax))

        ax3.coastlines()
        if cmpgridtype == 'll':
            plot3 = ax3.imshow(absdiff, extent=(cmpminlon, cmpmaxlon,
                                                cmpminlat, cmpmaxlat),
                               cmap=cmap_gray, vmin=vmin, vmax=vmax)
        else:
            for i in range(6):
                plot3 = ax3.pcolormesh(cmpgrid['lon_b'][i,:,:],
                                       cmpgrid['lat_b'][i,:,:],
                                       absdiff[i,:,:],
                                       cmap=cmap_nongray,
                                       vmin=vmin, vmax=vmax)
        if regridany:
            ax3.set_title('Difference ({})\nDev - Ref, Restricted Range [5%,95%]'.format(cmpres))
        else:
            ax3.set_title('Difference\nDev - Ref, Restricted Range [5%,95%]')

        # Define the colorbar for the plot.
        # If all values of absdiff = 0, then manually set a tick at 0.
        cb = plt.colorbar(plot3, ax=ax3, orientation='horizontal', pad=0.10)
        if np.all(absdiff==0):
            cb.locator = mpl.ticker.FixedLocator(locs=[0.0])
        else:
            if (vmax-vmin) < 0.1 or (vmax-vmin) > 100:
                cb.locator = mpl.ticker.MaxNLocator(nbins=4)
        cb.update_ticks()
        cb.set_label(units)

        ################################################################
        # Calculate fractional difference, set divides by zero to Nan
        ################################################################
        
        if cmpgridtype == 'll':
            fracdiff = (np.array(ds_dev_cmp) - np.array(ds_ref_cmp)) / np.array(ds_ref_cmp)
        else:
            fracdiff_raw = (ds_dev_cmp_reshaped - ds_ref_cmp_reshaped) / ds_ref_cmp_reshaped
            fracdiff = np.ma.masked_where(np.abs(cmpgrid['lon'] - 180) < 2, fracdiff_raw)

        fracdiff = np.where(fracdiff==np.inf, np.nan, fracdiff)
        fracdiffabsmax = max([np.abs(np.nanmin(fracdiff)), np.abs(np.nanmax(fracdiff))])

        ################################################################
        # Subplot (2,0): Fractional Difference, full dynamic range
        ################################################################
        
        [vmin, vmax] = [-fracdiffabsmax, fracdiffabsmax]
        if verbose: print('Subplot (2,0) vmin, vmax: {}, {}'.format(vmin, vmax))
  
        ax4.coastlines()
        if cmpgridtype == 'll':
            plot4 = ax4.imshow(fracdiff, extent=(cmpminlon, cmpmaxlon,
                                                 cmpminlat, cmpmaxlat),
                               vmin=vmin, vmax=vmax, cmap=cmap_gray)
        else:
            for i in range(6):
                plot4 = ax4.pcolormesh(cmpgrid['lon_b'][i,:,:],
                                       cmpgrid['lat_b'][i,:,:],
                                       fracdiff[i,:,:],
                                       cmap=cmap_nongray,
                                       vmin=vmin, vmax=vmax)
        if regridany:
            ax4.set_title('Fractional Difference ({})\n(Dev-Ref)/Ref, Dynamic Range'.format(cmpres)) 
        else:
            ax4.set_title('Fractional Difference\n(Dev-Ref)/Ref, Dynamic Range')

        # Define the colorbar for the plot.
        # If all values of absdiff = 0, then manually set a tick at 0.
        # Use absdiff in the test for all zeroes, because fracdiff may
        # have NaN's due to div by zero when both Dev and Ref are zero.
        cb = plt.colorbar(plot4, ax=ax4, orientation='horizontal', pad=0.10)
        if np.all(absdiff == 0):
            cb.locator = mpl.ticker.FixedLocator(locs=[0.0])
        else:
            if (vmax-vmin) < 0.1 or (vmax-vmin) > 100:
                cb.locator = mpl.ticker.MaxNLocator(nbins=4)
        cb.update_ticks()
        cb.set_label('unitless')

        ################################################################
        # Subplot (2,1): Fractional Difference, restricted range
        ################################################################
        
        [vmin, vmax] = [-2, 2]
        if verbose: print('Subplot (2,1) vmin, vmax: {}, {}'.format(vmin, vmax))

        ax5.coastlines()
        if cmpgridtype == 'll':
            plot5 = ax5.imshow(fracdiff, extent=(cmpminlon, cmpmaxlon,
                                                 cmpminlat, cmpmaxlat),
                               cmap=cmap_gray, vmin=vmin, vmax=vmax)
        else:
            for i in range(6):
                plot5 = ax5.pcolormesh(cmpgrid['lon_b'][i,:,:],
                                       cmpgrid['lat_b'][i,:,:],
                                       fracdiff[i,:,:],
                                       cmap=cmap_nongray,
                                       vmin=vmin, vmax=vmax)
        if regridany:
            ax5.set_title('Fractional Difference ({})\n(Dev-Ref)/Ref, Fixed Range'.format(cmpres))
        else:
            ax5.set_title('Fractional Difference\n(Dev-Ref)/Ref, Fixed Range') 

        # Define the colorbar for the plot.
        # If all values of absdiff = 0, then manually set a tick at 0.
        # Use absdiff in the test for all zeroes, because fracdiff may
        # have NaN's due to div by zero when both Dev and Ref are zero.
        cb = plt.colorbar(plot5, ax=ax5, orientation='horizontal', pad=0.10)
        if np.all(absdiff == 0):
            cb.locator = mpl.ticker.FixedLocator(locs=[0.0])
        cb.update_ticks()
        cb.set_label('unitless')

        ################################################################
        # Update the list of variables with significant differences.
        # Criterion: abs(max(fracdiff)) > 0.1
        # Do not include NaNs in the criterion, because these indicate
        # places where fracdiff could not be computed (div-by-zero).
        ################################################################
        if np.abs(np.nanmax(fracdiff)) > 0.1:
            sigdiff_list.append(varname)

        ################################################################
        # Add this page of 6-panel plots to the PDF file
        ################################################################
        if savepdf:    
            pdf.savefig(figs)
            plt.close(figs)

    ####################################################################
    # Finish
    ####################################################################
    if savepdf:
        pdf.close()

def compare_zonal_mean(refdata, refstr, devdata, devstr, varlist=None,
                       itime=0, weightsdir=None, pdfname='', cmpres=None,
                       match_cbar=True, pres_range=[0,2000],
                       normalize_by_area=False, enforce_units=True,
                       flip_ref=False, flip_dev=False, use_cmap_RdBu=False,
                       verbose=False, log_color_scale=False,
                       sigdiff_list=[]):

    '''
    Create single-level 3x2 comparison map plots for variables common in two xarray
    datasets. Optionally save to PDF. 

    Args:
        refdata : xarray dataset
            Dataset used as reference in comparison

        refstr  : str
            String description for reference data to be used in plots
     
        devdata : xarray dataset
            Dataset used as development in comparison

        devstr  : str
            String description for development data to be used in plots
 
    Keyword Args (optional):
        varlist : list of strings
            List of xarray dataset variable names to make plots for
            Default value: None (will compare all common 3D variables)

        itime : integer
            Dataset time dimension index using 0-based system
            Default value: 0

        weightsdir : str
            Directory path for storing regridding weights
            Default value: None (will create/store weights in current directory)

        pdfname : str
            File path to save plots as PDF
            Default value: Empty string (will not create PDF)

        cmpres : str
            String description of grid resolution at which to compare datasets
            Default value: None (will compare at highest resolution of ref and dev)

        match_cbar : boolean
            Logical indicating whether to use same colorbar bounds across plots
            Default value: True

        pres_range : list of two integers
            Pressure range of levels to plot [hPa]. Vertical axis will span outer
            pressure edges of levels that contain pres_range endpoints.
            Default value: [0,2000]

        normalize_by_area : boolean
            Logical indicating whether to normalize raw data by grid area
            Default value: False

        enforce_units : boolean
            Logical to force an error if reference and development variables units differ
            Default value: True

        flip_ref : boolean
            Logical to flip the vertical dimension of reference dataset 3D variables
            Default value: False

        flip_dev : boolean
            Logical to flip the vertical dimension of development dataset 3D variables
            Default value: False

        use_cmap_RdBu : boolean
            Logical to used a blue-white-red colormap for plotting raw reference and
            development datasets.
            Default value: False

        verbose : logical
            Logical to enable informative prints
            Default value: False

        log_color_scale: boolean         
            Logical to enable plotting data (not diffs) on a log color scale.
            Default value: False

        sigdiff_list: list of str
            Returns a list of all quantities having significant differences.
            The criteria is: |max(fractional difference)| > 0.1
            Default value: []

    Returns:
        Nothing

    Example:
        >>> import matplotlib.pyplot as plt
        >>> import xarray as xr
        >>> from gcpy import benchmark
        >>> refds = xr.open_dataset('path/to/ref.nc4')
        >>> devds = xr.open_dataset('path/to/dev.nc4')
        >>> varlist = ['SpeciesConc_O3', 'SpeciesConc_NO']
        >>> benchmark.compare_zonal_mean( refds, '12.3.2', devds, 'bug fix', varlist=varlist )
        >>> plt.show()
    '''

    # TODO: refactor this function and single level plot function. There is a lot of overlap and
    # repeated code that could be abstracted.

    if not isinstance(refdata, xr.Dataset):
        raise ValueError('The refdata argument must be an xarray Dataset!')

    if not isinstance(devdata, xr.Dataset):
        raise ValueError('The devdata argument must be an xarray Dataset!')

    # If no varlist is passed, plot all 3D variables in the dataset
    if varlist == None:
        [commonvars, commonvarsOther, commonvars2D, commonvars3D] = core.compare_varnames(refdata, devdata)
        varlist = commonvars3D
        print('Plotting all 3D variables')
    n_var = len(varlist)

    # Exit out if there are no 3D variables
    if not n_var: 
        print('WARNING: no 3D variables to plot zonal mean for!')
        return

    # If no weightsdir is passed, set to current directory in case it is needed
    if weightsdir == None:
        weightsdir = '.'

    # If no pdf name passed, then do not save to PDF
    savepdf = True
    if pdfname == '':
        savepdf = False

    # Get mid-point pressure and edge pressures for this grid (assume 72-level)
    pmid = GEOS_72L_grid.p_mid()
    pedge = GEOS_72L_grid.p_edge()

    # Get indexes of pressure subrange (full range is default)
    pedge_ind = np.where((pedge <= np.max(pres_range)) & (pedge >= np.min(pres_range)))
    pedge_ind = pedge_ind[0]
    # Pad edges if subset does not include surface or TOA so data spans entire subrange
    if min(pedge_ind) != 0:
        pedge_ind = np.append(min(pedge_ind)-1, pedge_ind)
    if max(pedge_ind) != 72:
        pedge_ind = np.append(pedge_ind, max(pedge_ind)+1)
    # pmid indexes do not include last pedge index
    pmid_ind = pedge_ind[:-1]
    nlev = len(pmid_ind)
        
    # Convert levels to pressures in ref and dev data
    if refdata.sizes['lev'] == 72:
        refdata['lev'] = pmid
    elif refdata.sizes['lev'] == 73:
        refdata['lev'] = pedge
    else:
        print('ERROR: compare_zonal_mean implemented for 72 or 73 levels only. Other values found in ref.')
        return
    refdata['lev'].attrs['units'] = 'hPa'
    refdata['lev'].attrs['long_name'] = 'level pressure'

    if devdata.sizes['lev'] == 72:
        devdata['lev'] = pmid
    elif devdata.sizes['lev'] == 73:
        devdata['lev'] = pedge
    else:
        print('ERROR: compare_zonal_mean implemented for 72 or 73 levels only. Other value found in dev.')
        return
    devdata['lev'].attrs['units'] = 'hPa'
    devdata['lev'].attrs['long_name'] = 'level pressure'

    ####################################################################
    # Reduce pressure range if reduced range passed as input
    ####################################################################

    refdata = refdata.isel(lev=pmid_ind)
    devdata = devdata.isel(lev=pmid_ind)
 
    ####################################################################
    # Determine input grid resolutions and types
    ####################################################################
    # GCC output and GCHP output using pre-v1.0.0 MAPL have lat and lon dims

    # ref
    vdims = refdata.dims
    if 'lat' in vdims and 'lon' in vdims:
        refnlat = refdata.sizes['lat']
        refnlon = refdata.sizes['lon']
        if refnlat == 46 and refnlon == 72:
            refres = '4x5'
            refgridtype = 'll'
        elif refnlat == 91 and refnlon == 144:
            refres = '2x2.5'
            refgridtype = 'll'
        elif refnlat/6 == refnlon:
            refres = refnlon
            refgridtype = 'cs'
        else:
            print('ERROR: ref {}x{} grid not defined in gcpy!'.format(refnlat,refnlon))
            return
    else:
        # GCHP data using MAPL v1.0.0+ has dims time, lev, nf, Ydim, and Xdim
        refres = refdata.dims['Xdim']
        refgridtype = 'cs'

    # dev
    vdims = devdata.dims
    if 'lat' in vdims and 'lon' in vdims:
        devnlat = devdata.sizes['lat']
        devnlon = devdata.sizes['lon']
        if devnlat == 46 and devnlon == 72:
            devres = '4x5'
            devgridtype = 'll'
        elif devnlat == 91 and devnlon == 144:
            devres = '2x2.5'
            devgridtype = 'll'
        elif devnlat/6 == devnlon:
            devres = devnlon
            devgridtype = 'cs'
        else:
            print('ERROR: dev {}x{} grid not defined in gcpy!'.format(refnlat,refnlon))
            return
    else:
        devres = devdata.dims['Xdim']
        devgridtype = 'cs'
    
    ####################################################################
    # Determine comparison grid resolution (if not passed)
    ####################################################################

    # If no cmpres is passed then choose highest resolution between ref and dev.
    # If both datasets are cubed sphere then default to 1x1.25 for comparison.
    cmpgridtype = 'll'
    if cmpres == None:
        if refres == devres and refgridtype == 'll':
            cmpres = refres
        elif refgridtype == 'll' and devgridtype == 'll':
            cmpres = min([refres, devres])
        else:
            cmpres = '1x1.25'
        
    # Determine what, if any, need regridding.
    regridref = refres != cmpres
    regriddev = devres != cmpres
    regridany = regridref or regriddev

    ####################################################################
    # Make grids (ref, dev, and comparison)
    ####################################################################

    # Ref
    if refgridtype == 'll':
        refgrid = make_grid_LL(refres)
    else:
        [refgrid, regrid_list] = make_grid_CS(refres)

    # Dev
    if devgridtype == 'll':
        devgrid = make_grid_LL(devres)
    else:
        [devgrid, devgrid_list] = make_grid_CS(devres)

    # Comparison    
    cmpgrid = make_grid_LL(cmpres)

    ####################################################################
    # Make regridders, if applicable
    ####################################################################

    if regridref:
        if refgridtype == 'll':
            refregridder = make_regridder_L2L(refres, cmpres, weightsdir=weightsdir, reuse_weights=True)
        else:
            refregridder_list = make_regridder_C2L(refres, cmpres, weightsdir=weightsdir, reuse_weights=True)
    if regriddev:
        if devgridtype == 'll':
            devregridder = make_regridder_L2L(devres, cmpres, weightsdir=weightsdir, reuse_weights=True)
        else:
            devregridder_list = make_regridder_C2L(devres, cmpres, weightsdir=weightsdir, reuse_weights=True)
            
    ####################################################################
    # Create pdf, if savepdf is passed as True
    ####################################################################

    # Universal plot setup
    xtick_positions = np.arange(-90,91,30)
    xticklabels = ['{}$\degree$'.format(x) for x in xtick_positions]    

    if savepdf:
        print('\nCreating {} for {} variables'.format(pdfname,n_var))
        pdf = PdfPages(pdfname)

    ####################################################################
    # Loop over variables
    ####################################################################

    # Loop over variables
    print_units_warning = True
    for ivar in range(n_var):
        if savepdf: print('{} '.format(ivar), end='')
        varname = varlist[ivar]
        varndim_ref = refdata[varname].ndim
        varndim_dev = devdata[varname].ndim

        # If units are mol/mol then convert to ppb
        conc_units = ['mol mol-1 dry','mol/mol','mol mol-1']
        if refdata[varname].units.strip() in conc_units:
            refdata[varname].attrs['units'] = 'ppbv'
            refdata[varname].values = refdata[varname].values * 1e9
        if devdata[varname].units.strip() in conc_units:
            devdata[varname].attrs['units'] = 'ppbv'
            devdata[varname].values = devdata[varname].values * 1e9

        # Binary diagnostic concentrations have units ppbv. Change to ppb.
        if refdata[varname].units.strip() == 'ppbv':
            refdata[varname].attrs['units'] = 'ppb'
        if devdata[varname].units.strip() == 'ppbv':
            devdata[varname].attrs['units'] = 'ppb'

        # Check that units match
        units_ref = refdata[varname].units.strip()
        units_dev = devdata[varname].units.strip()
        if units_ref != units_dev:
            if print_units_warning:
                print('WARNING: ref and dev concentration units do not match!')
                print('Ref units: {}'.format(units_ref))
                print('Dev units: {}'.format(units_dev))
            if enforce_units:
            # if enforcing units, stop the program if units do not match
               assert units_ref == units_dev, 'Units do not match for {}!'.format(varname)
            else:
               # if not enforcing units, just keep going after only printing warning once 
               print_units_warning = False

        ###############################################################
        # Slice the data, allowing for the
        # possibility of no time dimension (bpch)
        ###############################################################

        # Ref
        vdims = refdata[varname].dims
        if 'time' in vdims:
            ds_ref = refdata[varname].isel(time=itime)
        else:
            ds_ref = refdata[varname]

        # Dev
        vdims = devdata[varname].dims      
        if 'time' in vdims:
            ds_dev = devdata[varname].isel(time=itime)
        else:
            ds_dev = devdata[varname]

        ################################################################
        # Reshape cubed sphere data if using MAPL v1.0.0+
        # TODO: update function to expect data in this format
        ################################################################

        # ref
        vdims = refdata[varname].dims
        if 'nf' in vdims and 'Xdim' in vdims and 'Ydim' in vdims:
            ds_ref = ds_ref.stack(lat=('nf', 'Ydim'))
            ds_ref = ds_ref.rename({'Xdim':'lon'})
            ds_ref = ds_ref.transpose('lev', 'lat', 'lon')

        # dev
        vdims = devdata[varname].dims
        if 'nf' in vdims and 'Xdim' in vdims and 'Ydim' in vdims:
            ds_dev = ds_dev.stack(lat=('nf', 'Ydim'))
            ds_dev = ds_dev.rename({'Xdim':'lon'})
            ds_dev = ds_dev.transpose('lev', 'lat', 'lon')

        ###############################################################
        # Area normalization, if any
        ###############################################################
        
        # if normalizing by area, adjust units to be per m2, and adjust title string
        units = units_ref
        varndim = varndim_ref
        subtitle_extra = ''
        
        # if normalizing by area, transform on the native grid and adjust units and subtitle string
        exclude_list = ['WetLossConvFrac','Prod_','Loss_']
        if normalize_by_area and not any(s in varname for s in exclude_list):
            ds_ref.values = ds_ref.values / refdata['AREAM2'].values[np.newaxis,:,:]
            ds_dev.values = ds_dev.values / devdata['AREAM2'].values[np.newaxis,:,:]
            units = '{} m-2'.format(units)
            subtitle_extra = ', Normalized by Area'

        ###############################################################
        # Get comparison data sets, regridding input slices if needed
        ###############################################################

        # Flip in the veritical if applicable
        if flip_ref:
            ds_ref.data = ds_ref.data[::-1,:,:]
        if flip_dev:
            ds_ref.data = ds_ref.data[::-1,:,:]
            
        # Reshape ref/dev cubed sphere data, if any
        if refgridtype == 'cs':
            ds_ref_reshaped = ds_ref.data.reshape(nlev,6,refres,refres).swapaxes(0,1)
        if devgridtype == 'cs':
            ds_dev_reshaped = ds_dev.data.reshape(nlev,6,devres,devres).swapaxes(0,1)

        # Ref
        if regridref:
            if refgridtype == 'll':
                # regrid ll to ll
                ds_ref_cmp = refregridder(ds_ref)
            else:
                # regrid cs to ll
                ds_ref_cmp = np.zeros([nlev, cmpgrid['lat'].size, cmpgrid['lon'].size])
                for i in range(6):
                    regridder = refregridder_list[i]
                    ds_ref_cmp += regridder(ds_ref_reshaped[i])
        else:
            ds_ref_cmp = ds_ref

        # Dev
        if regriddev:
            if devgridtype == 'll':
                # regrid ll to ll
                ds_dev_cmp = devregridder(ds_dev)
            else:
                # regrid cs to ll
                ds_dev_cmp = np.zeros([nlev, cmpgrid['lat'].size, cmpgrid['lon'].size])
                for i in range(6):
                    regridder = devregridder_list[i]
                    ds_dev_cmp += regridder(ds_dev_reshaped[i])
        else:
            ds_dev_cmp = ds_dev

        ###############################################################
        # Calculate zonal mean
        ###############################################################
        
        # Ref
        if refgridtype == 'll':
            zm_ref = ds_ref.mean(dim='lon')
        else:
            zm_ref = ds_ref_cmp.mean(axis=2)

        # Dev
        if devgridtype == 'll':
            zm_dev = ds_dev.mean(dim='lon')
        else:
            zm_dev = ds_dev_cmp.mean(axis=2)

        # Comparison
        zm_dev_cmp = ds_dev_cmp.mean(axis=2)
        zm_ref_cmp = ds_ref_cmp.mean(axis=2)

        ###############################################################
        # Get min and max values for use in the colorbars
        ###############################################################

        # Ref
        vmin_ref = zm_ref.min()
        vmax_ref = zm_ref.max()

        # Dev
        vmin_dev = zm_dev.min()
        vmax_dev = zm_dev.max()

        # Comparison
        vmin_cmp = np.min([zm_ref_cmp.min(), zm_dev_cmp.min()])
        vmax_cmp = np.max([zm_ref_cmp.max(), zm_dev_cmp.max()])

        # Take min/max across all
        vmin_abs = np.min([vmin_ref, vmin_dev, vmin_cmp])
        vmax_abs = np.max([vmax_ref, vmax_dev, vmax_cmp])
        if match_cbar: [vmin, vmax] = [vmin_abs, vmax_abs]

        if verbose:
            print('vmin_ref: {}'.format(vmin_ref))
            print('vmin_dev: {}'.format(vmin_dev))
            print('vmin_cmp: {}'.format(vmin_cmp))
            print('vmin_abs: {}'.format(vmin_abs))
            print('vmax_ref: {}'.format(vmax_ref))
            print('vmax_dev: {}'.format(vmax_dev))
            print('vmax_cmp: {}'.format(vmax_cmp))
            print('vmax_abs: {}'.format(vmax_abs))

        ###############################################################
        # Create 3x2 figure
        ###############################################################
        
        figs, ((ax0, ax1), (ax2, ax3), (ax4, ax5)) = plt.subplots(3, 2, figsize=[12,15.3])
        offset = 0.96
        fontsize=25
        figs.suptitle('{}, Zonal Mean'.format(varname), fontsize=fontsize, y=offset)

        ###############################################################
        # Set colormap for raw data plots (first row)
        ###############################################################

        if use_cmap_RdBu:
            cmap1 = 'RdBu_r' # do not need to worry about cmap.set_bad since always lat-lon
        else:
            cmap1 = WhGrYlRd

        ###############################################################
        # Subplot (0,0): Ref
        ###############################################################

        # Set colorbar min/max
        if use_cmap_RdBu:
            if vmin_ref == 0 and vmax_ref:
                [vmin0, vmax0] = [vmin_ref, vmax_ref]
            else:
                if not match_cbar:
                    absmax_ref = max([np.abs(vmin_ref), np.abs(vmax_ref)])
                    [vmin0, vmax0] = [-absmax_ref, absmax_dev]
                else:
                    absmax = max([np.abs(vmin_ref), np.abs(vmax_ref),
                                  np.abs(vmin_dev), np.abs(vmax_dev)])
                    [vmin0, vmax0] = [-absmax, absmax]
        else:
            if not match_cbar:
                [vmin0, vmax0] = [vmin_ref, vmax_ref]
            else:
                [vmin0, vmax0] = [vmin_abs, vmax_abs]
        if verbose: print('Subplot (0,0) vmin0, vmax0: {}, {}'.format(vmin0, vmax0))

        # Check if all elements of Ref are zero
        ref_is_all_zeroes = np.all(zm_ref == 0)

        # Normalize colors (put in range [0..1]) for the plot
        norm0 = normalize_colors(vmin0, vmax0,
                                 is_all_zeroes=ref_is_all_zeroes,
                                 log_color_scale=log_color_scale)

        # Plot data for either lat-lon or cubed-sphere grids
        if refgridtype == 'll':
            plot0 = ax0.pcolormesh(refgrid['lat_b'],
                                   pedge[pedge_ind],
                                   zm_ref,
                                   cmap=cmap1, norm=norm0)
            ax0.set_title('{} (Ref){}\n{}'.format(
                refstr, subtitle_extra, refres))
        else:
            plot0 = ax0.pcolormesh(cmpgrid['lat_b'],
                                   pedge[pedge_ind],
                                   zm_ref,
                                   cmap=cmap1, norm=norm0)
            ax0.set_title('{} (Ref){}\n{} regridded from c{}'.format(
                refstr, subtitle_extra, cmpres, refres))
        ax0.invert_yaxis()
        ax0.set_ylabel('Pressure (hPa)')
        ax0.set_aspect('auto')
        ax0.set_xticks(xtick_positions)
        ax0.set_xticklabels(xticklabels)

        # Define the colorbar for log or linear color scales.
        # If all values of zm_ref = 0, then manually set a tick at zero.
        cb = plt.colorbar(plot0, ax=ax0, orientation='horizontal', pad=0.10)
        if ref_is_all_zeroes:
            cb.locator = mpl.ticker.FixedLocator(locs=[0.0])
        else:
            if log_color_scale:
                cb.formatter = mpl.ticker.LogFormatter(base=10)
            else:
                if (vmax-vmin) < 0.001 or (vmax-vmin) > 1000:
                    cb.locator = mpl.ticker.MaxNLocator(nbins=4)
        cb.update_ticks()
        cb.set_label(units)

        ###############################################################
        # Subplot (0,1): Dev
        ###############################################################

        # Set colorbar min/max
        if use_cmap_RdBu:
            if vmin_dev == 0 and vmax_dev == 0:
                [vmin1, vmax1] = [vmin_dev, vmax_dev]
            else:
                if not match_cbar:
                    absmax_dev = max([np.abs(vmin_dev), np.abs(vmax_dev)])
                    [vmin1, vmax1] = [-absmax_dev, absmax_dev]
                else:
                    absmax = max([np.abs(vmin_ref), np.abs(vmax_ref),
                                  np.abs(vmin_dev), np.abs(vmax_dev)])
                    [vmin1, vmax1] = [-absmax, absmax]
        else:
            if not match_cbar:
                [vmin1, vmax1] = [vmin_dev, vmax_dev]
            else:
                [vmin1, vmax1] = [vmin_abs, vmax_abs]
        if verbose: print('Subplot (0,1) vmin1, vmax1: {}, {}'.format(vmin1, vmax1))

        # Check if all elements of Ref are zero
        dev_is_all_zeroes = np.all(zm_dev == 0)

        # Normalize colors (put in range [0..1]) for the plot
        norm1 = normalize_colors(vmin1, vmax1,
                                 is_all_zeroes=dev_is_all_zeroes,
                                 log_color_scale=log_color_scale)

        # Plot data for either lat-lon or cubed-sphere grids.
        if devgridtype == 'll':
            plot1 = ax1.pcolormesh(devgrid['lat_b'], pedge[pedge_ind], 
                                   zm_dev, cmap=cmap1, norm=norm1)
            ax1.set_title('{} (Dev){}\n{}'.format(
                devstr, subtitle_extra, devres ))
        else:
            plot1 = ax1.pcolormesh(cmpgrid['lat_b'], pedge[pedge_ind], 
                                   zm_dev, cmap=cmap1, norm=norm1)
            ax1.set_title('{} (Dev){}\n{} regridded from c{}'.format(
                devstr, subtitle_extra, cmpres, devres))
        ax1.invert_yaxis()
        ax1.set_ylabel('Pressure (hPa)')
        ax1.set_aspect('auto')
        ax1.set_xticks(xtick_positions)
        ax1.set_xticklabels(xticklabels)

        # Define the colorbar for log or linear color scales.
        # If all values of zm_dev = 0, then manually set a tick at 0
        cb = plt.colorbar(plot1, ax=ax1, orientation='horizontal', pad=0.10)
        if dev_is_all_zeroes:
            cb.locator = mpl.ticker.FixedLocator(locs=[0.0])
        else:
            if log_color_scale:
                cb.formatter = mpl.ticker.LogFormatter(base=10)
            else:
                if (vmax-vmin) < 0.001 or (vmax-vmin) > 1000:
                    cb.locator = mpl.ticker.MaxNLocator(nbins=4)
        cb.update_ticks()
        cb.set_label(units)

        ################################################################
        # Configure colorbar for difference plots, use gray for NaNs
        ################################################################

        cmap_plot = copy.copy(mpl.cm.RdBu_r)
        cmap_plot.set_bad(color='gray')

        ################################################################
        # Calculate zonal mean difference
        ################################################################

        zm_diff = np.array(zm_dev_cmp) - np.array(zm_ref_cmp)

        ################################################################
        # Subplot (1,0): Difference, dynamic range
        ################################################################

        diffabsmax = max([np.abs(zm_diff.min()), np.abs(zm_diff.max())])
        [vmin, vmax] = [-diffabsmax, diffabsmax]
        if verbose: print('Subplot (1,0) vmin, vmax: {}, {}'.format(vmin, vmax))

        plot2 = ax2.pcolormesh(cmpgrid['lat_b'], pedge[pedge_ind],
                               zm_diff, cmap=cmap_plot,
                               vmin=vmin, vmax=vmax)
        ax2.invert_yaxis()
        if regridany:
            ax2.set_title('Difference ({})\nDev - Ref, Dynamic Range'.format(
                cmpres))
        else:
            ax2.set_title('Difference\nDev - Ref, Dynamic Range')
        ax2.set_ylabel('Pressure (hPa)')
        ax2.set_aspect('auto')
        ax2.set_xticks(xtick_positions)
        ax2.set_xticklabels(xticklabels)

        # Define the colorbar for the plot.
        # If all values of zm_diff = 0, then manually set a tick at 0
        cb = plt.colorbar(plot2, ax=ax2, orientation='horizontal', pad=0.10)
        if np.all(zm_diff==0): 
            cb.locator = mpl.ticker.FixedLocator(locs=[0.0])
        else:
            if (vmax-vmin) < 0.001 or (vmax-vmin) > 1000:
                cb.locator = mpl.ticker.MaxNLocator(nbins=4)
        cb.update_ticks()
        cb.set_label(units)

        ################################################################
        # Subplot (1,1): Difference, restricted range
        ################################################################

         # placeholder: use 5 and 95 percentiles as bounds
        [pct5, pct95] = [np.percentile(zm_diff,5), np.percentile(zm_diff, 95)]
        abspctmax = np.max([np.abs(pct5),np.abs(pct95)])
        [vmin,vmax] = [-abspctmax, abspctmax]
        if verbose: print('Subplot (1,1) vmin, vmax: {}, {}'.format(vmin, vmax))
  
        plot3 = ax3.pcolormesh(cmpgrid['lat_b'], pedge[pedge_ind],
                               zm_diff, cmap=cmap_plot, vmin=vmin, vmax=vmax)
        ax3.invert_yaxis()
        if regridany:
            ax3.set_title('Difference ({})\nDev - Ref, Restricted Range [5%,95%]'.format(cmpres))
        else:
            ax3.set_title('Difference\nDev - Ref, Restriced Range [5%,95%]')
        ax3.set_ylabel('Pressure (hPa)')
        ax3.set_aspect('auto')
        ax3.set_xticks(xtick_positions)
        ax3.set_xticklabels(xticklabels)

        # Define the colorbar for the plot.
        # If all values of zm_diff = 0, then manually set a tick at 0
        cb = plt.colorbar(plot3, ax=ax3, orientation='horizontal', pad=0.10)
        if np.all(zm_diff==0): 
            cb.locator = mpl.ticker.FixedLocator(locs=[0.0])
        else:
            if (vmax-vmin) < 0.001 or (vmax-vmin) > 1000:
                cb.locator = mpl.ticker.MaxNLocator(nbins=4)
        cb.update_ticks()
        cb.set_label(units)

        ################################################################
        # Calculate fractional difference, set divides by zero to Nan
        ################################################################

        zm_fracdiff = (np.array(zm_dev_cmp) - np.array(zm_ref_cmp)) / np.array(zm_ref_cmp)
        zm_fracdiff = np.where(zm_fracdiff==np.inf, np.nan, zm_fracdiff)

        ################################################################
        # Subplot (2,0): Fractional Difference, dynamic range
        ################################################################
        
        # Exclude NaN's in the fracdiffabsmax
        fracdiffabsmax = np.max([np.abs(np.nanmin(zm_fracdiff)), np.abs(np.nanmax(zm_fracdiff))])
        if np.all(zm_fracdiff == 0 ):
            [vmin, vmax] = [-2, 2]
        else:
            [vmin, vmax] = [-fracdiffabsmax, fracdiffabsmax]
        if verbose: print('Subplot (2,0) vmin, vmax: {}, {}'.format(vmin, vmax))

        plot4 = ax4.pcolormesh(cmpgrid['lat_b'], pedge[pedge_ind],
                               zm_fracdiff, cmap=cmap_plot, 
                               vmin=vmin, vmax=vmax)
        ax4.invert_yaxis()
        if regridany:
            ax4.set_title('Fractional Difference ({})\n(Dev-Ref)/Ref, Dynamic Range'.format(cmpres))
        else:
            ax4.set_title('Fractional Difference\n(Dev-Ref)/Ref, Dynamic Range')
        ax4.set_ylabel('Pressure (hPa)')
        ax4.set_aspect('auto')
        ax4.set_xticks(xtick_positions)
        ax4.set_xticklabels(xticklabels)

        # Define the colorbar for the plot.
        # If all values of absdiff = 0, then manually set a tick at 0.
        # Use zm_diff in the test for all zeroes, because fracdiff may
        # have NaN's due to div by zero when both Dev and Ref are zero.
        cb = plt.colorbar(plot4, ax=ax4, orientation='horizontal', pad=0.10)
        if np.all(zm_diff == 0):
            cb.locator = mpl.ticker.FixedLocator(locs=[0.0])
        else:
            if (vmax-vmin) < 0.1 or (vmax-vmin) > 100:
                cb.locator = mpl.ticker.MaxNLocator(nbins=4)
        cb.update_ticks()
        cb.set_clim(vmin=vmin, vmax=vmax)
        cb.set_label('unitless')   

        ################################################################
        # Subplot (2,1): Fractional Difference, restricted range
        ################################################################

        [vmin, vmax] = [-2, 2]
        if verbose: print('Subplot (2,1) vmin, vmax: {}, {}'.format(vmin, vmax))
  
        plot5 = ax5.pcolormesh(cmpgrid['lat_b'], pedge[pedge_ind],
                               zm_fracdiff, cmap=cmap_plot, 
                               vmin=vmin, vmax=vmax)
        ax5.invert_yaxis()
        if regridany:
            ax5.set_title('Fractional Difference ({})\n(Dev-Ref)/Ref, Fixed Range'.format(cmpres))
        else:
            ax5.set_title('Fractional Difference\n(Dev-Ref)/Ref, Fixed Range')
        ax5.set_ylabel('Pressure (hPa)')
        ax5.set_aspect('auto')
        ax5.set_xticks(xtick_positions)
        ax5.set_xticklabels(xticklabels)

        # Define the colorbar for the plot.
        # If all values of zm_diff = 0, then manually set a tick at 0.
        # Use zm_diff in the test for all zeroes, because fracdiff may
        # have NaN's due to div by zero when both Dev and Ref are zero.
        cb = plt.colorbar(plot5, ax=ax5, orientation='horizontal', pad=0.10)
        if np.all(zm_diff == 0): 
            cb.locator = mpl.ticker.FixedLocator(locs=[0.0])
        else:
            if (vmax-vmin) < 0.1 or (vmax-vmin) > 100:
                cb.locator = mpl.ticker.MaxNLocator(nbins=4)
        cb.update_ticks()
        cb.set_clim(vmin=vmin, vmax=vmax)
        cb.set_label('unitless') 

        ################################################################
        # Update the list of variables with significant differences.
        # Criterion: abs(max(fracdiff)) > 0.1
        # Do not include NaNs in the criterion, because these indicate
        # places where fracdiff could not be computed (div-by-zero).
        ################################################################
        if np.abs(np.nanmax(zm_fracdiff)) > 0.1:
            sigdiff_list.append(varname)

        ################################################################
        # Add this page of 6-panel plots to the PDF file
        ################################################################
        if savepdf:    
            pdf.savefig(figs)
            plt.close(figs)

    ####################################################################
    # Finish
    ####################################################################
    if savepdf:
        pdf.close()


def get_emissions_varnames(commonvars, template=None):
    '''
    Will return a list of emissions diagnostic variable names that
    contain a particular search string.

    Args:
        commonvars : list of strs
            A list of commmon variable names from two data sets.
            (This can be obtained with method gcpy.core.compare_varnames)

        template : str
            String template for matching variable names corresponding
            to emission diagnostics by sector.

    Returns:
        varnames : list of strs
           A list of variable names corresponding to emission
           diagnostics for a given species and sector.

    Example:
        >>> import gcpy
        >>> commonvars = ['EmisCO_Anthro', 'EmisNO_Anthro', 'AREA']
        >>> varnames = gcpy.get_emissions_varnames(commonvars, "Emis")
        >>> print(varnames)
        ['EmisCO_Anthro', 'EmisNO_Anthro']
    '''

    # Make sure the commonvars list has at least one element
    if len(commonvars) == 0:
        raise ValueError('No valid variable names were passed!')

    # Define template for emission diagnostics by sector
    if template is None:
        raise ValueError('The template argument was not passed!')

    # Find all emission diagnostics for the given species
    varnames = core.filter_names(commonvars, template)

    # Return list
    return varnames


def create_emission_display_name(diagnostic_name):
    '''
    Converts an emissions diagnostic name to a more easily digestible name
    that can be used as a plot title or in a table of emissions totals.

    Args:
        diagnostic_name : str
            Name of the diagnostic to be formatted

    Returns:
        display_name : str
            Formatted name that can be used as plot titles or in tables
            of emissions totals.

    Remarks:
        Assumes that diagnostic names will start with either "Emis"
        (for emissions by category) or "Inv" (for emissions by inventory).
        This should be an OK assumption to make since this routine is
        specifically geared towards model benchmarking.

    Example:
        >>> import gcpy
        >>> diag_name = "EmisCO_Anthro"
        >>> display_name = gcpy.create_emission_display_name(diag_name)
        >>> print(display_name)
        CO Anthro
    '''

    # Initialize
    display_name = diagnostic_name

    # Special handling for Inventory totals
    if 'INV' in display_name.upper():
        display_name = display_name.replace('_', ' ')

    # Replace text
    for v in ['Emis', 'EMIS', 'emis', 'Inv', 'INV', 'inv']:
        display_name = display_name.replace(v, '')

    # Replace underscores
    display_name = display_name.replace('_', ' ')

    # Return
    return display_name


def print_emission_totals(ref, refstr, dev, devstr, f):
    '''
    Computes and prints emission totals for two xarray DataArray objects.

    Args:
        ref : xarray DataArray
            The first DataArray to be compared (aka "Reference")

        refstr : str
            A string that can be used to identify refdata.
            (e.g. a model version number or other identifier).

        dev : xarray DatArray
            The second DataArray to be compared (aka "Development")

        devstr : str
            A string that can be used to identify devdata
            (e.g. a model version number or other identifier).

         f : file
            File object denoting a text file where output will be directed.

    Remarks:
        This is an internal method.  It is meant to be called from method
        create_total_emissions_table instead of being called directly.
    '''

    # Throw an error if the two DataArray objects have different units
    if ref.units != dev.units:
        msg = 'Ref has units "{}", but Dev array has units "{}"'.format(
            ref.units, dev.units)
        raise ValueError(msg)

    # Create the emissions display name from the diagnostic name
    diagnostic_name = dev.name
    display_name = create_emission_display_name(diagnostic_name)

    # Special handling for totals
    if '_TOTAL' in diagnostic_name.upper():
        print('-' * 79, file=f)

    # Compute sums and difference
    total_ref = np.sum(ref.values)
    total_dev = np.sum(dev.values)
    diff = total_dev - total_ref

    # Write output
    print('{} : {:13.6f}  {:13.6f}  {:13.6f} {}'.format(
        display_name.ljust(25), total_ref,
        total_dev, diff, dev.units), file=f)


def create_total_emissions_table(refdata, refstr, devdata, devstr,
                                 species, outfilename,
                                 interval=2678400.0, template='Emis{}_',
                                 ref_area_varname='AREA',
                                 dev_area_varname='AREA'):
    '''
    Creates a table of emissions totals (by sector and by inventory)
    for a list of species in contained in two data sets.  The data sets,
    which typically represent output from two differnet model versions,
    are usually contained in netCDF data files.

    Args:
        refdata : xarray Dataset
            The first data set to be compared (aka "Reference").

        refstr : str
            A string that can be used to identify refdata
            (e.g. a model version number or other identifier).

        devdata : xarray Dataset
            The second data set to be compared (aka "Development").

        devstr: str
            A string that can be used to identify the data set specified
            by devfile (e.g. a model version number or other identifier).

        species : dict
            Dictionary containing the name of each species and the target
            unit that emissions will be converted to. The format of
            species is as follows:

                { species_name : target_unit", etc. }

            where "species_name" and "target_unit" are strs.

        outfilename : str
            Name of the text file which will contain the table of
            emissions totals.

    Keyword Args (optional):
        interval : float
            The length of the data interval in seconds. By default, interval
            is set to the number of seconds in a 31-day month (86400 * 31),
            which corresponds to typical benchmark simulation output.

        template : str
            Template for the diagnostic names that are contained both
            "Reference" and "Development" data sets.  If not specified,
            template will be set to "Emis{}", where {} will be replaced
            by the species name.

        ref_area_varname : str
            Name of the variable containing the grid box surface areas
            (in m2) in the ref dataset.
            Default value: 'AREA'

        dev_area_varname : str
            Name of the variable containing the grid box surface areas
            (in m2) in the dev dataset.
            Default value: 'AREA'

    Returns:
        None.  Writes to a text file specified by the "outfilename" argument.

    Remarks:
        This method is mainly intended for model benchmarking purposes,
        rather than as a general-purpose tool.

        Species properties (such as molecular weights) are read from a
        JSON file called "species_database.json".

    Example:
        Print the total of CO and ACET emissions in two different
        data sets, which represent different model versions:

        >>> include gcpy
        >>> include xarray as xr
        >>> reffile = '~/output/12.1.1/HEMCO_diagnostics.201607010000.nc'
        >>> refstr = '12.1.1'
        >>> refdata = xr.open_dataset(reffile)
        >>> devfile = '~/output/12.2.0/HEMCO_sa.diagnostics.201607010000.nc'
        >>> devstr = '12.2.0'
        >>> devdata = xr.open_dataset(devfile)
        >>> outfilename = '12.2.0_emission_totals.txt'
        >>> species = { 'CO' : 'Tg', 'ACET', 'Tg C', 'ALK4', 'Gg C' }
        >>> create_total_emissions_table(refdata, refstr, devdata, devstr,
            species, outfilename)
    '''

    # Error check arguments
    if not isinstance(refdata, xr.Dataset):
        raise ValueError('The refdata argument must be an xarray Dataset!')

    if not isinstance(devdata, xr.Dataset):
        raise ValueError('The devdata argument must be an xarray Dataset!')

    if ref_area_varname not in refdata.data_vars.keys():
        raise ValueError('Area variable {} is not in the ref Dataset!'.format(
            ref_area_varname))

    if dev_area_varname not in devdata.data_vars.keys():
        raise ValueError('Area variable {} is not in the dev Dataset!'.format(
            dev_area_varname))

    # Load a JSON file containing species properties (such as
    # molecular weights), which we will need for unit conversions.
    # This is located in the "data" subfolder of this current directory.2
    properties_path = os.path.join(os.path.dirname(__file__),
                                   'species_database.json')
    properties = json_load_file(open(properties_path))

    # Find all common variables between the two datasets
    [cvars, cvarsOther, cvars2D, cvars3D] = core.compare_varnames(refdata,
                                                               devdata,
                                                               quiet=True)

    # Open file for output
    f = open(outfilename, 'w')

    # Loop through all of the species are in species_dict
    for species_name, target_units in species.items():

        # Get a list of emission variable names for each species
        diagnostic_template = template.format(species_name)
        varnames = get_emissions_varnames(cvars, diagnostic_template)

        # If no emissions are found, then skip to next species
        if len(varnames) == 0:
            print('No emissions found for {} ... skippping'.format(
                species_name))
            continue

        # Check if there is a total emissions variable in the list
        vartot = [v for v in varnames if '_TOTAL' in v.upper()]

        # Push the total variable to the last list element
        # so that it will be printed last of all
        if len(vartot) == 1:
            varnames.append(varnames.pop(varnames.index(vartot[0])))   

        # Title strings
        if 'Inv' in template:
            print('Computing inventory totals for {}'.format(species_name))
            title1 = '### Inventory totals for species {}'.format(species_name)
        else:
            print('Computing emissions totals for {}'.format(species_name))
            title1 = '### Emissions totals for species {}'.format(species_name)

        title2 = '### Ref = {}; Dev = {}'.format(refstr, devstr)

        # Print header to file
        print('#'*79, file=f)
        print('{}{}'.format(title1.ljust(76), '###'), file=f)
        print('{}{}'.format(title2.ljust(76), '###'), file=f)
        print('#'*79, file=f)
        print('{}{}{}{}'.format(' '.ljust(33), 'Ref'.ljust(15),
                                'Dev'.ljust(15), 'Dev - Ref'), file=f)

        # Loop over all emissions variable names
        for v in varnames:

            if 'Inv' in template:
                spc_name = v.split('_')[1]
            else:
                spc_name = species_name

            # Get a list of properties for the given species
            species_properties = properties.get(spc_name)

            # Convert units of Ref, and save to DataArray
            refarray = convert_units(refdata[v], spc_name,
                                     species_properties, target_units,
                                     interval, refdata[ref_area_varname])

            # Convert units of Dev, and save to DataArray
            devarray = convert_units(devdata[v], spc_name,
                                     species_properties, target_units,
                                     interval, devdata[dev_area_varname])

            # Print emission totals for Ref and Dev
            print_emission_totals(refarray, refstr, devarray, devstr, f)

        # Add newlines before going to the next species
        print(file=f)
        print(file=f)

    # Close file
    f.close()


def get_species_categories():
    '''
    Returns the list of benchmark categories that each species
    belongs to.  This determines which PDF files will contain the
    plots for the various species.

    Returns:
        spc_cat_dict : dict
            A nested dictionary of categories (and sub-categories)
            and the species belonging to each.

    NOTE: The benchmark categories are specified in JSON file
    benchmark_species.json.
    '''
    jsonfile = os.path.join(os.path.dirname(__file__), spc_categories)
    with open(jsonfile, 'r') as f:
        spc_cat_dict = json.loads(f.read())
    return spc_cat_dict


def archive_species_categories(dst):
    '''
    Writes the list of benchmark categories to a JSON file
    named "benchmark_species.json".

    Args:
        dst : str
            Name of the folder where the JSON file containing
            benchmark categories ("benchmark_species.json")
            will be written.
    '''

    src = os.path.join(os.path.dirname(__file__), spc_categories)
    print('Archiving {} in {}'.format(spc_categories, dst))
    shutil.copyfile(src, os.path.join(dst, spc_categories))


def make_benchmark_conc_plots(ref, refstr, dev, devstr, dst='./1mo_benchmark',
                              overwrite=False, verbose=False, restrict_cats=[],
                              plots=['sfc', '500hpa', 'zonalmean'], 
                              use_cmap_RdBu=False, log_color_scale=False,
                              sigdiff_files=None):
    '''
    Creates PDF files containing plots of species concentration
    for model benchmarking purposes.

    Args:
        ref: str
             Path name for the "Ref" (aka "Reference") data set.

        refstr : str
             A string to describe ref (e.g. version number)

        dev : str
             Path name for the "Dev" (aka "Development") data set.
             This data set will be compared against the "Reference"
             data set.

        devstr : str
             A string to describe dev (e.g. version number)

    Keyword Args (optional):
        dst : str
             A string denoting the destination folder where a PDF
             file containing plots will be written.
             Default value: ./1mo_benchmark

        overwrite : boolean
             Set this flag to True to overwrite files in the
             destination folder (specified by the dst argument).
             Default value: False.

        verbose : boolean
             Set this flag to True to print extra informational output.
             Default value: False.

        restrict_cats : list of strings
             List of benchmark categories in benchmark_categories.json to make
             plots for. If empty, plots are made for all categories.
             Default value: empty

        plots : list of strings
             List of plot types to create.
             Default value: ['sfc', '500hpa', 'zonalmean']

        log_color_scale: boolean         
            Logical to enable plotting data (not diffs) on a log color scale.
            Default value: False

        sigdiff_files : list of str
            Filenames that will contain the lists of species having
            significant differences in the 'sfc', '500hpa', and
            'zonalmean' plots.  These lists are needed in order to
            fill out the benchmark approval forms.
            Default value: None
    '''

    # NOTE: this function could use some refactoring; 
    # abstract processing per category?

    # ==================================================================
    # Initialization and data read
    # ==================================================================
    if os.path.isdir(dst) and not overwrite:
        print('Directory {} exists. Pass overwrite=True to overwrite files in that directory, if any.'.format(dst))
        return
    elif not os.path.isdir(dst):
        os.mkdir(dst)

    # Ref dataset
    try:
        refds = xr.open_dataset(ref)
    except FileNotFoundError:
        print('Could not find Ref file: {}'.format(ref))
        raise
    refds = core.add_lumped_species_to_dataset(refds, verbose=verbose)
    
    # Dev dataset
    try:
        devds = xr.open_dataset(dev)
    except FileNotFoundError:
        print('Could not find Dev file: {}!'.format(dev))
        raise
    devds = core.add_lumped_species_to_dataset(devds, verbose=verbose)

    catdict = get_species_categories()
    
    archive_species_categories(dst)
    core.archive_lumped_species_definitions(dst)

    # ==================================================================
    # Create the plots!
    # ==================================================================
    for i, filecat in enumerate(catdict):

        # If restrict_cats list is passed,
        # skip all categories except those in the list
        if restrict_cats and filecat not in restrict_cats:
            continue

        catdir = os.path.join(dst,filecat)
        if not os.path.isdir(catdir):
            os.mkdir(catdir)
        varlist = []
        warninglist = []
        for subcat in catdict[filecat]:
            for spc in catdict[filecat][subcat]:
                varname = 'SpeciesConc_'+spc
                if varname not in refds.data_vars or varname not in devds.data_vars:
                    warninglist.append(varname)
                    continue
                varlist.append(varname)
        if warninglist != []:
            print('\n\nWarning: variables in {} category not in dataset: {}'.format(filecat,warninglist))

        # Surface plots
        if 'sfc' in plots:
            pdfname = os.path.join(catdir,'{}_Surface.pdf'.format(filecat))
            diff_sfc = []
            compare_single_level(refds, refstr, devds, devstr, 
                                 varlist=varlist, ilev=0,
                                 pdfname=pdfname,
                                 use_cmap_RdBu=use_cmap_RdBu,
                                 log_color_scale=log_color_scale,
                                 sigdiff_list=diff_sfc)
            diff_sfc[:] = [v.replace('SpeciesConc_', '') for v in diff_sfc]
            add_nested_bookmarks_to_pdf(pdfname, filecat, 
                                        catdict, warninglist, 
                                        remove_prefix='SpeciesConc_')

        if '500hpa' in plots:
            pdfname = os.path.join(catdir,
                                   '{}_500hPa.pdf'.format(filecat))        
            diff_500 = []
            compare_single_level(refds, refstr, devds, devstr, 
                                 varlist=varlist, ilev=22,
                                 pdfname=pdfname, 
                                 use_cmap_RdBu=use_cmap_RdBu,
                                 log_color_scale=log_color_scale,
                                 sigdiff_list=diff_500)
            diff_500[:] = [v.replace('SpeciesConc_', '') for v in diff_500]
            add_nested_bookmarks_to_pdf(pdfname, filecat, 
                                        catdict, warninglist, 
                                        remove_prefix='SpeciesConc_')

        if 'zonalmean' in plots or 'zm' in plots:
            pdfname = os.path.join(catdir,'{}_FullColumn_ZonalMean.pdf'.format(
                filecat))
            diff_zm = []
            compare_zonal_mean(refds, refstr, devds, devstr, varlist=varlist,
                               pdfname=pdfname, use_cmap_RdBu=use_cmap_RdBu,
                               log_color_scale=log_color_scale,
                               sigdiff_list=diff_zm)
            diff_zm[:] = [v.replace('SpeciesConc_', '') for v in diff_zm]
            add_nested_bookmarks_to_pdf(pdfname, filecat, 
                                        catdict, warninglist, 
                                        remove_prefix='SpeciesConc_')

            pdfname = os.path.join(catdir,'{}_Strat_ZonalMean.pdf'.format(
                filecat))        
            compare_zonal_mean(refds, refstr, devds, devstr, varlist=varlist,
                               pdfname=pdfname, pres_range=[0,100], 
                               use_cmap_RdBu=use_cmap_RdBu,
                               log_color_scale=log_color_scale)
            add_nested_bookmarks_to_pdf(pdfname, filecat, 
                                        catdict, warninglist, 
                                        remove_prefix='SpeciesConc_')

        # ==============================================================
        # Write the list of species having significant differences,
        # which we need to fill out the benchmark approval forms.
        # ==============================================================
        if sigdiff_files != None:
            for filename in sigdiff_files:
                if 'sfc' in filename:
                    with open(filename, 'a+') as f:
                        print('* {}: '.format(filecat), file=f, end='')
                        for v in diff_sfc:
                            print('{} '.format(v), file=f, end='')
                        print(file=f)
                        f.close()

                if '500' in filename:
                    with open(filename, 'a+') as f:
                        print('* {}: '.format(filecat), file=f, end='')
                        for v in diff_500:
                            print('{} '.format(v), file=f, end='')
                        print(file=f)
                        f.close()

                if 'zonalmean' in filename or 'zm' in filename:
                    with open(filename, 'a+') as f:
                        print('* {}: '.format(filecat), file=f, end='')
                        for v in diff_zm:
                            print('{} '.format(v), file=f, end='')
                        print(file=f)
                        f.close()


def make_benchmark_emis_plots(ref, refstr, dev, devstr,
                              dst='./1mo_benchmark',
                              plot_by_benchmark_cat=False,
                              plot_by_hco_cat=False,
                              overwrite=False, verbose=False,
                              flip_ref=False, flip_dev=False,
                              log_color_scale=False,
                              sigdiff_files=None):
    '''
    Creates PDF files containing plots of emissions for model
    benchmarking purposes. This function is compatiblity with benchmark 
    simulation output only. It is not compatible with transport tracers
    emissions diagnostics.

    Args:
        ref: str
             Path name for the "Ref" (aka "Reference") data set.

        refstr : str
             A string to describe ref (e.g. version number)

        dev : str
             Path name for the "Dev" (aka "Development") data set.
             This data set will be compared against the "Reference"
             data set.

        devstr : str
             A string to describe dev (e.g. version number)

    Keyword Args (optional):
        dst : str
             A string denoting the destination folder where a
             PDF file containing plots will be written.
             Default value: './1mo_benchmark

        plot_by_benchmark_cat : boolean
             Set this flag to True to separate plots into PDF files
             according to the benchmark categories (e.g. Primary,
             Aerosols, Nitrogen, etc.)  These categories are specified
             in the JSON file benchmark_species.json.
             Default value: False

        plot_by_hco_cat : boolean
             Set this flag to True to separate plots into PDF files
             according to HEMCO emissions categories (e.g. Anthro,
             Aircraft, Bioburn, etc.)
             Default value: False

        overwrite : boolean
             Set this flag to True to overwrite files in the
             destination folder (specified by the dst argument).
             Default value: False

        verbose : boolean
             Set this flag to True to print extra informational output.
             Default value: False

        flip_ref : boolean
             Set this flag to True to reverse the vertical level
             ordering in the "Ref" dataset (in case "Ref" starts
             from the top of atmosphere instead of the surface).
             Default value: False

        flip_dev : boolean
             Set this flag to True to reverse the vertical level
             ordering in the "Dev" dataset (in case "Dev" starts
             from the top of atmosphere instead of the surface).
             Default value: False

        log_color_scale: boolean         
            Logical to enable plotting data (not diffs) on a log color scale.
            Default value: False

         sigdiff_files : list of str
            Filenames that will contain the lists of species having
            significant differences in the 'sfc', '500hpa', and
            'zonalmean' plots.  These lists are needed in order to
            fill out the benchmark approval forms.
            Default value: None

    Remarks:
        (1) If both plot_by_benchmark_cat and plot_by_hco_cat are
            False, then all emission plots will be placed into the
            same PDF file.

        (2) Emissions that are 3-dimensional will be plotted as
            column sums.
    '''
    # =================================================================
    # Initialization and data read
    # =================================================================
    if os.path.isdir(dst) and not overwrite:
        print('Directory {} exists. Pass overwrite=True to overwrite files in that directory, if any.'.format(dst))
        return
    elif not os.path.isdir(dst):
        os.mkdir(dst)
    emisdir = os.path.join(dst,'Emissions')
    if not os.path.isdir(emisdir):
        os.mkdir(emisdir)   

    # Ref dataset
    try:
        refds = xr.open_dataset(ref)
    except FileNotFoundError:
        print('Could not find Ref file: {}'.format(ref))
        raise

    # Dev dataset
    try:
        devds = xr.open_dataset(dev)
    except FileNotFoundError:
        print('Could not find Dev file: {}'.format(dev))
        raise

    # Find common variables
    quiet = not verbose
    vars, varsOther, vars2D, vars3D = core.compare_varnames(refds, devds, quiet)

    # Combine 2D and 3D variables into an overall list
    varlist = vars2D + vars3D

    # =================================================================
    # Compute column sums for 3D emissions
    # Make sure not to clobber the DataArray attributes
    # =================================================================
    with xr.set_options(keep_attrs=True):
        for v in vars3D:
            if 'lev' in refds[v].dims:
                refds[v] = refds[v].sum(dim='lev')
            if 'lev' in devds[v].dims:
                devds[v] = devds[v].sum(dim='lev')

    # =================================================================
    # If inputs plot_by* are both false, plot all emissions in same file
    # =================================================================
    if not plot_by_benchmark_cat and not plot_by_hco_cat:
        pdfname = os.path.join(emisdir,'Emissions.pdf')
        compare_single_level(refds, refstr, devds, devstr,
                             varlist=varlist, pdfname=pdfname,
                             log_color_scale=log_color_scale)
        add_bookmarks_to_pdf(pdfname, varlist,
                             remove_prefix='Emis', verbose=verbose)
        return

    # Get emissions variables (non-inventory), categories, and species
    emis_vars = [v for v in varlist if v[:4] == 'Emis']
    emis_cats = sorted(set([v.split('_')[1] for v in emis_vars]))
    emis_spc = sorted(set([v.split('_')[0][4:] for v in emis_vars]))

# This is fixed in 12.3.2, comment out for now (bmy, 5/1/19)
#    # Handle Bioburn and BioBurn as same categories (temporary until 12.3.1)
#    emis_cats.remove('BioBurn')
    
    # Sort alphabetically (assume English characters)
    emis_vars.sort(key=str.lower)

    # =================================================================
    # if plot_by_hco_cat is true, make a file for each HEMCO emissions
    # category that is in the diagnostics file
    #
    # Also write the list of emission quantities that have significant
    # diffs.  We'll need that to fill out the benchmark forms.
    # =================================================================
    if plot_by_hco_cat:
        emisspcdir = os.path.join(dst,'Emissions')
        if not os.path.isdir(emisspcdir):
            os.mkdir(emisspcdir)   

        diff_dict = {}
        for c in emis_cats:
            # Handle cases of bioburn and bioBurn (temporary until 12.3.1)
            if c == 'Bioburn':
                varnames = [k for k in emis_vars if any(b in k for b in ['Bioburn','BioBurn'])]
            else:
                varnames = [k for k in emis_vars if c in k]
            pdfname = os.path.join(emisspcdir,'{}_Emissions.pdf'.format(c))
            diff_emis = []
            compare_single_level(refds, refstr, devds, devstr,
                                 varlist=varnames, ilev=0, pdfname=pdfname,
                                 log_color_scale=log_color_scale,
                                 sigdiff_list=diff_emis)
            add_bookmarks_to_pdf(pdfname, varnames,
                                 remove_prefix='Emis', verbose=verbose)

            # Save the list of quantities with significant differences for
            # this category into the diff_dict dictionary for use below
            diff_emis[:] = [v.replace('Emis', '') for v in diff_emis]
            diff_emis[:] = [v.replace('_' + c, '') for v in diff_emis]
            diff_dict[c] = diff_emis

        # -------------------------------------------------------------
        # Write the list of species having significant differences,
        # which we need to fill out the benchmark approval forms.
        # -------------------------------------------------------------
        if sigdiff_files != None:
            for filename in sigdiff_files:
                if 'emis' in filename:
                    with open(filename, 'w+') as f:
                        for c, diff_list in diff_dict.items():
                            print('* {}: '.format(c), file=f, end='')
                            for v in diff_list:
                                print('{} '.format(v), file=f, end='')
                            print(file=f)
                        f.close()

    # =================================================================
    # if plot_by_benchmark_cat is true, make a file for each benchmark
    # species category with emissions in the diagnostics file
    # =================================================================
    if plot_by_benchmark_cat:
        
        catdict = get_species_categories()
        warninglist = [] # in case any emissions are skipped (for use in nested pdf bookmarks)
        allcatspc = []   # for checking if emissions species not defined in benchmark category file
        emisdict = {}    # used for nested pdf bookmarks
        for i, filecat in enumerate(catdict):

            # Get emissions for species in this benchmark category
            varlist = []
            emisdict[filecat] = {}
            for subcat in catdict[filecat]:
                for spc in catdict[filecat][subcat]:
                    allcatspc.append(spc)
                    if spc in emis_spc:
                        emisdict[filecat][spc] = []
                        emisvars = [v for v in emis_vars if spc == v.split('_')[0][4:]]
                        for var in emisvars:
                            emisdict[filecat][spc].append(var.replace('Emis',''))
                            varlist.append(var)
            if not varlist:
                print('\nWarning: no emissions species in benchmark species category {}'.format(filecat))
                continue        

            # Use same directory structure as for concentration plots
            catdir = os.path.join(dst,filecat)
            if not os.path.isdir(catdir):
                os.mkdir(catdir)

            # Create emissions file for this benchmark species category
            pdfname = os.path.join(catdir,'{}_Emissions.pdf'.format(filecat))
            compare_single_level(refds, refstr, devds, devstr, 
                                 varlist=varlist, ilev=0, pdfname=pdfname, 
                                 flip_ref=flip_ref, flip_dev=flip_dev,
                                 log_color_scale=log_color_scale)
            add_nested_bookmarks_to_pdf(pdfname, filecat, emisdict, warninglist)

        # Give warning if emissions species is not assigned a benchmark category
        for spc in emis_spc:
            if spc not in allcatspc:
                print('Warning: species {} has emissions diagnostics but is not in benchmark_categories.json'.format(spc))


def make_benchmark_emis_tables(reflist, refstr, devlist, devstr,
                               dst='./1mo_benchmark', overwrite=False,
                               interval=None, ref_area_varname=None,
                               dev_area_varname=None):
    '''
    Creates a text file containing emission totals by species and
    category for benchmarking purposes.

    Args:
        reflist: list of str
             List with the path names of the emissions and/or met field
             files that will constitute the "Ref" (aka "Reference")
             data set.

        refstr : str
             A string to describe ref (e.g. version number)

        devlist : list of str
             List with the path names of the emissions and/or met field
             files that will constitute the "Dev" (aka "Development") 
             data set.  The "Dev" data set will be compared against the
             "Ref" data set.

        devstr : str
             A string to describe dev (e.g. version number)

    Keyword Args (optional):
        dst : str
             A string denoting the destination folder where the file
             containing emissions totals will be written.
             Default value: ./1mo_benchmark

        overwrite : boolean
             Set this flag to True to overwrite files in the
             destination folder (specified by the dst argument).
             Default value : False

        interval : float
             Specifies the averaging period in seconds, which is used
             to convert fluxes (e.g. kg/m2/s) to masses (e.g kg).
             Default value : None

        ref_area_varname : str
            Name of the variable containing the grid box surface areas
            (in m2) in the ref dataset.  If not specified, then this
            will be set to "Met_AREAM2" if the ref dataset is on the
            cubed-sphere grid, or "AREA" if ref is on the lat-lon grid.
            Default value: None

        dev_area_varname : str
            Name of the variable containing the grid box surface areas
            (in m2) in the dev dataset.  If not specified, then this
            will be set to "Met_AREAM2" if the dev dataset is on the
            cubed-sphere grid, or "AREA" if dev is on the lat-lon grid.
            Default value: None
    '''

    # ===============================================================
    # Define destination directory
    # ===============================================================
    if os.path.isdir(dst) and not overwrite:
        print('Directory {} exists. Pass overwrite=True to overwrite files in that directory, if any.'.format(dst))
        return
    elif not os.path.isdir(dst):
        os.mkdir(dst)
    emisdir = os.path.join(dst,'Emissions')
    if not os.path.isdir(emisdir):
        os.mkdir(emisdir)

    # ===============================================================
    # Read data from netCDF into Dataset objects
    # ===============================================================

    # Ref
    try:
        refds = xr.open_mfdataset(reflist)
    except FileNotFoundError:
        print('Could not find one of the Ref files: {}'.format(reflist))
        raise

    # Dev
    try:
        devds = xr.open_mfdataset(devlist)
    except FileNotFoundError:
        print('Could not find one of the Dev files: {}'.format(devlist))
        raise

    # ===============================================================
    # Get the surface area variable name
    # by default, for lat-lon grids this is "AREA"
    # and for cubed-sphere grids, this is "Met_AREAM2"
    # ===============================================================

    # Ref
    if ref_area_varname == None:
        nlat = refds.sizes['lat']
        nlon = refds.sizes['lon']
        if nlat/6 == nlon:
            ref_area_varname = 'Met_AREAM2'
        else:
            ref_area_varname = 'AREA'

    # Error-check Ref
    if ref_area_varname not in refds.data_vars.keys():
        raise ValueError('Area variable {} is not in the ref dataset!'.format(
            ref_area_varname))

    # Dev
    if dev_area_varname == None:
        nlat = devds.sizes['lat']
        nlon = devds.sizes['lon']
        if nlat/6 == nlon:
            dev_area_varname = 'Met_AREAM2'
        else:
            dev_area_varname = 'AREA'

    # Error-check Dev
    if dev_area_varname not in devds.data_vars.keys():
        raise ValueError('Area variable {} is not in the dev dataset!'.format(
            dev_area_varname))

    # ===============================================================
    # Create table of emissions
    # ===============================================================

    # Emissions species dictionary
    species = json_load_file(open(os.path.join(os.path.dirname(__file__),
                                               emission_spc)))
    inventories = json_load_file(open(os.path.join(os.path.dirname(__file__),
                                                   emission_inv)))

    # Destination files
    file_emis_totals = os.path.join(dst, emisdir, 'Emission_totals.txt')
    file_inv_totals = os.path.join(dst, emisdir, 'Inventory_totals.txt')

    # If the averaging interval (in seconds) is not specified,
    # then assume July 2016 = 86400 seconds * 31 days
    if interval == None:
        interval = 86400.0 * 31.0

    # Write to file
    create_total_emissions_table(refds, refstr, devds, devstr, 
                                 species, file_emis_totals,
                                 interval, template='Emis{}_',
                                 ref_area_varname=ref_area_varname,
                                 dev_area_varname=dev_area_varname)
    create_total_emissions_table(refds, refstr, devds, devstr, 
                                 inventories, file_inv_totals,
                                 interval, template='Inv{}_',
                                 ref_area_varname=ref_area_varname,
                                 dev_area_varname=dev_area_varname)


def make_benchmark_jvalue_plots(ref, refstr, dev, devstr,
                                varlist=None, 
                                dst='./1mo_benchmark',
                                local_noon_jvalues=False,
                                plots=['sfc', '500hpa', 'zonalmean'],
                                overwrite=False, verbose=False,
                                flip_ref=False, flip_dev=False,
                                log_color_scale=False,
                                sigdiff_files=None):
    '''
    Creates PDF files containing plots of J-values for model
    benchmarking purposes.

    Args:
        ref: str
             Path name for the "Ref" (aka "Reference") data set.

        refstr : str
             A string to describe ref (e.g. version number)
    
        dev : str
             Path name for the "Dev" (aka "Development") data set.
             This data set will be compared against the "Reference"
             data set.

        devstr : str
             A string to describe dev (e.g. version number)
    
    Keyword Args (optional):
        varlist : list
             List of J-value variables to plot.  If not passed, 
             then all J-value variables common to both dev 
             and ref will be plotted.  The varlist argument can be
             a useful way of restricting the number of variables
             plotted to the pdf file when debugging.
             Default value: None

        dst : str
             A string denoting the destination folder where a
             PDF file  containing plots will be written.
             Default value: ./1mo_benchmark.
    
        local_noon_jvalues : boolean
             Set this switch to plot local noon J-values.  This will
             divide all J-value variables by the JNoonFrac counter,
             which is the fraction of the time that it was local noon
             at each location.
             Default value : False

        plots : list of strings
             List of plot types to create.
             Default value: ['sfc', '500hpa', 'zonalmean']

        overwrite : boolean
             Set this flag to True to overwrite files in the
             destination folder (specified by the dst argument).
             Default value: False.

        verbose : boolean
             Set this flag to True to print extra informational output.
             Default value: False

        flip_ref : boolean
             Set this flag to True to reverse the vertical level
             ordering in the "Ref" dataset (in case "Ref" starts
             from the top of atmosphere instead of the surface).
             Default value: False

        flip_dev : boolean
             Set this flag to True to reverse the vertical level
             ordering in the "Dev" dataset (in case "Dev" starts
             from the top of atmosphere instead of the surface).
             Default value: False

        log_color_scale: boolean         
            Logical to enable plotting data (not diffs) on a log color scale.
            Default value: False

        sigdiff_files : list of str
            Filenames that will contain the lists of J-values having
            significant differences in the 'sfc', '500hpa', and
            'zonalmean' plots.  These lists are needed in order to
            fill out the benchmark approval forms.
            Default value: None

    Remarks:
         Will create 4 files containing J-value plots:
            (1 ) Surface values
            (2 ) 500 hPa values
            (3a) Full-column zonal mean values.
            (3b) Stratospheric zonal mean values
         These can be toggled on/off with the plots keyword argument.

         At present, we do not yet have the capability to split the
         plots up into separate files per category (e.g. Primary,
         Aerosols, etc.).  This is primarily due to the fact that 
         we archive J-values from GEOS-Chem for individual species
         but not family species.  We could attempt to add this 
         functionality later if there is sufficient demand. 
    '''
    
    if os.path.isdir(dst) and not overwrite:
        print('Directory {} exists. Pass overwrite=True to overwrite files in tht directory, if any.'.format(dst))
        return
    elif not os.path.isdir(dst):
        os.mkdir(dst)

    # Ref dataset
    try:
        refds = xr.open_dataset(ref)
    except FileNotFoundError:
        print('Could not find Ref file: {}'.format(ref))
        raise

    # Dev dataset
    try:
        devds = xr.open_dataset(dev)
    except FileNotFoundError:
        print('Could not find Dev file: {}'.format(dev))
        raise

    # Find common variables in both datasets
    if varlist == None:
        quiet = not verbose
        [cmn, cmnOther, cmn2D, cmn3D] = core.compare_varnames(refds, devds, quiet)

    # =================================================================
    # Local noon or continuously-averaged J-values?
    # =================================================================
    if local_noon_jvalues:

        # Get a list of local noon J-value variables
        # (or use the varlist passed via tha argument list)
        prefix = 'JNoon_'
        if varlist == None:
            varlist = [v for v in cmn if prefix in v]

        # Make sure JNoonFrac (fraction of times it was local noon
        # in each column) is present in both Ref and Dev datasets
        if not 'JNoonFrac' in cmn:
            msg = 'JNoonFrac is not common to Ref and Dev datasets!'
            raise ValueError(msg)

        # JNoon_* are cumulative sums of local noon J-values; we need
        # to divide these by JNoonFrac to get the average value
        refds = core.divide_dataset_by_dataarray(refds,
                                                 refds['JNoonFrac'],
                                                 varlist)
        devds = core.divide_dataset_by_dataarray(devds,
                                                 devds['JNoonFrac'],
                                                 varlist)

        # Subfolder of dst where PDF files will be printed
        subdir= 'JValuesLocalNoon'

    else:

        # Get a list of continuously averaged J-value variables
        # (or use the varlist passed via tha argument list)
        prefix = 'Jval_'
        if varlist == None:
            varlist = [v for v in cmn if prefix in v]

        # Subfolder of dst where PDF files will be printed
        subdir = 'JValues'

    # =================================================================
    # Create the plots
    # =================================================================

    # Make the folder to contain plots if it doesn't exist
    jvdir = os.path.join(dst, subdir)
    if not os.path.isdir(jvdir):
        os.mkdir(jvdir)
    
    # Surface plots
    if 'sfc' in plots:
        pdfname = os.path.join(jvdir, '{}Surface.pdf'.format(prefix))
        diff_sfc = []
        compare_single_level(refds, refstr, devds, devstr,
                             varlist=varlist, ilev=0, pdfname=pdfname,
                             flip_ref=flip_ref, flip_dev=flip_dev,
                             log_color_scale=log_color_scale,
                             sigdiff_list=diff_sfc)
        diff_sfc[:] = [v.replace(prefix, '') for v in diff_sfc]
        add_bookmarks_to_pdf(pdfname, varlist,
                             remove_prefix=prefix, verbose=verbose)

    # 500hPa plots
    if '500hpa' in plots:
        pdfname = os.path.join(jvdir, '{}500hPa.pdf'.format(prefix))
        diff_500 = []
        compare_single_level(refds, refstr, devds, devstr,
                             varlist=varlist, ilev=22, pdfname=pdfname,
                             flip_ref=flip_ref, flip_dev=flip_dev,
                             log_color_scale=log_color_scale,
                             sigdiff_list=diff_500)
        diff_500[:] = [v.replace(prefix, '') for v in diff_500]
        add_bookmarks_to_pdf(pdfname, varlist,
                             remove_prefix=prefix, verbose=verbose)

    # Full-column zonal mean plots
    if 'zonalmean' in plots:
        pdfname = os.path.join(jvdir,
                               '{}FullColumn_ZonalMean.pdf'.format(prefix))
        diff_zm = []
        compare_zonal_mean(refds, refstr, devds, devstr,
                           varlist=varlist, pdfname=pdfname,
                           flip_ref=flip_ref, flip_dev=flip_dev,
                           log_color_scale=log_color_scale,
                           sigdiff_list=diff_zm)
        diff_zm[:] = [v.replace(prefix, '') for v in diff_zm]
        add_bookmarks_to_pdf(pdfname, varlist,
                             remove_prefix=prefix, verbose=verbose)

        # Stratospheric zonal mean plots
        pdfname = os.path.join(jvdir,'{}Strat_ZonalMean.pdf'.format(prefix))
        compare_zonal_mean(refds, refstr, devds, devstr,
                           varlist=varlist, pdfname=pdfname, pres_range=[0,100],
                           flip_ref=flip_ref, flip_dev=flip_dev,
                           log_color_scale=log_color_scale)
        add_bookmarks_to_pdf(pdfname, varlist,
                             remove_prefix=prefix, verbose=verbose)

        # =================================================================
        # Write the lists of J-values that have significant differences,
        # which we need to fill out the benchmark approval forms.
        # =================================================================
        if sigdiff_files != None:
            for filename in sigdiff_files:
                if 'sfc' in filename:
                    with open(filename, 'a+') as f:
                        print('* J-Values: ', file=f, end='')
                        for v in diff_sfc:
                            print('{} '.format(v), file=f, end='')
                        print(file=f)
                        f.close()

                if '500' in filename:
                    with open(filename, 'a+') as f:
                        print('* J-Values: ', file=f, end='')
                        for v in diff_500:
                            print('{} '.format(v), file=f, end='')
                        print(file=f)
                        f.close()

                if 'zonalmean' in filename or 'zm' in filename:
                    with open(filename, 'a+') as f:
                        print('* J-Values: ', file=f, end='')
                        for v in diff_zm:
                            print('{} '.format(v), file=f, end='')
                        print(file=f)
                        f.close()


def make_benchmark_aod_plots(ref, refstr, dev, devstr,
                             varlist=None, dst='./1mo_benchmark',
                             overwrite=False, verbose=False,
                             log_color_scale=False,
                             sigdiff_files=None):
    '''
    Creates PDF files containing plots of column aerosol optical
    depths (AODs) for model benchmarking purposes.

    Args:
        ref: str
             Path name for the "Ref" (aka "Reference") data set.

        refstr : str
             A string to describe ref (e.g. version number)

        dev : str
             Path name for the "Dev" (aka "Development") data set.
             This data set will be compared against the "Reference"
             data set.

        devstr : str
             A string to describe dev (e.g. version number)

    Keyword Args (optional):
        varlist : list
             List of AOD variables to plot.  If not passed, 
             then all AOD variables common to both dev 
             and ref will be plotted.  The varlist argument can be
             a useful way of restricting the number of variables
             plotted to the pdf file when debugging.
             Default value: None

        dst : str
             A string denoting the destination folder where a
             PDF file  containing plots will be written.
             Default value: ./1mo_benchmark.

        overwrite : boolean
             Set this flag to True to overwrite files in the
             destination folder (specified by the dst argument).
             Default value: False.

        verbose : boolean
             Set this flag to True to print extra informational output.
             Default value: False

        log_color_scale: boolean         
            Logical to enable plotting data (not diffs) on a log color scale.
            Default value: False

        sigdiff_files : list of str
            Filenames that will contain the list of quantities having
            having significant differences in the column AOD plots.
            These lists are needed in order to fill out the benchmark
            approval forms.
            Default value: None
    '''
    # =================================================================
    # Initialization and also read data
    # =================================================================
    if os.path.isdir(dst) and not overwrite:
        print('Directory {} exists. Pass overwrite=True to overwrite files in tht directory, if any.'.format(dst))
        return
    elif not os.path.isdir(dst):
        os.mkdir(dst)
    aoddir = os.path.join(dst, 'Aerosols')
    if not os.path.isdir(aoddir):
        os.mkdir(aoddir)

    # Read the Ref dataset
    try:
        refds = xr.open_dataset(ref)
    except FileNotFoundError:
        print('Could not find Ref file: {}'.format(ref))
        raise

    # Read the Dev dataset
    try:
        devds = xr.open_dataset(dev)
    except FileNotFoundError:
        print('Could not find Dev file: {}'.format(dev))
        raise

    # NOTE: GCHP diagnostic variable exports are defined before the
    # input.geos file is read.  This means "WL1" will not have been
    # replaced with "550nm" in the variable names.  Do this string
    # replace operation here, so that we can compare GCC and GCHP
    # data directly. (bmy, 4/29/19)
    with xr.set_options(keep_attrs=True):

        # Rename variables in the Ref dataset
        old2new = {}
        for v in refds.data_vars.keys():
            if 'WL1' in v:
                newname = v.replace('WL1', '550nm')
                old2new[v] = newname
        refds = refds.rename(old2new)

        # Rename variables in the Dev dataset
        old2new = {}
        for v in devds.data_vars.keys():
            if 'WL1' in v:
                newname = v.replace('WL1', '550nm')
                old2new[v] = newname
        devds = devds.rename(old2new)

    # Find common AOD variables in both datasets
    # (or use the varlist passed via keyword argument)
    if varlist == None:
        quiet = not verbose
        [cmn, cmnOther, cmn2D, cmn3D] = core.compare_varnames(refds, 
                                                           devds, quiet)
        varlist = [v for v in cmn3D if 'AOD' in v and '_bin' not in v]

    # Dictionary and list for new display names
    newvars = json_load_file(open(os.path.join(os.path.dirname(__file__),
                                               aod_spc)))
    newvarlist = []

    # =================================================================
    # Compute the total AOD by summing over the constituent members
    # =================================================================

    # Take one of the variables so we can use its dims, coords,
    # attrs to create the DataArray object for total AOD
    v = varlist[0]

    # Create a DataArray to hold total column AOD
    # This is the same shape as the DataArray objects in refds
    reftot = xr.DataArray(np.zeros(refds[v].values.shape),
                          name='AODTotal',
                          dims=refds[v].dims,
                          coords=refds[v].coords,
                          attrs=refds[v].attrs)

    # Create a DataArray to hold total column AOD
    # This is the same shape as the DataArray objects in devds
    devtot = xr.DataArray(np.zeros(devds[v].values.shape),
                          name='AODTotal',
                          dims=devds[v].dims,
                          coords=devds[v].coords,
                          attrs=devds[v].attrs)

    # Save the variable attributes so that we can reattach them
    refattrs = reftot.attrs
    devattrs = devtot.attrs

    # Compute the sum of all AOD variables
    for v in varlist:
        reftot = reftot + refds[v]
        devtot = devtot + devds[v]

    # Reattach the variable attributes
    reftot.name = 'AODTotal'
    reftot.attrs = refattrs
    reftot.attrs['long_name'] = 'Total aerosol optical depth'
    devtot.name = 'AODTotal'
    devtot.attrs = devattrs
    devtot.attrs['long_name'] = 'Total aerosol optical depth'

    # Merge these variables back into the dataset
    refds = xr.merge([refds, reftot])
    devds = xr.merge([devds, devtot])

    # Also add AODTotal to the list
    varlist.append('AODTotal')

    # =================================================================
    # Compute column AODs
    # Create a new DataArray for each column AOD variable,
    # using the new display name, and preserving attributes.
    # Merge the new DataArrays back into the DataSets.
    # =================================================================
    for v in varlist:

        # Get the new name for each AOD variable (it's easier to display)
        if v in newvars:
            newname = newvars[v]
            newvarlist.append(newname)
        else:
            raise ValueError('Could not find a display name for'.format(v))

        # Don't clobber existing DataArray and Dataset attributes
        with xr.set_options(keep_attrs=True):

            # Add column AOD of newname to Ref
            array = refds[v].sum(dim='lev')
            array.name = newname
            refds = xr.merge([refds, array])

            # Add column AOD of newname to Dev
            array = devds[v].sum(dim='lev')
            array.name = newname
            devds = xr.merge([devds, array])

    # =================================================================
    # Create the plots
    # =================================================================
    pdfname = os.path.join(aoddir, 'Aerosols_ColumnOptDepth.pdf')
    diff_aod = []
    compare_single_level(refds, refstr, devds, devstr,
                         varlist=newvarlist, ilev=0, pdfname=pdfname,
                         log_color_scale=log_color_scale,
                         sigdiff_list=diff_aod)
    diff_aod[:] = [v.replace('Column_AOD_', '') for v in diff_aod]
    add_bookmarks_to_pdf(pdfname, newvarlist,
                         remove_prefix='Column_AOD_', verbose=verbose)

    # =================================================================
    # Write the list of AOD quantities having significant differences,
    # which we will need to fill out the benchmark forms.
    # =================================================================
    for filename in sigdiff_files:
        if 'sfc' in filename:
            with open(filename, 'a+') as f:
                print('* Column AOD: ', file=f, end='')
                for v in diff_aod:
                    print('{} '.format(v), file=f, end='')
                print(file=f)
                f.close()


def create_budget_table(devdata, devstr, region, species, varnames,
                        outfilename, interval=2678400.0, template='Budget_{}'):
    '''
    Creates a table of budgets by species and component for a data set.

    Args:
        devdata : xarray Dataset
            The second data set to be compared (aka "Development").

        devstr: str
            A string that can be used to identify the data set specified
            by devfile (e.g. a model version number or other identifier).

        region : str
            Name of region for which budget will be computed.

        species : List of strings
            List of species  to include in budget tables.

        varnames : List of strings
            List of variable names in the budget diagnostics.

        outfilename : str
            Name of the text file which will contain the table of
            emissions totals.

    Keyword Args (optional):
        interval : float
            The length of the data interval in seconds. By default, interval
            is set to the number of seconds in a 31-day month (86400 * 31),
            which corresponds to typical benchmark simulation output.

        template : str
            Template for the diagnostic names that are contained in the
            data set. If not specified, template will be set to "Budget_{}",
            where {} will be replaced by the species name.

    Returns:
        None.  Writes to a text file specified by the "outfilename" argument.

    Remarks:
        This method is mainly intended for model benchmarking purposes,
        rather than as a general-purpose tool.
    '''

    # Error check arguments
    if not isinstance(devdata, xr.Dataset):
        raise ValueError('The devdata argument must be an xarray Dataset!')

    # Open file for output
    f = open(outfilename, 'w')

    for spc_name in species:

        # Title string
        title = '### {} budget totals for species {}'.format(devstr,spc_name)
    
        # Write header to file
        print('#'*79, file=f)
        print('{}{}'.format(title.ljust(76), '###'), file=f)
        print('#'*79, file=f)

        # Get variable names for this species
        spc_vars = [ v for v in varnames if v.endswith('_'+spc_name) ]

        for v in spc_vars:

            # Component name
            comp_name = v.replace('Budget', '')
            comp_name = comp_name.replace('_'+spc_name, '')
            comp_name = comp_name.replace(region, '')
        
            # Convert from kg/s to Tg
            devarray = devdata[v] * interval * 1e-9
            units    = 'Tg'
        
            # Compute sum
            total_dev = np.sum(devarray.values)

            # Write output
            print('{} : {:13.6e} {}'.format(comp_name.ljust(12),
                                            total_dev, units), file=f)

        # Add new lines before going to the next species
        print(file=f)
        print(file=f)

    # Close file
    f.close()

    
def make_benchmark_budget_tables(devlist, devstr, dst='./1mo_benchmark',
                                 overwrite=False, interval=None):
    '''
    Creates a text file containing budgets by species for benchmarking
    purposes.

    Args:
        devlist : list of str
             List with the path names of the emissions and/or met field
             files that will constitute the "Dev" (aka "Development") 
             data set.  The "Dev" data set will be compared against the
             "Ref" data set.

        devstr : str
             A string to describe dev (e.g. version number)

    Keyword Args (optional):
        dst : str
             A string denoting the destination folder where the file
             containing emissions totals will be written.
             Default value: ./1mo_benchmark

        overwrite : boolean
             Set this flag to True to overwrite files in the
             destination folder (specified by the dst argument).
             Default value : False

        interval : float
             Specifies the averaging period in seconds, which is used
             to convert fluxes (e.g. kg/m2/s) to masses (e.g kg).
             Default value : None
    '''

    # ===============================================================
    # Define destination directory
    # ===============================================================
    if os.path.isdir(dst) and not overwrite:
        print('Directory {} exists. Pass overwrite=True to overwrite files in that directory, if any.'.format(dst))
        return
    elif not os.path.isdir(dst):
        os.mkdir(dst)
    budgetdir = os.path.join(dst,'Budget')
    if not os.path.isdir(budgetdir):
        os.mkdir(budgetdir)

    # ===============================================================
    # Read data from netCDF into Dataset objects
    # ===============================================================

    # Dev
    try:
        devds = xr.open_mfdataset(devlist)
    except FileNotFoundError:
        print('Could not find one of the Dev files: {}'.format(devlist))
        raise

    # ===============================================================
    # Create budget table
    # ===============================================================

    # If the averaging interval (in seconds) is not specified,
    # then assume July 2016 = 86400 seconds * 31 days
    if interval == None:
        interval = 86400.0 * 31.0
    
    # Get budget variable and regions
    budget_vars    = [ k for k in devds.data_vars.keys() if k[:6] == 'Budget' ]
    budget_regions = sorted(set([v.split('_')[0][-4:] for v in budget_vars]))

    for region in budget_regions:
    
        # Destination file
        file_budget = os.path.join( dst, budgetdir,
                                        'Budget_'+region+'.txt')

        # Get variable names and species for this region
        region_vars = [ k for k in budget_vars if region in k ]
        region_spc  = sorted(set([v.split('_')[1] for v in region_vars]))
        
        # Write to file
        create_budget_table(devds, devstr, region, region_spc, region_vars,
                            file_budget, interval, template='Budget_{}')
                

def add_bookmarks_to_pdf(pdfname, varlist, remove_prefix='', verbose=False ):
    '''
    Adds bookmarks to an existing PDF file.

    Args:
        pdfname : str
            Name of an existing PDF file of species or emission plots
            to which bookmarks will be attached.

        varlist : list
            List of variables, which will be used to create the
            PDF bookmark names.

    Keyword Args (optional):
        remove_prefix : str
            Specifies a prefix to remove from each entry in varlist
            when creating bookmarks.  For example, if varlist has
            a variable name "SpeciesConc_NO", and you specify
            remove_prefix="SpeciesConc_", then the bookmark for
            that variable will be just "NO", etc.

         verbose : boolean
             Set this flag to True to print extra informational output.
             Default value: False
    '''

    # Setup
    pdfobj = open(pdfname,'rb')
    input = PdfFileReader(pdfobj)
    output = PdfFileWriter()
    
    for i, varname in enumerate(varlist):
        bookmarkname = varname.replace(remove_prefix,'')
        if verbose: print('Adding bookmark for {} with name {}'.format(varname, bookmarkname))
        output.addPage(input.getPage(i))
        output.addBookmark(bookmarkname,i)
        output.setPageMode('/UseOutlines')
        
    # Write to temp file
    pdfname_tmp = pdfname+'_with_bookmarks.pdf'
    outputstream = open(pdfname_tmp,'wb')
    output.write(outputstream) 
    outputstream.close()
    
    # Rename temp file with the target name
    os.rename(pdfname_tmp, pdfname)

      
def add_nested_bookmarks_to_pdf( pdfname, category, catdict, warninglist, remove_prefix=''):

    '''
    Add nested bookmarks to PDF.

    Args:
        pdfname : str
             Path of PDF to add bookmarks to

        category : str
             Top-level key name in catdict that maps to contents of PDF

        catdict : dictionary
             Dictionary containing key-value pairs where one top-level key matches
             category and has value fully describing pages in PDF. The value is a 
             dictionary where keys are level 1 bookmark names, and values are
             lists of level 2 bookmark names, with one level 2 name per PDF page. 
             Level 2 names must appear in catdict in the same order as in the PDF.

        warninglist : list of strings
             Level 2 bookmark names to skip since not present in PDF.

    Keyword Args (optional):
        remove_prefix : str
             Prefix to be remove from warninglist names before comparing with 
             level 2 bookmark names in catdict.
             Default value: empty string (warninglist names match names in catdict) 
    '''
    
    # Setup
    pdfobj = open(pdfname,'rb')
    input = PdfFileReader(pdfobj)
    output = PdfFileWriter()
    warninglist = [k.replace(remove_prefix,'') for k in warninglist]
        
    # Loop over the subcategories in this category; make parent bookmark
    i = -1
    for subcat in catdict[category]:

        # First check that there are actual variables for this subcategory; otherwise skip
        numvars = 0
        if catdict[category][subcat]:
            for varname in catdict[category][subcat]:
                if varname in warninglist:
                    continue
                else:
                    numvars += 1
        else:
                continue
        if numvars == 0:
            continue

        # There are non-zero variables to plot in this subcategory
        i = i+1
        output.addPage(input.getPage(i))
        parent = output.addBookmark(subcat,i)
        output.setPageMode('/UseOutlines')
        first = True
        
        # Loop over variables in this subcategory; make children bookmarks
        for varname in catdict[category][subcat]:
            if varname in warninglist:
                print('Warning: skipping {}'.format(varname))
                continue
            if first:
                output.addBookmark(varname, i, parent)
                first = False
            else:
                i = i+1
                output.addPage(input.getPage(i))
                output.addBookmark(varname, i, parent)
                output.setPageMode('/UseOutlines')
        
    # Write to temp file
    pdfname_tmp = pdfname+'_with_bookmarks.pdf'
    outputstream = open(pdfname_tmp,'wb')
    output.write(outputstream) 
    outputstream.close()

    # Rename temp file with the target name
    os.rename(pdfname_tmp, pdfname)


def normalize_colors(vmin, vmax,
                     is_all_zeroes=False,
                     log_color_scale=False):
    '''
    Normalizes colors to the range of 0..1 for input to
    matplotlib-based plotting functions, given the max & min values.

    For log-color scales, special handling is done to prevent
    taking the log of data that is all zeroes.

    Args:
    -----
    vmin, vmax : float
        Minimum and maximum values of a data array.

    Keyword Args:
    -------------
    is_all_zeroes : boolean
        Logical flag to denote if the data array is all zeroes.
        Default value: False

    log_color_scale : boolean
        Logical flag to denote that we are using a logarithmic
        color scale instead of a linear color scale.
        Default value: False:

    Returns:
    --------
    norm : matplotlib Norm
        The normalized colors, in a matplotlib Norm object.

    Remarks:
    --------
    (1) This is an internal routine, called from compare_single_level,
        and compare_zonal_mean.

    (2) For log scale, the min value will be 3 orders of magnitude
        less than the max.  But if the data is all zeroes, then we
        will revert to a linear scale (to avoid a math error of
        taking a log of zero).
    '''
    if is_all_zeroes:
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    else:
        if log_color_scale:
            norm = mpl.colors.LogNorm(vmin=vmax/1e3, vmax=vmax)
        else:
            norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

    return norm


# =============================================================================
# The rest of this file contains legacy code that may be deleted in the future
# =============================================================================

def plot_layer(dr, ax, title='', unit='', diff=False, vmin=None, vmax=None):
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

    NOTE: This is deprecated and will be removed in the future.
    '''

    assert isinstance(ax, GeoAxes), (
           'Input axis must be cartopy GeoAxes! '
           'Can be created by: \n'
           'plt.axes(projection=ccrs.PlateCarree()) \n or \n'
           'plt.subplots(n, m, subplot_kw={"projection": ccrs.PlateCarree()})'
           )
    assert ax.projection == ccrs.PlateCarree(), (
           'must use PlateCarree projection'
           )

    fig = ax.figure  # get parent figure

    if vmax == None and vmin == None:
        if diff:
            vmax = np.max(np.abs(dr.values))
            vmin = -vmax
            cmap = cmap_diff
        else:
            vmax = np.max(dr.values)
            vmin = 0
            cmap = cmap_abs

    if diff:
        cmap = cmap_diff
    else:
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

    NOTE: This is deprecated and will be removed in the future.
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

    # NOTE: The surface has a hybrid eta coordinate of 1.0 and the
    # atmosphere top has a hybrid eta coordinate of 0.0.  If we don't
    # flip the Y-axis, then the surface will be plotted at the top
    # of the plot. (bmy, 3/7/18)
    ax.invert_yaxis()

    # can also pass cbar_kwargs to dr.plot() to add colorbar
    # but it is easier to tweak colorbar afterwards
    cb = fig.colorbar(im, ax=ax, shrink=0.6, orientation='horizontal', pad=0.1)
    cb.set_label(unit)

    ax.set_title(title)


def make_pdf(ds1, ds2, filename, on_map=True, diff=False,
             title1='DataSet 1', title2='DataSet 2', unit='',
             matchcbar=False):
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

    NOTE: This is deprecated and will be removed in the future.
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
            unitmatch = False
            if matchcbar:
                vmin = min([ds1[varname].data.min(), ds2[varname].data.min()])
                vmax = max([ds1[varname].data.max(), ds2[varname].data.max()])
                unitmatch = ds1[varname].units == ds2[varname].units

            for j, ds in enumerate([ds1, ds2]):
                if on_map:
                    if matchcbar and unitmatch:
                        plot_func(ds[varname], axes[i][j], 
                                  unit=ds[varname].units, diff=diff,
                                  vmin=vmin, vmax=vmax)
                    else:
                        plot_func(ds[varname], axes[i][j], 
                                  unit=ds[varname].units, diff=diff)
                else:
                    # For now, assume zonal mean if plotting zonal (ewl)
                    if matchcbar and unitmatch:
                        plot_func(ds[varname].mean(axis=2), axes[i][j], 
                                  unit=ds[varname].units, diff=diff,
                                  vmin=vmin, vmax=vmax)
                    else:
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

