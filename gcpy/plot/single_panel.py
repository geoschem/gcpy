"""
Creates a single panel plot (geographic map or zonal mean).
"""
import copy
from matplotlib import ticker
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from dask.array import Array as DaskArray
import xarray as xr
import cartopy.crs as ccrs
from gcpy.grid import get_vert_grid, get_pressure_indices, \
    pad_pressure_edges, convert_lev_to_pres, get_grid_extents, \
    call_make_grid, get_input_res
from gcpy.regrid import regrid_comparison_data, create_regridders
from gcpy.util import reshape_MAPL_CS, all_zero_or_nan, verify_variable_type
from gcpy.plot.core  import gcpy_style, normalize_colors, WhGrYlRd

# Suppress numpy divide by zero warnings to prevent output spam
np.seterr(divide="ignore", invalid="ignore")

# Use a style sheet to control plot attributes
plt.style.use(gcpy_style)


def single_panel(
        plot_vals,
        ax=None,
        plot_type="single_level",
        grid=None,
        gridtype="",
        title="fill",
        comap=WhGrYlRd,
        norm=None,
        unit="",
        extent=None,
        masked_data=None,
        use_cmap_RdBu=False,
        log_color_scale=False,
        add_cb=True,
        pres_range=None,
        pedge=np.full((1, 1), -1),
        pedge_ind=np.full((1, 1), -1),
        log_yaxis=False,
        xtick_positions=None,
        xticklabels=None,
        proj=ccrs.PlateCarree(),
        sg_path='',
        ll_plot_func="imshow",
        vert_params=None,
        pdfname="",
        weightsdir='.',
        vmin=None,
        vmax=None,
        return_list_of_plots=False,
        **extra_plot_args
):
    """
    Core plotting routine -- creates a single plot panel.

    Args:
        plot_vals: xarray.DataArray, numpy.ndarray, or dask.array.Array
            Single data variable GEOS-Chem output to plot

    Keyword Args (Optional):
        ax: matplotlib axes
            Axes object to plot information
            Default value: None (Will create a new axes)
        plot_type: str
            Either "single_level" or "zonal_mean"
            Default value: "single_level"
        grid: dict
            Dictionary mapping plot_vals to plottable coordinates
            Default value: {} (will attempt to read grid from plot_vals)
        gridtype: str
            "ll" for lat/lon or "cs" for cubed-sphere
            Default value: "" (will automatically determine from grid)
        title: str
            Title to put at top of plot
            Default value: "fill" (will use name attribute of plot_vals
            if available)
        comap: matplotlib Colormap
            Colormap for plotting data values
            Default value: WhGrYlRd
        norm: list
            List with range [0..1] normalizing color range for matplotlib
            methods. Default value: None (will determine from plot_vals)
        unit: str
            Units of plotted data
            Default value: "" (will use units attribute of plot_vals
            if available)
        extent: tuple (minlon, maxlon, minlat, maxlat)
            Describes minimum and maximum latitude and longitude of input
            data.  Default value: None (Will use full extent of plot_vals
            if plot is single level).
        masked_data: numpy array
            Masked area for avoiding near-dateline cubed-sphere plotting
            issues  Default value: None (will attempt to determine from
            plot_vals)
        use_cmap_RdBu: bool
            Set this flag to True to use a blue-white-red colormap
            Default value: False
        log_color_scale: bool
            Set this flag to True to use a log-scale colormap
            Default value: False
        add_cb: bool
            Set this flag to True to add a colorbar to the plot
            Default value: True
        pres_range: list(int)
            Range from minimum to maximum pressure for zonal mean
            plotting. Default value: [0, 2000] (will plot entire
            atmosphere)
        pedge: numpy array
            Edge pressures of vertical grid cells in plot_vals
            for zonal mean plotting.  Default value: np.full((1, 1), -1)
            (will determine automatically)
        pedge_ind: numpy array
            Index of edge pressure values within pressure range in
            plot_vals for zonal mean plotting.
            Default value: np.full((1, 1), -1) (will determine
            automatically)
        log_yaxis: bool
            Set this flag to True to enable log scaling of pressure in
            zonal mean plots.  Default value: False
        xtick_positions: list(float)
            Locations of lat/lon or lon ticks on plot
            Default value: None (will place automatically for
            zonal mean plots)
        xticklabels: list(str)
            Labels for lat/lon ticks
            Default value: None (will determine automatically from
            xtick_positions)
        proj: cartopy projection
            Projection for plotting data
            Default value: ccrs.PlateCarree()
        sg_path: str
            Path to NetCDF file containing stretched-grid info
            (in attributes) for plot_vals.
            Default value: '' (will not be read in)
        ll_plot_func: str
            Function to use for lat/lon single level plotting with
            possible values 'imshow' and 'pcolormesh'. imshow is much
            faster but is slightly displaced when plotting from dateline
            to dateline and/or pole to pole.  Default value: 'imshow'
        vert_params: list(AP, BP) of list-like types
            Hybrid grid parameter A in hPa and B (unitless). Needed if
            grid is not 47 or 72 levels.  Default value: None
        pdfname: str
            File path to save plots as PDF
            Default value: "" (will not create PDF)
        weightsdir: str
            Directory path for storing regridding weights
            Default value: "." (will store regridding files in
            current directory)
        vmin: float
            minimum for colorbars
            Default value: None (will use plot value minimum)
        vmax: float
            maximum for colorbars
            Default value: None (will use plot value maximum)
        return_list_of_plots: bool
            Return plots as a list. This is helpful if you are using
            a cubedsphere grid and would like access to all 6 plots
            Default value: False
        extra_plot_args: various
            Any extra keyword arguments are passed to calls to
            pcolormesh() (CS) or imshow() (Lat/Lon).

    Returns:
        plot: matplotlib plot
            Plot object created from input
    """
    verify_variable_type(plot_vals, (xr.DataArray, np.ndarray, DaskArray))

    # Create empty lists for keyword arguments
    if pres_range is None:
        pres_range = [0, 2000]
    if vert_params is None:
        vert_params = [[], []]

    # Eliminate 1D level or time dimensions
    plot_vals = plot_vals.squeeze()
    data_is_xr = isinstance(plot_vals, xr.DataArray)
    if xtick_positions is None:
        xtick_positions = []
        if plot_type == "zonal_mean":
            xtick_positions = np.arange(-90, 90, 30)

    if xticklabels is None:
        xticklabels = [rf"{x}$\degree$" for x in xtick_positions]

    if unit == "" and data_is_xr:
        try:
            unit = plot_vals.units.strip()
        except BaseException:
            pass

    if title == "fill" and data_is_xr:
        try:
            title = plot_vals.name
        except BaseException:
            pass
    # Generate grid if not passed
    if grid is None:
        res, gridtype = get_input_res(plot_vals)
        sg_params = [1, 170, -90]
        if sg_path != '':
            sg_attrs = xr.open_dataset(sg_path).attrs
            sg_params = [
                sg_attrs['stretch_factor'],
                sg_attrs['target_longitude'],
                sg_attrs['target_latitude']]

        if plot_type == 'single_level':
            grid_extent = get_grid_extents(plot_vals)
            [grid, _] = call_make_grid(
                res,
                gridtype,
                in_extent=grid_extent,
                sg_params=sg_params
            )

        else:  # zonal mean
            if np.all(pedge_ind == -1) or np.all(pedge == -1):

                # Get mid-point pressure and edge pressures for this grid
                pedge, pmid, _ = get_vert_grid(plot_vals, *vert_params)

                # Get indexes of pressure subrange (full range is default)
                pedge_ind = get_pressure_indices(pedge, pres_range)

                # Pad edges if subset does not include surface or TOA so data spans
                # entire subrange
                pedge_ind = pad_pressure_edges(
                    pedge_ind, plot_vals.sizes["lev"], len(pmid))

                # pmid indexes do not include last pedge index
                pmid_ind = pedge_ind[:-1]
                # Convert levels to pressures in ref and dev data
                plot_vals = convert_lev_to_pres(plot_vals, pmid, pedge)
                # get proper levels
                plot_vals = plot_vals.isel(lev=pmid_ind)

            [input_res, input_gridtype, _, _,
             _, new_gridtype, regrid, _, _, _, _,
             grid, regridder, _, regridder_list, _] = create_regridders(
                plot_vals,
                plot_vals,
                weightsdir=weightsdir,
                cmpres=None,
                zm=True,
                sg_ref_params=sg_params
            )
            if gridtype == 'cs':
                plot_vals = reshape_MAPL_CS(plot_vals)
                nlev = len(plot_vals['lev'])
                # Ref
                plot_vals = regrid_comparison_data(
                    plot_vals,
                    input_res,
                    regrid,
                    regridder,
                    regridder_list,
                    grid,
                    input_gridtype,
                    new_gridtype,
                    nlev=nlev
                )
            # average across longitude bands
            # assume lon dim is index 2 (no time dim) if a numpy array is passed
            lon_ind = 2
            if isinstance(plot_vals, xr.DataArray):
                lon_ind = plot_vals.dims.index('lon')
            # calculate zonal means
            plot_vals = plot_vals.mean(axis=lon_ind)
    if gridtype == "":
        _, gridtype = get_input_res(plot_vals)
    if extent is None or extent == (None, None, None, None):
        extent = get_grid_extents(grid)
        # convert to -180 to 180 grid if needed (necessary if going
        # cross-dateline later)
        if extent[0] > 180 or extent[1] > 180:
            #extent = [((extent[0]+180)%360)-180, ((extent[1]+180)%360)-180, extent[2], extent[3]]
            extent = [extent[0] - 180, extent[1] - 180, extent[2], extent[3]]
        #'''
        #if extent[0] < -180 and 'x' in res:
        #    lon_res = float(res.split('x')[1])
        #    extent = [180,
        #if extent[1] > 180 and 'x' in res:
        #    extent[1] = 180
        #'''
    # Account for cross-dateline extent
    if extent[0] > extent[1]:
        if gridtype == "ll":
            # rearrange data with dateline in the middle instead of prime meridian
            # change extent / grid to where dateline is 0, prime meridian is -180 / 180
            # needed for numpy arrays if doing pcolormesh / imshow, and xarray DataArrays
            # if using imshow
            proj = ccrs.PlateCarree(central_longitude=180)
            if ll_plot_func == "imshow" or \
               not isinstance(plot_vals, xr.DataArray):
                i = 0
                while grid['lon_b'][i] < 0:
                    i = i+1
                plot_vals_holder = copy.deepcopy(plot_vals)
                if not isinstance(plot_vals, xr.DataArray):
                    plot_vals_holder[:,:-i] = plot_vals[:,i:]
                    plot_vals_holder[:,-i:] = plot_vals[:,:i]
                else:
                    plot_vals_holder.values[:,:-i] = plot_vals.values[:,i:]
                    plot_vals_holder.values[:,-i:] = plot_vals.values[:,:i]
                plot_vals = plot_vals_holder
            extent[0] = extent[0] % 360 - 180
            extent[1] = extent[1] % 360 - 180
            grid["lon_b"] = grid["lon_b"] % 360 - 180
            grid["lon"] = grid["lon"] % 360 - 180
            if isinstance(plot_vals, xr.DataArray):
                plot_vals['lon'] = plot_vals['lon'] % 360 - 180
            # realign grid also if doing imshow or using numpy arrays
            if ll_plot_func == "imshow" or \
               not isinstance(plot_vals, xr.DataArray):
                temp_grid = copy.deepcopy(grid)
                temp_grid['lon_b'][:-i] = grid['lon_b'][i:]
                temp_grid['lon_b'][-i:] = grid['lon_b'][:i]
                temp_grid['lon'][:-i] = grid['lon'][i:]
                temp_grid['lon'][-i:] = grid['lon'][:i]
                grid = temp_grid
                if isinstance(plot_vals, xr.DataArray):
                    plot_vals = plot_vals.assign_coords({'lon' : grid['lon']})
        if gridtype == "cs":
            proj = ccrs.PlateCarree(central_longitude=180)
            extent[0] = extent[0] % 360 - 180
            extent[1] = extent[1] % 360 - 180
            grid["lon_b"] = grid["lon_b"] % 360 - 180
            grid["lon"] = grid["lon"] % 360 - 180

    if ax is None:
        if plot_type == "zonal_mean":
            ax = plt.axes()
        if plot_type == "single_level":
            ax = plt.axes(projection=proj)

    fig = plt.gcf()
    data_is_xr = isinstance(plot_vals, xr.DataArray)
    # Normalize colors (put into range [0..1] for matplotlib methods)
    if norm is None:
        if data_is_xr:
            vmin = plot_vals.data.min() if vmin is None else vmin
            vmax = plot_vals.data.max() if vmax is None else vmax
        elif isinstance(plot_vals, np.ndarray):
            vmin = np.min(plot_vals) if vmin is None else vmin
            vmax = np.max(plot_vals) if vmax is None else vmax
        norm = normalize_colors(
            vmin,
            vmax,
            is_difference=use_cmap_RdBu,
            log_color_scale=log_color_scale)

    # Create plot
    ax.set_title(title)
    if plot_type == "zonal_mean":
        # Zonal mean plot
        plot = ax.pcolormesh(
            grid["lat_b"],
            pedge[pedge_ind],
            plot_vals,
            cmap=comap,
            norm=norm,
            **extra_plot_args)
        ax.set_aspect("auto")
        ax.set_ylabel("Pressure (hPa)")
        if log_yaxis:
            ax.set_yscale("log")
            ax.yaxis.set_major_formatter(
                ticker.FuncFormatter(lambda y, _: f"{y:g}")
            )
        ax.invert_yaxis()
        ax.set_xticks(xtick_positions)
        ax.set_xticklabels(xticklabels)

    elif gridtype == "ll":
        if ll_plot_func == 'imshow':
            # Lat/Lon single level
            [minlon, maxlon, minlat, maxlat] = extent
            # expand extent to minimize imshow distortion
            #[dlat,dlon] = list(map(float, res.split('x')))
            dlon = grid['lon'][2] - grid['lon'][1]
            dlat = grid['lat'][2] - grid['lat'][1]

            def get_nearest_extent(val, array, direction, spacing):
                # choose nearest values in grid to desired extent to minimize distortion
                grid_vals = np.asarray(array)
                diff = grid_vals - val
                if direction == 'greater':
                    diff[diff < 0] = np.inf
                    i = diff.argmin()
                    if diff[i] == np.inf:
                        # expand extent to value beyond grid limits if extent
                        # is already > max grid value
                        return grid_vals[(np.abs(grid_vals - val)).argmin()]
                    return grid_vals[i]
                # if direction is not "greater":
                diff[diff > 0] = -np.inf
                i = diff.argmax()
                if diff[i] == -np.inf:
                    # expand extent to value beyond grid limits if
                    # extent is already < min grid value
                    # plot will be distorted if full global to avoid
                    # cartopy issues
                    return grid_vals[(
                        np.abs(grid_vals - val)).argmin()] - spacing
                return max(grid_vals[i], -180)
            closest_minlon = get_nearest_extent(
                minlon, grid['lon_b'], 'less', dlon)
            closest_maxlon = get_nearest_extent(
                maxlon, grid['lon_b'], 'greater', dlon)
            # don't adjust if extent includes poles where points are not evenly
            # spaced anyway
            if np.abs(
                    grid['lat_b'][0] -
                    grid['lat_b'][1]) != np.abs(
                    grid['lat_b'][1] -
                    grid['lat_b'][2]) and minlat < grid['lat_b'][1]:
                closest_minlat = grid['lat_b'][0]
            else:
                closest_minlat = get_nearest_extent(
                    minlat, grid['lat_b'], 'less', dlat)

            if np.abs(grid['lat_b'][-1] - grid['lat_b'][-2]) != \
               np.abs(grid['lat_b'][-2] - grid['lat_b'][-3]) and \
               maxlat > grid['lat_b'][-2]:
                closest_maxlat = grid['lat_b'][-1]
            else:
                closest_maxlat = get_nearest_extent(
                    maxlat, grid['lat_b'], 'greater', dlat)

            extent = [
                closest_minlon,
                closest_maxlon,
                closest_minlat,
                closest_maxlat]
            if isinstance(plot_vals, xr.DataArray):
                # filter data by bounds of extent
                plot_vals = plot_vals.where(
                    plot_vals.lon > closest_minlon,
                    drop=True).where(
                    plot_vals.lon < closest_maxlon,
                    drop=True).where(
                    plot_vals.lat > minlat,
                    drop=True).where(
                    plot_vals.lat < maxlat,
                    drop=True)
            else:
                # filter data by indices of grid
                minlon_i = np.where(grid['lon_b']==closest_minlon)[0]
                if len(minlon_i) == 0:
                    minlon_i = 0
                else:
                    minlon_i = int(minlon_i)
                maxlon_i = np.where(grid['lon_b']==closest_maxlon)[0]
                if len(maxlon_i) == 0:
                    maxlon_i = -1
                else:
                    maxlon_i = int(maxlon_i)
                minlat_i = np.where(grid['lat_b']==closest_minlat)[0]
                if len(minlat_i) == 0:
                    minlat_i = 0
                else:
                    minlat_i = int(minlat_i)
                maxlat_i = np.where(grid['lat_b']==closest_maxlat)[0]
                if len(maxlat_i) == 0:
                    maxlat_i = -1
                else:
                    maxlat_i = int(maxlat_i)
                plot_vals = plot_vals[minlat_i:maxlat_i+1,
                                      minlon_i:maxlon_i+1]
            # Create a lon/lat plot
            plot = ax.imshow(
                plot_vals,
                extent=extent,
                transform=proj,
                cmap=comap,
                norm=norm,
                origin='lower',
                interpolation='nearest',
                **extra_plot_args
            )
        else:
            plot = ax.pcolormesh(
                grid["lon_b"],
                grid["lat_b"],
                plot_vals,
                transform=proj,
                cmap=comap,
                norm=norm,
                **extra_plot_args
            )
        ax.set_extent(extent, crs=proj)
        ax.coastlines()
        ax.set_xticks(xtick_positions)
        ax.set_xticklabels(xticklabels)

    else:
        # Cubed-sphere single level
        try:
            if masked_data is None:
                masked_data = np.ma.masked_where(
                    np.abs(
                        grid["lon"] -
                        180) < 2,
                    plot_vals.data.reshape(
                        6,
                        res,
                        res))
        except ValueError:
            # Comparison of numpy arrays throws errors
            pass
        [minlon, maxlon, minlat, maxlat] = extent
        # Catch issue with plots extending into both the western and eastern
        # hemisphere
        if np.max(grid["lon_b"] > 180):
            grid["lon_b"] = (((grid["lon_b"] + 180) % 360) - 180)

        plots = []
        for j in range(6):
            plot = ax.pcolormesh(
                grid["lon_b"][j, :, :],
                grid["lat_b"][j, :, :],
                masked_data[j, :, :],
                transform=proj,
                cmap=comap,
                norm=norm,
                **extra_plot_args
            )
            plots.append(plot)
        ax.set_extent(extent, crs=proj)
        ax.coastlines()
        ax.set_xticks(xtick_positions)
        ax.set_xticklabels(xticklabels)

    if add_cb:
        cbar = plt.colorbar(plot, ax=ax, orientation="horizontal", pad=0.10)
        cbar.mappable.set_norm(norm)
        if data_is_xr:
            all_zero, all_nan = all_zero_or_nan(plot_vals.values)
        else:
            all_zero, all_nan = all_zero_or_nan(plot_vals)
        if all_zero or all_nan:
            if use_cmap_RdBu:
                cbar.set_ticks([0.0])
            else:
                cbar.set_ticks([0.5])
            if all_nan:
                cbar.set_ticklabels(["Undefined throughout domain"])
            else:
                cbar.set_ticklabels(["Zero throughout domain"])
        else:
            if log_color_scale:
                cbar.formatter = ticker.LogFormatter(base=10)
            else:
                if (vmax - vmin) < 0.1 or (vmax - vmin) > 100:
                    cbar.locator = ticker.MaxNLocator(nbins=4)

        try:
            cbar.formatter.set_useOffset(False)
        except BaseException:
            # not all automatically chosen colorbar formatters properly handle
            # the above method
            pass
        cbar.update_ticks()
        cbar.set_label(unit)

    if pdfname != "":
        pdf = PdfPages(pdfname)
        pdf.savefig(fig)
        pdf.close()

    # in some cases users may wish to get a list of all associated plots
    # eg. cubedsphere grids have six plots associated with them
    if return_list_of_plots:
        return plots if 'plots' in locals() else [plot]
    return plot
