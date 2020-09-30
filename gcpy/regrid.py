''' Functions for creating xesmf regridder objects '''

import os
import xesmf as xe
from .grid import make_grid_LL, make_grid_CS, make_grid_SG, get_input_res, call_make_grid, get_grid_extents
import hashlib
import numpy as np
import xarray as xr
import pandas as pd

def make_regridder_L2L( llres_in, llres_out, weightsdir='.', reuse_weights=False,
                        in_extent=[-180,180,-90,90], out_extent=[-180,180,-90,90]):
    [in_minlon, in_maxlon, in_minlat, in_maxlat] = in_extent
    [out_minlon, out_maxlon, out_minlat, out_maxlat] = out_extent
    llgrid_in = make_grid_LL(llres_in, in_extent, out_extent)
    llgrid_out = make_grid_LL(llres_out, out_extent)
    if in_extent == [-180,180,-90,90] and out_extent == [-180,180,-90,90]:
        weightsfile = os.path.join(weightsdir,'conservative_{}_{}.nc'.format(llres_in, llres_out))
    else:
        in_extent_str = str(in_extent).replace('[', '').replace(']','').replace(', ', 'x')
        out_extent_str = str(out_extent).replace('[', '').replace(']','').replace(', ', 'x')
        weightsfile = os.path.join(weightsdir,'conservative_{}_{}_{}_{}.nc'.format(llres_in, llres_out, 
                                                                                   in_extent_str, out_extent_str))
    try:
        regridder = xe.Regridder(llgrid_in, llgrid_out, method='conservative', filename=weightsfile, reuse_weights=reuse_weights)
    except:
        regridder = xe.Regridder(llgrid_in, llgrid_out, method='conservative', filename=weightsfile, reuse_weights=reuse_weights)
    return regridder

def make_regridder_C2L( csres_in, llres_out, weightsdir='.', reuse_weights=True ):
    csgrid, csgrid_list = make_grid_CS(csres_in)
    llgrid = make_grid_LL(llres_out)
    regridder_list = []
    for i in range(6):
        weightsfile = os.path.join(weightsdir, 'conservative_c{}_{}_{}.nc'.format(str(csres_in), llres_out, str(i)))
        try:
            regridder = xe.Regridder(csgrid_list[i], llgrid, method='conservative', filename=weightsfile, reuse_weights=reuse_weights)
        except:
            regridder = xe.Regridder(csgrid_list[i], llgrid, method='conservative', filename=weightsfile, reuse_weights=reuse_weights)
        regridder_list.append(regridder)
    return regridder_list

def make_regridder_S2S(csres_in, csres_out, sf_in=1, tlat_in=-90, tlon_in=170, 
                       sf_out=1, tlat_out=-90, tlon_out=170, weightsdir='.', verbose=True):
    print('sf_in', sf_in)
    print('sf_out', sf_out)
    print('tlat_in', tlat_in)
    print('tlat_out', tlat_out)
    print('tlon_in', tlon_in)
    print('tlon_out', tlon_out)
    igrid, igrid_list = make_grid_SG(csres_in, stretch_factor=sf_in, target_lat=tlat_in, target_lon=tlon_in)
    ogrid, ogrid_list = make_grid_SG(csres_out, stretch_factor=sf_out, target_lat=tlat_out, target_lon=tlon_out)
    regridder_list = []
    for o_face in range(6):
        regridder_list.append({})
        for i_face in range(6):
            weights_fname = f'conservative_sg{sg_hash(csres_in, sf_in, tlat_in, tlon_in)}_F{i_face}_sg{sg_hash(csres_out, sf_out, tlat_out, tlon_out)}_F{o_face}.nc'
            weights_file = os.path.join(weightsdir, weights_fname)
            reuse_weights = os.path.exists(weights_file)
            try:
                regridder = xe.Regridder(igrid_list[i_face],
                                         ogrid_list[o_face],
                                         method='conservative',
                                         filename=weights_file,
                                         reuse_weights=reuse_weights)
                regridder_list[-1][i_face] = regridder
            except ValueError:
                if verbose:
                    print(f"iface {i_face} doesn't intersect oface {o_face}")
            
    return regridder_list

def create_regridders(refds, devds, weightsdir='.', reuse_weights=True, cmpres=None, zm=False, sg_ref_params=[1, -90, 170], sg_dev_params=[1, -90, 170]):
    #Take two lat/lon or cubed-sphere xarray datasets and regrid them if needed
    refres, refgridtype = get_input_res(refds)
    devres, devgridtype = get_input_res(devds)
    
    refminlon, refmaxlon, refminlat, refmaxlat = get_grid_extents(refds)
    devminlon, devmaxlon, devminlat, devmaxlat = get_grid_extents(devds)
    #choose smallest overlapping area for comparison 
    cmpminlon = max(x for x in [refminlon, devminlon] if x is not None)
    cmpmaxlon = min(x for x in [refmaxlon, devmaxlon] if x is not None)
    cmpminlat = max(x for x in [refminlat, devminlat] if x is not None)
    cmpmaxlat = min(x for x in [refmaxlat, devmaxlat] if x is not None)
    
    ref_extent = [refminlon, refmaxlon, refminlat, refmaxlat]
    cmp_extent = [cmpminlon, cmpmaxlon, cmpminlat, cmpmaxlat]
    dev_extent = [devminlon, devmaxlon, devminlat, devmaxlat]    
    # ==================================================================
    # Determine comparison grid resolution and type (if not passed)
    # ==================================================================
    
    # If no cmpres is passed then choose highest resolution between ref and dev.
    # If one dataset is lat-lon and the other is cubed-sphere, and no comparison
    # grid resolution is passed, then default to 1x1.25. If both cubed-sphere and
    # plotting zonal mean, over-ride to be 1x1.25 lat-lon with a warning
    sg_cmp_params=[1, -90, 170]
    if cmpres == None:
        if refres == devres and refgridtype == "ll":
            cmpres = refres
            cmpgridtype = refgridtype
        elif refgridtype == "ll" and devgridtype == "ll":
            cmpres = min([refres, devres])
            cmpgridtype = refgridtype
        elif refgridtype == "cs" and devgridtype == "cs":
            if zm:
                print("Warning: zonal mean comparison must be lat-lon. Defaulting to 1x1.25")
                cmpres='1x1.25'
                cmpgridtype = "ll"
            elif sg_ref_params!=[] or sg_dev_params!=[]:
                #pick ref grid when a stretched-grid and non-stretched-grid are passed
                cmpres=refres
                cmpgridtype="cs"
                sg_cmp_params=sg_ref_params
            else:
                #pick higher resolution CS grid out of two standard cubed-sphere grids
                cmpres = max([refres, devres])
                cmpgridtype="cs"
        elif refgridtype == "ll" and float(refres.split('x')[0])<1 and float(refres.split('x')[1])<1.25:
            cmpres = refres
            cmpgridtype = "ll"
        elif devgridtype == "ll" and float(devres.split('x')[0])<1 and float(devres.split('x')[1])<1.25:
            cmpres = devres
            cmpgridtype = "ll"
        else:
            cmpres = "1x1.25"
            cmpgridtype = "ll"
    elif "x" in cmpres:
        cmpgridtype = "ll"
    elif zm:
        print("Warning: zonal mean comparison must be lat-lon. Defaulting to 1x1.25")
        cmpres='1x1.25'
        cmpgridtype = "ll"
    elif refgridtype == "cs" and devgridtype == "cs":
        if type(cmpres) is list:
            #stretched-grid resolution
            #first element is cubed-sphere resolution, rest are sg params
            cmpres=int(cmpres[0])
            sg_cmp_params=cmpres[1:]
        else:
            cmpres = int(cmpres)  # must cast to integer for cubed-sphere
        cmpgridtype = "cs"
    else:
        print("Warning: lat/lon to CS regridding not currently implemented. Defaulting to 1x1.25")
        cmpres='1x1.25'
        cmpgridtype = "ll"

    # Determine what, if any, need regridding.
    regridref = refres != cmpres or sg_ref_params!=sg_cmp_params
    regriddev = devres != cmpres or sg_dev_params!=sg_cmp_params
    regridany = regridref or regriddev
    # ==================================================================
    # Make grids (ref, dev, and comparison)
    # ==================================================================
    [refgrid, refgrid_list] = call_make_grid(refres, refgridtype, ref_extent, cmp_extent, sg_ref_params)

    [devgrid, devgrid_list] = call_make_grid(devres, devgridtype, dev_extent, cmp_extent, sg_dev_params)

    [cmpgrid, cmpgrid_list] = call_make_grid(cmpres, cmpgridtype, cmp_extent, cmp_extent, sg_cmp_params)
    
    # =================================================================
    # Make regridders, if applicable
    # TODO: Make CS to CS regridders
    # =================================================================
    refregridder = None
    refregridder_list = None
    devregridder = None
    devregridder_list = None
    if regridref:
        if refgridtype == "ll":
            refregridder = make_regridder_L2L(
                refres, cmpres, weightsdir=weightsdir, reuse_weights=reuse_weights,
                in_extent = ref_extent, out_extent = cmp_extent
            )
        else:
            if cmpgridtype == "cs":
                refregridder_list = make_regridder_S2S(refres, cmpres, *sg_ref_params, *sg_cmp_params,
                                                       weightsdir=weightsdir, verbose=False)
            else:
                refregridder_list = make_regridder_C2L(
                    refres, cmpres, weightsdir=weightsdir, reuse_weights=reuse_weights
                )
    if regriddev:
        if devgridtype == "ll":
            devregridder = make_regridder_L2L(
                devres, cmpres, weightsdir=weightsdir, reuse_weights=reuse_weights,
                in_extent = dev_extent, out_extent = cmp_extent
            )
        else:
            if cmpgridtype == "cs":
                devregridder_list = make_regridder_S2S(devres, cmpres, *sg_dev_params, *sg_cmp_params,
                                                       weightsdir=weightsdir, verbose=False)
            else:
                devregridder_list = make_regridder_C2L(
                    devres, cmpres, weightsdir=weightsdir, reuse_weights=reuse_weights
                )

    print('refres', refres)
    print('refgridtype', refgridtype)
    print('devres', devres)
    print('devgridtype', devgridtype)
    print('cmpres', cmpres)
    print('cmpgridtype', cmpgridtype)
    print('regridref', regridref)
    print('regriddev', regriddev)
    print('sg_ref_params', sg_ref_params)
    print('sg_dev_params', sg_dev_params)
    print('sg_cmp_params', sg_cmp_params)
    return [refres, refgridtype, devres, devgridtype, cmpres, cmpgridtype,
    regridref, regriddev, regridany, refgrid, devgrid, cmpgrid, refregridder, 
    devregridder, refregridder_list, devregridder_list]

def dict_72_to_47():
    """
    Get 72L-to-47L conversion dict which stores weights from 72 levels to 47 levels
    Returns:
    --------
        conv_dict : dict {72L (int) : (47L (int), weight (int))}
              Mapping of 72L to 47L
    """

    conv_dict = {L72 : (xmat_72to47.col[L72], xmat_72to47.data[L72]) for L72 in xmat_72to47.row}
    return conv_dict

def reduce_72_to_47(DataArray, conv_dict, pmid_ind_72, pmid_ind_47):
    """
    Reduce 72 level DataArray to 47 level-equivalent for full or restricted
    pressure ranges.
    Args:
    -----
        DataArray : xarray DataArray
            72 level DataArray
    
        conv_dict : dict
            Mapping of 72L to 47L
        pmid_ind_72 : list(int)
            List of midpoint indices for 72L grid
        pmid_ind_47 : list(int)
            List of midpoint indices for 47L grid
        
    Returns:
    --------
         xarray DataArray
            DataArray now reduced to 47L grid
    """

    #reduce 72 level DataArray to 47 level-equivalent
    #This function works for both full and restricted pressure ranges
    new_shape = list(DataArray.data.shape)
    #assumes first dim is level
    new_shape[0] = len(pmid_ind_47)
    reduced_offset = min(pmid_ind_47)
    reduced_data = np.zeros(new_shape)
    for i in range(0, len(pmid_ind_72)):
        lev = pmid_ind_72[i]
        reduced_data[conv_dict[lev][0]-reduced_offset] = reduced_data[conv_dict[lev][0]-reduced_offset] + DataArray.data[i]*conv_dict[lev][1]
    new_coords = {coord : DataArray.coords[coord].data for coord in DataArray.coords if coord != 'lev'}
    new_coords['lev'] = np.array(pmid_ind_47)
    #GCHP-specific
    if 'lats' in DataArray.coords:
        new_coords['lats'] = (('lon', 'lat'), DataArray.coords['lats'].data)
    if 'lons' in DataArray.coords:
        new_coords['lons'] = (('lon', 'lat'), DataArray.coords['lons'].data)
    return xr.DataArray(reduced_data, dims=tuple([dim for dim in DataArray.dims]),
                        coords = new_coords, attrs = DataArray.attrs)

def regrid_comparison_data(data, res, regrid, regridder, regridder_list, global_cmp_grid, gridtype, cmpgridtype,
                           cmpminlat_ind=0, cmpmaxlat_ind=-2, cmpminlon_ind=0, cmpmaxlon_ind=-2, nlev=1):
    """
    Regrid comparison datasets to cubed-sphere or lat/lon format.
    Args:
    -----
        data : xarray DataArray
            DataArray containing a GEOS-Chem data variable
        res : int
            Cubed-sphere resolution for comparison grid        
        regrid : boolean
            Set to true to regrid dataset
        regridder : xESMF regridder
            Regridder between the original data grid and the comparison grid     
        regridder_list : list(xESMF regridder)
            List of regridders for cubed-sphere data
        global_cmp_grid : xarray DataArray
            Comparison grid    
        gridtype : str
            Type of input data grid (either 'll' or 'cs')       
        cmpgridtype : str
            Type of input data grid (either 'll' or 'cs')       
        cmpminlat_ind : int
            Index of minimum latitude extent for comparison grid
        cmpmaxlat_ind : int
            Index (minus 1) of maximum latitude extent for comparison grid
        cmpminlon_ind : int
            Index of minimum longitude extent for comparison grid
        cmpmaxlon_ind : int
            Index (minus 1) of maximum longitude extent for comparison grid
        nlev : int
            Number of levels of input grid and comparison grid
    Returns:
    --------
        data : xarray DataArray
            Original DataArray regridded to comparison grid (including resolution and extent changes)
    """

    if regrid:
        if gridtype == "ll":
            # regrid ll to ll
            return regridder(data)
        elif cmpgridtype == "ll":
            #CS to LL
            if nlev is 1:
                new_data = np.zeros([global_cmp_grid['lat'].size,
                                     global_cmp_grid['lon'].size])
                data_reshaped = data.data.reshape(6, res, res)
            else:
                new_data = np.zeros([nlev, global_cmp_grid['lat'].size,
                                     global_cmp_grid['lon'].size])
                data_reshaped = data.data.reshape(nlev, 6, res, res).swapaxes(0, 1)
            for j in range(6):
                regridder = regridder_list[j]
                new_data += regridder(data_reshaped[j])
            if nlev is 1:
                #limit to extent of cmpgrid
                return new_data[cmpminlat_ind:cmpmaxlat_ind+1,cmpminlon_ind:cmpmaxlon_ind+1].squeeze()
            else:
                return new_data
        elif cmpgridtype == "cs":

            # Reformat dimensions to T, Z, F, Y, X
            if 'Xdim' in data.dims:
                data_format='diagnostic'
            else:
                data_format='checkpoint'
            new_data = reformat_dims(data, format=data_format, towards_common=True)
            # Drop variables that don't look like fields
            #non_fields = [v for v in new_data.variables.keys() if len(set(new_data[v].dims) - {'T', 'Z', 'F', 'Y', 'X'})>0]
            #new_data = new_data.drop(non_fields)
            # Transpose to T, Z, F, Y, X
            if len(new_data.dims) == 5:
                new_data = new_data.transpose('T', 'Z', 'F', 'Y', 'X')
            elif len(new_data.dims) == 4:
                #no time
                new_data = new_data.transpose('Z', 'F', 'Y', 'X')
            elif len(new_data.dims) == 3:
                #no time or vertical
                new_data = new_data.transpose('F', 'Y', 'X')
            # For each output face, sum regridded input faces            
            print('new_data', new_data)
            oface_datasets = []
            for oface in range(6):
                oface_regridded = []
                for iface, regridder in regridder_list[oface].items():
                    ds_iface = new_data.isel(F=iface)
                    if 'nf' in ds_iface.coords:
                        ds_iface = ds_iface.drop('F')
                    oface_regridded.append(regridder(ds_iface, keep_attrs=True))
                oface_regridded = xr.concat(oface_regridded, dim='intersecting_ifaces').sum('intersecting_ifaces',
                                                                                            keep_attrs=True)
                oface_datasets.append(oface_regridded)
            print('oface  datasets', oface_datasets)
            new_data = xr.concat(oface_datasets, dim='F')

            new_data = new_data.rename({
                'y': 'Y',
                'x': 'X',
            })
            new_data = new_data.drop(['lat', 'lon'])  # lat, lon are from xESMF which we don't want

            #reformat dimensions to previous format
            new_data = reformat_dims(new_data, format=data_format, towards_common=False)
            return new_data
    else:
        return data


def reformat_dims(ds, format, towards_common):

    def unravel_checkpoint_lat(ds_in):
        if type(ds_in) is xr.Dataset:
            cs_res = ds_in.dims['lon']
            assert cs_res == ds_in.dims['lat'] // 6
        else:
            cs_res = ds_in['lon'].size
            assert cs_res == ds_in['lat'].size // 6            
        mi = pd.MultiIndex.from_product([
            np.linspace(1, 6, 6),
            np.linspace(1, cs_res, cs_res)
        ])
        ds_in = ds_in.assign_coords({'lat': mi})
        ds_in = ds_in.unstack('lat')
        return ds_in

    def ravel_checkpoint_lat(ds_out):
        if type(ds) is xr.Dataset:
            cs_res = ds_out.dims['lon']
        else:
            cs_res = ds_out['lon'].size
        ds_out = ds_out.stack(lat=['lat_level_0', 'lat_level_1'])
        ds_out = ds_out.assign_coords({
            'lat': np.linspace(1, 6*cs_res, 6*cs_res)
        })
        return ds_out


    dim_formats = {
        'checkpoint': {
            'unravel': [unravel_checkpoint_lat],
            'ravel': [ravel_checkpoint_lat],
            'rename': {
                'lon': 'X',
                'lat_level_0': 'F',
                'lat_level_1': 'Y',
                'time': 'T',
                'lev': 'Z',
            },
            'transpose': ('time', 'lev', 'lat', 'lon')
        },
        'diagnostic': {
            'rename': {
                'nf': 'F',
                'lev': 'Z',
                'Xdim': 'X',
                'Ydim': 'Y',
                'time': 'T',
            },
            'transpose': ('time', 'lev', 'nf', 'Ydim', 'Xdim')
        }
    }
    if towards_common:
        # Unravel dimensions
        for unravel_callback in dim_formats[format].get('unravel', []):
            ds = unravel_callback(ds)

        # Rename dimensions
        ds = ds.rename(dim_formats[format].get('rename', {}))

        return ds
    else:
        # Reverse rename
        ds = ds.rename({v: k for k, v in dim_formats[format].get('rename', {}).items()})

        # Ravel dimensions
        for ravel_callback in dim_formats[format].get('ravel', []):
            ds = ravel_callback(ds)

        # Transpose            
        if len(ds.dims)==5 or (len(ds.dims)==4 and 'lev' in list(ds.dims) and 'time' in list(ds.dims)):
            #full dim dataset
            ds = ds.transpose(*dim_formats[format].get('transpose', []))
        elif len(ds.dims)==4:
            #single time
            ds = ds.transpose(*dim_formats[format].get('transpose', [])[1:])
        elif len(ds.dims)==3:
            #single level / time
            ds = ds.transpose(*dim_formats[format].get('transpose', [])[2:])
        return ds

def sg_hash(cs_res, stretch_factor: float, target_lat: float, target_lon: float):
    return hashlib.sha1('cs={cs_res},sf={stretch_factor:.5f},tx={target_lon:.5f},ty={target_lat:.5f}'.format(
        stretch_factor=stretch_factor,
        target_lat=target_lat,
        target_lon=target_lon,
        cs_res=cs_res
    ).encode()).hexdigest()[:7]
