''' Core utilities for handling GEOS-Chem data '''

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os

import xarray as xr
import xbpch
import numpy as np
import json
import shutil
import matplotlib.colors as mcolors

# JSON files to read
lumped_spc = 'lumped_species.json'
bpch_to_nc_names = 'bpch_to_nc_names.json'


def open_dataset(filename, **kwargs):
    '''
    Load and decode a dataset from an output file generated by GEOS-Chem.

    This method inspects a GEOS-Chem output file and chooses a way to
    load it into memory as an xarray Dataset. Because two different
    libraries to support BPCH and netCDF outputs, you may need to pass
    additional keyword arguments to the function. See the Examples below.

    Args:
    -----
        filename : str
            Path to a GEOS-Chem output file (netCDF or BPCH format)
            which can be loaded through either xarray or xbpch.
            Note that xarray conventions for netCDF files apply.

    Keyword Args (optional):
    ------------------------
        Additional keyword arguments to be passed directly to
        `xarray.open_dataset` or `xbpch.open_bpchdataset`.

    Returns
    -------
        dataset : xarray.Dataset
            The dataset loaded from the referenced filename.

    See Also
    --------
        xarray.open_dataset
        xbpch.open_bpchdataset
        open_mfdataset

    Examples
    --------
        Open a legacy BPCH output file:

        >>> ds = open_dataset("my_GEOS-Chem_run.bpch",
        ...                   tracerinfo_file='tracerinfo.dat',
        ...                   diaginfo_file='diaginfo.dat')

        Open a netCDF output file, but disable metadata decoding:
        >>> ds = open_dataset("my_GCHP_data.nc",
        ...                   decode_times=False, decode_cf=False)
    '''

    basename, file_extension = os.path.splitext(filename)

    if file_extension == '.bpch':
        _opener = xbpch.open_bpchdataset
    elif file_extension == '.nc':
        _opener = xr.open_dataset
    else:
        raise ValueError('Found unknown file extension ({}); please '
                         'pass a BPCH or netCDF file with extension '
                         '"bpch" or "nc"!'.format(file_extension))

    return _opener(filename, **kwargs)


def open_mfdataset(filenames, concat_dim='time', compat='no_conflicts',
                   preprocess=None, lock=None, **kwargs):
    '''
    Load and decode multiple GEOS-Chem output files as a single Dataset.

    Parameters
    ----------

        filenames : list of str
            Paths to GEOS-Chem output files to load. Must have the same
            extension and be able to be concatenated along some common axis.

        concat_dim : str, default='time'
            Dimension to concatenate Datasets over. We default to "time"
            since this is how GEOS-Chem splits output files.

        compat : {'identical', 'equals', 'broadcast_equals',
                  'no_conflicts'}, optional
            String indicating how to compare variables of the same name for
            potential conflicts when merging:
                - 'broadcast_equals': all values must be equal when
                   variables are broadcast against each other to ensure
                   common dimensions.
                - 'equals': all values and dimensions must be the same.
                - 'identical': all values, dimensions and attributes
                   must be the same.
                - 'no_conflicts': only values which are not null in
                   both datasets must be equal. The returned dataset
                   then contains the combination of all non-null values.

    Keyword Arguments (optional):
    -----------------------------
        preprocess : callable (optional)
            A pre-processing function to apply to each Dataset prior to
            concatenation

        lock : False, True, or threading.Lock (optional)
            Passed to :py:func:`dask.array.from_array`. By default,
            xarray employs a per-variable lock when reading data from
            NetCDF files, but this model has not yet been extended or
            implemented for bpch files and so this is not actually used.
            However, it is likely necessary before dask's multi-threaded
            backend can be used.

        **kwargs
            Additional keyword arguments to be passed directly to
            `xbpch.open_mfbpchdataset` or `xarray.open_mfdataset`.

    Returns
    -------
        dataset : xarray.Dataset
            A dataset containing the data in the specified input filenames.

    See Also
    --------
        xarray.open_mfdataset
        xbpch.open_mfbpchdataset
        open_dataset
    '''

    try:
        test_fn = filenames[0]
    except:
        raise ValueError('Must pass a list with at least one filename')

    basename, file_extension = os.path.splitext(test_fn)

    if file_extension == '.bpch':
        _opener = xbpch.open_mfbpchdataset
    elif file_extension == '.nc':
        _opener = xr.open_mfdataset
    elif file_extension == '.nc4':
        _opener = xr.open_mfdataset
    else:
        raise ValueError('Found unknown file extension ({}); please ' \
                         'pass a BPCH or netCDF file with extension ' \
                         '"bpch" or "nc" or "nc4"'.format(file_extension))
        
    return _opener(filenames, concat_dim=concat_dim, compat=compat,
                   preprocess=preprocess, lock=lock, **kwargs)


def get_gcc_filepath(outputdir, collection, day, time):
    if collection == 'Emissions':
        filepath = os.path.join(outputdir, 'HEMCO_diagnostics.{}{}.nc'.format(day,time))
    else:
        filepath = os.path.join(outputdir, 'GEOSChem.{}.{}_{}z.nc4'.format(collection,day,time))
    return filepath


def get_gchp_filepath(outputdir, collection, day, time):
    filepath = os.path.join(outputdir, 'GCHP.{}.{}_{}z.nc4'.format(collection,day,time))
    return filepath


def check_paths(refpath, devpath):
    '''
    Checks to see if paths to data files exist.

    Args:
    -----
        refpath : str
            Path to the "Reference" data.

        devpath : str
            Path to the "Development" data.
    '''

    if not os.path.exists(refpath):
        print('ERROR! Path 1 does not exist: {}'.format(refpath))
    else:
        print('Path 1 exists: {}'.format(refpath))
    if not os.path.exists(devpath):
        print('ERROR! Path 2 does not exist: {}'.format(devpath))
    else:
        print('Path 2 exists: {}'.format(devpath))

        
def compare_varnames(refdata, devdata, refonly=[], devonly=[], quiet=False):
    '''
    Finds variables that are common to two xarray Dataset objects.

    Args:
    -----
        refdata : xarray Dataset
            The first Dataset to be compared.
            (This is often referred to as the "Reference" Dataset.)

        devdata : xarray Dataset
            The second Dataset to be compared.
            (This is often referred to as the "Development" Dataset.)

    Keyword Args (optional):
    ------------------------
        quiet : boolean
            Set this flag to True if you wish to suppress printing
            informational output to stdout.
            Default value: False

    Returns:
    --------
        commonvars: list of strs
            Variables that are common to both refdata and devdata,
            regardless of dimension.

        commonvarsOther: list of strs
            Variables that are common to refdata and devdata,
            and that do not have lat, lon, and level dimensions.

        commonvars2D: list of strs
            Variables that are common to refdata and devdata,
            and that have lat and lon dimensions, but not level.

        commonvars3D: list of strs
            Variables that are common to refdata and devdata,
            and that have lat, lon, and level dimensions.

        refonly: list of strs
            Plottable variables (i.e. 2D or 3D) that are only
            present in the Ref dataset.

        devonly: list of strs
            Plottable variables (i.e. 2D or 3D) that are only
            present in the Dev dataset.

    Examples:
    ---------
        >>> import gcpy
        >>> import xarray as xr
        >>> refdata = xr.open_dataset("ref_data_file.nc")
        >>> devdata = xr.open_dataset("dev_data_file.nc")
        >>> [commonvars, commonvarsOther, commonvars2D, commonvars3D, refonly, devonly ] = gcpy.compare_varnames(refdata, devdata)
    '''
    refvars = [k for k in refdata.data_vars.keys()]
    devvars= [k for k in devdata.data_vars.keys()]
    commonvars = sorted(list(set(refvars).intersection(set(devvars))))
    refonly = [v for v in refvars if v not in devvars]
    devonly = [v for v in devvars if v not in refvars]
    dimmismatch = [v for v in commonvars if refdata[v].ndim != devdata[v].ndim]
    commonvarsOther = [v for v in commonvars if \
                       (('lat' not in refdata[v].dims or
                         'Xdim' not in refdata[v].dims) and
                        ('lon' not in refdata[v].dims or
                         'Ydim' not in refdata[v].dims) and
                        ('lev' not in refdata[v].dims))]
    commonvars2D = [v for v in commonvars if \
                    (('lat' in refdata[v].dims or
                      'Xdim' in refdata[v].dims) and
                     ('lon' in refdata[v].dims or
                      'Ydim' in refdata[v].dims) and
                     ('lev' not in refdata[v].dims))]
    commonvars3D = [v for v in commonvars if \
                    (('lat' in refdata[v].dims or
                      'Xdim' in refdata[v].dims) and
                     ('lon' in refdata[v].dims or
                      'Ydim' in refdata[v].dims) and
                     ('lev' in refdata[v].dims))]

    # Print information on common and mismatching variables,
    # as well as dimensions
    if quiet == False:
        print('\nComparing variable names in compare_varnames')
        print('{} common variables'.format(len(commonvars)))
        if len(refonly) > 0:
            print('{} variables in ref only (skip)'.format(len(refonly)))
            print('   Variable names: {}'.format(refonly))
        else:
            print('0 variables in ref only')
            if len(devonly) > 0:
                print('{} variables in dev only (skip)'.format(len(devonly)))
                print('   Variable names: {}'.format(devonly))
            else:
                print('0 variables in dev only')
                if len(dimmismatch) > 0:
                    print('{} common variables have different dimensions'.format(len(dimmismatch)))
                    print('   Variable names: {}'.format(dimmismatch))
                else:
                    print('All variables have same dimensions in ref and dev')

    # For safety's sake, remove the 0-D and 1-D variables from
    # refonly and devonly.  This will ensure that refonly and
    # devonly will only contain variables that can be plotted.
    refonly = [v for v in refonly if v not in commonvarsOther]
    devonly = [v for v in devonly if v not in commonvarsOther]

    return [commonvars, commonvarsOther, commonvars2D,
            commonvars3D, refonly, devonly]


def compare_stats(refdata, refstr, devdata, devstr, varname):
    '''
    Prints out global statistics (array sizes, mean, min, max, sum)
    from two xarray Dataset objects.

    Args:
    ----
        refdata : xarray Dataset
            The first Dataset to be compared.
            (This is often referred to as the "Reference" Dataset.)

        refstr : str
            Label for refdata to be used in the printout

        devdata : xarray Dataset
            The second Dataset to be compared.
            (This is often referred to as the "Development" Dataset.)

        devstr : str
            Label for devdata to be used in the printout

        varname : str
            Variable name for which global statistics will be printed out.

    Examples:
    ---------
        >>> import gcpy
        >>> import xarray as xr
        >>> refdata = xr.open_dataset("ref_data_file.nc")
        >>> devdata = xr.open_dataset("dev_data_file.nc")
        >>> gcpy.compare_stats(ds_ref, "Ref", ds_dev, "Dev", "EmisNO2_Anthro")

        Data units:
            Ref:  molec/cm2/s
            Dev:  molec/cm2/s
        Array sizes:
            Ref:  (1, 47, 46, 72)
            Dev:  (1, 47, 46, 72)
        Global stats:
          Mean:
            Ref:  1770774.125
            Dev:  1770774.125
          Min:
            Ref:  0.0
            Dev:  0.0
          Max:
            Ref:  11548288000.0
            Dev:  11548288000.0
          Sum:
            Ref:  275645792256.0
            Dev:  275645792256.0
    '''

    refvar = refdata[varname]
    devvar = devdata[varname]
    units = refdata[varname].units
    print('Data units:')
    print('    {}:  {}'.format(refstr,units))
    print('    {}:  {}'.format(devstr,units))
    print('Array sizes:')
    print('    {}:  {}'.format(refstr,refvar.shape))
    print('    {}:  {}'.format(devstr,devvar.shape))
    print('Global stats:')
    print('  Mean:')
    print('    {}:  {}'.format(refstr,np.round(refvar.values.mean(),20)))
    print('    {}:  {}'.format(devstr,np.round(devvar.values.mean(),20)))
    print('  Min:')
    print('    {}:  {}'.format(refstr,np.round(refvar.values.min(),20)))
    print('    {}:  {}'.format(devstr,np.round(devvar.values.min(),20)))
    print('  Max:')
    print('    {}:  {}'.format(refstr,np.round(refvar.values.max(),20)))
    print('    {}:  {}'.format(devstr,np.round(devvar.values.max(),20)))
    print('  Sum:')
    print('    {}:  {}'.format(refstr,np.round(refvar.values.sum(),20)))
    print('    {}:  {}'.format(devstr,np.round(devvar.values.sum(),20)))

    
def get_collection_data(datadir, collection, day, time):
    datafile = get_gcc_filepath(datadir, collection, day, time)
    if not os.path.exists(datafile):
        print('ERROR! File does not exist: {}'.format(datafile))
    data_ds = xr.open_dataset(datafile)
    return data_ds


def get_gchp_collection_data(datadir, collection, day, time):
    datafile = get_gchp_filepath(datadir, collection, day, time)
    data_ds = xr.open_dataset(datafile)
    return data_ds


def convert_bpch_names_to_netcdf_names(ds, verbose=False):

    '''
    Function to convert the non-standard bpch diagnostic names
    to names used in the GEOS-Chem netCDF diagnostic outputs.
    
    Args:
    -----
        ds : xarray Dataset
            The xarray Dataset object whose names are to be replaced.

    Keyword Args (optional):
    ------------------------
        verbose : boolean
            Set this flag to True to print informational output.
            Default value: False


    Returns:
    --------
        ds_new : xarray Dataset
            A new xarray Dataset object all of the bpch-style
            diagnostic names replaced by GEOS-Chem netCDF names.

    Remarks:
    --------
        To add more diagnostic names, edit the dictionary contained
        in the bpch_to_nc_names.json.

    Examples:
    -----------------
       >>> import gcpy
       >>> ds_new = gcpy.convert_bpch_names_to_netcdf_names(ds)
    '''

    # Names dictionary (key = bpch id, value[0] = netcdf id,
    # value[1] = action to create full name using id)
    # Now read from JSON file (bmy, 4/5/19)
    jsonfile = os.path.join(os.path.dirname(__file__), bpch_to_nc_names)
    names = json.load(open(jsonfile))

    # define some special variable to overwrite above
    special_vars = {'Met_AIRNUMDE' : 'Met_AIRNUMDEN',
                    'Met_UWND'     : 'Met_U',
                    'Met_VWND'     : 'Met_V',
                    'Met_CLDTOP'   : 'Met_CLDTOPS',
                    'Met_GWET'     : 'Met_GWETTOP',
                    'Met_PRECON'   : 'Met_PRECCON',
                    'Met_PREACC'   : 'Met_PRECTOT',
                    'Met_PBL'      : 'Met_PBLH' }

    # Python dictionary for variable name replacement
    old_to_new = {}

    # Loop over all variable names in the data set
    for variable_name in ds.data_vars.keys():

        # Save the original variable name, since this is the name
        # that we actually need to replace in the dataset.
        original_variable_name = variable_name

        # Replace "__" with "_", in variable name (which will get tested
        # against the name sin the JSON file.  This will allow us to
        # replace variable names in files created with BPCH2COARDS.
        if '__' in variable_name:
            variable_name = variable_name.replace('__', '_')

        # Check if name matches anything in dictionary. Give warning if not.
        oldid = ''
        newid = ''
        idaction = ''
        for key in names:
            if key in variable_name:
                if names[key][1] == 'skip':
                    # Verbose output
                    if verbose:
                        print('WARNING: skipping {}'.format(key))
                else:
                    oldid = key
                    newid = names[key][0]
                    idaction = names[key][1]
                break

        # Go to the next line if no definition was found
        if oldid == '' or newid == '' or idaction == '':
            continue

        # If fullname replacement:
        if idaction == 'replace':
            oldvar = oldid
            newvar = newid

            # Update the dictionary of names with this pair
            # Use the original variable name.
            old_to_new.update({original_variable_name : newvar})

        # For all the rest:
        else:
            linearr = variable_name.split('_')
            varstr = linearr[-1]
            oldvar = oldid + varstr

            # These categories use append
            if oldid in ['IJ_AVG_S_', 'RN_DECAY_', 'WETDCV_S_',
                         'WETDLS_S_', 'BXHGHT_S_', 'DAO_3D_S_',
                         'PL_SUL_',   'CV_FLX_S_', 'EW_FLX_S_',
                         'NS_FLX_S_', 'UP_FLX_S_', 'MC_FRC_S_']:
                newvar = newid + '_' +varstr

            # DAO_FLDS
            # Skip certain fields that will cause conflicts w/ netCDF
            elif oldid in 'DAO_FLDS_':
                if oldid in [ 'DAO_FLDS_PS_PBL', 'DAO_FLDS_TROPPRAW' ]:

                    # Verbose output
                    if verbose:
                        print( 'Skipping: {}'.format(oldid) )
                else:
                    newvar = newid + '_' +varstr

            # Special handling for J-values: The bpch variable names all
            # begin with "J" (e.g. JNO, JACET), so we need to strip the first
            # character of the variable name manually (bmy, 4/8/19)
            elif oldid == 'JV_MAP_S_':
                newvar = newid + '_' + varstr[1:]

            # IJ_SOA_S_
            elif oldid == 'IJ_SOA_S_':
               newvar = newid + varstr

            # DRYD_FLX_, DRYD_VEL_
            elif 'DRYD_' in oldid:
                newvar = newid + '_' + varstr[:-2]

            # BIOBSRCE_, BIOFSRCE_, BIOGSRCE_. ANTHSRCE_
            elif oldid in ['BIOBSRCE_', 'BIOFSRCE_',
                           'BIOGSRCE_', 'ANTHSRCE_']:
                newvar = 'Emis' + varstr +'_' + newid
            
            # If nothing found...
            else:
                
                # Verbose output
                if verbose:
                    print('WARNING: Nothing defined for: {}'.
                          format(variable_name))
                continue

            # Overwrite certain variable names
            if newvar in special_vars:
                newvar = special_vars[newvar]

            # Update the dictionary of names with this pair
            old_to_new.update({original_variable_name : newvar})

    # Verbose output
    if verbose:
        print('\nList of bpch names and netCDF names')
        for key in old_to_new:
            print('{} ==> {}'.format(key.ljust(25),old_to_new[key].ljust(40)))

    # Rename the variables in the dataset
    if verbose:
        print( '\nRenaming variables in the data...')
    with xr.set_options(keep_attrs=True):
        ds = ds.rename(name_dict=old_to_new)
    
    # Return the dataset
    return ds


def get_lumped_species_definitions():
    jsonfile = os.path.join(os.path.dirname(__file__), lumped_spc)
    with open(jsonfile, 'r') as f:
        lumped_spc_dict = json.loads(f.read())
    return lumped_spc_dict


def archive_lumped_species_definitions(dst):
    src = os.path.join(os.path.dirname(__file__), lumped_spc)
    print('Archiving {} in {}'.format(lumped_spc, dst))
    shutil.copyfile(src, os.path.join(dst, lumped_spc))

    
def add_lumped_species_to_dataset(ds, lspc_dict={}, lspc_json='',
                                  verbose=False, overwrite=False,
                                  prefix='SpeciesConc_'):

    # Default is to add all benchmark lumped species.
    # Can overwrite by passing a dictionary
    # or a json file path containing one
    assert not (lspc_dict != {} and lspc_json != ''), \
        'Cannot pass both lspc_dict and lspc_json. Choose one only.'
    if lspc_dict == {} and lspc_json == '':
        lspc_dict = get_lumped_species_definitions()
    elif lspc_dict == {} and lspc_json != '':
        with open(lspc_json, 'r') as f:
            lspc_dict = json.loads(f.read())

    for lspc in lspc_dict:
        varname_new = prefix+lspc
        if varname_new in ds.data_vars and overwrite:
            ds.drop(varname_new)
        else:
            assert varname_new not in ds.data_vars, '{} already in dataset. To overwrite pass overwrite=True.'.format(varname_new)
        if verbose:
            print('Creating {}'.format(varname_new))
        for i, spc in enumerate(lspc_dict[lspc]):
            varname = prefix+spc
            if varname not in ds.data_vars:
                print('Warning: {} needed for {} not in dataset.'.format(spc,lspc))
                continue
            if verbose:
                print(' -> adding {} with scale {}'.format(spc,lspc_dict[lspc][spc]))
            if i == 0:
                darr = ds[varname] * lspc_dict[lspc][spc]
                units = ds[varname].units
            else:
                darr = darr + ds[varname] * lspc_dict[lspc][spc]
        darr.name = varname_new
        darr.attrs['units'] = units
        ds = xr.merge([ds,darr])
    return ds


def filter_names(names, text=''):
    '''
    Returns elements in a list that match a given substring.
    Can be used in conjnction with compare_varnames to return a subset
    of variable names pertaining to a given diagnostic type or species.
    
    Args:
    -----
        names: list of str
            Input list of names.

        text: str
            Target text string for restricting the search.
    
    Returns:
    --------
        filtered_names: list of str
            Returns all elements of names that contains the substring
            specified by the "text" argument.  If "text" is omitted,
            then the original contents of names will be returned.
        
    Examples:
    ---------
        Obtain a list of variable names that contain the substrings
        "CO", "NO", and "O3":
        
        >>> import gcpy
        >>> import xarray as xr
        >>> refdata = xr.open_dataset("ref_data_file.nc")
        >>> devdata = xr.open_dataset("dev_data_file.nc")
        >>> [var, varOther, var2D, var3D] = gcpy.compare_varnames(refdata, devdata)
        >>> var_CO = gcpy.filter_names(var, "CO")
        >>> var_NO = gcpy.filter_names(var, "NO")
        >>> var_O3 = gcpy.filter_names(var, "O3")
    '''

    if text != '':
        filtered_names = [k for k in names if text in k]
    else:
        filtered_names = [k for k in names if k]

    return filtered_names
    

def divide_dataset_by_dataarray(ds, dr, varlist=None):
    '''
    Divides variables in an xarray Dataset object by a single DataArray
    object.  Will also make sure that the Dataset variable attributes
    are preserved.

    This method can be useful for certain types of model diagnostics
    that have to be divided by a counter array.  For example, local
    noontime J-value variables in a Dataset can be divided by the
    fraction of time it was local noon in each grid box, etc.

    Args:
    -----
        ds: xarray Dataset
            The Dataset object containing variables to be divided.

        dr: xarray DataArray
            The DataArray object that will be used to divide the
            variables of ds.

    Keyword Args (optional):
    ------------------------
        varlist: list of str
            If passed, then only those variables of ds that are listed
            in varlist will be divided by dr.  Otherwise, all variables
            of ds will be divided by dr.

    Returns:
    --------
        ds_new : xarray Dataset
            A new xarray Dataset object with its variables divided by dr.
    '''

    # -----------------------------
    # Check arguments
    # -----------------------------
    if not isinstance(ds, xr.Dataset):
        raise TypeError('The ds argument must be of type xarray.Dataset!')

    if not isinstance(dr, xr.DataArray):
        raise TypeError('The dr argument must be of type xarray.DataArray!')

    if varlist == None:
        varlist = ds.data_vars.keys()

    # -----------------------------
    # Do the division
    # -----------------------------

    # Keep all Dataset attributes
    with xr.set_options(keep_attrs=True):

        # Loop over variables
        for v in varlist:

            # Divide each variable of ds by dr
            ds[v] = ds[v] / dr

    return ds


def get_dataarray_shape(dr):
    '''
    Convenience routine to return the shape of an xarray DataArray
    object.  Returns the size of each dimension in the preferred
    netCDF ordering (time, lev|ilev, lat, lon).

    Args:
    -----
    dr : xarray DataArray
        The input DataArray object.

    Returns:
    --------
    sizes : tuple of int
        Tuple containing the sizes of each dimension of dr in order:
        (time, lev|ilev, lat|YDim, lon|XDim).
    '''

    # Make sure that dr is an xarray Dataarray
    if not isinstance(dr, xr.DataArray):
        raise ValueError('dr must be of type xarray DataArray!')
    
    # Initialize
    sizes = ()
    dims = ['time', 'lev', 'ilev', 'lat', 'Ydim', 'lon', 'Xdim']

    # Return a tuple with the size of each found dimension
    for d in dims:
        if d in dr.sizes:
            sizes += (dr.sizes[d],)

    return sizes


def get_area_from_dataset(ds):
    '''
    Convenience routine to return the area variable (which is
    usually called "AREA" for GEOS-Chem "Classic" or "Met_AREAM2"
    for GCHP) from an xarray Dataset object.

    Args:
    -----
        ds : xarray Dataset
            The input dataset.

    Returns:
    --------
        area_m2 : xarray DataArray
            The surface area in m2, as found in ds.
    '''
    
    if 'Met_AREAM2' in ds.data_vars.keys():
        return ds['Met_AREAM2']
    elif 'AREA' in ds.data_vars.keys():
        return ds['AREA']
    else:
        msg = 'An area variable ("AREA" or "Met_AREAM2" is missing' + \
              ' from this dataset!'
        raise ValueError(msg)


def get_variables_from_dataset(ds, varlist):
    '''
    Convenience routine to return multiple selected DataArray
    variables from an xarray Dataset.  All variables must be
    found in the Dataset, or else an error will be raised.

    Args:
    -----
        ds : xarray Dataset
            The input dataset.

        varlist : list of str
            List of DataArray variables to extract from ds.

    Returns:
    --------
        ds_subset : xarray Dataset
            A new data set containing only the variables
            that were requested.

    Remarks:
    -------
    Use this routine if you absolutely need all of the requested
    variables to be returned.  Otherwise
    '''

    ds_subset = xr.Dataset()
    for v in varlist:
        if v in ds.data_vars.keys():
            ds_subset = xr.merge([ds_subset, ds[v]])
        else:
            msg = '{} was not found in this dataset!'.format(v)
            raise ValueError(msg)

    return ds_subset


def normalize_colors(vmin, vmax, is_difference=False, log_color_scale=False):
    '''
    Normalizes colors to the range of 0..1 for input to matplotlib-based
    plotting functions, given the max & min values of a data range.

    For log-color scales, special handling is done to prevent
    taking the log of data that is all zeroes.

    Args:
    -----
        vmin : float
            Minimum value of the data range.

        vmax : float
            Maximum value of the data range.

    Keyword Args:
    -------------
        is_difference : boolean
            Set this switch to denote that we are using a difference
            color scale (i.e. with zero in the middle of the range).
            Default value: False

        log_color_scale : boolean
            Logical flag to denote that we are using a logarithmic
            color scale instead of a linear color scale.
            Default value: False

    Returns:
    --------
        norm : matplotlib Norm
            The normalized colors (with range 0..1), stored in
            a matplotlib Norm object.

    Remarks:
    --------
         For log color scales, we will use a range of 3 orders of
         magnitude (i.e. from vmax/1e3 to vmax).
    '''
    if (vmin == 0 and vmax == 0) or (np.isnan(vmin) and np.isnan(vmax)):

        # If the min and max of the data are both zero, then normalize
        # the data range so that the color corresponding to zero (i.e.
        # white) is placed at 0.5 for difference colormaps (like RdBu)
        # and at 0.0 for non-difference colormaps (like WhGrYlRd).
        if is_difference:
            return mcolors.Normalize(vmin=-1.0, vmax=1.0)
        else:
            return mcolors.Normalize(vmin=0.0, vmax=1.0)
            
    else:

        # For log color scales, assume a range 3 orders of magnitude
        # below the maximum value.  Otherwise use a linear scale.
        if log_color_scale:
            return mcolors.LogNorm(vmin=vmax/1e3, vmax=vmax)
        else:
            return mcolors.Normalize(vmin=vmin, vmax=vmax)
