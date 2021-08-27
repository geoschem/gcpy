"""
Internal utilities for helping to manage xarray and numpy
objects used throughout GCPy
"""

import os
import warnings
import shutil
import yaml
import numpy as np
import xarray as xr
from PyPDF2 import PdfFileWriter, PdfFileReader

def convert_lon(
        data,
        dim='lon',
        format='atlantic',
        neg_dateline=True
):
    """
    Convert longitudes from -180..180 to 0..360, or vice-versa.

    Args:
        data: DataArray or Dataset
             The container holding the data to be converted; the dimension
             indicated by 'dim' must be associated with this container

    Keyword Args (optional):
        dim: str
             Name of dimension holding the longitude coordinates
             Default value: 'lon'
        format: str
             Control whether or not to shift from -180..180 to 0..360 (
             ('pacific') or from 0..360 to -180..180 ('atlantic')
             Default value: 'atlantic'
        neg_dateline: logical
             If True, then the international dateline is set to -180
             instead of 180.
             Default value: True

    Returns:
        data, with dimension 'dim' altered according to conversion rule
    """

    data_copy = data.copy()

    lon = data_copy[dim].values
    new_lon = np.empty_like(lon)

    # Tweak offset for rolling the longitudes later
    offset = 0 if neg_dateline else 1

    if format not in ['atlantic', 'pacific']:
        raise ValueError("Cannot convert longitudes for format '{}'; "
                         "please choose one of 'atlantic' or 'pacific'"
                         .format(format))

    # Create a mask to decide how to mutate the longitude values
    if format == 'atlantic':
        mask = lon >= 180 if neg_dateline else lon > 180

        new_lon[mask] = -(360. - lon[mask])
        new_lon[~mask] = lon[~mask]

        roll_len = len(data[dim]) // 2 - offset

    elif format == 'pacific':
        mask = lon < 0.

        new_lon[mask] = lon[mask] + 360.
        new_lon[~mask] = lon[~mask]

        roll_len = -len(data[dim]) // 2 - offset

    # Copy mutated longitude values into copied data container
    data_copy[dim].values = new_lon
    data_copy = data_copy.roll(**{dim: roll_len})

    return data_copy


def get_emissions_varnames(
        commonvars,
        template=None
):
    """
    Will return a list of emissions diagnostic variable names that
    contain a particular search string.

    Args:
        commonvars: list of strs
            A list of commmon variable names from two data sets.
            (This can be obtained with method gcpy.util.compare_varnames)
        template: str
            String template for matching variable names corresponding
            to emission diagnostics by sector
            Default Value: None
    Returns:
        varnames: list of strs
            A list of variable names corresponding to emission
            diagnostics for a given species and sector
    """

    # Make sure the commonvars list has at least one element
    if len(commonvars) == 0:
        raise ValueError("No valid variable names were passed!")

    # Define template for emission diagnostics by sector
    if template is None:
        raise ValueError("The template argument was not passed!")

    # Find all emission diagnostics for the given species
    varnames = filter_names(commonvars, template)

    return varnames


def create_display_name(
        diagnostic_name
):
    """
    Converts a diagnostic name to a more easily digestible name
    that can be used as a plot title or in a table of totals.

    Args:
        diagnostic_name: str
            Name of the diagnostic to be formatted

    Returns:
        display_name: str
            Formatted name that can be used as plot titles or in tables
            of emissions totals.

    Remarks:
        Assumes that diagnostic names will start with either "Emis"
        (for emissions by category) or "Inv" (for emissions by inventory).
        This should be an OK assumption to make since this routine is
        specifically geared towards model benchmarking.
    """

    # Initialize
    display_name = diagnostic_name

    if "SpeciesRst" in display_name:
        display_name = display_name.split("_")[1]

    # Special handling for Inventory totals
    if "INV" in display_name.upper():
        display_name = display_name.replace("_", " ")

    # Replace text
    for v in ["Emis", "EMIS", "emis", "Inv", "INV", "inv"]:
        display_name = display_name.replace(v, "")

    # Replace underscores
    display_name = display_name.replace("_", " ")

    return display_name


def print_totals(
        ref,
        dev,
        f,
        masks=None
):
    """
    Computes and prints Ref and Dev totals (as well as the difference
    Dev - Ref) for two xarray DataArray objects.

    Args:
        ref: xarray DataArray
            The first DataArray to be compared (aka "Reference")
        dev: xarray DataArray
            The second DataArray to be compared (aka "Development")
        f: file
            File object denoting a text file where output will be directed.

    Keyword Args (optional):
        masks: dict of xarray DataArray
            Dictionary containing the tropospheric mask arrays
            for Ref and Dev.  If this keyword argument is passed,
            then print_totals will print tropospheric totals.
            Default value: None (i.e. print whole-atmosphere totals)

    Remarks:
        This is an internal method.  It is meant to be called from method
        create_total_emissions_table or create_global_mass_table instead of
        being called directly.
    """

    # ==================================================================
    # Initialization and error checks
    # ==================================================================

    # Make sure that both Ref and Dev are xarray DataArray objects
    if not isinstance(ref, xr.DataArray):
        raise TypeError("The ref argument must be an xarray DataArray!")
    if not isinstance(dev, xr.DataArray):
        raise TypeError("The dev argument must be an xarray DataArray!")

    # Determine if either Ref or Dev have all NaN values:
    ref_is_all_nan = np.isnan(ref.values).all()
    dev_is_all_nan = np.isnan(dev.values).all()

    # If Ref and Dev do not contain all NaNs, then make sure
    # that Ref and Dev have the same units before proceeding.
    if (not ref_is_all_nan) and (not dev_is_all_nan):
        if ref.units != dev.units:
            msg = 'Ref has units "{}", but Dev array has units "{}"'.format(
                ref.units, dev.units
            )
            raise ValueError(msg)

    # ==================================================================
    # Get the diagnostic name and units
    # ==================================================================
    if dev_is_all_nan:
        diagnostic_name = ref.name
        units = ref.units
    else:
        diagnostic_name = dev.name
        units = dev.units

    # Create the display name by editing the diagnostic name
    display_name = create_display_name(diagnostic_name)

    # Special handling for totals
    if "_TOTAL" in diagnostic_name.upper():
        print("-"*83, file=f)

    # ==================================================================
    # Sum the Ref array (or set to NaN if missing)
    # ==================================================================
    if ref_is_all_nan:
        total_ref = np.nan
    else:
        if masks is None:
            total_ref = np.sum(ref.values)
        else:
            arr = np.ma.masked_array(ref.values, masks["Ref_TropMask"])
            total_ref = np.sum(arr)

    # ==================================================================
    # Sum the Dev array (or set to NaN if missing)
    # ==================================================================
    if dev_is_all_nan:
        total_dev = np.nan
    else:
        if masks is None:
            total_dev = np.sum(dev.values)
        else:
            arr = np.ma.masked_array(dev.values, masks["Dev_TropMask"])
            total_dev = np.sum(arr)

    # ==================================================================
    # Compute differences (or set to NaN if missing)
    # ==================================================================
    if ref_is_all_nan or dev_is_all_nan:
        diff = np.nan
    else:
        diff = total_dev - total_ref

    # ==================================================================
    # Compute % differences (or set to NaN if missing)
    # If ref is very small, near zero, also set the % diff to NaN
    # ==================================================================
    if np.isnan(total_ref) or np.isnan(total_dev):
        pctdiff = np.nan
    else:
        pctdiff = ((total_dev - total_ref) / total_ref) * 100.0
        if total_ref < 1.0e-15:
            pctdiff = np.nan

    # ==================================================================
    # Write output to file
    # ==================================================================
    print(
        "{} : {:18.6f}  {:18.6f}  {:12.6f}  {:8.3f}".format(
            display_name.ljust(18), total_ref, total_dev, diff, pctdiff
        ),
        file=f,
    )


def get_species_categories(
        benchmark_type="FullChemBenchmark"
):
    """
    Returns the list of benchmark categories that each species
    belongs to.  This determines which PDF files will contain the
    plots for the various species.

    Args:
        benchmark_type: str
            Specifies the type of the benchmark (either
            FullChemBenchmark (default) or TransportTracersBenchmark).

    Returns:
        spc_cat_dict: dict
            A nested dictionary of categories (and sub-categories)
            and the species belonging to each.

    NOTE: The benchmark categories are specified in YAML file
    benchmark_species.yml.
    """
    spc_categories = "benchmark_categories.yml"
    yamlfile = os.path.join(os.path.dirname(__file__), spc_categories)
    with open(yamlfile, "r") as f:
        spc_cat_dict = yaml.load(f.read(), Loader=yaml.FullLoader)
    return spc_cat_dict[benchmark_type]


def archive_species_categories(
        dst
):
    """
    Writes the list of benchmark categories to a YAML file
    named "benchmark_species.yml".

    Args:
        dst: str
            Name of the folder where the YAML file containing
            benchmark categories ("benchmark_species.yml")
            will be written.
    """
    spc_categories = "benchmark_categories.yml"
    src = os.path.join(os.path.dirname(__file__), spc_categories)
    copy = os.path.join(dst, spc_categories)
    if not os.path.exists(copy):
        print("\nArchiving {} in {}".format(spc_categories, dst))
        shutil.copyfile(src, copy)


def add_bookmarks_to_pdf(
        pdfname,
        varlist,
        remove_prefix="",
        verbose=False
):
    """
    Adds bookmarks to an existing PDF file.

    Args:
        pdfname: str
            Name of an existing PDF file of species or emission plots
            to which bookmarks will be attached.
        varlist: list
            List of variables, which will be used to create the
            PDF bookmark names.

    Keyword Args (optional):
        remove_prefix: str
            Specifies a prefix to remove from each entry in varlist
            when creating bookmarks.  For example, if varlist has
            a variable name "SpeciesConc_NO", and you specify
            remove_prefix="SpeciesConc_", then the bookmark for
            that variable will be just "NO", etc.
         verbose: bool
            Set this flag to True to print extra informational output.
            Default value: False
    """

    # Setup
    pdfobj = open(pdfname, "rb")
    input_pdf = PdfFileReader(pdfobj, overwriteWarnings=False)
    output_pdf = PdfFileWriter()

    for i, varname in enumerate(varlist):
        bookmarkname = varname.replace(remove_prefix, "")
        if verbose:
            print(
                "Adding bookmark for {} with name {}".format(
                    varname, bookmarkname))
        output_pdf.addPage(input_pdf.getPage(i))
        output_pdf.addBookmark(bookmarkname, i)
        output_pdf.setPageMode("/UseOutlines")

    # Write to temp file
    pdfname_tmp = pdfname + "_with_bookmarks.pdf"
    outputstream = open(pdfname_tmp, "wb")
    output_pdf.write(outputstream)
    outputstream.close()

    # Rename temp file with the target name
    os.rename(pdfname_tmp, pdfname)
    pdfobj.close()


def add_nested_bookmarks_to_pdf(
        pdfname,
        category,
        catdict,
        warninglist,
        remove_prefix=""
):
    """
    Add nested bookmarks to PDF.

    Args:
        pdfname: str
            Path of PDF to add bookmarks to
        category: str
            Top-level key name in catdict that maps to contents of PDF
        catdict: dictionary
            Dictionary containing key-value pairs where one top-level
            key matches category and has value fully describing pages
            in PDF. The value is a dictionary where keys are level 1
            bookmark names, and values are lists of level 2 bookmark
            names, with one level 2 name per PDF page.  Level 2 names
            must appear in catdict in the same order as in the PDF.
        warninglist: list of strings
            Level 2 bookmark names to skip since not present in PDF.

    Keyword Args (optional):
        remove_prefix: str
            Prefix to be remove from warninglist names before comparing with
            level 2 bookmark names in catdict.
            Default value: empty string (warninglist names match names
            in catdict)
    """

    # ==================================================================
    # Setup
    # ==================================================================
    pdfobj = open(pdfname, "rb")
    input_pdf = PdfFileReader(pdfobj, overwriteWarnings=False)
    output_pdf = PdfFileWriter()
    warninglist = [k.replace(remove_prefix, "") for k in warninglist]

    # ==================================================================
    # Loop over the subcategories in this category; make parent bookmark
    # ==================================================================
    i = -1
    for subcat in catdict[category]:

        # First check that there are actual variables for
        # this subcategory; otherwise skip
        numvars = 0
        if catdict[category][subcat]:
            for varname in catdict[category][subcat]:
                if varname in warninglist:
                    continue
                numvars += 1
        else:
            continue
        if numvars == 0:
            continue

        # There are non-zero variables to plot in this subcategory
        i = i + 1
        output_pdf.addPage(input_pdf.getPage(i))
        parent = output_pdf.addBookmark(subcat, i)
        output_pdf.setPageMode("/UseOutlines")
        first = True

        # Loop over variables in this subcategory; make children bookmarks
        for varname in catdict[category][subcat]:
            if varname in warninglist:
                print("Warning: skipping {}".format(varname))
                continue
            if first:
                output_pdf.addBookmark(varname, i, parent)
                first = False
            else:
                i = i + 1
                output_pdf.addPage(input_pdf.getPage(i))
                output_pdf.addBookmark(varname, i, parent)
                output_pdf.setPageMode("/UseOutlines")

    # ==================================================================
    # Write to temp file
    # ==================================================================
    pdfname_tmp = pdfname + "_with_bookmarks.pdf"
    outputstream = open(pdfname_tmp, "wb")
    output_pdf.write(outputstream)
    outputstream.close()

    # Rename temp file with the target name
    os.rename(pdfname_tmp, pdfname)
    pdfobj.close()


def add_missing_variables(
        refdata,
        devdata,
        verbose=False,
        **kwargs
):
    """
    Compares two xarray Datasets, "Ref", and "Dev".  For each variable
    that is present  in "Ref" but not in "Dev", a DataArray of missing
    values (i.e. NaN) will be added to "Dev".  Similarly, for each
    variable that is present in "Dev" but not in "Ref", a DataArray
    of missing values will be added to "Ref".
    This routine is mostly intended for benchmark purposes, so that we
    can represent variables that were removed from a new GEOS-Chem
    version by missing values in the benchmark plots.
    NOTE: This function assuming incoming datasets have the same sizes and
    dimensions, which is not true if comparing datasets with different grid
    resolutions or types.

    Args:
        refdata: xarray Dataset
            The "Reference" (aka "Ref") dataset
        devdata: xarray Dataset
            The "Development" (aka "Dev") dataset

    Keyword Args (optional):
        verbose: bool
            Toggles extra debug print output
            Default value: False

    Returns:
        refdata, devdata: xarray Datasets
            The returned "Ref" and "Dev" datasets, with
            placeholder missing value variables added
    """
    # ==================================================================
    # Initialize
    # ==================================================================

    # Make sure that refdata and devdata are both xarray Dataset objects
    if not isinstance(refdata, xr.Dataset):
        raise TypeError("The refdata object must be an xarray Dataset!")
    if not isinstance(devdata, xr.Dataset):
        raise TypeError("The refdata object must be an xarray Dataset!")

    # Find common variables as well as variables only in one or the other
    vardict = compare_varnames(refdata, devdata, quiet=True)
    refonly = vardict["refonly"]
    devonly = vardict["devonly"]
    # Don't clobber any DataArray attributes
    with xr.set_options(keep_attrs=True):

        # ==============================================================
        # For each variable that is in refdata but not in devdata,
        # add a new DataArray to devdata with the same sizes but
        # containing all NaN's.  This will allow us to represent those
        # variables as missing values # when we plot against refdata.
        # ==============================================================
        devlist = [devdata]
        for v in refonly:
            if verbose:
                print("Creating array of NaN in devdata for: {}".format(v))
            dr = create_dataarray_of_nan(
                name=refdata[v].name,
                sizes=devdata.sizes,
                coords=devdata.coords,
                attrs=refdata[v].attrs,
                **kwargs
            )
            devlist.append(dr)
        devdata = xr.merge(devlist)

        # ==============================================================
        # For each variable that is in devdata but not in refdata,
        # add a new DataArray to refdata with the same sizes but
        # containing all NaN's.  This will allow us to represent those
        # variables as missing values # when we plot against devdata.
        # ==================================================================
        reflist = [refdata]
        for v in devonly:
            if verbose:
                print("Creating array of NaN in refdata for: {}".format(v))
            dr = create_dataarray_of_nan(
                name=devdata[v].name,
                sizes=refdata.sizes,
                coords=refdata.coords,
                attrs=devdata[v].attrs,
                **kwargs
            )
            reflist.append(dr)
        refdata = xr.merge(reflist)

    return refdata, devdata


def reshape_MAPL_CS(
        da
):
    """
    Reshapes data if contains dimensions indicate MAPL v1.0.0+ output
    Args:
        da: xarray DataArray
            Data array variable

    Returns:
        data: xarray DataArray
            Data with dimensions renamed and transposed to match old MAPL format
    """

    # Suppress annoying future warnings for now
    warnings.filterwarnings("ignore", category=FutureWarning)

    if type(da) != np.ndarray:
        vdims = da.dims
        if "nf" in vdims and "Xdim" in vdims and "Ydim" in vdims:
            da = da.stack(lat=("nf", "Ydim"))
            da = da.rename({"Xdim": "lon"})

        if "lev" in da.dims and "time" in da.dims:
            da = da.transpose("time", "lev", "lat", "lon")
        elif "lev" in da.dims:
            da = da.transpose("lev", "lat", "lon")
        elif "time" in da.dims:
            da = da.transpose("time", "lat", "lon")
        else:
            da = da.transpose("lat", "lon")
    return da


def get_diff_of_diffs(
        ref,
        dev
):
    """
    Generate datasets containing differences between two datasets

    Args:
        ref: xarray Dataset
            The "Reference" (aka "Ref") dataset.
        dev: xarray Dataset
            The "Development" (aka "Dev") dataset

    Returns:
         absdiffs: xarray Dataset
            Dataset containing dev-ref values
         fracdiffs: xarray Dataset
            Dataset containing dev/ref values
    """

    # get diff of diffs datasets for 2 datasets
    # limit each pair to be the same type of output (GEOS-Chem Classic or GCHP)
    # and same resolution / extent
    vardict = compare_varnames(ref, dev, quiet=True)
    varlist = vardict["commonvars"]
    # Select only common fields between the Ref and Dev datasets
    ref = ref[varlist]
    dev = dev[varlist]
    if 'nf' not in ref.dims and 'nf' not in dev.dims:
        with xr.set_options(keep_attrs=True):
            absdiffs = dev - ref
            fracdiffs = dev / ref
            for v in dev.data_vars.keys():
                # Ensure the diffs Dataset includes attributes
                absdiffs[v].attrs = dev[v].attrs
                fracdiffs[v].attrs = dev[v].attrs
    elif 'nf' in ref.dims and 'nf' in dev.dims:

        # Include special handling if cubed sphere grid dimension names are different
        # since they changed in MAPL v1.0.0.
        if "lat" in ref.dims and "Xdim" in dev.dims:
            ref_newdimnames = dev.copy()
            for v in dev.data_vars.keys():
                if "Xdim" in dev[v].dims:
                    ref_newdimnames[v].values = ref[v].values.reshape(
                        dev[v].values.shape)
                # NOTE: the reverse conversion is gchp_dev[v].stack(lat=("nf","Ydim")).transpose(
                # "time","lev","lat","Xdim").values

        with xr.set_options(keep_attrs=True):
            absdiffs = dev.copy()
            fracdiffs = dev.copy()
            for v in dev.data_vars.keys():
                if "Xdim" in dev[v].dims or "lat" in dev[v].dims:
                    absdiffs[v].values = dev[v].values - ref[v].values
                    fracdiffs[v].values = dev[v].values / ref[v].values
                    # NOTE: The diffs Datasets are created without variable
                    # attributes; we have to reattach them
                    absdiffs[v].attrs = dev[v].attrs
                    fracdiffs[v].attrs = dev[v].attrs
    else:
        print('Diff-of-diffs plot supports only identical grid types (lat/lon or cubed-sphere)' + \
              ' within each dataset pair')
        raise ValueError

    return absdiffs, fracdiffs


def slice_by_lev_and_time(
        ds,
        varname,
        itime,
        ilev,
        flip
):
    """
    Slice a DataArray by desired time and level.

    Args:
        ds: xarray Dataset
            Dataset containing GEOS-Chem data.
        varname: str
            Variable name for data variable to be sliced
        itime: int
            Index of time by which to slice
        ilev: int
            Index of level by which to slice
        flip: bool
            Whether to flip ilev to be indexed from ground or top of atmosphere

    Returns:
        ds[varname]: xarray DataArray
            DataArray of data variable sliced according to ilev and itime
    """
    # used in compare_single_level and compare_zonal_mean to get dataset slices
    vdims = ds[varname].dims
    if "time" in vdims and "lev" in vdims:
        if flip:
            maxlev_i = len(ds['lev'])-1
            return ds[varname].isel(time=itime, lev=maxlev_i - ilev)
        return ds[varname].isel(time=itime, lev=ilev)
    if ("time" not in vdims or itime == -1) and "lev" in vdims:
        if flip:
            maxlev_i = len(ds['lev'])-1
            return ds[varname].isel(lev=maxlev_i - ilev)
        return ds[varname].isel(lev=ilev)
    if "time" in vdims and "lev" not in vdims and itime != -1:
        return ds[varname].isel(time=itime)
    return ds[varname]


def rename_and_flip_gchp_rst_vars(
        ds
):
    '''
    Transforms a GCHP restart dataset to match GCC names and level convention

    Args:
        ds: xarray Dataset
            Dataset containing GCHP restart file data, such as variables
            SPC_{species}, BXHEIGHT, DELP_DRY, and TropLev, with level
            convention down (level 0 is top-of-atmosphere).

    Returns:
        ds: xarray Dataset
            Dataset containing GCHP restart file data with names and level
            convention matching GCC restart. Variables include
            SpeciesRst_{species}, Met_BXHEIGHT, Met_DELPDRY, and Met_TropLev,
            with level convention up (level 0 is surface).
    '''
    for v in ds.data_vars.keys():
        if v.startswith('SPC_'):
            spc = v.replace('SPC_', '')
            ds = ds.rename({v: 'SpeciesRst_' + spc})
        elif v == 'DELP_DRY':
            ds = ds.rename({"DELP_DRY": "Met_DELPDRY"})
        elif v == 'BXHEIGHT':
            ds = ds.rename({"BXHEIGHT": "Met_BXHEIGHT"})
        elif v == 'TropLev':
            ds = ds.rename({"TropLev": "Met_TropLev"})
    ds = ds.sortby('lev', ascending=False)
    return ds


def dict_diff(
        dict0,
        dict1
):
    """
    Function to take the difference of two dict objects.
    Assumes that both objects have the same keys.

    Args:
        dict0, dict1: dict
            Dictionaries to be subtracted (dict1 - dict0)

    Returns:
        result: dict
            Key-by-key difference of dict1 - dict0
    """
    result = {}
    for key, _ in dict0.items():
        result[key] = dict1[key] - dict0[key]

    return result


def compare_varnames(
        refdata,
        devdata,
        refonly=None,
        devonly=None,
        quiet=False):
    """
    Finds variables that are common to two xarray Dataset objects.

    Args:
        refdata: xarray Dataset
            The first Dataset to be compared.
            (This is often referred to as the "Reference" Dataset.)
        devdata: xarray Dataset
            The second Dataset to be compared.
            (This is often referred to as the "Development" Dataset.)

    Keyword Args (optional):
        quiet: bool
            Set this flag to True if you wish to suppress printing
            informational output to stdout.
            Default value: False

    Returns:
        vardict: dict of lists of str
            Dictionary containing several lists of variable names:
            Key              Value
            -----            -----
            commonvars       List of variables that are common to
                             both refdata and devdata
            commonvarsOther  List of variables that are common
                             to both refdata and devdata, but do
                             not have lat, lon, and/or level
                             dimensions (e.g. index variables).
            commonvars2D     List of variables that are common to
                             common to refdata and devdata, and that
                             have lat and lon dimensions, but not level.
            commonvars3D     List of variables that are common to
                             refdata and devdata, and that have lat,
                             lon, and level dimensions.
            refonly          List of 2D or 3D variables that are only
                             present in refdata.
            devonly          List of 2D or 3D variables that are only
                             present in devdata
    """
    refvars = [k for k in refdata.data_vars.keys()]
    devvars = [k for k in devdata.data_vars.keys()]
    commonvars = sorted(list(set(refvars).intersection(set(devvars))))
    refonly = [v for v in refvars if v not in devvars]
    devonly = [v for v in devvars if v not in refvars]
    dimmismatch = [v for v in commonvars if refdata[v].ndim != devdata[v].ndim]
    commonvarsOther = [
        v
        for v in commonvars
        if (
            ("lat" not in refdata[v].dims or "Xdim" not in refdata[v].dims)
            and ("lon" not in refdata[v].dims or "Ydim" not in refdata[v].dims)
            and ("lev" not in refdata[v].dims)
        )
    ]
    commonvars2D = [
        v
        for v in commonvars
        if (
            ("lat" in refdata[v].dims or "Xdim" in refdata[v].dims)
            and ("lon" in refdata[v].dims or "Ydim" in refdata[v].dims)
            and ("lev" not in refdata[v].dims)
        )
    ]
    commonvars3D = [
        v
        for v in commonvars
        if (
            ("lat" in refdata[v].dims or "Xdim" in refdata[v].dims)
            and ("lon" in refdata[v].dims or "Ydim" in refdata[v].dims)
            and ("lev" in refdata[v].dims)
        )
    ]

    # Print information on common and mismatching variables,
    # as well as dimensions
    if not quiet:
        print("\nComparing variable names in compare_varnames")
        print("{} common variables".format(len(commonvars)))
        if len(refonly) > 0:
            print("{} variables in ref only (skip)".format(len(refonly)))
            print("   Variable names: {}".format(refonly))
        else:
            print("0 variables in ref only")
            if len(devonly) > 0:
                print("{} variables in dev only (skip)".format(len(devonly)))
                print("   Variable names: {}".format(devonly))
            else:
                print("0 variables in dev only")
                if len(dimmismatch) > 0:
                    print(
                        "{} common variables have different dimensions".format(
                            len(dimmismatch)
                        )
                    )
                    print("   Variable names: {}".format(dimmismatch))
                else:
                    print("All variables have same dimensions in ref and dev")

    # For safety's sake, remove the 0-D and 1-D variables from
    # refonly and devonly.  This will ensure that refonly and
    # devonly will only contain variables that can be plotted.
    refonly = [v for v in refonly if v not in commonvarsOther]
    devonly = [v for v in devonly if v not in commonvarsOther]

    return {
        "commonvars": commonvars,
        "commonvarsOther": commonvarsOther,
        "commonvars2D": commonvars2D,
        "commonvars3D": commonvars3D,
        "refonly": refonly,
        "devonly": devonly,
    }


def compare_stats(refdata, refstr, devdata, devstr, varname):
    """
    Prints out global statistics (array sizes, mean, min, max, sum)
    from two xarray Dataset objects.

    Args:
        refdata: xarray Dataset
            The first Dataset to be compared.
            (This is often referred to as the "Reference" Dataset.)
        refstr: str
            Label for refdata to be used in the printout
        devdata: xarray Dataset
            The second Dataset to be compared.
            (This is often referred to as the "Development" Dataset.)
        devstr: str
            Label for devdata to be used in the printout
        varname: str
            Variable name for which global statistics will be printed out.
    """

    refvar = refdata[varname]
    devvar = devdata[varname]
    units = refdata[varname].units
    print("Data units:")
    print("    {}:  {}".format(refstr, units))
    print("    {}:  {}".format(devstr, units))
    print("Array sizes:")
    print("    {}:  {}".format(refstr, refvar.shape))
    print("    {}:  {}".format(devstr, devvar.shape))
    print("Global stats:")
    print("  Mean:")
    print("    {}:  {}".format(refstr, np.round(refvar.values.mean(), 20)))
    print("    {}:  {}".format(devstr, np.round(devvar.values.mean(), 20)))
    print("  Min:")
    print("    {}:  {}".format(refstr, np.round(refvar.values.min(), 20)))
    print("    {}:  {}".format(devstr, np.round(devvar.values.min(), 20)))
    print("  Max:")
    print("    {}:  {}".format(refstr, np.round(refvar.values.max(), 20)))
    print("    {}:  {}".format(devstr, np.round(devvar.values.max(), 20)))
    print("  Sum:")
    print("    {}:  {}".format(refstr, np.round(refvar.values.sum(), 20)))
    print("    {}:  {}".format(devstr, np.round(devvar.values.sum(), 20)))


def convert_bpch_names_to_netcdf_names(
        ds,
        verbose=False
):
    """
    Function to convert the non-standard bpch diagnostic names
    to names used in the GEOS-Chem netCDF diagnostic outputs.

    Args:
        ds: xarray Dataset
            The xarray Dataset object whose names are to be replaced.

    Keyword Args (optional):
        verbose: bool
            Set this flag to True to print informational output.
            Default value: False

    Returns:
        ds_new: xarray Dataset
            A new xarray Dataset object all of the bpch-style
            diagnostic names replaced by GEOS-Chem netCDF names.

    Remarks:
        To add more diagnostic names, edit the dictionary contained
        in the bpch_to_nc_names.yml.
    """

    # Names dictionary (key = bpch id, value[0] = netcdf id,
    # value[1] = action to create full name using id)
    # Now read from YAML file (bmy, 4/5/19)
    bpch_to_nc_names = "bpch_to_nc_names.yml"
    yamlfile = os.path.join(os.path.dirname(__file__), bpch_to_nc_names)
    names = yaml.load(open(yamlfile), Loader=yaml.FullLoader)

    # define some special variable to overwrite above
    special_vars = {
        "Met_AIRNUMDE": "Met_AIRNUMDEN",
        "Met_UWND": "Met_U",
        "Met_VWND": "Met_V",
        "Met_CLDTOP": "Met_CLDTOPS",
        "Met_GWET": "Met_GWETTOP",
        "Met_PRECON": "Met_PRECCON",
        "Met_PREACC": "Met_PRECTOT",
        "Met_PBL": "Met_PBLH",
    }

    # Tags for the UVFlux* diagnostics
    uvflux_tags = [
        "187nm",
        "191nm",
        "193nm",
        "196nm",
        "202nm",
        "208nm",
        "211nm",
        "214nm",
        "261nm",
        "267nm",
        "277nm",
        "295nm",
        "303nm",
        "310nm",
        "316nm",
        "333nm",
        "380nm",
        "574nm",
    ]

    # Python dictionary for variable name replacement
    old_to_new = {}

    # Loop over all variable names in the data set
    for variable_name in ds.data_vars.keys():

        # Save the original variable name, since this is the name
        # that we actually need to replace in the dataset.
        original_variable_name = variable_name

        # Replace "__" with "_", in variable name (which will get tested
        # against the name sin the YAML file.  This will allow us to
        # replace variable names in files created with BPCH2COARDS.
        if "__" in variable_name:
            variable_name = variable_name.replace("__", "_")

        # Check if name matches anything in dictionary. Give warning if not.
        oldid = ""
        newid = ""
        idaction = ""
        for key in names:
            if key in variable_name:
                if names[key][1] == "skip":
                    # Verbose output
                    if verbose:
                        print("WARNING: skipping {}".format(key))
                else:
                    oldid = key
                    newid = names[key][0]
                    idaction = names[key][1]
                break

        # Go to the next line if no definition was found
        if oldid == "" or newid == "" or idaction == "":
            continue

        # If fullname replacement:
        if idaction == "replace":
            newvar = newid

            # Update the dictionary of names with this pair
            # Use the original variable name.
            old_to_new.update({original_variable_name: newvar})

        # For all the rest:
        else:
            linearr = variable_name.split("_")
            varstr = linearr[-1]

            # These categories use append
            if oldid in [
                    "IJ_AVG_S_",
                    "RN_DECAY_",
                    "WETDCV_S_",
                    "WETDLS_S_",
                    "BXHGHT_S_",
                    "DAO_3D_S_",
                    "PL_SUL_",
                    "CV_FLX_S_",
                    "EW_FLX_S_",
                    "NS_FLX_S_",
                    "UP_FLX_S_",
                    "MC_FRC_S_",
            ]:
                newvar = newid + "_" + varstr

            # DAO_FLDS
            # Skip certain fields that will cause conflicts w/ netCDF
            elif oldid in "DAO_FLDS_":
                if oldid in ["DAO_FLDS_PS_PBL", "DAO_FLDS_TROPPRAW"]:

                    # Verbose output
                    if verbose:
                        print("Skipping: {}".format(oldid))
                else:
                    newvar = newid + "_" + varstr

            # Special handling for J-values: The bpch variable names all
            # begin with "J" (e.g. JNO, JACET), so we need to strip the first
            # character of the variable name manually (bmy, 4/8/19)
            elif oldid == "JV_MAP_S_":
                newvar = newid + "_" + varstr[1:]

            # IJ_SOA_S_
            elif oldid == "IJ_SOA_S_":
                newvar = newid + varstr

            # DRYD_FLX_, DRYD_VEL_
            elif "DRYD_" in oldid:
                newvar = newid + "_" + varstr[:-2]

            # BIOBSRCE_, BIOFSRCE_, BIOGSRCE_. ANTHSRCE_
            elif oldid in ["BIOBSRCE_", "BIOFSRCE_", "BIOGSRCE_", "ANTHSRCE_"]:
                newvar = "Emis" + varstr + "_" + newid

            # Special handling for UV radiative flux diagnostics:
            # We need to append the bin descriptor to the new name.
            elif "FJX_FLXS" in oldid:
                uvind = int(original_variable_name[-2:]) - 1
                newvar = newid + "_" + uvflux_tags[uvind]

            # If nothing found...
            else:

                # Verbose output
                if verbose:
                    print(
                        "WARNING: Nothing defined for: {}".format(
                            variable_name))
                continue

            # Overwrite certain variable names
            if newvar in special_vars:
                newvar = special_vars[newvar]

            # Update the dictionary of names with this pair
            old_to_new.update({original_variable_name: newvar})

    # Verbose output
    if verbose:
        print("\nList of bpch names and netCDF names")
        for key in old_to_new:
            print("{} ==> {}".format(key.ljust(25), old_to_new[key].ljust(40)))

    # Rename the variables in the dataset
    if verbose:
        print("\nRenaming variables in the data...")
    with xr.set_options(keep_attrs=True):
        ds = ds.rename(name_dict=old_to_new)

    # Return the dataset
    return ds


def get_lumped_species_definitions():
    """
    Returns lumped species definitions from a YAML file.

    Returns:
        lumped_spc_dict : dict of str
            Dictionary of lumped species
    """
    lumped_spc = "lumped_species.yml"
    yamlfile = os.path.join(os.path.dirname(__file__), lumped_spc)
    with open(yamlfile, "r") as f:
        lumped_spc_dict = yaml.load(f.read(), Loader=yaml.FullLoader)
    return lumped_spc_dict


def archive_lumped_species_definitions(
        dst
):
    """
    Archives lumped species definitions to a YAML file.

    Args:
        dst : str
            Name of the folder where the YAML file containing
            benchmark categories ("benchmark_species.yml")
            will be written.
    """
    lumped_spc = "lumped_species.yml"
    src = os.path.join(os.path.dirname(__file__), lumped_spc)
    copy = os.path.join(dst, lumped_spc)
    if not os.path.exists(copy):
        print("\nArchiving {} in {}".format(lumped_spc, dst))
        shutil.copyfile(src, copy)


def add_lumped_species_to_dataset(
        ds,
        lspc_dict={},
        lspc_yaml="",
        verbose=False,
        overwrite=False,
        prefix="SpeciesConc_",
):
    """
    Function to calculate lumped species concentrations and add
    them to an xarray Dataset. Lumped species definitions may be passed
    as a dictionary or a path to a yaml file. If neither is passed then
    the lumped species yaml file stored in gcpy is used. This file is
    customized for use with benchmark simuation SpeciesConc diagnostic
    collection output.

    Args:
        ds: xarray Dataset
            An xarray Dataset object prior to adding lumped species.

    Keyword Args (optional):
        lspc_dict: dictionary
            Dictionary containing list of constituent species and their
            integer scale factors per lumped species.
            Default value: False
        lspc_yaml: str
            Set this flag to True to print informational output.
            Default value: False
        verbose: bool
            Whether to print informational output.
            Default value: False
        overwrite: bool
            Whether to overwrite an existing species dataarray in a dataset
            if it has the same name as a new lumped species. If False and
            overlapping names are found then the function will raise an error.
            Default value: False
        prefix: str
            Prefix to prepend to new lumped species names. This argument is
            also used to extract an existing dataarray in the dataset with
            the correct size and dimensions to use during initialization of
            new lumped species dataarrays.
            Default value: "SpeciesConc_"

    Returns:
        ds_new: xarray Dataset
            A new xarray Dataset object containing all of the original
            species plus new lumped species.
    """

    # Default is to add all benchmark lumped species.
    # Can overwrite by passing a dictionary
    # or a yaml file path containing one
    assert not (
        lspc_dict != {} and lspc_yaml != ""
    ), "Cannot pass both lspc_dict and lspc_yaml. Choose one only."
    if lspc_dict == {} and lspc_yaml == "":
        lspc_dict = get_lumped_species_definitions()
    elif lspc_dict == {} and lspc_yaml != "":
        with open(lspc_yaml, "r") as f:
            lspc_dict = yaml.load(f.read(), Loader=yaml.FullLoader)

    # Make sure attributes are transferred when copying dataset / dataarrays
    with xr.set_options(keep_attrs=True):

        # Get a dummy array to use for initialization
        dummy_darr = None
        for var in ds.data_vars:
            if prefix in var:
                dummy_darr = ds[var]
                break
        if dummy_darr is None:
            msg = "Invalid prefix: " + prefix
            raise ValueError(msg)

        # Create a new dataset equivalent to the old
        ds_new = ds.copy(deep=True)

        # Free memory of original dataset
        ds = xr.Dataset()

        # Loop over lumped species
        for lspc in lspc_dict:

            # Assemble lumped species variable name
            varname_new = prefix + lspc

            # Check if overlap with existing species
            if varname_new in ds_new.data_vars and overwrite:
                ds_new.drop(varname_new)
            else:
                assert(varname_new not in ds_new.data_vars), \
                    "{} already in dataset. To overwrite pass overwrite=True.".\
                    format(varname_new)

            # Verbose prints
            if verbose:
                print("Creating {}".format(varname_new))

            # Initialize new dataarray for the lumped species
            darr = dummy_darr
            darr.name = varname_new
            darr.values = np.full(darr.shape, 0.0)

            # Loop over and sum constituent species values
            num_spc = 0
            for _, spc in enumerate(lspc_dict[lspc]):
                varname = prefix + spc
                if varname not in ds_new.data_vars:
                    if verbose:
                        print("Warning: {} needed for {} not in dataset.".
                              format(spc, lspc))
                    continue
                if verbose:
                    print(" -> adding {} with scale {}".
                          format(spc, lspc_dict[lspc][spc]))
                darr.values += ds_new[varname].values * lspc_dict[lspc][spc]
                num_spc += 1

            # Replace values with NaN is no species found in dataset
            if num_spc == 0:
                if verbose:
                    print("No constituent species found! Setting to NaN.")
                darr.values = np.full(darr.shape, np.nan)

            # Merge new variable into dataset
            ds_new = xr.merge([ds_new, darr])

    return ds_new


def filter_names(
        names,
        text=""):
    """
    Returns elements in a list that match a given substring.
    Can be used in conjnction with compare_varnames to return a subset
    of variable names pertaining to a given diagnostic type or species.

    Args:
        names: list of str
            Input list of names.
        text: str
            Target text string for restricting the search.

    Returns:
        filtered_names: list of str
            Returns all elements of names that contains the substring
            specified by the "text" argument.  If "text" is omitted,
            then the original contents of names will be returned.
    """

    if text != "":
        filtered_names = [k for k in names if text in k]
    else:
        filtered_names = [k for k in names if k]

    return filtered_names


def divide_dataset_by_dataarray(
        ds,
        dr,
        varlist=None
):
    """
    Divides variables in an xarray Dataset object by a single DataArray
    object.  Will also make sure that the Dataset variable attributes
    are preserved.
    This method can be useful for certain types of model diagnostics
    that have to be divided by a counter array.  For example, local
    noontime J-value variables in a Dataset can be divided by the
    fraction of time it was local noon in each grid box, etc.

    Args:
        ds: xarray Dataset
            The Dataset object containing variables to be divided.
        dr: xarray DataArray
            The DataArray object that will be used to divide the
            variables of ds.

    Keyword Args (optional):
        varlist: list of str
            If passed, then only those variables of ds that are listed
            in varlist will be divided by dr.  Otherwise, all variables
            of ds will be divided by dr.
            Default value: None
    Returns:
        ds_new: xarray Dataset
            A new xarray Dataset object with its variables divided by dr.
    """

    # -----------------------------
    # Check arguments
    # -----------------------------
    if not isinstance(ds, xr.Dataset):
        raise TypeError("The ds argument must be of type xarray.Dataset!")

    if not isinstance(dr, xr.DataArray):
        raise TypeError("The dr argument must be of type xarray.DataArray!")

    if varlist is None:
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


def get_shape_of_data(
        data,
        vertical_dim="lev",
        return_dims=False
):
    """
    Convenience routine to return a the shape (and dimensions, if
    requested) of an xarray Dataset, or xarray DataArray.  Can also
    also take as input a dictionary of sizes (i.e. {'time': 1,
    'lev': 72, ...} from an xarray Dataset or xarray Datarray object.

    Args:
        data: xarray Dataset, xarray DataArray, or dict
            The data for which the size is requested.

    Keyword Args (optional):
        vertical_dim: str
            Specify the vertical dimension that you wish to
            return: lev or ilev.
            Default value: 'lev'
        return_dims: bool
            Set this switch to True if you also wish to return a list of
            dimensions in the same order as the tuple of dimension sizes.
            Default value: False

    Returns:
        shape: tuple of int
            Tuple containing the sizes of each dimension of dr in order:
            (time, lev|ilev, nf, lat|YDim, lon|XDim).
        dims: list of str
            If return_dims is True, then dims will contain a list of
            dimension names in the same order as shape
            (['time', 'lev', 'lat', 'lon'] for GEOS-Chem "Classic",
             or ['time', 'lev', 'nf', 'Ydim', 'Xdim'] for GCHP.
    """

    # Validate the data argument
    if isinstance(data, (xr.Dataset, xr.DataArray)):
        sizelist = data.sizes
    elif isinstance(data, dict):
        sizelist = data
    else:
        msg = (
            'The "dataset" argument must be either an xarray Dataset, '
            + " xarray DataArray, or a dictionary!"
        )
        raise ValueError(msg)

    # Initialize
    dimlist = ["time", vertical_dim, "lat", "nf", "Ydim", "lon", "Xdim"]
    shape = ()
    dims = []

    # Return a tuple with the shape of each dimension (and also a
    # list of each dimension if return_dims is True).
    for d in dimlist:
        if d in sizelist:
            shape += (sizelist[d],)
            dims.append(d)

    if return_dims:
        return shape, dims
    return shape


def get_area_from_dataset(
        ds
):
    """
    Convenience routine to return the area variable (which is
    usually called "AREA" for GEOS-Chem "Classic" or "Met_AREAM2"
    for GCHP) from an xarray Dataset object.

    Args:
        ds: xarray Dataset
            The input dataset.
    Returns:
        area_m2: xarray DataArray
            The surface area in m2, as found in ds.
    """

    if "Met_AREAM2" in ds.data_vars.keys():
        return ds["Met_AREAM2"]
    if "AREA" in ds.data_vars.keys():
        return ds["AREA"]
    msg = (
        'An area variable ("AREA" or "Met_AREAM2" is missing'
        + " from this dataset!"
    )
    raise ValueError(msg)


def get_variables_from_dataset(
        ds,
        varlist
):
    """
    Convenience routine to return multiple selected DataArray
    variables from an xarray Dataset.  All variables must be
    found in the Dataset, or else an error will be raised.

    Args:
        ds: xarray Dataset
            The input dataset.
        varlist: list of str
            List of DataArray variables to extract from ds.

    Returns:
        ds_subset: xarray Dataset
            A new data set containing only the variables
            that were requested.

    Remarks:
    Use this routine if you absolutely need all of the requested
    variables to be returned.  Otherwise
    """

    ds_subset = xr.Dataset()
    for v in varlist:
        if v in ds.data_vars.keys():
            ds_subset = xr.merge([ds_subset, ds[v]])
        else:
            msg = "{} was not found in this dataset!".format(v)
            raise ValueError(msg)

    return ds_subset


def create_dataarray_of_nan(
        name,
        sizes,
        coords,
        attrs,
        vertical_dim="lev"
):
    """
    Given an xarray DataArray dr, returns a DataArray object with
    the same dimensions, coordinates, attributes, and name, but
    with its data set to missing values (NaN) everywhere.
    This is useful if you need to plot or compare two DataArray
    variables, and need to represent one as missing or undefined.

    Args:
    name: str
        The name for the DataArray object that will contain NaNs.
    sizes: dict of int
        Dictionary of the dimension names and their sizes (e.g.
        {'time': 1 ', 'lev': 72, ...} that will be used to create
        the DataArray of NaNs.  This can be obtained from an
        xarray Dataset as ds.sizes.
    coords: dict of lists of float
        Dictionary containing the coordinate variables that will
        be used to create the DataArray of NaNs.  This can be obtained
        from an xarray Dataset with ds.coords.
    attrs: dict of str
        Dictionary containing the DataArray variable attributes
        (such as "units", "long_name", etc.).  This can be obtained
        from an xarray Dataset with dr.attrs.

    Returns:
    dr: xarray DataArray
        The output DataArray object, which will contain NaN values
        everywhere.  This will denote missing data.
    """

    # Save dims and coords into local variables
    # NOTE: Cast to type dict so that we can delete keys and values
    new_sizes = dict(sizes)
    new_coords = dict(coords)

    # Only keep one of the vertical dimensions (lev or ilev)
    if vertical_dim == "lev":
        if "ilev" in new_sizes:
            del new_sizes["ilev"]
            del new_coords["ilev"]
    elif vertical_dim == "ilev":
        if "lev" in new_sizes:
            del new_sizes["lev"]
            del new_coords["lev"]
    else:
        msg = 'The "vertical_lev" argument must be either "lev" or "ilev"!'
        raise ValueError(msg)

    # Get the names and sizes of the dimensions
    # after discarding one of "lev" or "ilev"
    [new_shape, new_dims] = get_shape_of_data(new_sizes, return_dims=True)

    # Create an array full of NaNs of the required size
    nan_arr = np.empty(new_shape, np.float)
    nan_arr.fill(np.nan)

    # Create a DataArray of NaN's
    return xr.DataArray(
        nan_arr, name=name, dims=new_dims, coords=new_coords, attrs=attrs
    )


def check_for_area(
        ds,
        gcc_area_name="AREA",
        gchp_area_name="Met_AREAM2"
):
    """
    Makes sure that a dataset has a surface area variable contained
    within it.
    GEOS-Chem Classic files all contain surface area as variable AREA.
    GCHP files do not and area must be retrieved from the met-field
    collection from variable Met_AREAM2. To simplify comparisons,
    the GCHP area name will be appended to the dataset under the
    GEOS-Chem "Classic" area name if it is present.

    Args:
        ds: xarray Dataset
            The Dataset object that will be checked.

    Keyword Args (optional):
        gcc_area_name: str
            Specifies the name of the GEOS-Chem "Classic" surface
            area varaible
            Default value: "AREA"
        gchp_area_name: str
            Specifies the name of the GCHP surface area variable.
            Default value: "Met_AREAM2"

    Returns:
        ds: xarray Dataset
            The modified Dataset object
    """

    found_gcc = gcc_area_name in ds.data_vars.keys()
    found_gchp = gchp_area_name in ds.data_vars.keys()

    if (not found_gcc) and (not found_gchp):
        msg = "Could not find {} or {} in the dataset!".format(
            gcc_area_name, gchp_area_name
        )
        raise ValueError(msg)

    if found_gchp:
        ds[gcc_area_name] = ds[gchp_area_name]

    return ds


def get_filepath(
        datadir,
        col,
        date,
        is_gchp=False,
        gchp_format_is_legacy=False
):
    """
    Routine to return file path for a given GEOS-Chem "Classic"
    (aka "GCC") or GCHP diagnostic collection and date.

    Args:
        datadir: str
            Path name of the directory containing GCC or GCHP data files.
        col: str
            Name of collection (e.g. Emissions, SpeciesConc, etc.)
            for which file path will be returned.
        date: numpy.datetime64
            Date for which file paths are requested.

    Keyword Args (optional):
        is_gchp: bool
            Set this switch to True to obtain file pathnames to
            GCHP diagnostic data files. If False, assumes GEOS-Chem "Classic"

        gchp_format_is_legacy: bool
            Set this switch to True to obtain GCHP file pathnames of
            the legacy format for diagnostics, which do not match GC-Classic
            filenames. Set to False to use same format as GC-Classic.

    Returns:
        path: str
            Pathname for the specified collection and date.
    """

    # Set filename template, extension, separator, and date string from
    # the collection, date, and data directory arguments
    separator = "_"
    extension = "z.nc4"
    date_str = np.datetime_as_string(date, unit="m")
    if is_gchp:
        if "Restart" in col:
            file_tmpl = os.path.join(datadir,
                                     "gcchem_internal_checkpoint.restart.")
            extension = ".nc4"
            date_str = np.datetime_as_string(date, unit="s")
        else:
            if gchp_format_is_legacy:
                file_tmpl = os.path.join(datadir, "GCHP.{}.".format(col))
            else:
                file_tmpl = os.path.join(datadir, "GEOSChem.{}.".format(col))
    else:
        if "Emissions" in col:
            file_tmpl = os.path.join(datadir, "HEMCO_diagnostics.")
            extension = ".nc"
            separator = ""
        else:
            file_tmpl = os.path.join(datadir, "GEOSChem.{}.".format(col))
    if isinstance(date_str, np.str_):
        date_str = str(date_str)
    date_str = date_str.replace("T", separator)
    date_str = date_str.replace("-", "")
    date_str = date_str.replace(":", "")

    # Set file path to return
    path = file_tmpl + date_str + extension
    return path


def get_filepaths(
        datadir,
        collections,
        dates,
        is_gchp=False,
        gchp_format_is_legacy=False
):
    """
    Routine to return filepaths for a given GEOS-Chem "Classic"
    (aka "GCC") or GCHP diagnostic collection.

    Args:
        datadir: str
            Path name of the directory containing GCC or GCHP data files.
        collections: list of str
            Names of collections (e.g. Emissions, SpeciesConc, etc.)
            for which file paths will be returned.
        dates: array of numpy.datetime64
            Array of dates for which file paths are requested.

    Keyword Args (optional):
        is_gchp: bool
            Set this switch to True to obtain file pathnames to
            GCHP diagnostic data files. If False, assumes GEOS-Chem "Classic"

        gchp_format_is_legacy: bool
            Set this switch to True to obtain GCHP file pathnames of
            the legacy format for diagnostics, which do not match GC-Classic
            filenames. Set to False to use same format as GC-Classic.

    Returns:
        paths: 2D list of str
            A list of pathnames for each specified collection and date.
            First dimension is collection, and second is date.
    """

    # ==================================================================
    # Initialization
    # ==================================================================

    # If collections is passed as a scalar
    # make it a list so that we can iterate
    if not isinstance(collections, list):
        collections = [collections]

    # Create the return variable
    rows, cols = (len(collections), len(dates))
    paths = [[''] * cols] * rows

    # ==================================================================
    # Create the file list
    # ==================================================================
    for c, collection in enumerate(collections):

        separator = "_"
        extension = "z.nc4"
        if is_gchp:
            # ---------------------------------------
            # Get the file path template for GCHP
            # ---------------------------------------
            if "Restart" in collection:
                file_tmpl = os.path.join(datadir,
                                         "gcchem_internal_checkpoint.restart.")
                extension = ".nc4"
            else:
                if gchp_format_is_legacy:
                    file_tmpl = os.path.join(datadir,
                                             "GCHP.{}.".format(collection))
                else:
                    file_tmpl = os.path.join(datadir, "GEOSChem.{}.".format(collection))
        else:
            # ---------------------------------------
            # Get the file path template for GCC
            # ---------------------------------------
            if "Emissions" in collection:
                file_tmpl = os.path.join(datadir,
                                         "HEMCO_diagnostics.")
                separator = ""
                extension = ".nc"
            else:
                file_tmpl = os.path.join(datadir,
                                         "GEOSChem.{}.".format(collection))

        # --------------------------------------------
        # Create a list of files for each date/time
        # --------------------------------------------
        for d, date in enumerate(dates):
            if is_gchp and "Restart" in collection:
                date_time = str(np.datetime_as_string(date, unit="s"))
            else:
                date_time = str(np.datetime_as_string(date, unit="m"))
            date_time = date_time.replace("T", separator)
            date_time = date_time.replace("-", "")
            date_time = date_time.replace(":", "")
            paths[c][d] = file_tmpl + date_time + extension

    return paths


def extract_pathnames_from_log(
        filename,
        prefix_filter=""
):
    """
    Returns a list of pathnames from a GEOS-Chem log file.
    This can be used to get a list of files that should be
    downloaded from gcgrid or from Amazon S3.

    Args:
        filename: str
            GEOS-Chem standard log file
        prefix_filter (optional): str
            Restricts the output to file paths starting with
            this prefix (e.g. "/home/ubuntu/ExtData/HEMCO/")
            Default value: ''
    Returns:
        data list: list of str
            List of full pathnames of data files found in
            the log file.
    Author:
        Jiawei Zhuang (jiaweizhuang@g.harvard.edu)
    """

    # Initialization
    prefix_len = len(prefix_filter)
    data_list = set()  # only keep unique files

    # Open file (or die with error)
    try:
        f = open(filename, "r")
    except FileNotFoundError:
        raise FileNotFoundError("Could not find file {}".format(filename))

    # Read data from the file line by line.
    # Add file paths to the data_list set.
    line = f.readline()
    while line:
        upcaseline = line.upper()
        if (": OPENING" in upcaseline) or (": READING" in upcaseline):
            data_path = line.split()[-1]
            # remove common prefix
            if data_path.startswith(prefix_filter):
                trimmed_path = data_path[prefix_len:]
                data_list.add(trimmed_path)

        # Read next line
        line = f.readline()

    # Close file and return
    f.close()
    data_list = sorted(list(data_list))
    return data_list


def get_gcc_filepath(
        outputdir,
        collection,
        day,
        time
):
    '''
    Routine for getting filepath of GEOS-Chem Classic output

    Args:
        outputdir: str
             Path of the OutputDir directory
        collection: str
             Name of output collection, e.g. Emissions or SpeciesConc
        day: str
             Number day of output, e.g. 31
        time: str
             Z time of output, e.g. 1200z

    Returns:
        filepath: str
             Path of requested file
    '''
    if collection == "Emissions":
        filepath = os.path.join(
            outputdir, "HEMCO_diagnostics.{}{}.nc".format(day, time)
        )
    else:
        filepath = os.path.join(
            outputdir, "GEOSChem.{}.{}_{}z.nc4".format(collection, day, time)
        )
    return filepath


def get_gchp_filepath(
        outputdir,
        collection,
        day,
        time
):
    '''
    Routine for getting filepath of GCHP output

    Args:
        outputdir: str
             Path of the OutputDir directory
        collection: str
             Name of output collection, e.g. Emissions or SpeciesConc
        day: str
             Number day of output, e.g. 31
        time: str
             Z time of output, e.g. 1200z

    Returns:
        filepath: str
             Path of requested file
    '''

    filepath = os.path.join(
        outputdir, "GCHP.{}.{}_{}z.nc4".format(collection, day, time)
    )
    return filepath


def get_nan_mask(
        data
):
    """
    Create a mask with NaN values removed from an input array

    Args:
        data: numpy array
            Input array possibly containing NaNs

    Returns:
        new_data: numpy array
            Original array with NaN values removed
    """

    # remove NaNs
    fill = np.nanmax(data) + 100000
    new_data = np.where(np.isnan(data), fill, data)
    new_data = np.ma.masked_where(data == fill, data)
    return new_data


def all_zero_or_nan(
        ds
):
    """
    Return whether ds is all zeros, or all nans

    Args:
        ds: numpy array
            Input GEOS-Chem data
    Returns:
        all_zero, all_nan: bool, bool
            All_zero is whether ds is all zeros, all_nan is whether ds i s all NaNs
    """

    return not np.any(ds), np.isnan(ds).all()


def dataset_mean(
        ds,
        dim="time",
        skipna=True
):
    """
    Convenience wrapper for taking the mean of an xarray Dataset.

    Args:
       ds : xarray Dataset
           Input data

    Keyword Args:
       dim : str
           Dimension over which the mean will be taken.
           Default: "time"
       skipna : bool
           Flag to omit missing values from the mean.
           Default: True

    Returns:
       ds_mean : xarray Dataset or None
           Dataset containing mean values
           Will return None if ds is not defined
    """
    if ds is None:
        return ds

    if not isinstance(ds, xr.Dataset):
        raise ValueError("Argument ds must be None or xarray.Dataset!")

    with xr.set_options(keep_attrs=True):
        return ds.mean(dim=dim, skipna=skipna)


def dataset_reader(
        multi_files
):
    """
    Returns a function to read an xarray Dataset.

    Args:
        multi_files : bool
            Denotes whether we will be reading multiple files into
            an xarray Dataset.
            Default value: False

    Returns:
         reader : either xr.open_mfdataset or xr.open_dataset
    """
    if multi_files:
        reader = xr.open_mfdataset
    else:
        reader = xr.open_dataset

    return reader

def read_config_file(config_file):
    """
    Reads configuration information from a YAML file.
    """
    # Read the configuration file in YAML format
    try:
        print(f"Using configuration file {config_file}")
        config = yaml.safe_load(open(config_file))
    except Exception as err:
        msg = f"Error reading configuration in {config_file}: {err}"
        raise Exception(msg)

    return config

