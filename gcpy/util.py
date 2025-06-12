"""
Internal utilities for helping to manage xarray and numpy
objects used throughout GCPy
"""
import os
from shutil import copyfile
import warnings
from textwrap import wrap
from yaml import safe_load
import numpy as np
from pandas import Series
import xarray as xr
from pypdf import PdfWriter, PdfReader
from gcpy.constants import ENCODING, TABLE_WIDTH
from gcpy.cstools import is_cubed_sphere_rst_grid

# ======================================================================
# %%%%% METHODS %%%%%
# ======================================================================

def convert_lon(
        data,
        dim='lon',
        fmt='atlantic',
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
    verify_variable_type(data, (xr.DataArray, xr.Dataset))

    data_copy = data.copy()

    lon = data_copy[dim].values
    new_lon = np.empty_like(lon)

    # Tweak offset for rolling the longitudes later
    offset = 0 if neg_dateline else 1

    if fmt not in ['atlantic', 'pacific']:
        msg = f"Cannot convert longitudes for format '{fmt}'; "
        msg += "please choose one of 'atlantic' or 'pacific'"
        raise ValueError(msg)

    # Create a mask to decide how to mutate the longitude values
    if fmt == 'atlantic':
        mask = lon >= 180 if neg_dateline else lon > 180

        new_lon[mask] = -(360. - lon[mask])
        new_lon[~mask] = lon[~mask]

        roll_len = len(data[dim]) // 2 - offset

    elif fmt == 'pacific':
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

    # For restart files, just split at the first underscore and return
    # the text followiong the underscore.  This will preserve certain
    # species names, such as the TransportTracers species CO_25, etc.
    if "SpeciesRst" in display_name:
        return display_name.split("_", 1)[1]

    # Special handling for Inventory totals
    if "INV" in display_name.upper():
        display_name = display_name.replace("_", " ")

    # Replace text
    for var in ["Emis", "EMIS", "emis", "Inv", "INV", "inv"]:
        display_name = display_name.replace(var, "")

    # Replace only the first underscore with a space
    display_name = display_name.replace("_", " ", 1)

    return display_name


def format_number_for_table(
        number,
        max_thresh=1.0e8,
        min_thresh=1.0e-6,
        f_fmt="18.6f",
        e_fmt="18.8e"
):
    """
    Returns a format string for use in the "print_totals" routine.
    If the number is greater than a maximum threshold or smaller
    than a minimum threshold, then use scientific notation format.
    Otherwise use floating-piont format.

    Special case: do not convert 0.0 to exponential notation.

    Args:
    -----
    number : float
        Number to be printed

    max_thresh, min_thresh: float
        If |number| > max_thresh, use scientific notation.
        If |number| < min_thresh, use scientific notation

    f_fmt, e_fmt : str
        The default floating point string and default scientific
        notation string.
        Default values: 18.6f, 18.6e

    Returns:
    --------
    fmt_str : str
        Formatted string that can be inserted into the print
        statement in print_totals.
    """
    abs_number = np.abs(number)

    if not abs_number > 1e-60:
        return f"{number:{f_fmt}}"

    if abs_number > max_thresh or abs_number < min_thresh:
        return f"{number:{e_fmt}}"
    return f"{number:{f_fmt}}"


def print_totals(
        ref,
        dev,
        ofile,
        diff_list,
        masks=None,
):
    """
    Computes and prints Ref and Dev totals (as well as the difference
    Dev - Ref) for two xarray DataArray objects.

    Args:
        ref: xarray DataArray
            The first DataArray to be compared (aka "Reference")
        dev: xarray DataArray
            The second DataArray to be compared (aka "Development")
        ofile: file
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
    verify_variable_type(ref, xr.DataArray)
    verify_variable_type(dev, xr.DataArray)
    verify_variable_type(diff_list, list)

    # Determine if either Ref or Dev have all NaN values:
    ref_is_all_nan = np.isnan(ref.values).all()
    dev_is_all_nan = np.isnan(dev.values).all()

    # If Ref and Dev do not contain all NaNs, then make sure
    # that Ref and Dev have the same units before proceeding.
    if not ref_is_all_nan and not dev_is_all_nan:
        if ref.units != dev.units:
            msg = f"Ref has units {ref.units}, but Dev has units {dev.units}!"
            raise ValueError(msg)

    # ==================================================================
    # Get the diagnostic name and units
    # ==================================================================
    diagnostic_name = dev.name
    if dev_is_all_nan:
        diagnostic_name = ref.name

    # Create the display name for the table
    display_name = create_display_name(diagnostic_name)

    # Get the species name from the display name
    species_name = display_name
    cidx = species_name.find(" ")
    if cidx > 0:
        species_name = display_name[0:cidx]

    # Special handling for totals
    if "_TOTAL" in diagnostic_name.upper():
        print("-" * TABLE_WIDTH, file=ofile)

    # ==================================================================
    # Sum the Ref array (or set to NaN if missing)
    # ==================================================================
    refarr = ref.values
    if ref_is_all_nan:
        total_ref = np.nan
    else:
        if masks is not None:
            refarr = np.ma.masked_array(refarr, masks["Ref_TropMask"])
        total_ref = np.sum(refarr, dtype=np.float64)

    # ==================================================================
    # Sum the Dev array (or set to NaN if missing)
    # ==================================================================
    devarr = dev.values
    if dev_is_all_nan:
        total_dev = np.nan
    else:
        if masks is not None:
            devarr = np.ma.masked_array(devarr, masks["Dev_TropMask"])
        total_dev = np.sum(devarr, dtype=np.float64)

    # ==================================================================
    # Compute differences (or set to NaN if missing)
    # ==================================================================
    if ref_is_all_nan or dev_is_all_nan:
        diff = np.nan
    else:
        diff = total_dev - total_ref
    has_diffs = abs(diff) > np.float64(0.0)

    # Append to the list of differences.  If no differences then append
    # None.  Duplicates can be stripped out in the calling routine.
    if has_diffs:
        diff_str = " * "
        diff_list.append(species_name)
    else:
        diff_str = ""
        diff_list.append(None)

    # ==================================================================
    # Compute % differences (or set to NaN if missing)
    # If ref is very small, near zero, also set the % diff to NaN
    # ==================================================================
    if np.isnan(total_ref) or np.isnan(total_dev):
        pctdiff = np.nan
    else:
        pctdiff = ((total_dev - total_ref) / total_ref) * 100.0
        if np.abs(total_ref) < 1.0e-15:
            pctdiff = np.nan

    # ==================================================================
    # Write output to file and return
    # ==================================================================
    ref_fmt = format_number_for_table(total_ref)
    dev_fmt = format_number_for_table(total_dev)
    diff_fmt = format_number_for_table(diff)
    pctdiff_fmt = format_number_for_table(pctdiff)

    print(f"{display_name[0:19].ljust(19)}: {ref_fmt}  {dev_fmt}  {diff_fmt}  {pctdiff_fmt}  {diff_str}", file=ofile)

    return diff_list


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
            a variable name "SpeciesConcVV_NO", and you specify
            remove_prefix="SpeciesConcVV_", then the bookmark for
            that variable will be just "NO", etc.
         verbose: bool
            Set this flag to True to print extra informational output.
            Default value: False
    """

    # Setup
    with open(pdfname, "rb") as pdfobj:
        input_pdf = PdfReader(pdfobj) #, overwriteWarnings=False)
        output_pdf = PdfWriter()

        for i, varname in enumerate(varlist):
            bookmarkname = varname.replace(remove_prefix, "")
            if verbose:
                print(f"Adding bookmark for {varname} with name {bookmarkname}")
            output_pdf.add_page(input_pdf.pages[i])
            output_pdf.add_outline_item(bookmarkname, i)
            output_pdf.page_mode = "/UseOutlines"

        # Write to temp file
        pdfname_tmp = pdfname + "_with_bookmarks.pdf"

        with open(pdfname_tmp, "wb") as output_stream:
            output_pdf.write(output_stream)
            output_stream.close()

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
    with open(pdfname, "rb") as pdfobj:
        input_pdf = PdfReader(pdfobj)
        output_pdf = PdfWriter()
        warninglist = [k.replace(remove_prefix, "") for k in warninglist]

        # ===============================================================
        # Loop over the subcats in this category; make parent bookmark
        # ===============================================================
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
            output_pdf.add_page(input_pdf.pages[i])
            parent = output_pdf.add_outline_item(subcat, i)
            output_pdf.page_mode = "/UseOutlines"
            first = True

            # Loop over variables in this subcategory; make children bookmarks
            for varname in catdict[category][subcat]:
                if varname in warninglist:
                    print(f"Warning: skipping {varname}")
                    continue
                if first:
                    output_pdf.add_outline_item(varname, i, parent)
                    first = False
                else:
                    i = i + 1
                    output_pdf.add_page(input_pdf.pages[i])
                    output_pdf.add_outline_item(varname, i, parent)
                    output_pdf.page_mode = "/UseOutlines"

        # ==============================================================
        # Write to temp file
        # ==============================================================
        pdfname_tmp = pdfname + "_with_bookmarks.pdf"
        with open(pdfname_tmp, "wb") as output_stream:
            output_pdf.write(output_stream)
            output_stream.close()

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
    verify_variable_type(refdata, xr.Dataset)
    verify_variable_type(devdata, xr.Dataset)

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
        for var in refonly:
            if verbose:
                print(f"Creating array of NaN in devdata for: {var}")
            darr = create_blank_dataarray(
                name=refdata[var].name,
                sizes=devdata.sizes,
                coords=devdata.coords,
                attrs=refdata[var].attrs,
                **kwargs
            )
            devlist.append(darr)
        devdata = xr.merge(devlist)

        # ==============================================================
        # For each variable that is in devdata but not in refdata,
        # add a new DataArray to refdata with the same sizes but
        # containing all NaN's.  This will allow us to represent those
        # variables as missing values # when we plot against devdata.
        # ==================================================================
        reflist = [refdata]
        for var in devonly:
            if verbose:
                print(f"Creating array of NaN in refdata for: {var}")
            darr = create_blank_dataarray(
                name=devdata[var].name,
                sizes=refdata.sizes,
                coords=refdata.coords,
                attrs=devdata[var].attrs,
                **kwargs
            )
            reflist.append(darr)
        refdata = xr.merge(reflist)

    return refdata, devdata


def reshape_MAPL_CS(darr):
    """
    Reshapes data if contains dimensions indicate MAPL v1.0.0+ output
    (i.e. reshapes from "diagnostic" to "checkpoint" dimension format.)

    Args:
    -----
    darr: xarray DataArray
        The input data array.

    Returns:
    --------
    darr: xarray DataArray
        The modified data array (w/ dimensions renamed & transposed).

    Remarks:
    --------
    Currently only used for GCPy plotting code.
    """
    # Suppress annoying future warnings for now
    warnings.filterwarnings("ignore", category=FutureWarning)

    # Only do the following for DataArray objects
    # (otherwise just fall through and return the original argument as-is)
    if isinstance(darr, xr.DataArray):
        with xr.set_options(keep_attrs=True):
            if "nf" in darr.dims and \
               "Xdim" in darr.dims and "Ydim" in darr.dims:
                darr = darr.stack(lat=("nf", "Ydim"))
                darr = darr.rename({"Xdim": "lon"})
            if "lev" in darr.dims and "time" in darr.dims:
                darr = darr.transpose("time", "lev", "lat", "lon")
            elif "lev" in darr.dims:
                darr = darr.transpose("lev", "lat", "lon")
            elif "time" in darr.dims:
                darr = darr.transpose("time", "lat", "lon")
            else:
                darr = darr.transpose("lat", "lon")
    return darr


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
        # if the coords do not align then set time dimensions equal
        try:
            xr.align(dev, ref, join='exact')
        except BaseException:
            ref.coords["time"] = dev.coords["time"]
        with xr.set_options(keep_attrs=True):
            absdiffs = dev - ref
            fracdiffs = dev / ref
            for var in dev.data_vars.keys():
                # Ensure the diffs Dataset includes attributes
                absdiffs[var].attrs = dev[var].attrs
                fracdiffs[var].attrs = dev[var].attrs
    elif 'nf' in ref.dims and 'nf' in dev.dims:

        # Include special handling if cubed sphere grid dimension names are different
        # since they changed in MAPL v1.0.0.
        if "lat" in ref.dims and "Xdim" in dev.dims:
            ref_newdimnames = dev.copy()
            for var in dev.data_vars.keys():
                if "Xdim" in dev[var].dims:
                    ref_newdimnames[var].values = ref[var].values.reshape(
                        dev[var].values.shape)
                # NOTE: the reverse conversion is gchp_dev[v].stack(lat=("nf","Ydim")).transpose(
                # "time","lev","lat","Xdim").values

        with xr.set_options(keep_attrs=True):
            absdiffs = dev.copy()
            fracdiffs = dev.copy()
            for var in dev.data_vars.keys():
                if "Xdim" in dev[var].dims or "lat" in dev[var].dims:
                    absdiffs[var].values = dev[var].values - ref[var].values
                    fracdiffs[var].values = dev[var].values / ref[var].values
                    # NOTE: The diffs Datasets are created without variable
                    # attributes; we have to reattach them
                    absdiffs[var].attrs = dev[var].attrs
                    fracdiffs[var].attrs = dev[var].attrs
    else:
        print('Diff-of-diffs plot supports only identical grid types (lat/lon or cubed-sphere)' + \
              ' within each dataset pair')
        raise ValueError

    return absdiffs, fracdiffs


def slice_by_lev_and_time(
        dset,
        varname,
        itime,
        ilev,
        flip
):
    """
    Given a Dataset, returns a DataArray sliced by desired time and level.

    Args:
        dset: xarray Dataset
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
        darr: xarray DataArray
            DataArray of data variable sliced according to ilev and itime
    """
    # used in compare_single_level and compare_zonal_mean to get dataset slices
    verify_variable_type(dset, xr.Dataset)
    if not varname in dset.data_vars.keys():
        msg="Could not find 'varname' in ds!"
        raise ValueError(msg)

    # NOTE: isel no longer seems to work on a Dataset, so
    # first createthe DataArray object, then use isel on it.
    #  -- Bob Yantosca (19 Jan 2023)
    darr = dset[varname]
    vdims = darr.dims
    if ("time" in vdims and darr.time.size > 0) and "lev" in vdims:
        if flip:
            fliplev=len(darr['lev']) - 1 - ilev
            return darr.isel(time=itime, lev=fliplev)
        return darr.isel(time=itime, lev=ilev)
    if ("time" not in vdims or itime == -1) and "lev" in vdims:
        if flip:
            fliplev= len(darr['lev']) - 1 - ilev
            return darr.isel(lev=fliplev)
        return darr.isel(lev=ilev)
    if ("time" in vdims and darr.time.size > 0 and itime != -1) and \
       "lev" not in vdims:
        return darr.isel(time=itime)
    return darr


def rename_and_flip_gchp_rst_vars(
        dset
):
    '''
    Transforms a GCHP restart dataset to match GCClassic names
    and level conventions.

    Args:
        dset: xarray Dataset
            The input dataset.

    Returns:
        dset: xarray Dataset
            If the input dataset is from a GCHP restart file, then
            dset will contain the original data with variables renamed
            to match the GEOS-Chem Classic naming conventions, and
            with levels indexed as lev:positive="up".  Otherwise, the
            original data will be returned.
    '''
    verify_variable_type(dset, xr.Dataset)

    # Return if this dataset is not from a GCHP checkpoint/restart file
    if not is_cubed_sphere_rst_grid(dset):
        return dset

    # Create dictionary of variable name replacements
    old_to_new = {}
    for var in dset.data_vars.keys():
        # TODO: Think of better algorithm in case we ever change
        # the internal state to start with something else than "SPC_".
        if var.startswith("SPC_"):
            spc = var.replace('SPC_', '')
            old_to_new[var] = 'SpeciesRst_' + spc
        if var == "DELP_DRY":
            old_to_new["DELP_DRY"] = "Met_DELPDRY"
        if var == "DELPDRY":
            old_to_new["DELPDRY"] = "Met_DELPDRY"
        if var == "BXHEIGHT":
            old_to_new["BXHEIGHT"] = "Met_BXHEIGHT"
        if var == "TropLev":
            old_to_new["TropLev"] = "Met_TropLev"

    # Replace variable names in one operation
    dset = dset.rename(old_to_new)

    # Flip levels
    dset = dset.sortby('lev', ascending=False)
    dset.lev.attrs["positive"] = "up"

    return dset


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
    verify_variable_type(dict0, dict)
    verify_variable_type(dict1, dict)

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
            commonvarsData   List of all commmon 2D or 3D data variables,
                             excluding index variables.  This is the
                             list of "plottable" variables.
            refonly          List of 2D or 3D variables that are only
                             present in refdata.
            devonly          List of 2D or 3D variables that are only
                             present in devdata
    """
    verify_variable_type(refdata, xr.Dataset)
    verify_variable_type(devdata, xr.Dataset)

    refvars = list(refdata.data_vars.keys())
    devvars = list(devdata.data_vars.keys())
    commonvars = sorted(list(set(refvars).intersection(set(devvars))))
    refonly = [var for var in refvars if var not in devvars]
    devonly = [var for var in devvars if var not in refvars]
    dimmismatch = [v for v in commonvars if refdata[v].ndim != devdata[v].ndim]
    # Assume plottable data has (lon, lat) or (Xdim, Ydim)
    # This is OK for purposes of benchmarking
    #  -- Bob Yantosca (09 Feb 2023)
    commonvars_data = [
        var for var in commonvars if (
            ("lat" in refdata[var].dims or "Ydim" in refdata[var].dims) and
            ("lon" in refdata[var].dims or "Xdim" in refdata[var].dims) and
            ("lat" in devdata[var].dims or "Ydim" in devdata[var].dims) and
            ("lon" in devdata[var].dims or "Xdim" in devdata[var].dims)
        )
    ]
    commonvars_other = [
        var for var in commonvars if (
           var not in commonvars_data
        )
    ]
    commonvars_2d = [
        var for var in commonvars if (
            (var in commonvars_data) and ("lev" not in refdata[var].dims)
                                     and ("lev" not in devdata[var].dims)
        )
    ]
    commonvars_3d = [
        var for var in commonvars if (
            (var in commonvars_data) and ("lev" in refdata[var].dims)
                                     and ("lev" in devdata[var].dims)
        )
    ]

    # Print information on common and mismatching variables,
    # as well as dimensions
    if not quiet:
        print("\nComparing variable names in compare_varnames")
        print(f"{len(commonvars)} common variables")
        if len(refonly) > 0:
            print(f"{len(refonly)} variables in ref only (skip)")
            print(f"   Variable names: {refonly}")
        else:
            print("0 variables in ref only")
            if len(devonly) > 0:
                print(f"len({devonly} variables in dev only (skip)")
                print(f"   Variable names: {devonly}")
            else:
                print("0 variables in dev only")
                if len(dimmismatch) > 0:
                    print(f"{dimmismatch} common variables have different dimensions")
                    print(f"   Variable names: {dimmismatch}")
                else:
                    print("All variables have same dimensions in ref and dev")

    # For safety's sake, remove the 0-D and 1-D variables from
    # commonvarsData, refonly, and devonly.  This will ensure that
    # these lists will only contain variables that can be plotted.
    commonvars_data = [var for var in commonvars if var not in commonvars_other]
    refonly = [var for var in refonly if var not in commonvars_other]
    devonly = [var for var in devonly if var not in commonvars_other]

    return {
        "commonvars": commonvars,
        "commonvars2D": commonvars_2d,
        "commonvars3D": commonvars_3d,
        "commonvarsData": commonvars_data,
        "commonvarsOther": commonvars_other,
        "refonly": refonly,
        "devonly": devonly
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
    print(f"    {refstr}:  {units}")
    print(f"    {devstr}:  {units}")
    print("Array sizes:")
    print(f"    {refstr}:  {refvar.shape}")
    print(f"    {devstr}:  {devvar.shape}")
    print("Global stats:")
    print("  Mean:")
    print(f"    {refstr}:  {np.round(refvar.values.mean(), 20)}")
    print(f"    {devstr}:  {np.round(devvar.values.mean(), 20)}")
    print("  Min:")
    print(f"    {refstr}:  {np.round(refvar.values.min(), 20)}")
    print(f"    {devstr}:  {np.round(devvar.values.min(), 20)}")
    print("  Max:")
    print(f"    {refstr}:  {np.round(refvar.values.max(), 20)}")
    print(f"    {devstr}:  {np.round(devvar.values.max(), 20)}")
    print("  Sum:")
    print(f"    {refstr}:  {np.round(refvar.values.sum(), 20)}")
    print(f"    {devstr}:  {np.round(devvar.values.sum(), 20)}")


def convert_bpch_names_to_netcdf_names(
        dset,
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
    names = read_config_file(
        os.path.join(
            os.path.dirname(__file__),
            "bpch_to_nc_names.yml"
        ),
        quiet=True
    )

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
    for variable_name in dset.data_vars.keys():

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
                        print(f"WARNING: skipping {key}")
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
                        print(f"Skipping: {oldid}")
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
                    print(f"WARNING: Nothing defined for: {variable_name}")
                continue

            # Overwrite certain variable names
            if newvar in special_vars:
                newvar = special_vars.get(newvar)

            # Update the dictionary of names with this pair
            old_to_new.update({original_variable_name: newvar})

    # Verbose output
    if verbose:
        print("\nList of bpch names and netCDF names")
        for key in old_to_new:
            print(f"{key : <25} ==> {old_to_new[key] : <40}")

    # Rename the variables in the dataset
    if verbose:
        print("\nRenaming variables in the data...")
    with xr.set_options(keep_attrs=True):
        dset = dset.rename(name_dict=old_to_new)

    # Return the dataset
    return dset


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
        return [var for var in names if text in var]
    return [var for var in names if var]


def divide_dataset_by_dataarray(
        dset,
        darr,
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
        dset: xarray Dataset
            The Dataset object containing variables to be divided.
        darr: xarray DataArray
            The DataArray object that will be used to divide the
            variables of ds.

    Keyword Args (optional):
        varlist: list of str
            If passed, then only those variables of ds that are listed
            in varlist will be divided by dr.  Otherwise, all variables
            of ds will be divided by dr.
            Default value: None
    Returns:
        dset_new: xarray Dataset
            A new xarray Dataset object with its variables divided
            by darr.
    """

    # -----------------------------
    # Check arguments
    # -----------------------------
    verify_variable_type(dset, xr.Dataset)
    verify_variable_type(darr, xr.DataArray)
    if varlist is None:
        varlist = dset.data_vars.keys()

    # -----------------------------
    # Do the division
    # -----------------------------

    # Keep all Dataset attributes
    with xr.set_options(keep_attrs=True):

        # Loop over variables
        for var in varlist:

            # Divide each variable of ds by dr
            dset[var] = dset[var] / darr

    return dset


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
    for dim in dimlist:
        if dim in sizelist:
            shape += (sizelist[dim],)
            dims.append(dim)

    if return_dims:
        return shape, dims
    return shape


def get_area_from_dataset(
        dset
):
    """
    Convenience routine to return the area variable (which is
    usually called "AREA" for GEOS-Chem "Classic" or "Met_AREAM2"
    for GCHP) from an xarray Dataset object.

    Args:
        dset: xarray Dataset
            The input dataset.
    Returns:
        area_m2: xarray DataArray
            The surface area in m2, as found in ds.
    """
    verify_variable_type(dset, xr.Dataset)

    if "Met_AREAM2" in dset.data_vars.keys():
        return dset["Met_AREAM2"]
    if "AREA" in dset.data_vars.keys():
        return dset["AREA"]
    msg = (
        'An area variable ("AREA" or "Met_AREAM2" is missing'
        + " from this dataset!"
    )
    raise ValueError(msg)


def get_variables_from_dataset(
        dset,
        varlist
):
    """
    Convenience routine to return multiple selected DataArray
    variables from an xarray Dataset.  All variables must be
    found in the Dataset, or else an error will be raised.

    Args:
        dset: xarray Dataset
            The input dataset.
        varlist: list of str
            List of DataArray variables to extract from ds.

    Returns:
        dset_subset: xarray Dataset
            A new data set containing only the variables
            that were requested.

    Remarks:
    Use this routine if you absolutely need all of the requested
    variables to be returned.  Otherwise
    """
    verify_variable_type(dset, xr.Dataset)

    dset_subset = xr.Dataset()
    for var in varlist:
        if var in dset.data_vars.keys():
            dset_subset = xr.merge([dset_subset, dset[var]])
        else:
            msg = f"{var} was not found in this dataset!"
            raise ValueError(msg)

    return dset_subset


def create_blank_dataarray(
        name,
        sizes,
        coords,
        attrs,
        fill_value=np.nan,
        fill_type=np.float64,
        vertical_dim="lev"
):
    """
    Given an xarray DataArray dr, returns a DataArray object with
    the same dimensions, coordinates, attributes, and name, but
    with its data set to missing values (default=NaN) everywhere.
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
    fill_value: np.nan or numeric type
        Value with which the DataArray object will be filled.
        Default value: np.nan
    fill_type: numeric type
        Specifies the numeric type of the DataArray object.
        Default value: np.float64 (aka "double")
    vertical_dim: str
        Specifies the name of the vertical dimension (e.g. "lev", "ilev")
        Default: "lev"

    Returns:
    dr: xarray DataArray
        The output DataArray object, which will be set to the value
        specified by the fill_value argument everywhere.
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
    fill_arr = np.empty(new_shape, dtype=fill_type)
    fill_arr.fill(fill_value)

    # Create a DataArray of NaN's
    return xr.DataArray(
        fill_arr,
        name=name,
        dims=new_dims,
        coords=new_coords,
        attrs=attrs
    )


def check_for_area(
        dset,
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
        dset: xarray Dataset
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
    verify_variable_type(dset, xr.Dataset)

    found_gcc = gcc_area_name in dset.data_vars.keys()
    found_gchp = gchp_area_name in dset.data_vars.keys()

    if not found_gcc and not found_gchp:
        msg = f"Could not find {gcc_area_name} or {gchp_area_name} "
        msg += "in the dataset!"
        raise ValueError(msg)

    if found_gchp:
        dset[gcc_area_name] = dset[gchp_area_name]

    return dset


def get_filepath(
        datadir,
        col,
        date,
        is_gchp=False,
        gchp_res="c00",
        gchp_is_pre_14_0=False
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

        gchp_res: str
            Cubed-sphere resolution of GCHP data grid.
            Only needed for restart files.
            Default value: "c00".

        gchp_is_pre_14_0: bool
            Set this switch to True to obtain GCHP file pathnames used in
            versions before 14.0. Only needed for restart files.

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
            extension = ".nc4"
            date_str = np.datetime_as_string(date, unit="s")
            if gchp_is_pre_14_0:
                file_tmpl = os.path.join(
                    datadir,
                    "gcchem_internal_checkpoint.restart."
                )
            else:
                file_tmpl = os.path.join(
                    datadir,
                    "GEOSChem.Restart."
                )
        else:
            file_tmpl = os.path.join(datadir, f"GEOSChem.{col}.")
    else:
        if "Emissions" in col:
            file_tmpl = os.path.join(datadir, "HEMCO_diagnostics.")
            extension = ".nc"
            separator = ""
        elif "Restart" in col:
            file_tmpl = os.path.join(datadir, "GEOSChem.Restart.")
        else:
            file_tmpl = os.path.join(datadir, f"GEOSChem.{col}.")
    if isinstance(date_str, np.str_):
        date_str = str(date_str)
    date_str = date_str.replace("T", separator)
    date_str = date_str.replace("-", "")
    date_str = date_str.replace(":", "")

    # Set file path. Include grid resolution if GCHP restart file.
    path = file_tmpl + date_str + extension
    if is_gchp and "Restart" in col and not gchp_is_pre_14_0:
        path = file_tmpl + date_str[:len(date_str)-2] + \
            "z." + gchp_res + extension

    return path


def get_filepaths(
        datadir,
        collections,
        dates,
        is_gchp=False,
        gchp_res="c00",
        gchp_is_pre_14_0=False
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

        gchp_res: str
            Cubed-sphere resolution of GCHP data grid.
            Only needed for restart files.
            Default value: "c00".

        gchp_is_pre_14_0: bool
            Set this switch to True to obtain GCHP file pathnames used in
            versions before 14.0. Only needed for diagnostic files.

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
    for c_idx, collection in enumerate(collections):

        separator = "_"
        extension = "z.nc4"
        if is_gchp:
            # ---------------------------------------
            # Get the file path template for GCHP
            # ---------------------------------------
            if "Restart" in collection:
                extension = ".nc4"
                if gchp_is_pre_14_0:
                    file_tmpl = os.path.join(
                        datadir,
                        "gcchem_internal_checkpoint.restart."
                    )
                else:
                    file_tmpl = os.path.join(
                        datadir,
                        "GEOSChem.Restart."
                    )
            else:
                file_tmpl = os.path.join(
                    datadir,
                    f"GEOSChem.{collection}."
                )
        else:
            # ---------------------------------------
            # Get the file path template for GCC
            # ---------------------------------------
            if "Emissions" in collection:
                file_tmpl = os.path.join(
                    datadir,
                    "HEMCO_diagnostics."
                )
                separator = ""
                extension = ".nc"
            elif "Restart" in collection:
                file_tmpl = os.path.join(
                    datadir,
                    "GEOSChem.Restart."
                )

            else:
                file_tmpl = os.path.join(
                    datadir,
                    f"GEOSChem.{collection}."
                )

        # --------------------------------------------
        # Create a list of files for each date/time
        # --------------------------------------------
        for d_idx, date in enumerate(dates):
            if is_gchp and "Restart" in collection:
                date_time = str(np.datetime_as_string(date, unit="s"))
            else:
                date_time = str(np.datetime_as_string(date, unit="m"))
            date_time = date_time.replace("T", separator)
            date_time = date_time.replace("-", "")
            date_time = date_time.replace(":", "")

            # Set file path. Include grid resolution if GCHP restart file.
            paths[c_idx][d_idx] = file_tmpl + date_time + extension
            if is_gchp and "Restart" in collection and not gchp_is_pre_14_0:
                paths[c_idx][d_idx] = file_tmpl + \
                    date_time[:len(date_time)-2] + \
                    "z." + gchp_res + extension

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

    # Open file
    with open(filename, "r", encoding=ENCODING) as ifile:

        # Read data from the file line by line.
        # Add file paths to the data_list set.
        line = ifile.readline()
        while line:
            upcaseline = line.upper()
            if (": OPENING" in upcaseline) or (": READING" in upcaseline):
                data_path = line.split()[-1]
                # remove common prefix
                if data_path.startswith(prefix_filter):
                    trimmed_path = data_path[prefix_len:]
                    data_list.add(trimmed_path)

            # Read next line
            line = ifile.readline()

        # Close file and return
        ifile.close()

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
            outputdir,
            f"HEMCO_diagnostics.{day}{time}.nc"
        )
    else:
        filepath = os.path.join(
            outputdir,
            f"GEOSChem.{collection}.{day}_{time}z.nc4"
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
        outputdir,
        f"GCHP.{collection}.{day}_{time}z.nc4"
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
        dset
):
    """
    Return whether ds is all zeros, or all nans

    Args:
        dset: numpy array
            Input GEOS-Chem data
    Returns:
        all_zero, all_nan: bool, bool
            all_zero is whether ds is all zeros,
            all_nan  is whether ds is all NaNs
    """

    return not np.any(dset), np.isnan(dset).all()


def dataset_mean(
        dset,
        dim="time",
        skipna=True
):
    """
    Convenience wrapper for taking the mean of an xarray Dataset.

    Args:
       dset : xarray Dataset
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
    verify_variable_type(dset, (xr.Dataset, type(None)))

    if dset is None:
        return dset

    with xr.set_options(keep_attrs=True):
        return dset.mean(dim=dim, skipna=skipna)


def dataset_reader(
        multi_files,
        verbose=False
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
        if verbose:
            print('Reading data via xarray open_mfdataset\n')
    else:
        reader = xr.open_dataset
        if verbose:
            print('Reading data via xarray open_dataset\n')

    return reader


def read_config_file(config_file, quiet=False):
    """
    Reads configuration information from a YAML file.
    """
    # Read the configuration file in YAML format
    try:
        if not quiet:
            print(f"Using configuration file {config_file}")
        with open(config_file, encoding=ENCODING) as stream:
            return safe_load(stream)
    except Exception as err:
        msg = f"Error reading configuration in {config_file}: {err}"
        raise Exception(msg) from err

def unique_values(
        this_list,
        drop=None,
):
    """
    Given a list, returns a sorted list of unique values.

    Args:
    -----
    this_list : list
        Input list (may contain duplicate values)

    drop: list of str
        List of variable names to exclude

    Returns:
    --------
    unique: list
        List of unique values from this_list
    """
    verify_variable_type(this_list, list)
    verify_variable_type(drop, list)

    unique = list(set(this_list))

    if drop is not None:
        for var in drop:
            if var in unique:
                unique.remove(var)

    unique.sort()

    return unique


def wrap_text(
        text,
        width=80
):
    """
    Wraps text so that it fits within a certain line width.

    Args:
    -----
    text: str or list of str
        Input text to be word-wrapped.
    width: int
        Line width, in characters.
        Default value: 80

    Returns:
    --------
    Original text reformatted so that it fits within lines
    of 'width' characters or less.
    """
    if not isinstance(text, str):
        if isinstance(text, list):
            text = ' '.join(text)  # List -> str conversion
        else:
            raise ValueError("Argument 'text' must be either str or list!")

    text = wrap(text, width=width)
    text = '\n'.join(text)

    return text


def insert_text_into_file(
        filename,
        search_text,
        replace_text,
        width=80
):
    """
    Convenience routine to insert text into a file.  The best way
    to do this is to read the contents of the file, manipulate the
    text, and then overwrite the file.

    Args:
    -----
    filename: str
        The file with text to be replaced.
    search_text: str
        Text string in the file that will be replaced.
    replace_text: str or list of str
        Text that will replace 'search_text'
    width: int
        Will "word-wrap" the text in 'replace_text' to this width
    """
    verify_variable_type(filename, str)
    verify_variable_type(search_text, str)
    verify_variable_type(replace_text, (str, list))

    # Word-wrap the replacement text
    # (does list -> str conversion if necessary)
    replace_text = wrap_text(
        replace_text,
        width=width
    )

    with open(filename, "r", encoding=ENCODING) as ifile:
        filedata = ifile.read()
        ifile.close()

    filedata = filedata.replace(
        search_text,
        replace_text
    )

    with open(filename, "w", encoding=ENCODING) as ofile:
        ofile.write(filedata)
        ofile.close()


def array_equals(
        refdata,
        devdata,
        dtype=np.float64
):
    """
    Tests two arrays for equality.  Useful for checking which
    species have nonzero differences in benchmark output.

    Args:
    -----
    refdata: xarray DataArray or numpy ndarray
        The first array to be checked.
    devdata: xarray DataArray or numpy ndarray
        The second array to be checked.
    dtype : np.float32 or np.float64
        The precision that will be used to make the evaluation.
        Default: np.float64

    Returns:
    --------
    True if both arrays are equal; False if not
    """
    if not isinstance(refdata, np.ndarray):
        if isinstance(refdata, xr.DataArray):
            refdata = refdata.values
        else:
            raise ValueError(
            "Argument 'refdata' must be an xarray DataArray or numpy ndarray!"
            )
    if not isinstance(devdata, np.ndarray):
        if isinstance(devdata, xr.DataArray):
            devdata = devdata.values
        else:
            raise ValueError(
            "Argument 'devdata' must be an xarray DataArray or numpy ndarray!"
            )

    # This method will work if the arrays hve different dimensions
    # but an element-by-element search will not!
    refsum = np.nansum(refdata, dtype=dtype)
    devsum = np.nansum(devdata, dtype=dtype)
    return (not np.abs(devsum - refsum) > dtype(0.0))


def make_directory(
        dir_name,
        overwrite
):
    """
    Creates a directory where benchmark plots/tables will be placed.

    Args:
    -----
    dir_name : str
        Name of the directory to be created.
    overwrite : bool
        Set to True if you wish to overwrite prior contents in
        the directory 'dir_name'
    """
    verify_variable_type(dir_name, str)
    verify_variable_type(overwrite, bool)

    if os.path.isdir(dir_name) and not overwrite:
        msg = f"Directory {dir_name} exists!\n"
        msg += "Pass overwrite=True to overwrite files in that directory."
        raise ValueError(msg)

    if not os.path.isdir(dir_name):
        os.makedirs(dir_name)


def trim_cloud_benchmark_label(
        label
):
    """
    Removes the first part of the cloud benchmark label string
    (e.g. "gchp-c24-1Hr", "gcc-4x5-1Mon", etc) to avoid clutter.
    """
    verify_variable_type(label, str)

    for var in [
        "gcc-4x5-1Hr",
        "gchp-c24-1Hr",
        "gcc-4x5-1Mon",
        "gchp-c24-1Mon",
    ]:
        if var in label:
            label.replace(var, "")

    return label


def verify_variable_type(
        var,
        var_type
):
    """
    Convenience routine that will raise a TypeError if a variable's
    type does not match a list of expected types.

    Args:
    -----
    var : variable of any type
        The variable to check.

    var_type : type or tuple of types
        A single type definition (list, str, pandas.Series, etc.)
        or a tuple of type definitions.
    """
    if isinstance(var, var_type):
        return
    raise TypeError( f"{var} is not of type: {var_type}!")


def copy_file_to_dir(
        ifile,
        dest,
):
    """
    Convenience wrapper for shutil.copyfile, used to copy a file to
    a directory.

    Args
    ifile : str : Input file in original location
    dest  : str : Destination folder where ifile will be copied.
    """
    ifile = os.path.realpath(ifile)
    ofile = os.path.join(dest, os.path.basename(ifile))
    if not os.path.exists(ofile):
        copyfile(ifile, ofile)


def replace_whitespace(
        string,
        repl_char="_"
):
    """
    Replaces whitespace in a string with underscores.
    Useful for removing spaces in filename strings.

    Args
    string    : str : The input string
    repl_char : str : Replacement character (default is "_")

    Returns
    string    : str : String with whitespace replaced
    """
    verify_variable_type(string, str)
    verify_variable_type(repl_char, str)

    return repl_char.join(string.split())


def get_element_of_series(series, element):
    """
    Returns a specified element of a pd.Series object.

    Args
    serie   : pd.Series : A pd.Series object
    element : int       : Element of the pd.Series object to return

    Returns
    value   : various   : The returned element
    """
    verify_variable_type(series, Series)
    verify_variable_type(element, int)

    return list(series)[element]
