#!/usr/bin/env python3
"""
Creates species mass tables from GEOS-Chem benchmark simulation output.
"""
import os
import warnings
import gc
import numpy as np
import xarray as xr
from gcpy.grid import get_troposphere_mask
from gcpy.util import \
    compare_varnames, get_area_from_dataset, get_variables_from_dataset, \
    insert_text_into_file, print_totals, read_species_metadata, \
    rename_and_flip_gchp_rst_vars, replace_whitespace, verify_variable_type
from gcpy.units import convert_units
from gcpy.benchmark.modules.benchmark_utils import \
    diff_list_to_text
from gcpy.constants import \
    COL_WIDTH, ENCODING, SKIP_THESE_VARS, TABLE_WIDTH

# Suppress numpy divide by zero warnings to prevent output spam
np.seterr(divide="ignore", invalid="ignore")


def create_global_mass_table(
        refdata,
        refstr,
        devdata,
        devstr,
        varlist,
        met_and_masks,
        spcdb_files,
        ref_hdr_label="",
        dev_hdr_label="",
        trop_only=False,
        outfilename="GlobalMass_TropStrat.txt",
        verbose=False,
):
    """
    Creates a table of global masses for a list of species in contained in
    two data sets.  The data sets,  which typically represent output from two
    different model versions, are usually contained in netCDF data files.

    Args:
        refdata: xarray Dataset
            The first data set to be compared (aka "Reference").
        refstr: str
            A string that can be used to identify refdata
            (e.g. a model version number or other identifier).
        devdata: xarray Dataset
            The second data set to be compared (aka "Development").
        devstr: str
            A string that can be used to identify the data set specified
            by devfile (e.g. a model version number or other identifier).
        varlist: list of strings
            List of species concentation variable names to include
            in the list of global totals.
        met_and_masks: dict of xarray DataArray
            Dictionary containing the meterological variables and
            masks for the Ref and Dev datasets.
        spcdb_files : list
            Paths to species_database.yml files in Ref & Dev rundirs

    Keyword Args (optional):
        ref_hdr_label : str
            Label for Ref, placed after refstr in the file header
            Default value: ""
        dev_hdr_label : str
            Label for Dev, placed after devstr in the file header
            Default value: ""
        trop_only: bool
            Set this switch to True if you wish to print totals
            only for the troposphere.
            Default value: False (i.e. print whole-atmosphere totals).
        outfilename: str
            Name of the text file which will contain the table of
            emissions totals.
            Default value: "GlobalMass_TropStrat.txt"
        verbose: bool
            Set this switch to True if you wish to print out extra
            informational messages.
            Default value: False

    Remarks:
        This method is mainly intended for model benchmarking purposes,
        rather than as a general-purpose tool.

        Species metadata (such as molecular weights) are read from a
        YAML file called "species_database.yml".
    """

    # ==================================================================
    # Initialization
    # ==================================================================
    verify_variable_type(refdata, xr.Dataset)
    verify_variable_type(devdata, xr.Dataset)

    # Make sure required arguments are passed
    if varlist is None:
        raise ValueError('The "varlist" argument was not passed!')
    if met_and_masks is None:
        raise ValueError('The "met_and_masks" argument was not passed!')

    # Read the species database files in the Ref & Dev rundirs, and
    # return a dict containing metadata for each.
    ref_metadata, dev_metadata = read_species_metadata(
        spcdb_files,
        quiet=True
    )

    # Replace whitespace in the ref and dev labels
    refstr = replace_whitespace(refstr)
    devstr = replace_whitespace(devstr)

    # ==================================================================
    # Open file for output
    # ==================================================================

    # Create file
    try:
        f = open(outfilename, "w", encoding=ENCODING)
    except (IOError, OSError, FileNotFoundError) as e:
        raise e(f"Could not open {outfilename} for writing!") from e

    # Define a list for differences
    diff_list = []

    # Title strings
    title1 = "### Global mass (Gg) (Trop + Strat)"
    if trop_only:
        title1 = "### Global mass (Gg) (Trop only)"
    title2 = f"### Ref = {refstr} {ref_hdr_label}"
    title3 = f"### Dev = {devstr} {dev_hdr_label}"

    # Write a placeholder to the file that denotes where
    # the list of species with differences will be written
    placeholder = "@%% insert diff status here %%@"

    # Print header to file
    print("#" * TABLE_WIDTH, file=f)
    print(f"{title1 : <{TABLE_WIDTH-3}}{'###'}", file=f)
    print(f"{'###'  : <{TABLE_WIDTH-3}}{'###'}", file=f)
    print(f"{title2 : <{TABLE_WIDTH-3}}{'###'}", file=f)
    print(f"{title3 : <{TABLE_WIDTH-3}}{'###'}", file=f)
    print(f"{'###'  : <{TABLE_WIDTH-3}}{'###'}", file=f)
    print(f"{placeholder}", file=f)
    print("#" * TABLE_WIDTH, file=f)

    # Column headers
    print(f"{'' : <{COL_WIDTH-1}}{'Ref' : >{COL_WIDTH}}{'Dev' : >{COL_WIDTH}}{'Dev - Ref' : >{COL_WIDTH}}{'% diff' : >{COL_WIDTH}} {'diffs'}", file=f)

    # ==================================================================
    # Print global masses for all species
    #
    # NOTE: By this point, all secies will be in both Ref and Dev'
    # because we have added them in the calling routine
    # ==================================================================
    for v in varlist:

        # Get the species name
        spc_name = v.split("_")[1]

        # Get metadta for the given species
        ref_species_metadata = ref_metadata.get(spc_name)
        dev_species_metadata = dev_metadata.get(spc_name)
        if ref_species_metadata is None and dev_species_metadata is None:
            if verbose:
                msg = f"No metadata found for {spc_name} ... skippping"
                print(msg)
            continue

        # Specify target units
        target_units = "Gg"

        # ==============================================================
        # Convert units of Ref and save to a DataArray
        # (or skip if Ref contains NaNs everywhere)
        # ==============================================================
        refarray = refdata[v]
        if not np.isnan(refdata[v].values).all():
            refarray = convert_units(
                refarray,
                spc_name,
                ref_species_metadata,
                target_units,
                area_m2=met_and_masks["Ref_Area"],
                delta_p=met_and_masks["Ref_Delta_P"],
                box_height=met_and_masks["Ref_BxHeight"],
            )

        # ==============================================================
        # Convert units of Dev and save to a DataArray
        # (or skip if Dev contains NaNs everywhere)
        # ==============================================================
        devarray = devdata[v]
        if not np.isnan(devdata[v].values).all():
            devarray = convert_units(
                devarray,
                spc_name,
                dev_species_metadata,
                target_units,
                area_m2=met_and_masks["Dev_Area"],
                delta_p=met_and_masks["Dev_Delta_P"],
                box_height=met_and_masks["Dev_BxHeight"],
            )

        # ==============================================================
        # Print global masses for Ref and Dev
        # (we will mask out tropospheric boxes in print_totals)
        # ==============================================================
        if trop_only:
            print_totals(
                refarray,
                devarray,
                f,
                diff_list,
                masks=met_and_masks
            )
        else:
            print_totals(
                refarray,
                devarray,
                f,
                diff_list
            )

    # ==================================================================
    # Cleanup and quit
    # ==================================================================

    # Close file
    f.close()

    # Reopen file and replace placeholder text by diff_text
    insert_text_into_file(
        filename=outfilename,
        search_text=placeholder,
        replace_text=diff_list_to_text(
            refstr,
            devstr,
            diff_list,
            fancy_format=True
        ),
        width=TABLE_WIDTH
    )


def create_mass_accumulation_table(
        refdatastart,
        refdataend,
        refstr,
        refperiodstr,
        devdatastart,
        devdataend,
        devstr,
        devperiodstr,
        varlist,
        met_and_masks,
        label,
        spcdb_files,
        trop_only=False,
        outfilename="GlobalMassAccum_TropStrat.txt",
        verbose=False,
):
    """
    Creates a table of global mass accumulation for a list of species in
    two data sets.  The data sets, which typically represent output from two
    different model versions, are usually contained in netCDF data files.

    Args:
        refdatastart: xarray Dataset
            The first data set to be compared (aka "Reference").
        refdataend: xarray Dataset
            The first data set to be compared (aka "Reference").
        refstr: str
            A string that can be used to identify refdata
            (e.g. a model version number or other identifier).
        refperiodstr: str
            Ref simulation period start and end
        devdatastart: xarray Dataset
            The second data set to be compared (aka "Development").
        devdataend: xarray Dataset
            The second data set to be compared (aka "Development").
        devstr: str
            A string that can be used to identify the data set specified
            by devfile (e.g. a model version number or other identifier).
        devperiodstr: str
            Ref simulation period start and end
        varlist: list of strings
            List of species concentation variable names to include
            in the list of global totals.
        met_and_masks: dict of xarray DataArray
            Dictionary containing the meterological variables and
            masks for the Ref and Dev datasets.
        label: str
            Label to go in the header string.  Can be used to
            pass the month & year.
        spcdb_files : list
            Paths to species_database.yml files in Ref & Dev rundirs

    Keyword Args (optional):
        trop_only: bool
            Set this switch to True if you wish to print totals
            only for the troposphere.
            Default value: False (i.e. print whole-atmosphere totals).
        outfilename: str
            Name of the text file which will contain the table of
            emissions totals.
            Default value: "GlobalMass_TropStrat.txt"
        verbose: bool
            Set this switch to True if you wish to print out extra
            informational messages.
            Default value: False

    Remarks:
        This method is mainly intended for model benchmarking purposes,
        rather than as a general-purpose tool.

        Species metadata (such as molecular weights) are read from a
        YAML file called "species_database.yml".
    """

    # ==================================================================
    # Initialization
    # ==================================================================
    verify_variable_type(refdatastart, xr.Dataset)
    verify_variable_type(refdataend, xr.Dataset)
    verify_variable_type(devdatastart, xr.Dataset)
    verify_variable_type(devdataend, xr.Dataset)

    # Make sure required arguments are passed
    if varlist is None:
        raise ValueError('The "varlist" argument was not passed!')
    if met_and_masks is None:
        raise ValueError('The "met_and_masks" argument was not passed!')

    # Read the species database files in the Ref & Dev rundirs, and
    # return a dict containing metadata for each.
    ref_metadata, dev_metadata = read_species_metadata(
        spcdb_files,
        quiet=True
    )

    # Replace whitespace in the ref and dev labels
    refstr = replace_whitespace(refstr)
    devstr = replace_whitespace(devstr)

    # ==================================================================
    # Open file for output
    # ==================================================================

    # Create file
    try:
        f = open(outfilename, "w", encoding=ENCODING)
    except (IOError, OSError, FileNotFoundError) as e:
        raise e(f"Could not open {outfilename} for writing!") from e

    # Define a list for differences
    diff_list = []

    # Title strings
    title1 = f"### Global mass accumulation (Gg) {label} (Trop + Strat)"
    if trop_only:
        title1 = f"### Global mass accumulation (Gg) {label} (Trop only)"
    title2 = "### Computed as change in instantaneous mass across period"
    title3 = f"### Ref = {refstr}"
    title4 = f"### Dev = {devstr}"
    title5 = f"### Ref period: {refperiodstr}"
    title6 = f"### Dev period: {devperiodstr}"

    # Write a placeholder to the file that denotes where
    # the list of species with differences will be written
    placeholder = "@%% insert diff status here %%@"

    # Print header to file
    print("#" * TABLE_WIDTH, file=f)
    print(f"{title1 : <{TABLE_WIDTH-3}}{'###'}", file=f)
    print(f"{'###'  : <{TABLE_WIDTH-3}}{'###'}", file=f)
    print(f"{title2 : <{TABLE_WIDTH-3}}{'###'}", file=f)
    print(f"{'###'  : <{TABLE_WIDTH-3}}{'###'}", file=f)
    print(f"{title3 : <{TABLE_WIDTH-3}}{'###'}", file=f)
    print(f"{title4 : <{TABLE_WIDTH-3}}{'###'}", file=f)
    print(f"{'###'  : <{TABLE_WIDTH-3}}{'###'}", file=f)
    print(f"{title5 : <{TABLE_WIDTH-3}}{'###'}", file=f)
    print(f"{title6 : <{TABLE_WIDTH-3}}{'###'}", file=f)
    print(f"{'###'  : <{TABLE_WIDTH-3}}{'###'}", file=f)
    print(f"{placeholder}", file=f)
    print("#" * TABLE_WIDTH, file=f)

    # Column headers
    print(f"{'' : <{COL_WIDTH-1}}{'Ref' : >{COL_WIDTH}}{'Dev' : >{COL_WIDTH}}{'Dev - Ref' : >{COL_WIDTH}}{'% diff' : >{COL_WIDTH}} {'diffs'}", file=f)

    # ==================================================================
    # Print global masses for all species
    #
    # NOTE: By this point, all secies will be in both Ref and Dev'
    # because we have added them in the calling routine
    # ==================================================================
    for v in varlist:

        # Get the species name
        spc_name = v.split("_")[1]

        # Get a list of metadata for the given species
        ref_species_metadata = ref_metadata.get(spc_name)
        dev_species_metadata = dev_metadata.get(spc_name)
        if ref_species_metadata is None and dev_species_metadata is None:
            if verbose:
                msg = f"No properties found for {spc_name} ... skippping"
                print(msg)
            continue

        # Specify target units
        target_units = "Gg"

        # ==============================================================
        # Convert units of Ref and save to a DataArray
        # (or skip if Ref contains NaNs everywhere)
        # ==============================================================
        refarrays = refdatastart[v]
        if not np.isnan(refdatastart[v].values).all():
            refarrays = convert_units(
                refarrays,
                spc_name,
                ref_species_metadata,
                target_units,
                area_m2=met_and_masks["Refs_Area"],
                delta_p=met_and_masks["Refs_Delta_P"],
                box_height=met_and_masks["Refs_BxHeight"],
            )

        refarraye = refdataend[v]
        if not np.isnan(refdataend[v].values).all():
            refarraye = convert_units(
                refarraye,
                spc_name,
                ref_species_metadata,
                target_units,
                area_m2=met_and_masks["Refe_Area"],
                delta_p=met_and_masks["Refe_Delta_P"],
                box_height=met_and_masks["Refe_BxHeight"],
            )

        refarray = refarrays
        refarray.values = refarraye.values - refarrays.values

        # ==============================================================
        # Convert units of Dev and save to a DataArray
        # (or skip if Dev contains NaNs everywhere)
        # ==============================================================
        devarrays = devdatastart[v]
        if not np.isnan(devdatastart[v].values).all():
            devarrays = convert_units(
                devarrays,
                spc_name,
                dev_species_metadata,
                target_units,
                area_m2=met_and_masks["Devs_Area"],
                delta_p=met_and_masks["Devs_Delta_P"],
                box_height=met_and_masks["Devs_BxHeight"],
            )
        #print('devarrays: {}'.format(devarrays.values))

        devarraye = devdataend[v]
        if not np.isnan(devdataend[v].values).all():
            devarraye = convert_units(
                devarraye,
                spc_name,
                dev_species_metadata,
                target_units,
                area_m2=met_and_masks["Deve_Area"],
                delta_p=met_and_masks["Deve_Delta_P"],
                box_height=met_and_masks["Deve_BxHeight"],
            )

        devarray = devarrays
        devarray.values = devarraye.values - devarrays.values

        # ==============================================================
        # Print global masses for Ref and Dev
        # (we will mask out tropospheric boxes in print_totals)
        # ==============================================================
        # ewl: for now trop_only is always false for accumulation table
        if trop_only:
            print_totals(
                refarray,
                devarray,
                f,
                diff_list,
                masks=met_and_masks
            )
        else:
            print_totals(
                refarray,
                devarray,
                f,
                diff_list
            )

    # ==================================================================
    # Cleanup and quit
    # ==================================================================

    # Close file
    f.close()

    # Reopen file and replace placeholder text by diff_text
    insert_text_into_file(
        filename=outfilename,
        search_text=placeholder,
        replace_text=diff_list_to_text(
            refstr,
            devstr,
            diff_list,
            fancy_format=True
        ),
        width=TABLE_WIDTH
    )


def make_benchmark_mass_tables(
        ref,
        refstr,
        dev,
        devstr,
        spcdb_files,
        varlist=None,
        dst="./benchmark",
        subdst=None,
        overwrite=False,
        verbose=False,
        ref_hdr_label="",
        dev_hdr_label="",
        ref_met_extra=None,
        dev_met_extra=None,
):
    """
    Creates a text file containing global mass totals by species and
    category for benchmarking purposes.

    Args:
        reflist: str
            Pathname that will constitute
            the "Ref" (aka "Reference") data set.
        refstr: str
            A string to describe ref (e.g. version number)
        dev: list of str
            Pathname that will constitute
            the "Dev" (aka "Development") data set.  The "Dev"
            data set will be compared against the "Ref" data set.
        devstr: str
            A string to describe dev (e.g. version number)
        spcdb_files : list
            Paths to species_database.yml files in Ref & Dev rundirs

    Keyword Args (optional):
        varlist: list of str
            List of variables to include in the list of totals.
            If omitted, then all variables that are found in either
            "Ref" or "Dev" will be included.  The varlist argument
            can be a useful way of reducing the number of
            variables during debugging and testing.
            Default value: None
        dst: str
            A string denoting the destination folder where the file
            containing emissions totals will be written.
            Default value: ./benchmark
        subdst: str
            A string denoting the sub-directory of dst where PDF
            files containing plots will be written.  In practice,
            subdst is only needed for the 1-year benchmark output,
            and denotes a date string (such as "Jan2016") that
            corresponds to the month that is being plotted.
            Default value: None
        overwrite: bool
            Set this flag to True to overwrite files in the
            destination folder (specified by the dst argument).
            Default value: False
        verbose: bool
            Set this flag to True to print extra informational output.
            Default value: False.
        ref_hdr_label : str
            Label for Ref, placed after refstr in the file header
            Default value: ""
        dev_hdr_label : str
            Label for Dev, placed after devstr in the file header
            Default value: ""
        ref_met_extra: str
            Path to ref Met file containing area data for use with restart files
            which do not contain the Area variable.
            Default value: ''
        dev_met_extra: str
            Path to dev Met file containing area data for use with restart files
            which do not contain the Area variable.
            Default value: ''
    """
    # ==================================================================
    # Initialization
    # ==================================================================

    # Replace whitespace in the ref and dev labels
    refstr = replace_whitespace(refstr)
    devstr = replace_whitespace(devstr)

    # ==================================================================
    # Define destination directory
    # ==================================================================
    if os.path.isdir(dst) and not overwrite:
        msg = "Directory {} exists. Pass overwrite=True to overwrite " \
            + "files in that directory, if any."
        msg = msg.format(dst)
        raise ValueError(msg)
    if not os.path.isdir(dst):
        try:
            os.makedirs(dst)
        except FileExistsError:
            pass

    # ==================================================================
    # Read data from netCDF into Dataset objects
    # ==================================================================

    # Read data
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=xr.SerializationWarning)
        refds = xr.open_dataset(ref, drop_variables=SKIP_THESE_VARS)
        devds = xr.open_dataset(dev, drop_variables=SKIP_THESE_VARS)

    # ==================================================================
    # Update GCHP restart dataset (if any)
    # ==================================================================

    # If the data is from a GCHP restart file, rename variables and
    # flip levels to match the GEOS-Chem Classic naming and level
    # conventions.  Otherwise no changes will be made.
    refds = rename_and_flip_gchp_rst_vars(refds)
    devds = rename_and_flip_gchp_rst_vars(devds)

    # ==================================================================
    # Make sure that all necessary meteorological variables are found
    # ==================================================================
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=xr.SerializationWarning)

        # Find the area variable in Ref
        if ref_met_extra is None:
            ref_area = get_area_from_dataset(refds)
        else:
            ref_area = get_area_from_dataset(
                xr.open_dataset(
                    ref_met_extra,
                    drop_variables=SKIP_THESE_VARS
                )
            )

        # Find the area variable in Dev
        if dev_met_extra is None:
            dev_area = get_area_from_dataset(devds)
        else:
            dev_area = get_area_from_dataset(
                xr.open_dataset(
                    dev_met_extra,
                    drop_variables=SKIP_THESE_VARS
                )
            )

    # Find required meteorological variables in Ref
    # (or exit with an error if we can't find them)
    metvar_list = ["Met_DELPDRY", "Met_BXHEIGHT", "Met_TropLev"]
    refmet = get_variables_from_dataset(refds, metvar_list)
    devmet = get_variables_from_dataset(devds, metvar_list)

    # ==================================================================
    # Make sure that all necessary species are found
    # ==================================================================

    # Get lists of variables names in datasets
    vardict = compare_varnames(refds, devds, quiet=not verbose)
    commonvars = vardict["commonvars3D"]
    refonly = vardict['refonly']
    devonly = vardict['devonly']

    # Narrow down the lists to only include species
    commonspc = [v for v in commonvars if "SpeciesRst_" in v]
    refonlyspc = [v for v in refonly if v.startswith('SpeciesRst_')]
    devonlyspc = [v for v in devonly if v.startswith('SpeciesRst_')]

    # Add ref only species to dev dataset with all nan values
    if refonlyspc:
        for v in refonlyspc:
            devds[v] = devds[commonspc[0]]
            devds[v].data = np.full(devds[v].shape, np.nan)
            devds[v].attrs['units'] = refds[v].units
            commonspc.append(v)

    # Add dev only species to ref dataset with all nan values
    if devonlyspc:
        for v in devonlyspc:
            refds[v] = refds[commonspc[0]]
            refds[v].data = np.full(refds[v].shape, np.nan)
            devds[v].attrs['units'] = refds[v].units
            commonspc.append(v)

    # Set list of variables to print in mass table. If this list was passed
    # as argument, check that all the vars are now in commonspc to ensure
    # in both datasets.
    if varlist:
        for v in varlist:
            if v not in commonspc:
                raise ValueError(
                    f"{dst} folder error: Variable {v} in varlist passed to make_benchmark_mass_tables is not present in Ref and Dev datasets"
                )
    else:
        varlist = commonspc

    # Sort the list of species to be printed alphabetically
    varlist.sort()

    # ==================================================================
    # Create the mask arrays for the troposphere for Ref and Dev
    # ==================================================================
    ref_tropmask = get_troposphere_mask(refmet)
    dev_tropmask = get_troposphere_mask(devmet)

    # ==================================================================
    # Create a dictionary to hold all of the meterological
    # variables and mask variables that we need to pass down
    # ==================================================================
    met_and_masks = {
        "Ref_Area": ref_area,
        "Dev_Area": dev_area,
        "Ref_Delta_P": refmet["Met_DELPDRY"],
        "Dev_Delta_P": devmet["Met_DELPDRY"],
        "Ref_BxHeight": refmet["Met_BXHEIGHT"],
        "Dev_BxHeight": devmet["Met_BXHEIGHT"],
        "Ref_TropMask": ref_tropmask,
        "Dev_TropMask": dev_tropmask,
    }

    # ==================================================================
    # Create global mass table
    # ==================================================================
    if subdst is not None:
        mass_filename = f"GlobalMass_TropStrat_{subdst}.txt"
    else:
        mass_filename = "GlobalMass_TropStrat.txt"
    mass_file = os.path.join(dst, mass_filename)
    create_global_mass_table(
        refds,
        refstr,
        devds,
        devstr,
        varlist,
        met_and_masks,
        spcdb_files,
        ref_hdr_label=ref_hdr_label,
        dev_hdr_label=dev_hdr_label,
        outfilename=mass_file,
        verbose=verbose,
    )

    # ==================================================================
    # Create tropospheric mass table
    # ==================================================================

    # If a file name has not been specified, then use the "filename"
    # keyword argument.  Otherwise generate a default filename.
    if subdst is not None:
        mass_filename = f"GlobalMass_Trop_{subdst}.txt"
    else:
        mass_filename = 'GlobalMass_Trop.txt'
    mass_file = os.path.join(dst, mass_filename)
    create_global_mass_table(
        refds,
        refstr,
        devds,
        devstr,
        varlist,
        met_and_masks,
        spcdb_files,
        ref_hdr_label=ref_hdr_label,
        dev_hdr_label=dev_hdr_label,
        outfilename=mass_file,
        trop_only=True,
        verbose=verbose,
    )

    # -------------------------------------------
    # Clean up
    # -------------------------------------------
    del refds
    del devds
    gc.collect()


def make_benchmark_mass_accumulation_tables(
        ref_start,
        ref_end,
        refstr,
        refperiodstr,
        dev_start,
        dev_end,
        devstr,
        devperiodstr,
        spcdb_files,
        varlist=None,
        dst="./benchmark",
        subdst=None,
        overwrite=False,
        verbose=False,
        label="at end of simulation",
):
    """
    Creates a text file containing global mass totals by species and
    category for benchmarking purposes.

    Args:
        ref_start: list of str
            Pathname that will constitute
            the "Ref" (aka "Reference") data set.
        ref_end: list of str
            Pathname that will constitute
            the "Ref" (aka "Reference") data set.
        refstr: str
            A string to describe ref (e.g. version number)
        refperiodstr: str
            Ref simulation period start and end
        dev_start: list of str
            Pathname that will constitute
            the "Dev" (aka "Development") data set.  The "Dev"
            data set will be compared against the "Ref" data set.
        dev_end: list of str
            Pathname that will constitute
            the "Dev" (aka "Development") data set.  The "Dev"
            data set will be compared against the "Ref" data set.
        devstr: str
            A string to describe dev (e.g. version number)
        devperiodstr: str
            Dev simulation period start and end
        spcdb_files : list
            Paths to species_database.yml files in Ref & Dev rundirs

    Keyword Args (optional):
        varlist: list of str
            List of variables to include in the list of totals.
            If omitted, then all variables that are found in either
            "Ref" or "Dev" will be included.  The varlist argument
            can be a useful way of reducing the number of
            variables during debugging and testing.
            Default value: None
        dst: str
            A string denoting the destination folder where the file
            containing emissions totals will be written.
            Default value: ./benchmark
        subdst: str
            A string denoting the sub-directory of dst where PDF
            files containing plots will be written.  In practice,
            subdst is only needed for the 1-year benchmark output,
            and denotes a date string (such as "Jan2016") that
            corresponds to the month that is being plotted.
            Default value: None
        overwrite: bool
            Set this flag to True to overwrite files in the
            destination folder (specified by the dst argument).
            Default value: False
        verbose: bool
            Set this flag to True to print extra informational output.
            Default value: False.
    """
    # ==================================================================
    # Initialization
    # ==================================================================

    # Replace whitespace in the ref and dev labels
    refstr = replace_whitespace(refstr)
    devstr = replace_whitespace(devstr)

    # ==================================================================
    # Define destination directory
    # ==================================================================
    if os.path.isdir(dst) and not overwrite:
        msg = "Directory {} exists. Pass overwrite=True to overwrite " \
            + "files in that directory, if any."
        msg = msg.format(dst)
        raise ValueError(msg)
    if not os.path.isdir(dst):
        try:
            os.makedirs(dst)
        except FileExistsError:
            pass

    # ==================================================================
    # Read data from netCDF into Dataset objects
    # ==================================================================

    print('Creating mass accumulation tables from four restart files:')
    print(f'   Ref start: {ref_start}')
    print(f'   Ref end:   {ref_end}')
    print(f'   Dev start: {dev_start}')
    print(f'   Dev end:   {dev_end}')

    # Read data
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=xr.SerializationWarning)
        refSds = xr.open_dataset(ref_start, drop_variables=SKIP_THESE_VARS)
        refEds = xr.open_dataset(ref_end, drop_variables=SKIP_THESE_VARS)
        devSds = xr.open_dataset(dev_start, drop_variables=SKIP_THESE_VARS)
        devEds = xr.open_dataset(dev_end, drop_variables=SKIP_THESE_VARS)

    # ==================================================================
    # Update GCHP restart dataset if needed
    # ==================================================================

    # If the data is from a GCHP restart file, rename variables and
    # flip levels to match the GEOS-Chem Classic naming and level
    # conventions.  Otherwise no changes will be made.
    refSds = rename_and_flip_gchp_rst_vars(refSds)
    refEds = rename_and_flip_gchp_rst_vars(refEds)
    devSds = rename_and_flip_gchp_rst_vars(devSds)
    devEds = rename_and_flip_gchp_rst_vars(devEds)

    # Add area to start restart dataset if area in end but not start
    # Need to consider area variable names used in both GC-Classic and GCHP
    # Should put this in a function (todo)
    refSkeys = refSds.data_vars.keys()
    refEkeys = refEds.data_vars.keys()
    devSkeys = devSds.data_vars.keys()
    devEkeys = devEds.data_vars.keys()
    areaVars = ["Met_AREAM2", "AREA"]
    for areaVar in areaVars:
        if areaVar in refEkeys and areaVar not in refSkeys:
            refSds[areaVar] = refEds[areaVar]
        if areaVar in devEkeys and areaVar not in devSkeys:
            devSds[areaVar] = devEds[areaVar]

    # ==================================================================
    # Make sure that all necessary meteorological variables are found
    # ==================================================================
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=xr.SerializationWarning)

        # Find the area variable in Ref
        refs_area = get_area_from_dataset(refSds)
        refe_area = get_area_from_dataset(refEds)

        # Find the area variable in Dev
        devs_area = get_area_from_dataset(devSds)
        deve_area = get_area_from_dataset(devEds)

    # Find required meteorological variables in Ref
    # (or exit with an error if we can't find them)
    metvar_list = ["Met_DELPDRY", "Met_BXHEIGHT", "Met_TropLev"]
    refsmet = get_variables_from_dataset(refSds, metvar_list)
    refemet = get_variables_from_dataset(refEds, metvar_list)
    devsmet = get_variables_from_dataset(devSds, metvar_list)
    devemet = get_variables_from_dataset(devEds, metvar_list)

    # ==================================================================
    # Make sure that all necessary species are found
    # ==================================================================

    # Get lists of variables names in datasets
    vardict = compare_varnames(refSds, devSds, quiet=not verbose)
    commonvars = vardict["commonvars3D"]
    refonly = vardict['refonly']
    devonly = vardict['devonly']

    # Narrow down the lists to only include species
    commonspc = [v for v in commonvars if "SpeciesRst_" in v]
    refonlyspc = [v for v in refonly if v.startswith('SpeciesRst_')]
    devonlyspc = [v for v in devonly if v.startswith('SpeciesRst_')]

    # Add ref only species to dev dataset with all nan values
    if refonlyspc:
        for v in refonlyspc:
            devSds[v] = devSds[commonspc[0]]
            devSds[v].data = np.full(devSds[v].shape, np.nan)
            devSds[v].attrs['units'] = refSds[v].units
            devEds[v] = devEds[commonspc[0]]
            devEds[v].data = np.full(devEds[v].shape, np.nan)
            devEds[v].attrs['units'] = refEds[v].units
            commonspc.append(v)

    # Add dev only species to ref dataset with all nan values
    if devonlyspc:
        for v in devonlyspc:
            refSds[v] = refSds[commonspc[0]]
            refSds[v].data = np.full(refSds[v].shape, np.nan)
            devSds[v].attrs['units'] = refSds[v].units
            refEds[v] = refEds[commonspc[0]]
            refEds[v].data = np.full(refEds[v].shape, np.nan)
            devEds[v].attrs['units'] = refEds[v].units
            commonspc.append(v)

    # Set list of variables to print in mass table. If this list was passed
    # as argument, check that all the vars are now in commonspc to ensure
    # in both datasets.
    if varlist:
        for v in varlist:
            if v not in commonspc:
                raise ValueError(
                    f"{dst} folder error: Variable {v} in varlist passed to make_benchmark_mass_tables is not present in Ref and Dev datasets"
                )
    else:
        varlist = commonspc

    # Sort the list of species to be printed alphabetically
    varlist.sort()

    # ==================================================================
    # Create the mask arrays for the troposphere for Ref and Dev
    # ==================================================================
    refs_tropmask = get_troposphere_mask(refsmet)
    refe_tropmask = get_troposphere_mask(refemet)
    devs_tropmask = get_troposphere_mask(devsmet)
    deve_tropmask = get_troposphere_mask(devemet)

    # ==================================================================
    # Create a dictionary to hold all of the meterological
    # variables and mask variables that we need to pass down
    # ==================================================================
    met_and_masks = {
        "Refs_Area": refs_area,
        "Refe_Area": refe_area,
        "Devs_Area": devs_area,
        "Deve_Area": deve_area,
        "Refs_Delta_P": refsmet["Met_DELPDRY"],
        "Refe_Delta_P": refemet["Met_DELPDRY"],
        "Devs_Delta_P": devsmet["Met_DELPDRY"],
        "Deve_Delta_P": devemet["Met_DELPDRY"],
        "Refs_BxHeight": refsmet["Met_BXHEIGHT"],
        "Refe_BxHeight": refemet["Met_BXHEIGHT"],
        "Devs_BxHeight": devsmet["Met_BXHEIGHT"],
        "Deve_BxHeight": devemet["Met_BXHEIGHT"],
        "Refs_TropMask": refs_tropmask,
        "Refe_TropMask": refe_tropmask,
        "Devs_TropMask": devs_tropmask,
        "Deve_TropMask": deve_tropmask,
    }

    # ==================================================================
    # Create global mass accumulation table
    # ==================================================================
    if subdst is not None:
        mass_filename = f"GlobalMassAccumulation_TropStrat_{subdst}.txt"
    else:
        mass_filename = "GlobalMassAccumulation_TropStrat.txt"
    mass_file = os.path.join(dst, mass_filename)
    create_mass_accumulation_table(
        refSds,
        refEds,
        refstr,
        refperiodstr,
        devSds,
        devEds,
        devstr,
        devperiodstr,
        varlist,
        met_and_masks,
        label,
        outfilename=mass_file,
        verbose=verbose,
        spcdb_files=spcdb_files,
    )

    ## ==================================================================
    ## Create tropospheric mass table
    ## ==================================================================
    #if subdst is not None:
    #    mass_filename = f"GlobalMassAccumulation_Trop_{subdst}.txt"
    #else:
    #    mass_filename = 'GlobalMassAccumulation_Trop.txt'
    #mass_file = os.path.join(dst, mass_filename)
    #create_mass_accumulation_table(
    #    refSds,
    #    refEds,
    #    refstr,
    #    devSds,
    #    devEds,
    #    devstr,
    #    varlist,
    #    met_and_masks,
    #    label,
    #    outfilename=mass_file,
    #    trop_only=True,
    #    verbose=verbose,
    #    spcdb_files=spcdb_files,
    #)

    # -------------------------------------------
    # Clean up
    # -------------------------------------------
    del refSds
    del refEds
    del devSds
    del devEds
    gc.collect()
