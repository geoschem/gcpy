#!/usr/bin/env python3
"""
Generates aerosol plots and tables from GEOS-Chem benchmark
simulation output.
"""
import os
import gc
import numpy as np
import xarray as xr
from gcpy.regrid import create_regridders
from gcpy.grid import get_troposphere_mask
from gcpy.util import \
    add_bookmarks_to_pdf, add_missing_variables, \
    compare_varnames, dataset_reader, \
    dataset_mean, make_directory, read_config_file, \
    read_species_metadata, replace_whitespace
from gcpy.constants import \
    ENCODING, MW_AIR_g, SKIP_THESE_VARS
from gcpy.plot.compare_single_level import compare_single_level
from gcpy.benchmark.modules.benchmark_utils import \
    AOD_SPC, rename_speciesconc_to_speciesconcvv

# Suppress numpy divide by zero warnings to prevent output spam
np.seterr(divide="ignore", invalid="ignore")


def make_benchmark_aod_plots(
        ref,
        refstr,
        dev,
        devstr,
        spcdb_files,
        varlist=None,
        dst="./benchmark",
        subdst=None,
        cmpres=None,
        overwrite=False,
        verbose=False,
        log_color_scale=False,
        sigdiff_files=None,
        weightsdir='.',
        n_job=-1,
        time_mean=False,
):
    """
    Creates PDF files containing plots of column aerosol optical
    depths (AODs) for model benchmarking purposes.

    Args:
        ref: str
            Path name for the "Ref" (aka "Reference") data set.
        refstr: str
            A string to describe ref (e.g. version number)
        dev: str
            Path name for the "Dev" (aka "Development") data set.
            This data set will be compared against the "Reference"
            data set.
        devstr: str
            A string to describe dev (e.g. version number)
        spcdb_files : list
            Paths to species_database.yml files in Ref & Dev rundirs

    Keyword Args (optional):
        varlist: list of str
            List of AOD variables to plot.  If not passed, then all
            AOD variables common to both Dev and Ref will be plotted.
            Use the varlist argument to restrict the number of
            variables plotted to the pdf file when debugging.
            Default value: None
        dst: str
            A string denoting the destination folder where a
            PDF file  containing plots will be written.
            Default value: ./benchmark.
        subdst: str
            A string denoting the sub-directory of dst where PDF
            files containing plots will be written.  In practice,
            subdst is only needed for the 1-year benchmark output,
            and denotes a date string (such as "Jan2016") that
            corresponds to the month that is being plotted.
            Default value: None
        cmpres: string
            Grid resolution at which to compare ref and dev data, e.g. '1x1.25'
        overwrite: bool
            Set this flag to True to overwrite files in the
            destination folder (specified by the dst argument).
            Default value: False.
        verbose: bool
            Set this flag to True to print extra informational output.
            Default value: False
        log_color_scale: bool
            Set this flag to True to enable plotting data (not diffs)
            on a log color scale.
            Default value: False
        sigdiff_files: list of str
            Filenames that will contain the list of quantities having
            having significant differences in the column AOD plots.
            These lists are needed in order to fill out the benchmark
            approval forms.
            Default value: None
        weightsdir: str
            Directory in which to place (and possibly reuse) xESMF regridder
            netCDF files.
            Default value: '.'
        n_job: int
            Defines the number of simultaneous workers for parallel plotting.
            Set to 1 to disable parallel plotting. Value of -1 allows the
            application to decide.
            Default value: -1
        time_mean : bool
            Determines if we should average the datasets over time
            Default value: False
    """
    # ==================================================================
    # Initialization and also read data
    # ==================================================================

    # Create destination plots directory
    make_directory(dst, overwrite)

    # Create the "Aerosols" directory as a subfolder of dst.
    # If subdst is passed, then create a subdirectory of the "Aerosols"
    # directory (e.g. which is needed for the 1-year benchmarks).
    aoddir = os.path.join(dst, "Aerosols")
    if not os.path.isdir(aoddir):
        os.mkdir(aoddir)
    if subdst is not None:
        aoddir = os.path.join(aoddir, subdst)
        if not os.path.isdir(aoddir):
            os.mkdir(aoddir)
        extra_title_txt = subdst
    else:
        extra_title_txt = None

    # Replace whitespace in the ref and dev labels
    refstr = replace_whitespace(refstr)
    devstr = replace_whitespace(devstr)

    # Get the function that will read file(s) into a dataset
    reader = dataset_reader(time_mean, verbose=verbose)

    # Read the Ref dataset
    try:
        refds = reader(ref, drop_variables=SKIP_THESE_VARS)
    except (OSError, IOError, FileNotFoundError) as e:
        raise e(f"Could not find Ref file: {ref}") from e

    # Read the Dev dataset
    try:
        devds = reader(dev, drop_variables=SKIP_THESE_VARS)
    except (OSError, IOError, FileNotFoundError) as e:
        raise e(f"Could not find Ref file: {dev}") from e

    # Compute mean of data over the time dimension (if time_mean=True)
    if time_mean:
        refds = dataset_mean(refds)
        devds = dataset_mean(devds)

    # Create regridding files if necessary
    [_ for _ in create_regridders(refds, devds, weightsdir=weightsdir)]

    # NOTE: GCHP diagnostic variable exports are defined before the
    # input.geos file is read.  This means "WL1" will not have been
    # replaced with "550nm" in the variable names.  Do this string
    # replace operation here, so that we can compare GCC and GCHP
    # data directly. (bmy, 4/29/19)
    with xr.set_options(keep_attrs=True):

        # Rename variables in the Ref dataset
        old2new = {}
        for v in refds.data_vars.keys():
            if "WL1" in v:
                newname = v.replace("WL1", "550nm")
                old2new[v] = newname
        refds = refds.rename(old2new)

        # Rename variables in the Dev dataset
        old2new = {}
        for v in devds.data_vars.keys():
            if "WL1" in v:
                newname = v.replace("WL1", "550nm")
                old2new[v] = newname
        devds = devds.rename(old2new)

    # Make sure that Ref and Dev datasets have the same variables.
    # Variables that are in Ref but not in Dev will be added to Dev
    # with all missing values (NaNs). And vice-versa.
    [refds, devds] = add_missing_variables(refds, devds)

    # Find common AOD variables in both datasets
    # (or use the varlist passed via keyword argument)
    if varlist is None:
        vardict = compare_varnames(refds, devds, quiet=not verbose)
        cmn3D = vardict["commonvars3D"]
        varlist = [v for v in cmn3D if "AOD" in v]

    # Dictionary and list for new display names
    newvars = read_config_file(
        os.path.join(
            os.path.dirname(__file__),
            AOD_SPC
        ),
        quiet=True
    )
    newvarlist = []

    # ==================================================================
    # Compute the total AOD by summing over the constituent members
    # ==================================================================

    # Take one of the variables so we can use its dims, coords,
    # attrs to create the DataArray object for total AOD
    v = varlist[0]

    # Create a DataArray to hold total column AOD
    # This is the same shape as the DataArray objects in refds
    reftot = xr.DataArray(
        np.zeros(refds[v].values.shape),
        name="AODTotal",
        dims=refds[v].dims,
        coords=refds[v].coords,
        attrs=refds[v].attrs,
    )

    # Create a DataArray to hold total column AOD
    # This is the same shape as the DataArray objects in devds
    devtot = xr.DataArray(
        np.zeros(devds[v].values.shape),
        name="AODTotal",
        dims=devds[v].dims,
        coords=devds[v].coords,
        attrs=devds[v].attrs,
    )

    # Save the variable attributes so that we can reattach them
    refattrs = reftot.attrs
    devattrs = devtot.attrs

    # Compute the sum of all AOD variables
    # Avoid double-counting the following:
    # (1) Individual dust AODs, use the column AOD total instead
    # (2) SOA from aqueous isoprene, which is already accounted
    #     for in AODHyg550nm_OCPI.  Also see Github issue:
    #     https://github.com/geoschem/gcpy/issues/65
    for v in varlist:
        if "_bin" in v or "AODSOAfromAqIsoprene550nm" in v:
            continue
        reftot = reftot + refds[v]
        devtot = devtot + devds[v]

    # Reattach the variable attributes
    reftot.name = "AODTotal"
    reftot.attrs = refattrs
    reftot.attrs["long_name"] = "Total aerosol optical depth"
    devtot.name = "AODTotal"
    devtot.attrs = devattrs
    devtot.attrs["long_name"] = "Total aerosol optical depth"

    # Merge these variables back into the dataset
    refds = xr.merge([refds, reftot])
    devds = xr.merge([devds, devtot])

    # Also add AODTotal to the list
    varlist.append("AODTotal")

    # ==================================================================
    # Compute column AODs
    # Create a new DataArray for each column AOD variable,
    # using the new display name, and preserving attributes.
    # Merge the new DataArrays back into the DataSets.
    # ==================================================================
    for v in varlist:

        # Get the new name for each AOD variable (it's easier to display)
        if v in newvars:
            newname = newvars[v]
            newvarlist.append(newname)
        else:
            raise ValueError(f"Could not find a display name for {v}")

        # Don't clobber existing DataArray and Dataset attributes
        with xr.set_options(keep_attrs=True):

            # Add column AOD of newname to Ref
            array = refds[v].sum(dim="lev")
            array.name = newname
            array.attrs["units"] = "1"
            refds = xr.merge([refds, array])

            # Add column AOD of newname to Dev
            array = devds[v].sum(dim="lev")
            array.name = newname
            array.attrs["units"] = "1"
            devds = xr.merge([devds, array])

    # ==================================================================
    # Create the plots
    # ==================================================================
    if subdst is not None:
        pdfname = os.path.join(
            aoddir,
            f"Aerosols_ColumnOptDepth_{subdst}.pdf"
        )
    else:
        pdfname = os.path.join(
            aoddir,
            "Aerosols_ColumnOptDepth.pdf"
        )

    diff_aod = []
    compare_single_level(
        refds,
        refstr,
        devds,
        devstr,
        varlist=newvarlist,
        cmpres=cmpres,
        ilev=0,
        pdfname=pdfname,
        log_color_scale=log_color_scale,
        extra_title_txt=extra_title_txt,
        sigdiff_list=diff_aod,
        weightsdir=weightsdir,
        n_job=n_job,
        spcdb_files=spcdb_files
    )
    diff_aod[:] = [v.replace("Column_AOD_", "") for v in diff_aod]
    add_bookmarks_to_pdf(
        pdfname,
        newvarlist,
        remove_prefix="Column_AOD_",
        verbose=verbose
    )

    # ==================================================================
    # Write the list of AOD quantities having significant differences,
    # which we will need to fill out the benchmark forms.
    # ==================================================================
    if sigdiff_files is not None:
        for filename in sigdiff_files:
            if "sfc" in filename:
                with open(filename, "a+", encoding=ENCODING) as f:
                    print("* Column AOD: ", file=f, end="")
                    for v in diff_aod:
                        print(f"{v} ", file=f, end="")
                    print(file=f)
                    f.close()

    # -------------------------------------------
    # Clean up
    # -------------------------------------------
    del refds
    del devds
    gc.collect()


def make_benchmark_aerosol_tables(
        devdir,
        devlist_aero,
        devlist_spc,
        devlist_met,
        devstr,
        year,
        days_per_mon,
        spcdb_files,
        dst='./benchmark',
        overwrite=False,
        is_gchp=False,
):
    """
    Compute FullChemBenchmark aerosol budgets & burdens

    Args:
        devdir: str
            Path to development ("Dev") data directory
        devlist_aero: list of str
            List of Aerosols collection files (different months)
        devlist_spc: list of str
            List of SpeciesConc collection files (different months)
        devlist_met: list of str
            List of meteorology collection files (different months)
        devstr: str
            Descriptive string for datasets (e.g. version number)
        year: str
            The year of the benchmark simulation (e.g. '2016').
        days_per_month: list of int
            List of number of days per month for all months
        spcdb_files : list
            Paths to species_database.yml files in Ref & Dev rundirs

    Keyword Args (optional):
        dst: str
            Directory where budget tables will be created.
            Default value: './benchmark'
        overwrite: bool
            Overwrite burden & budget tables? (default=True)
            Default value: False
        is_gchp: bool
            Whether datasets are for GCHP
            Default value: False

    """
    # Create destination directory
    make_directory(dst, overwrite)

    # List of species (and subsets for the trop & strat)
    species_list = [
        "BCPI", "OCPI", "SO4", "DST1", "DST2", "DST3", "DST4",
        "DSTbin1", "DSTbin2", "DSTbin3", "DSTbin4", "DSTbin5",
        "DSTbin6", "DSTbin7", "SALA", "SALC"
    ]

    # Read the species database files in the Ref & Dev rundirs, and
    # return a dict containing metadata for the union of species.
    # We'll need properties such as mol. wt. for unit conversions, etc.
    _, dev_metadata = read_species_metadata(spcdb_files, quiet=True)

    # Get the list of relevant AOD diagnostics from a YAML file
    ifile= AOD_SPC
    aod = read_config_file(
        os.path.join(
            os.path.dirname(__file__),
            ifile,
        ),
        quiet=True
    )
    ds_aer = xr.open_mfdataset(
        devlist_aero,
        drop_variables=SKIP_THESE_VARS
    )
    ds_spc = xr.open_mfdataset(
        devlist_spc,
        drop_variables=SKIP_THESE_VARS
    )
    ds_met = xr.open_mfdataset(
        devlist_met,
        drop_variables=SKIP_THESE_VARS
    )

    # Rename SpeciesConc_ to SpeciesConcVV_ for consistency with new
    # naming introduced in GEOS-Chem 14.1.0
    ds_spc = rename_speciesconc_to_speciesconcvv(ds_spc)

    # Trim species_list to the variables that are only in the dataset
    varlist = [
        var for var in species_list if f"SpeciesConcVV_{var}" in ds_spc.data_vars
    ]

    # Molecular weights [g mol-1], as taken from the species database
    mw = {}
    full_names = {}
    for var in varlist:
        full_names[var] = dev_metadata[var]["FullName"].strip()
        mw[var] = dev_metadata[var]["MW_g"]
    mw["Air"] = MW_AIR_g

    # Get troposphere mask
    tropmask = get_troposphere_mask(ds_met)

    # Get number of months
    n_mon = len(days_per_mon)

    # ------------------------------------------------------------------
    # Surface area
    # (kludgey but it works - revisit this)
    # ------------------------------------------------------------------

    # Get number of vertical levels
    N_LEVS = ds_spc.dims["lev"]

    if is_gchp:
        area = ds_met["Met_AREAM2"].values
        a = area.shape
        area_m2 = np.zeros([a[0], N_LEVS, a[1], a[2], a[3]])
        for t in range(n_mon):
            for k in range(N_LEVS):
                area_m2[t, k, :, :, :] = area[t, :, :, :]
        total_area_m2 = np.sum(area_m2[0, 0, :, :, :])
    else:
        area = ds_met["AREA"].values
        a = area.shape
        area_m2 = np.zeros([a[0], N_LEVS, a[1], a[2]])
        for t in range(n_mon):
            for k in range(N_LEVS):
                area_m2[t, k, :, :] = area[t, :, :]
        total_area_m2 = np.sum(area_m2[0, 0, :, :])

    # ------------------------------------------------------------------
    # Conversion factors and time increments
    # ------------------------------------------------------------------
    # v/v dry --> Tg
    vv_to_Tg = {}
    for var in varlist:
        if f"SpeciesConcVV_{var}" in ds_spc.data_vars:
            vv_to_Tg[var] = \
                ds_met["Met_AD"].values * (mw[var] / mw["Air"]) * 1e-9

    # Days in the benchmark duration
    days_per_yr = np.sum(days_per_mon)

    # ------------------------------------------------------------------
    # Define function to print tables
    # ------------------------------------------------------------------
    def print_aerosol_metrics(data, varlist, namelist, filename, title, label):

        with open(filename, "w+", encoding=ENCODING) as ofile:

            # Print top header
            print("%" * 79, file=ofile)
            print(f" {title} for {year} in {devstr}", file=ofile)
            print(" (weighted by the number of days per month)", file=ofile)
            print("%" * 79, file=ofile)
            line = "\n" + " " *67 + "Strat         Trop         Strat+Trop\n"
            line += " " * 67 + "-----------   ----------   ----------"
            print(line, file=ofile)

            # Print data
            for var in varlist:
                line = f"{namelist[var] : <41} ({var : <7}) {label} :  {data[var + '_s']:11.9f}   {data[var + '_t']:10.8f}   {data[var + '_f']:10.8f}\n"
                print(line, file=ofile)

    # ------------------------------------------------------------------
    # Compute aerosol burdens [Tg] and print
    # ------------------------------------------------------------------

    # Table info
    filename = f"{dst}/Aerosol_Burdens.txt"
    if n_mon == 12:
        title = "Annual average global aerosol burdens"
    else:
        title = f"Average global aerosol burdens across {n_mon} months"
    label = "burden [Tg]"

    # Initialize
    q = {}
    q_sum_f = np.zeros(n_mon)
    q_sum_t = np.zeros(n_mon)
    q_sum_s = np.zeros(n_mon)
    burdens = {}

    # Define the axes we need to sum over to make monthly sums
    if is_gchp:
        sum_axes = (1, 2, 3, 4)
    else:
        sum_axes = (1, 2, 3)

    # Loop over species
    for var in varlist:

        # Whole-atmosphere and trop-only quantities [g]
        # NOTE: DryDep is by nature trop-only
        varname = f"SpeciesConcVV_{var}"
        q[var + "_f"] = ds_spc[varname].values * vv_to_Tg[var]
        q[var + "_t"] = np.ma.masked_array(q[var + "_f"], tropmask)

        # Compute monthly sums, weighted by the number of days per month
        q_sum_f = np.sum(q[var + "_f"], axis=sum_axes) * days_per_mon
        q_sum_t = np.sum(q[var + "_t"], axis=sum_axes) * days_per_mon
        q_sum_s = q_sum_f - q_sum_t

        # Compute annual averages
        burdens[var + "_f"] = np.sum(q_sum_f) / days_per_yr
        burdens[var + "_t"] = np.sum(q_sum_t) / days_per_yr
        burdens[var + "_s"] = np.sum(q_sum_s) / days_per_yr

    print_aerosol_metrics(burdens, varlist, full_names, filename, title, label)

    # ------------------------------------------------------------------
    # Define function to print tables
    # ------------------------------------------------------------------
    def print_aods(data, varlist, filename, title):

        with open(filename, "w+", encoding=ENCODING) as ofile:

            # Print top header
            print("%" * 79, file=ofile)
            print(f" {title} for {year} in {devstr}", file=ofile)
            print(" (weighted by the number of days per month)", file=ofile)
            print("%" * 79, file=ofile)
            line = "\n" + " " *32 + "Strat         Trop         Strat+Trop\n"
            line += " " * 32 + "-----------   ----------   ----------"
            print(line, file=ofile)

            # Print data
            for var in varlist:
                line = f"{var : <4} column optical depth [1]:  {data[var + '_s']:11.9f}   {data[var + '_t']:10.8f}   {data[var + '_f']:10.8f}\n"
                print(line, file=ofile)

    # ------------------------------------------------------------------
    # Compute average AOD's [Tg] and print
    # ------------------------------------------------------------------

    # Table info
    filename = f"{dst}/Global_Mean_AOD.txt"
    if n_mon == 12:
        title = "Annual average global AODs"
    else:
        title = f"Average global AODs across {n_mon} months"
    label = "mean AOD [1]"

    # Initialize
    q = {}
    q_sum_f = np.zeros(n_mon)
    q_sum_t = np.zeros(n_mon)
    q_sum_s = np.zeros(n_mon)
    aods = {}

    # Define axes to sum over, and total surface area
    if is_gchp:
        sum_axes = (1, 2, 3, 4)
    else:
        sum_axes = (1, 2, 3)

    # List of AOD species to plot
    varlist = [var for var in aod.keys() if "Dust" in var or "Hyg" in var]
    if is_gchp:
        varlist = [var.replace('550nm', 'WL1') for var in varlist]

    # Loop over AOD variables
    aod_list = []
    for varname in varlist:

        # Extract the s
        if "Dust" in varname:
            var = "Dust"
        elif "Total" in varname:
            continue
        else:
            var = varname.split("_")[1]
        aod_list.append(var)

        # Whole-atmosphere AOD [1]
        q[var + "_f"] = ds_aer[varname].values

        # Tropospheric-only AOD [1]
        q[var + "_t"] = np.ma.masked_array(q[var + "_f"], tropmask)

        # Create monthly sums, weighted by the number of days per month
        q_sum_f = np.sum(q[var + "_f"] * area_m2, axis=sum_axes) * days_per_mon
        q_sum_t = np.sum(q[var + "_t"] * area_m2, axis=sum_axes) * days_per_mon
        q_sum_s = q_sum_f - q_sum_t

        # Take annual averages
        aods[var + "_f"] = np.sum(q_sum_f) / total_area_m2 / days_per_yr
        aods[var + "_t"] = np.sum(q_sum_t) / total_area_m2 / days_per_yr
        aods[var + "_s"] = np.sum(q_sum_s) / total_area_m2 / days_per_yr

    print_aods(aods, aod_list, filename, title)

    # ------------------------------------------------------------------
    # Clean up
    # ------------------------------------------------------------------
    del ds_aer
    del ds_spc
    del ds_met
    gc.collect()
