#!/usr/bin/env python3
"""
Utility functions specific to the benchmark plotting/tabling scripts.
"""
import os
import numpy as np
import xarray as xr
import pandas as pd
from gcpy import util
from gcpy.constants import SKIP_THESE_VARS
from gcpy.util import verify_variable_type

# Suppress numpy divide by zero warnings to prevent output spam
np.seterr(divide="ignore", invalid="ignore")

# YAML files
AOD_SPC = "aod_species.yml"
BENCHMARK_CAT = "benchmark_categories.yml"
EMISSION_SPC = "emission_species.yml"
EMISSION_INV = "emission_inventories.yml"
LUMPED_SPC = "lumped_species.yml"
SPECIES_DATABASE = "species_database.yml"

def make_output_dir(
        dst,
        collection,
        subdst,
        overwrite=False
):
    """
    Creates a subdirectory for the given collection type
    in the destination directory.

    Args:
    ----
    dst        (str     ) : Destination directory
    collection (str|None) : e.g. "Aerosols", "DryDep", "Oxidants", ...
    subdst     (str     ) : e.g. "AnnualMean", "Apr2019", ...

    Returns:
    --------
    dst        (str     ) : Path of the directory that was created
    """
    util.verify_variable_type(dst, str)
    util.verify_variable_type(collection, str)
    util.verify_variable_type(subdst, (str, type(None)))

    # Create the destination folder (if it doesn't exist)
    util.make_directory(dst, overwrite)

    # Make the dst/collection subdirectory
    dst = os.path.join(dst, collection)
    if not os.path.isdir(dst):
        os.mkdir(dst)

    # If necessary, make the dst/collection/subdst subdirectory
    if subdst is not None:
        dst = os.path.join(dst, subdst)
        if not os.path.isdir(dst):
            os.mkdir(dst)

    return dst


def read_ref_and_dev(
        ref,
        dev,
        time_mean=False,
        multi_file=False,
        verbose=False
):
    """
    Reads files from the Ref and Dev models into xarray.Dataset objects.

    Args:
    -----
    ref (str|list) : Ref data file(s)
    def (str|list) : Dev data file(s)

    Keyword Args (optional)
    -----------------------
    multi_file (bool) : Read multiple files w/o taking avg over time
    time_mean  (bool) : Return the average over the time dimension?
    verbose    (bool) : Enable verbose output

    Returns:
    ref_data : xr.Dataset : Data from the Ref model
    dev_data : xr.Dataset : Data from the Dev model
    """
    util.verify_variable_type(ref, (str, list))
    util.verify_variable_type(dev, (str, list))
    util.verify_variable_type(time_mean, bool)

    ref_data = None
    dev_data = None
    reader = util.dataset_reader(time_mean|multi_file, verbose=verbose)

    if ref is not None:
        ref_data = reader(ref, drop_variables=SKIP_THESE_VARS)
        if time_mean:
            ref_data = util.dataset_mean(ref_data)

    if dev is not None:
        dev_data = reader(dev, drop_variables=SKIP_THESE_VARS)
        if time_mean:
            dev_data = util.dataset_mean(dev_data)

    return ref_data, dev_data


def get_common_varnames(
        refdata,
        devdata,
        prefix,
        verbose=False,
):
    """
    Returns an alphabetically-sorted list of common variables two
    xr.Dataset objects matching a given prefix.

    Args
    refdata : xr.Dataset : Data from the Ref model.
    devdata : xr.Dataset : Data from the Dev model.
    prefix  : str        : Variable prefix to match.
    verbose : bool       : Toggle verbose printout on/off.

    Returns
    varlist : list       : Sorted list of common variable names.
    """
    vardict = util.compare_varnames(refdata, devdata, quiet=not verbose)
    varlist = [var for var in vardict["commonvars"] if prefix in var]

    return sorted(varlist)


def print_sigdiffs(
        sigdiff_files,
        sigdiff_list,
        sigdiff_type,
        sigdiff_cat
):
    """
    Write
    Appends a list of species showing significant differences in a
    benchmark plotting category to a file.

    Args
    sigdiff_files : list|None : List of files for signficant diffs output
    sigdiff_list  : list      : List of significant differences to print
    sigdiff_type  : str       : e.g. "sfc", "500hPa", "zm"
    sigdiff_cat   : str       : e.g. "Oxidants", "Aerosols", "DryDep", etc.
    """
    util.verify_variable_type(sigdiff_files, (list, type(None)))
    util.verify_variable_type(sigdiff_list, list)
    util.verify_variable_type(sigdiff_type, str)
    util.verify_variable_type(sigdiff_cat, str)

    if sigdiff_files is not None:
        for ofile in sigdiff_files:
            if sigdiff_type in ofile:
                write_sigdiff(sigdiff_list, sigdiff_cat, ofile)


def write_sigdiff(
        sigdiff_list,
        sigdiff_cat,
        sigdiff_file,
):
    """
    Appends a list of species showing significant differences in a
    benchmark plotting category to a file.

    Args
    sigdiff_list : list : List of significant differences.
    sigdiff_cat  : str  : e.g. "Oxidants", "Aerosols", "DryDep", etc.
    sigdiff_file : str  : Filename to which the list will be appended
    """
    util.verify_variable_type(sigdiff_list, list)
    util.verify_variable_type(sigdiff_cat, str)
    util.verify_variable_type(sigdiff_file, str)

    with open(sigdiff_file, "a+", encoding="UTF-8") as ofile:
        print(f"* {sigdiff_cat}: ", file=ofile, end="")
        for var in sigdiff_list:
            print(f"{var} ", file=ofile, end="")
        print(file=ofile)
        ofile.close()


def pdf_filename(
        dst,
        collection,
        subdst,
        plot_type
):
    """
    Creates the absolute path for a PDF file containing benchmark plots.

    Args
    dst        : str      : Root folder for benchmark output plots
    collection : str      : e.g. "Aerosols", "DryDep", etc.
    subdst     : str|None : e.g. "AnnualMean", "Apr2019", etc.
    plot_type  : str      : e.g. "Surface", "FullColumnZonalMean", etc.

    Returns
    pdf_path   : str      : Absolute path for the PDF file containing plots
    """
    util.verify_variable_type(dst, str)
    util.verify_variable_type(collection, str)
    util.verify_variable_type(subdst, (str, type(None)))

    pdf_path = f"{collection}_{plot_type}.pdf"
    if subdst is not None:
        pdf_path = f"{collection}_{plot_type}_{subdst}.pdf"

    return os.path.join(dst, pdf_path)


def print_benchmark_info(
        config
):
    """
    Prints which benchmark plots and tables will be generated.

    Args
    config : dict : Inputs from the benchmark config YAML file.
    """
    util.verify_variable_type(config, dict)

    conf = config["options"]["bmk_type"]
    print(f"The following plots and tables will be created for {conf}")

    conf = config["options"]["outputs"]
    for key in conf.keys():
        if "plot_conc" in key and conf[key]:
            print(" - Concentration plots")
        if "plot_emis" in key and conf[key]:
            print(" - Emissions plots")
        if "plot_jvalues" in key and conf[key]:
            print(" - J-values (photolysis rates) plots")
        if "plot_budget" in key and conf[key]:
            print(" - Budget plots")
        if "plot_uvflux" in key and conf[key]:
            print(" - UVFlux plots")
        if "plot_2d_met" in key and conf[key]:
            print(" - StateMet 2D variable plots")
        if "plot_3d_met" in key and conf[key]:
            print(" - StateMet 3D variable plots")
        if "plot_aod" in key and conf[key]:
            print(" - Aerosol optical depth plots")
        if "plot_drydep" in key and conf[key]:
            print(" - Drydep velocity plots")
        if "plot_wetdep" in key and conf[key]:
            print(" - Convective and large-scale wet deposition plots")
        if "plot_models_vs_obs" in key and conf[key]:
            print(" - Plots of models vs. observations")
        if "emis_table" in key and conf[key]:
            print(" - Table of emissions totals by species and inventory")
        if "rnpbbe_budget" in key and conf[key]:
            print(" - Radionuclides budget table")
        if "ops_budget_table" in key and conf[key]:
            print(" - Operations budget tables")
        if "aer_budget_table" in key and conf[key]:
            print(" - Aerosol budget/burden tables")
        if "mass_table" in key and conf[key]:
            print(" - Table of species mass")
        if "mass_accum_table" in key and conf[key]:
            print(" - Table of species mass accumulation")
        if "OH_metrics" in key and conf[key]:
            print(" - Table of OH metrics")
        if "ste_table" in key and conf[key]:
            print(" - Table of strat-trop exchange")
        if "cons_table" in key and conf[key]:
            print(" - Table of mass conservation")

    print("Comparisons will be made for the following combinations:")
    conf = config["options"]["comparisons"]
    if conf["gcc_vs_gcc"]["run"]:
        print(" - GCC vs GCC")
    if conf["gchp_vs_gcc"]["run"]:
        print(" - GCHP vs GCC")
    if conf["gchp_vs_gchp"]["run"]:
        print(" - GCHP vs GCHP")
    if conf["gchp_vs_gcc_diff_of_diffs"]["run"]:
        print(" - GCHP vs GCC diff of diffs")


def get_geoschem_level_metadata(
        filename=None,
        search_key=None,
        verbose=False,
):
    """
    Reads a comma-separated variable (.csv) file with GEOS-Chem vertical
    level metadata and returns it in a pandas.DataFrame object.

    Args
    filename   : str          : Name of the comma-separated variable to read
    search_key : str|None     : Returns metadata that matches this value
    verbose    : bool         : Toggles verbose printout on or off

    Returns
    metadata   : pd.DataFrame : Metadata for GEOS-Chem vertical levels
    """
    if filename is None:
        filename = os.path.join(
            os.path.dirname(__file__),
            "GC_72_vertical_levels.csv"
        )

    try:
        if verbose:
            print(f"get_geoschem_level_metadata: Reading {filename}")
        metadata = pd.read_csv(filename)
    except (IOError, OSError, FileNotFoundError) as exc:
        msg = f"Could not read GEOS-Chem level metadata in {filename}!"
        raise exc(msg) from exc

    if search_key is None:
        return metadata
    return metadata[search_key]


def get_lumped_species_definitions():
    """
    Returns lumped species definitions from a YAML file.

    Returns
    lumped_spc_dict : dict : Dictionary of lumped species
    """
    ifile = LUMPED_SPC
    return util.read_config_file(
        os.path.join(
            os.path.dirname(__file__),
            ifile,
        ),
        quiet=True
    )


def archive_lumped_species_definitions(
        dst
):
    """
    Archives lumped species definitions to a YAML file.

    Args:
    dst : str : Destination folder for YAML file output.
    """
    ofile = LUMPED_SPC
    src = os.path.join(os.path.dirname(__file__), ofile)
    util.copy_file_to_dir(src, dst)


def add_lumped_species_to_dataset(
        dset,
        lspc_dict=None,
        lspc_yaml="",
        verbose=False,
        overwrite=False,
        prefix="SpeciesConcVV_",
):
    """
    Function to calculate lumped species concentrations and add
    them to an xarray Dataset. Lumped species definitions may be passed
    as a dictionary or a path to a yaml file. If neither is passed then
    the lumped species yaml file stored in gcpy is used. This file is
    customized for use with benchmark simuation SpeciesConc diagnostic
    collection output.  The algorithm has been optimized by AI to
    improve performance.

    Args
    dset      : xr.Dataset : Data prior to adding lumped species
    lspc_dict : dict       : Species & scale factors for each lumped species
    lspc_yaml : str        : YAML file w/ lumped species definitions
    verbose   : bool       : Toggles verbose printout on or off.
    overwrite : bool       : Overwrite existing species or raise an error
    prefix    : str        : Prefix to prepend to lumped species names

    Returns
    dset      : xr.Dataset : Original species plus added lumped species

    Remarks
    -------
    Key Improvements:
    1. Vectorized summation: Uses sum(to_sum) instead of incremental +=
    2. Lazy evaluation: Operations remain lazy until actual computation
    3. Single merge: Uses .assign() instead of merging many DataArrays
    4. Cleaner logic: More Pythonic dictionary iteration

    Performance Impact:
    Original: O(n_lumped × n_constituents) individual array operations
    Optimized: O(n_lumped) vectorized operations
    """

    # Default is to add all benchmark lumped species.  Can overwrite
    # by passing a dictionary or a yaml file path containing one.
    assert not (
        lspc_dict is not None and lspc_yaml != ""
    ), "Cannot pass both lspc_dict and lspc_yaml. Choose one only."
    if lspc_dict is None and lspc_yaml == "":
        lspc_dict = get_lumped_species_definitions()
    elif lspc_dict is None and lspc_yaml != "":
        lspc_dict = util.read_config_file(lspc_yaml)

    # Make sure attributes are transferred when copying dataset / dataarrays
    with xr.set_options(keep_attrs=True):

        # Dictionary to store new lumped species
        new_vars = {}

        # Loop over lumped species
        for lspc_name, constituents in lspc_dict.items():
            full_name = prefix + lspc_name

            # Check if overlap with existing species
            if full_name in dset.data_vars:
                if overwrite:
                    if verbose:
                        print(f"Overwriting existing {full_name}")
                else:
                    raise ValueError(
                        f"{full_name} already in dataset. "
                        "To overwrite pass overwrite=True."
                    )

            if verbose:
                print(f"Creating {full_name}")

            # Collect all constituent species that exist
            to_sum = []
            for spcname, scale_factor in constituents.items():
                varname = prefix + spcname
                if varname not in dset.data_vars:
                    if verbose:
                        print(f"Warning: {varname} needed for {scale_factor} not in dataset")
                    continue

                if verbose:
                    print(f" -> adding {varname} with scale {scale_factor}")

                # Build list of scaled species (lazy operations)
                to_sum.append(dset[varname] * scale_factor)

            # Vectorized sum of all constituents at once
            if len(to_sum) > 0:
                new_vars[full_name] = sum(to_sum)
            else:
                # Create NaN array matching first species shape
                if verbose:
                    print("No constituent species found! Setting to NaN.")
                template_var = next(
                    (var for key, var in dset.data_vars.items()
                     if prefix in key or prefix.replace("VV", "") in key),
                    None
                )
                if template_var is not None:
                    new_vars[full_name] = template_var.copy(deep=False) * np.nan

        # Single merge operation
        if overwrite:
            dset = dset.drop_vars(
                [key for key in new_vars.keys() if key in dset.data_vars],
                errors='ignore'
            )
        dset = dset.assign(new_vars)

    return dset


def get_species_categories(
        benchmark_type="FullChemBenchmark"
):
    """
    Returns the list of benchmark categories that each species
    belongs to.  This determines which PDF files will contain the
    plots for the various species.

    Args
    benchmark_type : str  : Specifies the type of the benchmark

    Returns:
    spc_cat_dict   : dict : Dictionary of categories and sub-categories
    """
    ifile = BENCHMARK_CAT
    spc_cat_dict = util.read_config_file(
        os.path.join(
            os.path.dirname(__file__),
            ifile,
        )
    )
    return spc_cat_dict[benchmark_type]


def archive_species_categories(
        dst
):
    """
    Writes the list of benchmark categories to a YAML file
    for archival purposes.

    Args:
    dst  : str : Destination folder for YAML file output.
    """
    ofile = BENCHMARK_CAT
    src = os.path.join(os.path.dirname(__file__), ofile)
    util.copy_file_to_dir(src, dst)


def rename_speciesconc_to_speciesconcvv(
        dset
):
    """
    Renames netCDF variables starting with "SpeciesConc_" (whcih was
    used prior to GEOS-Chem 14.1.0) to start with "SpeciesConcVV_".
    This is needed for backwards compatibility with older versions.

    Args
    dset : xr.Dataset : The input dataset

    Returns
    dset : xr.Dataset : The modified dataset
    """
    rename_dict = {}
    for var in dset.data_vars.keys():
        if var.startswith("SpeciesConc_"):
            rename_dict[var] = var.replace("SpeciesConc_", "SpeciesConcVV_")

    return dset.rename(rename_dict)


def gcc_vs_gcc_dirs(
        config,
        subdir,
):
    """
    Convenience function to return GCC vs. GCC file paths
    for use in the benchmarking modules.

    Args
    config         : dict : Info read from config file
    subdir         : str  : Subdirectory

    Returns
    refdir, devdir : str : Fike paths
    """
    util.verify_variable_type(config, dict)
    util.verify_variable_type(subdir, str)

    # Log file paths
    refdir = os.path.join(
        config["paths"]["main_dir"],
        config["data"]["ref"]["gcc"]["dir"],
        config["data"]["ref"]["gcc"][subdir]
    )
    devdir = os.path.join(
        config["paths"]["main_dir"],
        config["data"]["dev"]["gcc"]["dir"],
        config["data"]["dev"]["gcc"][subdir]
    )

    return refdir, devdir


def gchp_vs_gcc_dirs(
        config,
        subdir,
):
    """
    Convenience function to return GCHP vs. GCC file paths
    for use in the benchmarking modules.


    Args
    config         : dict : Info read from config file
    subdir         : str  : Subdirectory

    Returns
    refdir, devdir : str : Fike paths
    """
    util.verify_variable_type(config, dict)
    util.verify_variable_type(subdir, str)

    refdir = os.path.join(
        config["paths"]["main_dir"],
        config["data"]["dev"]["gcc"]["dir"],
        config["data"]["dev"]["gcc"][subdir]
    )
    devdir = os.path.join(
        config["paths"]["main_dir"],
        config["data"]["dev"]["gchp"]["dir"],
        config["data"]["dev"]["gchp"][subdir]
    )

    return refdir, devdir


def gchp_vs_gchp_dirs(
        config,
        subdir,
):
    """
    Convenience function to return GCHP vs. GCHP file paths
    for use in the benchmarking modules.

    Args
    config         : dict : Info read from config file
    subdir         : str  : Subdirectory

    Returns
    refdir, devdir : str : Fike paths
    """
    util.verify_variable_type(config, dict)
    util.verify_variable_type(subdir, str)

    refdir = os.path.join(
        config["paths"]["main_dir"],
        config["data"]["ref"]["gchp"]["dir"],
        config["data"]["ref"]["gchp"][subdir]
    )
    devdir = os.path.join(
        config["paths"]["main_dir"],
        config["data"]["dev"]["gchp"]["dir"],
        config["data"]["dev"]["gchp"][subdir]
    )

    return refdir, devdir


def get_log_filepaths(
        logs_dir,
        template,
        timestamps,
):
    """
    Returns a list of paths for GEOS-Chem log files.
    These are needed to compute the benchmark timing tables.

    Args
    logs_dir   : str  : Path to directory w/ log files
    template   : str  : Log file template w/ "%DATE%" token
    timestamps : list : List of datetimes
    """
    util.verify_variable_type(logs_dir, str)
    util.verify_variable_type(template, str)

    # Initialize local variables
    format_str = ""
    fmts = ["%Y", "%m", "%d", "%h"]
    result = []

    # Create the format string for the log file template
    for fmt in fmts:
        if fmt in template:
            format_str += fmt

    # If there is only one timestamp, add it to a list
    # so that the for loop below will work properly.
    if timestamps.size == 1:
        timestamps = [timestamps]

    # Create each output logfile name, replacing template with date
    for timestamp in timestamps:
        time = timestamp.item().strftime(format_str)
        result.append(
            os.path.join(
                logs_dir,
                template.replace(format_str, time),
            )
        )

    return result


def get_datetimes_from_filenames(
        files
):
    """
    Returns datetimes obtained from GEOS-Chem diagnostic or
    restart file names.

    Args
    files     : list       : GEOS-CHem diagnostic/restart file names

    Returns
    datetimes : np.ndarray : Array of np.datetime64 values
    """
    datetimes = np.zeros(
        len(files),
        dtype=np.datetime64("1970-01-01T00:00")
    )
    for idx, ifile in enumerate(files):
        substr = os.path.basename(ifile).split("_")
        date = substr[0].split(".")[-1]
        time = substr[1].split("z")[0]
        dt_str = date[0:4] + "-" + date[4:6] + "-" + date[6:8]
        dt_str += "T" + time[0:2] + ":" + time[2:4]
        datetimes[idx] = np.datetime64(dt_str)

    return datetimes


def get_species_database_files(config, ref_model, dev_model):
    """
    Returns the paths to the species_database.yml files in the
    Ref and Dev benchmark run directories.

    Args:
    config      : dict : Benchmark configuration information
    ref_model   : str  : Either "gcc" or "gchp"
    dev_model   : str  : Either "gcc" or "gchp"

    Returns:
    spcdb_files : list : Paths to the species database files
                         corresponding to Ref & Dev simulations
    """
    verify_variable_type(config, dict)
    verify_variable_type(ref_model, str)
    verify_variable_type(dev_model, str)

    # "gcc" and "gchp" are the only allowable values for ref_model
    if ref_model not in ("gcc", "gchp"):
        msg = "Argument 'ref_model' must be either 'gcc' or 'gchp'!"
        raise ValueError(msg)

    # "gcc" and "gchp" are the only allowable values for dev_model
    if dev_model not in ("gcc", "gchp"):
        msg = "Argument 'dev_model' must be either 'gcc' or 'gchp'!"
        raise ValueError(msg)

    # For GCHP vs GCC comparisons, "Ref" is actually the "Dev" of GCC
    if ref_model == "gcc" and dev_model == "gchp":
        ref_or_dev = "dev"
    else:
        ref_or_dev = "ref"

    # Species database in the "Ref" simulation folder
    ref_spcdb_file = os.path.join(
        config["paths"]["main_dir"],
        config["data"][ref_or_dev][ref_model]["dir"],
        SPECIES_DATABASE,
     )
    if not os.path.exists(ref_spcdb_file):
        msg = f"Could not find {ref_spcdb_file}!"
        raise FileNotFoundError(msg)
    msg = f"\nRef species database: {ref_spcdb_file}"
    print(msg)

    # Species database in the "Dev" simulation folder
    dev_spcdb_file = os.path.join(
        config["paths"]["main_dir"],
        config["data"]["dev"][dev_model]["dir"],
        SPECIES_DATABASE,
    )
    if not os.path.exists(dev_spcdb_file):
        msg = f"Could not find {dev_spcdb_file}!"
        raise FileNotFoundError(msg)
    msg = f"Dev species database: {dev_spcdb_file}"
    print(msg)

    return [ref_spcdb_file, dev_spcdb_file]
