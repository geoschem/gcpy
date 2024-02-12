"""
Utility functions specific to the benchmark plotting/tabling scripts.
TODO: Migrate other benchmark-specific utilities from gcpy/benchmark.py to here.
"""
import os
import numpy as np
import xarray as xr
import pandas as pd
from gcpy import util
from gcpy.constants import skip_these_vars

# Suppress numpy divide by zero warnings to prevent output spam
np.seterr(divide="ignore", invalid="ignore")

# YAML files
AOD_SPC = "aod_species.yml"
BENCHMARK_CAT = "benchmark_categories.yml"
EMISSION_SPC = "emission_species.yml"
EMISSION_INV = "emission_inventories.yml"
LUMPED_SPC = "lumped_species.yml"


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
        ref_data = reader(ref, drop_variables=skip_these_vars)
        if time_mean:
            ref_data = util.dataset_mean(ref_data)

    if dev is not None:
        dev_data = reader(dev, drop_variables=skip_these_vars)
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
    copy = os.path.join(dst, ofile)
    if not os.path.exists(copy):
        print(f"\nArchiving {ofile} in {dst}")
        shutil.copyfile(src, copy)


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
    collection output.

    Args
    dset      : xr.Dataset : Data prior to adding lumped species
    lspc_dict : dict       : Species & scale factors for each lumped species
    lspc_yaml : str        : YAML file w/ lumped species definitions
    verbose   : bool       : Toggles verbose printout on or off.
    overwrite : bool       : Overwrite existing species or raise an error
    prefix    : str        : Prefix to prepend to lumped species names

    Returns
    dset      : xr.Dataset : Original species plus added lumped species
    """

    # Default is to add all benchmark lumped species.
    # Can overwrite by passing a dictionary
    # or a yaml file path containing one
    assert not (
        lspc_dict is not None and lspc_yaml != ""
    ), "Cannot pass both lspc_dict and lspc_yaml. Choose one only."
    if lspc_dict is None and lspc_yaml == "":
        lspc_dict = get_lumped_species_definitions()
    elif lspc_dict is None and lspc_yaml != "":
        lspc_dict = util.read_config_file(lspc_yaml)

    # Make sure attributes are transferred when copying dataset / dataarrays
    with xr.set_options(keep_attrs=True):

        # Get a dummy DataArray to use for initialization
        dummy_darr = None
        for var in dset.data_vars:
            if prefix in var or prefix.replace("VV", "") in var:
                dummy_darr = dset[var]
                dummy_type = dummy_darr.dtype
                dummy_shape = dummy_darr.shape
                break
        if dummy_darr is None:
            msg = "Invalid prefix: " + prefix
            raise ValueError(msg)

        # Create a list with a copy of the dummy DataArray object
        n_lumped_spc = len(lspc_dict)
        lumped_spc = [None] * n_lumped_spc
        for var, spcname in enumerate(lspc_dict):
            lumped_spc[var] = dummy_darr.copy(deep=False)
            lumped_spc[var].name = prefix + spcname
            lumped_spc[var].values = np.full(dummy_shape, 0.0, dtype=dummy_type)

        # Loop over lumped species list
        for var, lspc in enumerate(lumped_spc):

            # Search key for lspc_dict is lspc.name minus the prefix
            cidx = lspc.name.find("_")
            key = lspc.name[cidx+1:]

            # Check if overlap with existing species
            if lspc.name in dset.data_vars and overwrite:
                dset.drop(lspc.name)
            else:
                assert(lspc.name not in dset.data_vars), \
                    f"{lspc.name} already in dataset. To overwrite pass overwrite=True."

            # Verbose prints
            if verbose:
                print(f"Creating {lspc.name}")

            # Loop over and sum constituent species values
            num_spc = 0
            for _, spcname in enumerate(lspc_dict[key]):
                varname = prefix + spcname
                if varname not in dset.data_vars:
                    if verbose:
                        print(f"Warning: {varname} needed for {lspc_dict[key][spcname]} not in dataset")
                    continue
                if verbose:
                    print(f" -> adding {varname} with scale {lspc_dict[key][spcname]}")
                lspc.values += dset[varname].values * lspc_dict[key][spcname]
                num_spc += 1

            # Replace values with NaN if no species found in dataset
            if num_spc == 0:
                if verbose:
                    print("No constituent species found! Setting to NaN.")
                lspc.values = np.full(lspc.shape, np.nan)

        # Insert the DataSet into the list of DataArrays
        # so that we can only do the merge operation once
        lumped_spc.insert(0, dset)
        dset = xr.merge(lumped_spc)

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
    copy = os.path.join(dst, ofile)
    if not os.path.exists(copy):
        print(f"\nArchiving {ofile} in {dst}")
        shutil.copyfile(src, copy)
