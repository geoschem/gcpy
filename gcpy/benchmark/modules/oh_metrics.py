#!/usr/bin/env python3
"""
Prints key metrics (e.g. global mean OH, MCF lifetime, and CH4 lifetimes)
for a GEOS-Chem full-chemistry simulation or methane simulation.
"""
# =====================================================================
# %%% IMPORTS ETC. %%%
# =====================================================================
import os
import warnings
import numpy as np
import xarray as xr
from gcpy.constants import \
    AVOGADRO, ENCODING, MW_AIR_kg, SKIP_THESE_VARS
from gcpy.util import \
    make_directory, read_species_metadata, replace_whitespace

# =====================================================================
# %%% METHODS %%%
# =====================================================================


def combine_dataset(file_list=None):
    """
    Wrapper for xarray.open_mfdataset

    Args:
    file_list : list       : List of files to read

    REturns
    data      : xr.Dataset : Object w/ "Metrics" collection data
    """

    # Return a single Dataset containing data from all MeanOH files.
    try:
        data = xr.open_mfdataset(
            file_list,
            drop_variables=SKIP_THESE_VARS,
        )
    except FileNotFoundError as exc:
        msg = f"Could not find one or more files in {file_list}"
        raise FileNotFoundError(msg) from exc

    return data


def validate_metrics_collection(ds):
    """
    Determines if a Dataset contains variables for computing
    metrics from a CH4 simulation or a fullchem simulation.

    Args:
        ds: xarray Dataset

    Returns:
        is_ch4_sim: bool
    """

    # CH4 and fullchem simulations have these variables
    must_have = [
        "AirMassColumnFull",
        "LossOHbyCH4columnTrop",
        "LossOHbyMCFcolumnTrop",
        "OHwgtByAirMassColumnFull",
    ]

    # Keep a count
    count = 0

    # Look for the common variables in the dataset
    for v in must_have:
        if v in ds.data_vars.keys():
            count += 1

    return count == len(must_have)


def read_metrics_collection(files):
    """
    Reads data from all "Metrics" collection netCDF files
    into a single xarray Dataset.

    Args:
    files : list       : List of "Metrics" collection netCDF files

    Returns
    data  : xr.Dataset : Object containing "Metrics" collection data
    """

    # If files a scalar, promote it to a list
    # so that we can use open_mfdataset
    if len(files) == 0:
        files = [files]

    # Combine data into a single dataset
    data = combine_dataset(files)

    # Exit if we do not have all necessary metrics variables
    if not validate_metrics_collection(data):
        msg = "Dataset does not have enough variables for computing metrics!"
        raise ValueError(msg)

    return data


def total_airmass(data):
    """
    Computes the total airmass (in both kg and molec).

    Args
    data       : xr.DataSet : Object w/ "Metrics" collection data

    Returns
    airmass_kg : np.float64 : Total atmospheric air mass [kg]
    airmass_m  : np.float64 : Total atmospheric air mass [molecules]
    """

    airmass_kg = np.nansum(data["AirMassColumnFull"].values)
    airmass_m = airmass_kg * (AVOGADRO / MW_AIR_kg)

    return airmass_kg, airmass_m


def global_mean_oh(data, airmass_kg, mw_oh_kg):
    """
    Computes the global mean OH concentration (1e5 molec cm-3)

    Args
    data        : xr.DataSet : Object w/ "Metrics" collection data
    airmass_kg  : np.float64 : Total atmospheric air mass [kg]
    mw_oh_kg    : np.flaot64 : Mol. wt. of OH [kg]

    Returns
    sum_mean_oh : np.float64 : Sum of Mean OH [1e5 molec/cm3]
    """
    # Divide out total airmass to get total mean OH concentration [kg m-3]
    # Then convert mean OH from [kg m-3] to [1e5 molec cm-3]
    mean_oh = np.nansum(data["OHwgtByAirMassColumnFull"].values)
    mean_oh = mean_oh / airmass_kg
    mean_oh *= (AVOGADRO / (mw_oh_kg * 1.0e6)) * 1.0e-5

    return mean_oh


def lifetimes_wrt_oh(data, airmass_m):
    """
    Computes the lifetimes (in years) of CH4 and CH3CCl3 (aka MCF)
    against tropospheric OH.

    Args
    data            : xr.DataSet : Object w/ "Metrics" collection data
    airmass_m       : np.float64 : Total airmass [molecules]
    s_per_yr        : np.float64 : Seconds in this year

    Returns
    ch4_life_wrt_oh : np.float64 : CH4 lifetime w/r/t OH [years]
    mcf_life_wrt_oh : np.float64 : MCF lifetime w/r/t OH [years]
    """

    # Seconds per year
    s_per_yr = np.float64(86400.0) * np.float64(365.25)

    # Loss of OH by CH4+OH and MCF+OH reactions [molec]
    oh_loss_by_ch4 = np.nansum(data["LossOHbyCH4columnTrop"].values)
    oh_loss_by_mcf = np.nansum(data["LossOHbyMCFcolumnTrop"].values)

    # CH4 and MCF lifetimes against OH [years]
    ch4_life_wrt_oh = (airmass_m / oh_loss_by_ch4) / s_per_yr
    mcf_life_wrt_oh = (airmass_m / oh_loss_by_mcf) / s_per_yr

    return ch4_life_wrt_oh, mcf_life_wrt_oh


def init_common_vars(ref, refstr, dev, devstr, spcdb_files):
    """
    Returns a dictionary containing various quantities that
    need to be passed between methods.

    Args
    ref         : str  : Path name of "Ref" (aka "Reference") data file
    refstr      : str  : Label to describe Ref
    dev         : str  : Path name of "Dev" (aka "Development") data file
    devstr      : str  : Label to describe Dev
    spcdb_files : list : Paths to Ref & Dev species_database.yml files

    Returns
    common_vars : dict : OH Metrics data
    """
    # Read the species database files in the Ref & Dev rundirs, and
    # return a dict containing metadata for the union of species.
    # We'll need properties such as mol. wt. for unit conversions, etc.
    _, spcdb = read_species_metadata(spcdb_files, quiet=True)

    # Define common_vars dictionary
    common_vars = {
        #
        # Molecular weights
        "mw_ch4_kg": spcdb["CH4"]["MW_g"] * 1.0e-3,
        "mw_oh_kg": spcdb["OH"]["MW_g"] * 1.0e-3,
        #
        # Conversion factors
        "kg_to_m_ch4": (spcdb["CH4"]["MW_g"] * 1.0e-3) / AVOGADRO,
        #
        # Ref data
        "refdata": read_metrics_collection(ref),
        "refstr": refstr,
        "mean_oh_ref": 0.0,
        "mcf_life_ref": 0.0,
        "ch4_life_ref": 0.0,
        #
        # Dev data
        "devdata": read_metrics_collection(dev),
        "devstr": devstr,
        "mean_oh_dev": 0.0,
        "mcf_life_dev": 0.0,
        "ch4_life_dev": 0.0,
    }

    return common_vars


def compute_oh_metrics(common_vars):
    """
    Computes the mass-weighted mean OH concentration, CH3CCl3 (aka MCF)
    lifetime w/r/t OH, and CH4 lifetime w/r/t OH.

    Args
    common_vars : dict : OH Metrics data

    Returns
    common_vars : dict : Updated OH Metrics data
    """

    # ==================================================================
    # Ref dataset
    # ==================================================================

    # Get total airmasses ([kg] and [molec])
    airmass_kg, airmass_m = total_airmass(
        common_vars["refdata"]
    )

    # Get mean OH [1e-5 molec cm-3]
    common_vars["mean_oh_ref"] = global_mean_oh(
        common_vars["refdata"],
        airmass_kg,
        common_vars["mw_oh_kg"]
    )

    # Get lifetimes of CH4 and MCF against tropospheric OH [years]
    common_vars["ch4_life_ref"], common_vars["mcf_life_ref"] = \
        lifetimes_wrt_oh(
            common_vars["refdata"],
            airmass_m
    )

    # ==================================================================
    # Dev dataset
    # ==================================================================

    # Get total airmasses ([kg] and [molec])
    airmass_kg, airmass_m = total_airmass(
        common_vars["devdata"]
    )

    # Get mean OH [1e-5 molec cm-3]
    common_vars["mean_oh_dev"] = global_mean_oh(
        common_vars["devdata"],
        airmass_kg,
        common_vars["mw_oh_kg"]
    )

    # Get lifetimes of CH4 and MCF against tropospheric OH [years]
    common_vars["ch4_life_dev"], common_vars["mcf_life_dev"] = \
        lifetimes_wrt_oh(
            common_vars["devdata"],
            airmass_m
    )

    return common_vars


def write_to_file(f, title, ref, dev, absdiff, pctdiff, is_mean_oh=False):
    """
    Internal routine used by print_metrics to write a specific
    quantity (mean OH, MCF lifetime, CH4 lifetime) to a file.

    Args
    f          : file       : File object
    title      : str        : Title for the data
    ref        : np.float64 : Ref data value
    dev        : np.float64 : Dev data value
    absdiff    : np.float64 : Absolute difference
    pctdiff    : np.float64 : Percent difference

    Keyword Args
    is_mean_oh : bool       : Denotes if this data is Mean OH or not
    """
    print(file=f)
    print("-" * 60, file=f)
    print(title, file=f)
    print("-" * 60, file=f)

    if is_mean_oh:
        print(f"Ref      : {ref:14.11f}", file=f)
        print(f"Dev      : {dev:14.11f}", file=f)
        print(f"Abs diff : {absdiff:14.11f}", file=f)
        print(f"%   diff : {pctdiff:9.6f}", file=f)
    else:
        print(f"Ref      : {ref:9.6f}", file=f)
        print(f"Dev      : {dev:9.6f}", file=f)
        print(f"Abs diff : {absdiff:9.6f}", file=f)
        print(f"%   diff : {pctdiff:9.6f}", file=f)


def print_metrics(common_vars, dst):
    """
    Prints the mass-weighted mean OH (full atmospheric column)
    from a GEOS-Chem simulation.

    Args
    common_vars : dict : Data containing OH Metrics data
    dst         : str  : Folder where OH Metrics output will be written
    """

    # Create file
    outfilename = os.path.join(dst, "OH_metrics.txt")

    # Open filename for output
    with open(outfilename, "w", encoding=ENCODING) as f:

        # ==============================================================
        # Write header
        # ==============================================================
        print("#" * 79, file=f)
        print("### OH Metrics", file=f)
        print(f"### Ref = {common_vars['refstr']}", file=f)
        print(f"### Dev = {common_vars['devstr']}", file=f)
        print("#" * 79, file=f)

        # ==============================================================
        # Mean OH concentration [1e5 molec/cm3]
        # ==============================================================
        title = "Global mass-weighted OH concentration [10^5 molec cm^-3]"
        absdiff = common_vars["mean_oh_dev"] - common_vars["mean_oh_ref"]
        pctdiff = (absdiff / common_vars["mean_oh_ref"]) * 100.0
        write_to_file(
            f,
            title,
            common_vars["mean_oh_ref"],
            common_vars["mean_oh_dev"],
            absdiff,
            pctdiff,
            is_mean_oh=True
        )

        # ==============================================================
        # Write MCF lifetime [years]
        # ==============================================================
        title = "CH3CCl3 (aka MCF) lifetime w/r/t tropospheric OH [years]"
        absdiff = common_vars["mcf_life_dev"] - common_vars["mcf_life_ref"]
        pctdiff = (absdiff / common_vars["mcf_life_ref"]) * 100.0
        write_to_file(
            f,
            title,
            common_vars["mcf_life_ref"],
            common_vars["mcf_life_dev"],
            absdiff,
            pctdiff
        )

        # ==============================================================
        # Write CH4 lifetime [years]
        # ==============================================================
        title = "CH4 lifetime w/r/t tropospheric OH [years]"
        absdiff = common_vars["ch4_life_dev"] - common_vars["ch4_life_ref"]
        pctdiff = (absdiff / common_vars["ch4_life_ref"]) * 100.0
        write_to_file(
            f,
            title,
            common_vars["ch4_life_ref"],
            common_vars["ch4_life_dev"],
            absdiff,
            pctdiff
        )

        # Close file
        f.close()


def make_benchmark_oh_metrics(
        ref,
        refstr,
        dev,
        devstr,
        spcdb_files,
        dst="./benchmark",
        overwrite=True,
):
    """
    Creates a text file containing metrics of global mean OH, MCF lifetime,
    and CH4 lifetime for benchmarking purposes.

    Args
    ref         : str  : Path name of "Ref" (aka "Reference") data file
    refstr      : str  : Label to describe Ref
    dev         : str  : Path name of "Dev" (aka "Development") data file
    devstr      : str  : Label to describe Dev
    spcdb_files : list : Paths to Ref & Dev species_database.yml files

    Keyword Args
    dst         : str  : Folder where OH Metrics output will be written
    overwrite   : bool : Overwrite previously-generated files? (T/F)
    """
    # Replace whitespace in the ref and dev labels
    refstr = replace_whitespace(refstr)
    devstr = replace_whitespace(devstr)

    # Tell matplotlib not to look for an X-window
    os.environ["QT_QPA_PLATFORM"] = "offscreen"

    # Suppress harmless run-time warnings (mostly about underflow in division)
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    warnings.filterwarnings("ignore", category=UserWarning)

    # Make sure that the destination directory exists
    # (or create it if it does not)
    make_directory(dst, overwrite)

    # Initialize a dictionary containing common variables
    common_vars = init_common_vars(
        ref,
        refstr,
        dev,
        devstr,
        spcdb_files
    )

    # Compute the OH metrics
    common_vars = compute_oh_metrics(
        common_vars
    )

    # Print the OH metrics
    print_metrics(
        common_vars,
        dst
    )

    # Free Dataset memory
    common_vars["refdata"] = xr.Dataset()
    common_vars["devdata"] = xr.Dataset()
