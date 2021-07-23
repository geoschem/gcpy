#!/usr/bin/env python3
"""
Prints key metrics (e.g. global mean OH, MCF lifetime, and CH4 lifetimes)
for a GEOS-Chem full-chemistry simulation or methane simulation.
Requires Python3.

Calling sequence:
./oh_metrics.py
"""
# =====================================================================
# %%% IMPORTS ETC. %%%
# =====================================================================
import os
import warnings
import numpy as np
import xarray as xr
import yaml
import gcpy.constants as const

# =====================================================================
# %%% METHODS %%%
# =====================================================================


def combine_dataset(file_list=None):
    """
    Wrapper for xarray.open_mfdataset, taking into account the
    extra arguments needed in xarray 0.15 and later.

    Args:
        file_list: list of str

    Returns:
        ds: xarray Dataset
    """

    # Return a single Dataset containing data from all MeanOH files.
    # NOTE: Need to add combine="nested" and concat_dim="time"
    # for xarray 0.15 and higher!!!
    v = xr.__version__.split(".")
    if int(v[0]) == 0 and int(v[1]) >= 15:
        try:
            ds = xr.open_mfdataset(
                file_list,
                drop_variables=const.skip_these_vars,
                combine="nested",
                concat_dim="time"
            )
        except FileNotFoundError:
            msg = "Could not find one or more files in {}".format(file_list)
            raise FileNotFoundError(msg)
    else:
        try:
            ds = xr.open_mfdataset(
                file_list,
                drop_variables=const.skip_these_vars
            )
        except FileNotFoundError:
            msg = "Could not find one or more files in {}".format(file_list)
            raise FileNotFoundError(msg)

    return ds


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
        data_dir: str
            Directory containing data files.
            Default: "./OutputDir".

    Returns:
        ds: xarray Dataset
    """

    # If files a scalar, promote it to a list
    # so that we can use open_mfdataset
    if len(files) == 0:
        files = [files]

    # Combine data into a single dataset
    ds = combine_dataset(files)

    # Exit if we do not have all necessary metrics variables
    if not validate_metrics_collection(ds):
        msg = "Dataset does not have enough variables for computing metrics!"
        raise ValueError(msg)

    return ds


def total_airmass(ds):
    """
    Computes the total airmass (in both kg and molec).

    Args:
        ds: xarray Dataset

    Returns:
        airmass_kg, airmass_m: numpy float64
            Total atmospheric air mass in [kg] and [molec]
    """

    airmass_kg = np.nansum(ds["AirMassColumnFull"].values)
    airmass_m = airmass_kg * (const.AVOGADRO / const.MW_AIR_kg)

    return airmass_kg, airmass_m


def global_mean_oh(ds, airmass_kg, mw_oh_kg):
    """
    Computes the global mean OH concentration (1e5 molec cm-3)

    Args:
        sum_airmass_kg: numpy float64
        ds: xarray Dataset

    Returns:
        sum_mean_oh: numpy float64
    """
    # Divide out total airmass to get total mean OH concentration [kg m-3]
    # Then convert mean OH from [kg m-3] to [1e5 molec cm-3]
    mean_oh = np.nansum(ds["OHwgtByAirMassColumnFull"].values)
    mean_oh = (mean_oh / airmass_kg)
    mean_oh *= (const.AVOGADRO / (mw_oh_kg * 1.0e6)) * 1.0e-5

    return mean_oh


def lifetimes_wrt_oh(ds, airmass_m):
    """
    Computes the lifetimes (in years) of CH4 and CH3CCl3 (aka MCF)
    against tropospheric OH.

    Args:
        ds: xarray Dataset

        airmass_m: numpy float64
           Total airmass [molecules]

        s_per_yr: numpy float64
           Conversion factor: seconds to year.

    Returns:
        ch4_life_wrt_oh, mcf_life_wrt_oh: numpy float64
    """

    # Seconds per year
    s_per_yr = np.float64(86400.0) * np.float64(365.25)

    # Loss of OH by CH4+OH and MCF+OH reactions [molec]
    oh_loss_by_ch4 = np.nansum(ds["LossOHbyCH4columnTrop"].values)
    oh_loss_by_mcf = np.nansum(ds["LossOHbyMCFcolumnTrop"].values)

    # CH4 and MCF lifetimes against OH [years]
    ch4_life_wrt_oh = (airmass_m / oh_loss_by_ch4) / s_per_yr
    mcf_life_wrt_oh = (airmass_m / oh_loss_by_mcf) / s_per_yr

    return ch4_life_wrt_oh, mcf_life_wrt_oh


def init_common_vars(ref, refstr, dev, devstr, spcdb_dir):
    """
    Returns a dictionary containing various quantities that
    need to be passed between methods.

    Args:
        ref: str
            Path name of "Ref" (aka "Reference") data set file.

        refstr: str
            A string to describe ref (e.g. version number)

        dev: str
            Path name of "Dev" (aka "Development") data set file.
            The "Dev" data set will be compared against the "Ref" data set.

        devstr: str
            A string to describe dev (e.g. version number)

        spcdb_dir: str
            Directory of species_datbase.yml file
            Default value: Directory of GCPy code repository

    Returns:
        common_vars: dict
    """

    # Get species database
    spcdb_file = os.path.join(spcdb_dir, "species_database.yml")
    spcdb = yaml.load(open(spcdb_file), Loader=yaml.FullLoader)

    # Define common_vars dictionary
    common_vars = {
        #
        # Molecular weights
        "mw_ch4_kg": spcdb["CH4"]["MW_g"] * 1.0e-3,
        "mw_oh_kg": spcdb["OH"]["MW_g"] * 1.0e-3,
        #
        # Conversion factors
        "kg_to_m_ch4": (spcdb["CH4"]["MW_g"] * 1.0e-3) / const.AVOGADRO,
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

    Args:
        common_vars: dict

    Returns:
        common_vars: dict
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

    Args:
       f: file

       title: str

       ref, dev, absdiff, pctdiff: numpy float64

       is_mean_oh: bool
    """
    print(file=f)
    print("-" * 60, file=f)
    print(title, file=f)
    print("-" * 60, file=f)

    if is_mean_oh:
        print("Ref      : {:14.11f}".format(ref), file=f)
        print("Dev      : {:14.11f}".format(dev), file=f)
        print("Abs diff : {:14.11f}".format(absdiff), file=f)
        print("%   diff : {:9.6f}".format(pctdiff), file=f)
    else:
        print("Ref      : {:9.6f}".format(ref), file=f)
        print("Dev      : {:9.6f}".format(dev), file=f)
        print("Abs diff : {:9.6f}".format(absdiff), file=f)
        print("%   diff : {:9.6f}".format(pctdiff), file=f)


def print_metrics(common_vars, dst):
    """
    Prints the mass-weighted mean OH (full atmospheric column)
    from a GEOS-Chem simulation.

    Args:
        ds: xarray Dataset
        is_ch4_sim: bool
    """

    # Create file
    outfilename = os.path.join(dst, "OH_metrics.txt")

    # Open filename for output
    with open(outfilename, "w") as f:

        # ==============================================================
        # Write header
        # ==============================================================
        print("#" * 79, file=f)
        print("### OH Metrics", file=f)
        print("### Ref = {}; Dev = {}".format(
            common_vars["refstr"],
            common_vars["devstr"]
        ), file=f)
        print("#" * 79, file=f)
        print("\n")

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
        dst="./benchmark",
        overwrite=True,
        spcdb_dir=os.path.dirname(__file__)
):
    """
    Creates a text file containing metrics of global mean OH, MCF lifetime,
    and CH4 lifetime for benchmarking purposes.

    Args:
        ref: str
            Path name of "Ref" (aka "Reference") data set file.

        refstr: str
            A string to describe ref (e.g. version number)

        dev: str
            Path name of "Dev" (aka "Development") data set file.
            The "Dev" data set will be compared against the "Ref" data set.

        devstr: str
            A string to describe dev (e.g. version number)

    Keyword Args (optional):
        dst: str
            A string denoting the destination folder where the file
            containing emissions totals will be written.
            Default value: ./benchmark

        overwrite: bool
            Set this flag to True to overwrite files in the
            destination folder (specified by the dst argument).
            Default value: False

        spcdb_dir: str
            Directory of species_datbase.yml file
            Default value: Directory of GCPy code repository
    """
    # Tell matplotlib not to look for an X-window
    os.environ["QT_QPA_PLATFORM"] = "offscreen"

    # Suppress harmless run-time warnings (mostly about underflow in division)
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    warnings.filterwarnings("ignore", category=UserWarning)

    # Make sure that the destination directory exists
    # (or create it if it does not)
    if os.path.isdir(dst):
        if not overwrite:
            msg = "Directory {} exists. Pass overwrite=True to overwrite " \
                + "files in that directory, if any."
            msg = msg.format(dst)
            raise ValueError(msg)
    else:
        os.makedirs(dst)

    # Initialize a dictionary containing common variables
    common_vars = init_common_vars(
        ref,
        refstr,
        dev,
        devstr,
        spcdb_dir
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
