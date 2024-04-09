"""
Creates mass conservation tables from passive tracer concentrations
stored in GEOS-Chem Classic and/or GCHP restart files.
"""
import os
import warnings
import numpy as np
import xarray as xr
from gcpy.constants import skip_these_vars
from gcpy.units import convert_units
from gcpy.util import dataset_reader, get_area_from_dataset, \
    make_directory, read_config_file, verify_variable_type


# Constants
SPC_NAME = "PassiveTracer"
TARGET_UNITS = "Tg"


def get_area(
        area_path,
        dset
):
    """
    Returns the area variable from a dataset (if present),
    or reads it from the supplied file path.

    Args
    area_path : str|None    : Full file path of area data
    dset      : xr.Dataset  : Input data

    Returns
    area      : xr.DataArray : Grid box areas [m2]
    """
    verify_variable_type(area_path, (str, type(None)))
    verify_variable_type(dset, xr.Dataset)

    # If the area variable is present in the data set, return it
    if area_path is None:
        return get_area_from_dataset(dset)

    # Otherwise read the data from the supplied area_path)
    reader = dataset_reader(multi_files=False, verbose=False)
    return get_area_from_dataset(
        reader(area_path, drop_variables=skip_these_vars).load()
    )


def get_delta_pressure(
        dset
):
    """
    Returns the delta-pressure variable from GEOS-Chem Classic
    or GCHP data files.

    Args:
    dset : xr.Dataset|xr.DataArray : Input data
    """
    verify_variable_type(dset, (xr.Dataset, xr.DataArray))

    # GEOS-Chem Classic
    if 'Met_DELPDRY' in list(dset.data_vars):
        return dset['Met_DELPDRY']

    # GCHP
    return dset['DELP_DRY']


def get_passive_tracer_metadata(
        spcdb_dir
):
    """
    Returns a dictionary with metadata for the passive tracer.

    Args
    spcdb_dir  : str  : Directory containing species_database.yml

    Returns
    properties : dict : Dictionary with species metadata
    """
    verify_variable_type(spcdb_dir, str)

    spc_name = SPC_NAME
    properties = read_config_file(
        os.path.join(
            spcdb_dir,
            "species_database.yml"
        ),
        quiet=True
    )

    return properties.get(spc_name)


def get_passive_tracer_varname(
        dset
):
    """
    Returns the variable name under which the passive tracer
    is stored GEOS-Chem Classic or GCHP restart files.

    Args
    dset    : xr.Dataset : The input data

    Returns
    varname : str        : Variable name for passive tracer
    """
    verify_variable_type(dset, xr.Dataset)

    # Name of species (it's more efficient to copy to local variable!)
    name = SPC_NAME

    # GEOS-Chem Classic
    if f"SpeciesRst_{name}" in dset.data_vars:
        return f"SpeciesRst_{name}"

    # GCHP
    return f"SPC_{name}"


def compute_total_mass(
        t_idx,
        dset,
        area,
        delta_p,
        metadata,
):
    """
    Computes the total mass (in Tg) for the passive tracer.

    Args
    t_idx      : int          : Time index
    dset       : xr.Dataset   : Data [mol/mol dry air]
    area       : xr.DataArray : Grid box areas [m2]
    delta_p    : xr.Dataset   : Pressure thicknesses [hPa]
    metadata   : dict         : Dictionary w/ species metdata

    Returns
    total_mass : np.float64   : Total mass [Tg] of species.
    """
    with xr.set_options(keep_attrs=True):

        # Local variables
        units = TARGET_UNITS
        varname = get_passive_tracer_varname(dset)

        # If area has multiple time slices, take the first one
        if "time" in area.dims:
            area = area.isel(time=0)

        # Compute mass in Tg
        darr = convert_units(
            dset[varname].astype(np.float64).isel(time=t_idx),
            varname,
            metadata,
            units,
            area_m2=area,
            delta_p=delta_p.isel(time=t_idx),
        )

        return np.sum(darr)


def compute_statistics(masses):
    """
    Returns a dictionary with statistics for total masses.

    Args
    masses     : np.ndarray : Total masses in Tg

    Returns
    statistics : dict       : Dictionary with statistics
    """
    verify_variable_type(masses, (np.ndarray, list))

    max_mass = np.max(masses)
    min_mass = np.min(masses)
    start_mass = masses[0]
    end_mass = masses[-1]

    return {
        "min_mass"           : min_mass,
        "max_mass"           : max_mass,
        "minmax_absdiff_g"   : (max_mass - min_mass) * 1.0e12,
        "minmax_pctdiff"     : (max_mass - min_mass)/min_mass * 100.0,
        "start_mass"         : start_mass,
        "end_mass"           : end_mass,
        "startend_absdiff_g" : (end_mass - start_mass) * 1.0e12,
        "startend_pctdiff"   : (end_mass - start_mass)/start_mass * 100.0,
        "mean_mass"          : np.mean(masses, dtype=np.float64),
        "variance"           : np.var(masses, dtype=np.float64),
    }


def compute_diff(
        key,
        ref,
        dev
):
    """
    Computes the difference in two dictionaries (Dev - Ref) for
    a given search key.

    key    : str   : Search key
    ref    : dict  : Dictionary of values from Ref model
    dev    : dict  : Dictionary of values from Dev model

    Returns
    diffs  : dict : Absolute & percent differences btw Dev & Ref for key
    """
    verify_variable_type(key, str)
    verify_variable_type(ref, dict)
    verify_variable_type(dev, dict)

    return {
        "absdiff": dev[key] - ref[key],
        "pctdiff": ((dev[key] - ref[key]) / ref[key]) * 100.0
    }


def compute_diff_statistics(
        ref,
        dev
):
    """
    Computes difference statistics between the Ref and Dev versions.

    Args
    ref_masses : dict : Statistics for Ref model
    dev_masses : dict : Statistics for Dev model

    Returns
    diff_stats : dict : Difference statistics between Dev and Ref
    """
    verify_variable_type(ref, dict)
    verify_variable_type(dev, dict)

    min_mass           = compute_diff("min_mass",           ref, dev)
    max_mass           = compute_diff("max_mass",           ref, dev)
    minmax_absdiff_g   = compute_diff("minmax_absdiff_g",   ref, dev)
    minmax_pctdiff     = compute_diff("minmax_pctdiff",     ref, dev)
    start_mass         = compute_diff("start_mass",         ref, dev)
    end_mass           = compute_diff("start_mass",         ref, dev)
    startend_absdiff_g = compute_diff("startend_absdiff_g", ref, dev)
    startend_pctdiff   = compute_diff("startend_pctdiff",   ref, dev)
    mean_mass          = compute_diff("mean_mass",          ref, dev)
    variance           = compute_diff("variance",           ref, dev)

    return {
        "min_mass__absdiff"           : min_mass["absdiff"],
        "min_mass__pctdiff"           : min_mass["pctdiff"],
        "max_mass__absdiff"           : max_mass["absdiff"],
        "max_mass__pctdiff"           : max_mass["pctdiff"],
        "minmax_absdiff_g__absdiff"   : minmax_absdiff_g["absdiff"],
        "minmax_absdiff_g__pctdiff"   : minmax_absdiff_g["pctdiff"],
        "minmax_pctdiff__absdiff"     : minmax_pctdiff["absdiff"],
        "minmax_pctdiff__pctdiff"     : minmax_pctdiff["pctdiff"],

        "start_mass__absdiff"         : start_mass["absdiff"],
        "start_mass__pctdiff"         : start_mass["pctdiff"],
        "end_mass__absdiff"           : end_mass["absdiff"],
        "end_mass__pctdiff"           : end_mass["pctdiff"],
        "startend_absdiff_g__absdiff" : startend_absdiff_g["absdiff"],
        "startend_absdiff_g__pctdiff" : startend_absdiff_g["pctdiff"],
        "startend_pctdiff__absdiff"   : startend_pctdiff["absdiff"],
        "startend_pctdiff__pctdiff"   : startend_pctdiff["pctdiff"],

        "mean_mass__absdiff"          : mean_mass["absdiff"],
        "mean_mass__pctdiff"          : mean_mass["pctdiff"],
        "variance__absdiff"           : variance["absdiff"],
        "variance__pctdiff"           : variance["pctdiff"],
    }


def make_benchmark_mass_conservation_table(
        ref_files,
        ref_label,
        dev_files,
        dev_label,
        dst="./benchmark",
        overwrite=False,
        ref_areapath=None,
        dev_areapath=None,
        spcdb_dir=os.path.dirname(__file__)
):
    """
    Creates a text file containing global mass of passive species
    contained in GEOS-Chem Classic and/or GCHP restart files.

    Args
    ref_files    : list|str : List of files from the Ref model
    ref_label    : str      : Ref version label
    dev_files    : list|str : List of files from the Dev model
    dev_label    : str      : Dev version label
    dst          : str      : Destination folder for file output
    overwrite    : bool     : Overwrite pre-existing files?
    ref_areapath : list|str : Path to file w/ Ref area data (optional)
    dev_areapath : list|str : Path to file w/ Dev area data (optional)
    spcdb_dir    : str      : Path to species database file
    """
    verify_variable_type(ref_files, (list, str))
    verify_variable_type(ref_label, str)
    verify_variable_type(dev_files, (list, str))
    verify_variable_type(dev_label, str)
    verify_variable_type(dst, (str, type(None)))
    verify_variable_type(overwrite, bool)
    verify_variable_type(ref_areapath, (str, type(None)))
    verify_variable_type(ref_areapath, (str, type(None)))
    verify_variable_type(spcdb_dir, str)

    # ==================================================================
    # Initialize
    # ==================================================================

    # Create the destination folder
    make_directory(dst, overwrite)

    # Get a list of properties for the given species
    metadata = get_passive_tracer_metadata(spcdb_dir)

    # Preserve xarray attributes
    with xr.set_options(keep_attrs=True):

        # ==============================================================
        # Read data and make sure time dimensions are consistent
        # ==============================================================
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=xr.SerializationWarning)

        # Pick the proper function to read the data
        reader = dataset_reader(multi_files=True, verbose=False)

        # Get data
        ref_data = reader(ref_files, drop_variables=skip_these_vars).load()
        dev_data = reader(dev_files, drop_variables=skip_these_vars).load()
        ref_area = get_area(ref_areapath, ref_data)
        dev_area = get_area(dev_areapath, dev_data)
        ref_delta_prs = get_delta_pressure(ref_data)
        dev_delta_prs = get_delta_pressure(dev_data)

        # Get datetime values
        ref_time = ref_data["time"].values
        dev_time = dev_data["time"].values

        # Throw an error if Ref & Dev have differing datetime values
        if not np.all(ref_time == dev_time):
            msg = "Ref and Dev have inconsistent time values!\n"
            raise ValueError(msg)

        # Lists for holding the sum of masses in Ref & Dev
        ref_masses = np.zeros(len(dev_time), dtype=np.float64)
        dev_masses = np.zeros(len(dev_time), dtype=np.float64)

        # List for holding the datetimes
        display_dates = []

        # ==================================================================
        # Calculate global mass for the tracer at all restart dates
        # ==================================================================
        for t_idx, time in enumerate(dev_time):

            # Save datetime string into display_dates list
            time = str(np.datetime_as_string(time, unit="m"))
            display_dates.append(time.replace("T", " "))

            # Compute total masses [Tg] for Ref & Dev
            ref_masses[t_idx] = compute_total_mass(
                t_idx,
                ref_data,
                ref_area,
                ref_delta_prs,
                metadata,
            )
            dev_masses[t_idx] = compute_total_mass(
                t_idx,
                dev_data,
                dev_area,
                dev_delta_prs,
                metadata,
            )

    # ==================================================================
    # Print masses and statistics to file
    # ==================================================================

    # Get min, max, absdiff, maxdiff for Ref & Dev
    ref_stats = compute_statistics(ref_masses)
    dev_stats = compute_statistics(dev_masses)
    diff_stats = compute_diff_statistics(ref_stats, dev_stats)

    # Create file
    outfilename = os.path.join(
        dst,
        f"Passive_mass.{ref_label}_vs_{dev_label}.txt"
    )
    with open(outfilename, 'w', encoding="utf-8") as ofile:

        # Title
        print("="*79, file=ofile)
        print("Global mass of PassiveTracer", file=ofile)
        print("", file=ofile)
        print(f"Ref = {ref_label}", file=ofile)
        print(f"Dev = {dev_label}", file=ofile)
        print("="*79, file=ofile)

        # Headers
        print("", file=ofile)
        template  = " Date & Time" + " "*18 + "Ref mass [Tg]" + " "*13
        template += "Dev mass [Tg]"+ " "*6 + "Abs Diff" + "    % Diff"
        print(template, file=ofile)
        template  = " " + "-"*17 + " "*5 + "-"*20 + " "*6 + "-"*20
        template += " " + "-"*13 + "  " + "-"*8
        print(template, file=ofile)

        # Total masses
        for t_idx, time in enumerate(display_dates):
            absdiff   = dev_masses[t_idx] - ref_masses[t_idx]
            pctdiff   = (absdiff / ref_masses[t_idx]) * 100.0
            template  = f" {time}      "
            template += f"{ref_masses[t_idx] : >20.15f}      "
            template += f"{dev_masses[t_idx] : >20.15f} "
            template += f"{absdiff : >13.4e}  "
            template += f"{pctdiff : >8.3f}"
            print(template, file=ofile)
        print(" ", file=ofile)

        # Statistics
        template  = " Summary" + " "*32+ "Ref" + " "*23 + "Dev"
        template += " "*6 + "Abs Diff" + "    % Diff"
        print(template, file=ofile)
        template  = " " + "-"*17 + " "*5 + "-"*20 + " "*6 + "-"*20
        template += " " + "-"*13 + "  " + "-"*8
        print(template, file=ofile)
        template  =  " Maximum mass [Tg]     "
        template += f"{ref_stats['max_mass'] : >20.15f}      "
        template += f"{dev_stats['max_mass'] : >20.15f} "
        template += f"{diff_stats['max_mass__absdiff'] : >13.4e}  "
        template += f"{diff_stats['max_mass__pctdiff'] : >8.3f}"
        print(template, file=ofile)
        template  =  " Minimum mass [Tg]     "
        template += f"{ref_stats['min_mass'] : >20.13f}      "
        template += f"{dev_stats['min_mass'] : >20.13f} "
        template += f"{diff_stats['min_mass__absdiff'] : >13.4e}  "
        template += f"{diff_stats['min_mass__pctdiff'] : >8.3f}"
        print(template, file=ofile)
        template  =  " Abs diff [g]          "
        template += f"{ref_stats['minmax_absdiff_g'] : >20.13f}      "
        template += f"{dev_stats['minmax_absdiff_g'] : >20.13f} "
        template += f"{diff_stats['minmax_absdiff_g__absdiff'] : >13.4e}  "
        template += f"{diff_stats['minmax_absdiff_g__pctdiff'] : >8.3f}"
        print(template, file=ofile)
        template  =  " % difference          "
        template += f"{ref_stats['minmax_pctdiff'] : >20.15f}      "
        template += f"{dev_stats['minmax_pctdiff'] : >20.15f} "
        template += f"{diff_stats['minmax_pctdiff__absdiff'] : >13.4e}  "
        template += f"{diff_stats['minmax_pctdiff__pctdiff'] : >8.3f}"
        print(template, file=ofile)
        print("", file=ofile)
        template  =  " Start mass [Tg]       "
        template += f"{ref_stats['start_mass'] : >20.15f}      "
        template += f"{dev_stats['start_mass'] : >20.15f} "
        template += f"{diff_stats['start_mass__absdiff'] : >13.4e}  "
        template += f"{diff_stats['start_mass__pctdiff'] : >8.3f}"
        print(template, file=ofile)
        template  =  " End mass [Tg]         "
        template += f"{ref_stats['end_mass'] : >20.15f}      "
        template += f"{dev_stats['end_mass'] : >20.15f} "
        template += f"{diff_stats['end_mass__absdiff'] : >13.4e}  "
        template += f"{diff_stats['end_mass__pctdiff'] : >8.3f}"
        print(template, file=ofile)
        template  =  " Abs diff [g]          "
        template += f"{ref_stats['startend_absdiff_g'] : >20.13f}      "
        template += f"{dev_stats['startend_absdiff_g'] : >20.13f} "
        template += f"{diff_stats['startend_absdiff_g__absdiff'] : >13.4e}  "
        template += f"{diff_stats['startend_absdiff_g__pctdiff'] : >8.3f}"
        print(template, file=ofile)
        template  =  " % difference          "
        template += f"{ref_stats['startend_pctdiff'] : >20.15f}      "
        template += f"{dev_stats['startend_pctdiff'] : >20.15f} "
        template += f"{diff_stats['startend_pctdiff__absdiff'] : >13.4e}  "
        template += f"{diff_stats['startend_pctdiff__pctdiff'] : >8.3f}"
        print(template, file=ofile)
        print("", file=ofile)
        template  =  " Mean mass [Tg]        "
        template += f"{ref_stats['mean_mass']:>20.15f}      "
        template += f"{dev_stats['mean_mass']:>20.15f} "
        template += f"{diff_stats['mean_mass__absdiff']:>13.4e}  "
        template += f"{diff_stats['mean_mass__pctdiff']:>8.3f}"
        print(template, file=ofile)
        template  = " Variance [Tg]         "
        template += f"{ref_stats['variance']:>20.13e}      "
        template += f"{dev_stats['variance']:>20.13e} "
        template += f"{diff_stats['variance__absdiff']:>13.4e}  "
        template += f"{diff_stats['variance__pctdiff']:>8.3f}"
        print(template, file=ofile)
