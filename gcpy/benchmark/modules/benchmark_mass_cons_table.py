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
        units = TARGET_UNITS
        varname = get_passive_tracer_varname(dset)
        darr = convert_units(
            dset[varname].astype(np.float64).isel(time=t_idx),
            varname,
            metadata,
            units,
            area_m2=area.isel(time=0),
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

    return {
        "max_mass": max_mass,
        "min_mass": min_mass,
        "absdiff_g": (max_mass - min_mass) * 10**12,
        "pctdiff":  (max_mass-min_mass)/min_mass * 100,
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

        # Number of points in the time dimension
        ref_time = ref_data["time"].values
        dev_time = dev_data["time"].values

        # Throw an error if Ref & Dev have differing time values
        if not np.all(ref_time == dev_time):
            msg = "Ref and Dev have inconsistent time values!\n"
            raise ValueError(msg)

        # Lists for holding the sum of masses in Ref & Dev
        ref_masses = np.zeros(len(dev_time), dtype=np.float64)
        dev_masses = np.zeros(len(dev_time), dtype=np.float64)

        # List for holding the dates & times
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
        template = " Date & Time" + " "*18 + "Ref mass [Tg]"
        template += " "*13 + "Dev mass [Tg]"
        print(template, file=ofile)
        template = " " + "-"*17 + " "*5 + "-"*20 + " "*6 + "-"*20
        print(template, file=ofile)

        # Total masses
        for t_idx, time in enumerate(display_dates):
            template = f" {time}      "
            template +=f"{ref_masses[t_idx] : >20.13f}      "
            template +=f"{dev_masses[t_idx] : >20.13f}"
            print(template, file=ofile)
        print(" ", file=ofile)

        # Statistics
        template = " Summary" + " "*32+ "Ref" + " "*23 + "Dev"
        print(template, file=ofile)
        template = " " + "-"*17 + " "*5 + "-"*20 + " "*6 + "-"*20
        print(template, file=ofile)
        template = f" Maximum mass [Tg]     {ref_stats['max_mass'] : >20.13f}"
        template+= f"      {dev_stats['max_mass'] : >20.13f}"
        print(template, file=ofile)
        template = f" Minimum mass [Tg]     {ref_stats['min_mass'] : >20.13f}"
        template+= f"      {dev_stats['min_mass'] : >20.13f}"
        print(template, file=ofile)
        template = f" Abs diff [g]          {ref_stats['absdiff_g'] : >20.13f}"
        template+= f"      {dev_stats['absdiff_g'] : >20.13f}"
        print(template, file=ofile)
        template = f" % difference          {ref_stats['pctdiff'] : >20.13f}"
        template+= f"      {dev_stats['pctdiff'] : >20.13f}"
        print(template, file=ofile)
