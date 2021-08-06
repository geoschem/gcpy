#!/usr/bin/env python

"""
Computes the budget of Pb, Be7, and Be10 from 1-year
TransportTracersBenchmark simulations.
"""

# ======================================================================
# Imports etc.
# ======================================================================

import os
from os.path import join
from glob import glob
import warnings
from calendar import monthrange
import numpy as np
import xarray as xr
from yaml import load as yaml_load_file
import gcpy.constants as constants
from gcpy.grid import get_troposphere_mask
from gcpy.util import rename_and_flip_gchp_rst_vars, dict_diff, reshape_MAPL_CS

# Suppress harmless run-time warnings (mostly about underflow in division)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# ======================================================================
# Define a class for passing global variables to the methods below
# ======================================================================

class _GlobVars:
    """
    Private class _GlobVars contains global data that needs to be
    shared among the methods in this module.
    """

    def __init__(self, devstr, devdir, devrstdir,
                 year, dst, is_gchp, overwrite, spcdb_dir):
        """
        Initializes the _GlobVars class.

        Args:
            devstr: str
                Label denoting the "Dev" version.
            devdir: str
                Directory where diagnostic files are found.
            devrstdir: str
                Directory where restart files are found.
            dst: str
                Directory where plots & tables will be created.
            year: int
                Year of the benchmark simulation.
            is_gchp: bool
                Denotes if this is GCHP (True) or GCC (False) data.
            overwrite: bool
                Denotes whether to ovewrite existing budget tables.
            spcdb_dir: str
                Directory where species_database.yml is stored.
        """
        # ------------------------------
        # Arguments from outside
        # ------------------------------
        self.devstr = devstr
        self.devdir = devdir
        self.devrstdir = devrstdir
        self.dst = dst
        self.is_gchp = is_gchp
        self.overwrite = overwrite

        # ------------------------------
        # Benchmark year
        # ------------------------------
        self.y0 = year
        self.y1 = self.y0 + 1
        self.y0_str = "{}".format(self.y0)
        self.y1_str = "{}".format(self.y1)

        # ------------------------------
        # Collection file lists
        # ------------------------------

        # Restarts
        if self.is_gchp:
            RstInit = join(
                self.devrstdir,
                "initial_GEOSChem_rst.c48_TransportTracers.nc"
            )
            RstFinal = join(
                self.devrstdir,
                "gcchem_internal_checkpoint.restart.{}*.nc4".format(
                    self.y1_str)
            )
        else:
            RstInit = join(
                self.devrstdir,
                "GEOSChem.Restart.{}*nc4".format(self.y0_str)
            )
            RstFinal = join(
                self.devrstdir,
                "GEOSChem.Restart.{}*.nc4".format(self.y1_str)
            )
        # Diagnostics
        HemcoDiag = join(self.devdir,
                         "HEMCO_diagnostics.{}*.nc".format(self.y0_str))
        DryDep = join(self.devdir,
                      "*.DryDep.{}*.nc4".format(self.y0_str))
        RadioNucl = join(self.devdir,
                         "*.RadioNuclide.{}*.nc4".format(self.y0_str))
        if is_gchp:
            StateMetAvg = join(self.devdir,
                               "*.StateMet_avg.{}*.nc4".format(self.y0_str))

            StateMet = join(self.devdir,
                            "*.StateMet.{}*.nc4".format(self.y0_str))

            # Set a logical if we need to read StateMet_avg or StateMet
            gchp_use_statemet_avg  = os.path.exists(StateMetAvg)

        else:
            StateMet = join(self.devdir,
                            "*.StateMet.{}*.nc4".format(self.y0_str))
        SpeciesConc = join(self.devdir,
                           "*.SpeciesConc.{}*.nc4".format(self.y0_str))
        WetLossConv = join(self.devdir,
                           "*.WetLossConv.{}*.nc4".format(self.y0_str))
        WetLossLS = join(self.devdir,
                         "*.WetLossLS.{}*.nc4".format(self.y0_str))
        GCHPEmiss = join(self.devdir,
                         "*.Emissions.{}*.nc4".format(self.y0_str))

        # ------------------------------
        # Read data collections
        # ------------------------------

        # Restarts
        skip_vars = constants.skip_these_vars
        extra_kwargs = {}

        self.ds_ini = xr.open_mfdataset(
            RstInit, drop_variables=skip_vars, **extra_kwargs)
        self.ds_end = xr.open_mfdataset(
            RstFinal, drop_variables=skip_vars, **extra_kwargs)

        # Change the restart datasets into format similar to GCC, and flip
        # vertical axis.  Also test if the restart files have the BXHEIGHT
        # variable contained within them.
        if is_gchp:
            self.ds_ini = rename_and_flip_gchp_rst_vars(self.ds_ini)
            self.ds_end = rename_and_flip_gchp_rst_vars(self.ds_end)

        # Diagnostics
        self.ds_dcy = xr.open_mfdataset(
            RadioNucl, drop_variables=skip_vars, **extra_kwargs)
        self.ds_dry = xr.open_mfdataset(
            DryDep, drop_variables=skip_vars, **extra_kwargs)
        self.ds_cnc = xr.open_mfdataset(
            SpeciesConc, drop_variables=skip_vars, **extra_kwargs)
        self.ds_wcv = xr.open_mfdataset(
            WetLossConv, drop_variables=skip_vars, **extra_kwargs)
        self.ds_wls = xr.open_mfdataset(
            WetLossLS, drop_variables=skip_vars, **extra_kwargs)

        # Met fields
        # For GCHP: Read from StateMet_avg if present (otherwise StateMet)
        if is_gchp and gchp_use_statemet_avg:
            self.ds_met = xr.open_mfdataset(
                StateMetAvg, drop_variables=skip_vars, **extra_kwargs)
        else:
            self.ds_met = xr.open_mfdataset(
                StateMet, drop_variables=skip_vars, **extra_kwargs)

        # Emissions
        if self.is_gchp:
            self.ds_hco = xr.open_mfdataset(
                GCHPEmiss, drop_variables=skip_vars, **extra_kwargs)
        else:
            self.ds_hco = xr.open_mfdataset(
                HemcoDiag, drop_variables=skip_vars, **extra_kwargs)

        # Area and troposphere mask
        # The area in m2 is on the restart file grid
        # The area in cm2 is on the History diagnostic grid
        # Both grids are identical in GCClassic but differ in GCHP
        if self.is_gchp:
            if 'Met_AREAM2' not in self.ds_met.data_vars.keys():
                msg = 'Could not find Met_AREAM2 in StateMet_avg collection!'
                raise ValueError(msg)
            area_m2 = self.ds_met["Met_AREAM2"].isel(time=0)
            area_m2 = reshape_MAPL_CS(area_m2)
            self.area_m2 = area_m2
            self.area_cm2 = self.ds_met["Met_AREAM2"] * 1.0e4
        else:
            self.area_m2 = self.ds_met["AREA"].isel(time=0)
            self.area_cm2 = self.area_m2 * 1.0e4
        self.tropmask = get_troposphere_mask(self.ds_met)

        # ------------------------------
        # Months and days
        # ------------------------------
        self.N_MONTHS = 12
        self.N_MONTHS_FLOAT = self.N_MONTHS * 1.0

        # Days per month in the benchmark year
        self.d_per_mon = np.zeros(self.N_MONTHS)
        for t in range(self.N_MONTHS):
            self.d_per_mon[t] = monthrange(self.y0, t + 1)[1] * 1.0

        # Days in the benchmark year
        self.d_per_yr = np.sum(self.d_per_mon)

        # Fraction of year occupied by each month
        self.frac_of_yr = np.zeros(self.N_MONTHS)
        for t in range(self.N_MONTHS):
            self.frac_of_yr[t] = self.d_per_mon[t] / self.d_per_yr

        # ------------------------------
        # Species info
        # ------------------------------

        # List of species (and subsets for the trop & strat)
        self.species_list = ["Pb210", "Be7", "Be10"]

        # Read the species database
        path = join(spcdb_dir, "species_database.yml")
        spcdb = yaml_load_file(open(path))

        # Molecular weights [g mol-1], as taken from the species database
        self.mw = {}
        for v in self.species_list:
            self.mw[v] = spcdb[v]["MW_g"]
        self.mw["Air"] = constants.MW_AIR_g

        # kg/s --> g/day
        self.kg_s_to_g_d_value = 86400.0 * 1000.0

        # ------------------------------
        # Conversion factors
        # ------------------------------
        self.kg_per_mol = {}
        self.vv_to_g = {}
        self.mcm2s_to_g_d = {}
        self.kg_s_to_g_d = {}
        for spc in self.species_list:

            # kg/s --> g/day (same for all species)
            self.kg_s_to_g_d[spc] = self.kg_s_to_g_d_value

            # kg/mole for each species
            self.kg_per_mol[spc] = constants.AVOGADRO / (self.mw[spc] * 1e-3)

            # v/v dry --> g
            self.vv_to_g[spc] = self.ds_met["Met_AD"].values  \
                * (self.mw[spc] / self.mw["Air"]) * 1000.0

            # molec/cm2/s --> g/day
            self.mcm2s_to_g_d[spc] = self.area_cm2.values \
                / self.kg_per_mol[spc] \
                * self.kg_s_to_g_d[spc]

# ======================================================================
# Functions
# ======================================================================

def total(globvars, dict_list):
    """
    Function to take the difference of two dict objects.
    Assumes that all objects have the same keys.

    Args:
        globvars: obj of type _GlobVars
            Global variables needed for budget computations.
        dict_list: list of dict
            Dictionaries to be summed.

    Returns:
        result: dict
            Key-by-key sum of all dicts in dict_list.
    """
    # Initialize
    result = {}
    for spc in globvars.species_list:
        result[spc + "_f"] = 0.0   # full-atmosphere
        result[spc + "_t"] = 0.0   # trop-only
        result[spc + "_s"] = 0.0   # strat-only

    # Sum over all dictionaries
    for d in dict_list:
        for k, v in d.items():
            result[k] += v

    return result


def mass_from_rst(globvars, ds, tropmask):
    """
    Computes global species mass from a restart file.

    Args:
        globvars: obj of type _GlobVars
            Global variables needed for budget computations.
        ds: xarray Dataset
            Data containing species mass to be summed.
        tropmask: numpy ndarray
            Mask to denote tropospheric grid boxes.

    Returns:
        result: dict
            Species mass in strat, trop, and strat+trop regimes.
    """
    # Initialize
    vv_to_g = {}
    rst_f = {}
    rst_t = {}
    result = {}

    # Conversion factors based on restart-file met fields
    # NOTE: Convert DataArray to numpy ndarray so that the
    # broadcasting will be done properly!!!
    g100 = 100.0 / constants.G
    if 'time' in ds["Met_DELPDRY"].dims:
        deltap = ds["Met_DELPDRY"].isel(time=0).values
    else:
        deltap = ds["Met_DELPDRY"].values
    area = globvars.area_m2.values
    airmass = deltap * area * g100

    # GCClassic: Restart variables begin with SpeciesRst
    prefix = "SpeciesRst_"

    # Determine if restart file variables have a time dimension
    # (if one does, they all do)
    varname = prefix + globvars.species_list[0]
    if 'time' in ds[varname].dims:
        has_time_dim = True
    else:
        has_time_dim = False

    # Loop over species
    for spc in globvars.species_list:

        # Add the prefix to the species name
        varname = prefix + spc

        # Conversion factor from mixing ratio to g, w/ met from rst file
        vv_to_g[spc] = airmass \
            * (globvars.mw[spc] / globvars.mw["Air"]) * 1000.0

        if has_time_dim:
            rst_f[spc] = ds[varname].isel(time=0).values * vv_to_g[spc]
        else:
            rst_f[spc] = ds[varname].values * vv_to_g[spc]

        # Troposphere-only mass
        rst_t[spc] = np.ma.masked_array(rst_f[spc], tropmask)

        # Sums
        result[spc + "_f"] = np.sum(rst_f[spc])
        result[spc + "_t"] = np.sum(rst_t[spc])
        result[spc + "_s"] = result[spc + "_f"] - result[spc + "_t"]

    return result


def annual_average(globvars, ds, collection, conv_factor):
    """
    Computes the annual average of budgets or fluxes.

    Args:
        globvars: obj of type _GlobVars
            Global variables needed for budget computations.
        ds: xarray Dataset
            Data to be averaged
        collection: str
            Name of the diagnostic collection.
        conv_factor: str
            Conversion factor to be applied.

     Returns:
        result: dict
            Annual-average budgets or fluxes in
            in strat, trop, and strat+trop regimes.
    """

    # Initialize
    q = {}
    q_sum_f = np.zeros(globvars.N_MONTHS)
    q_sum_t = np.zeros(globvars.N_MONTHS)
    q_sum_s = np.zeros(globvars.N_MONTHS)
    result = {}

    for spc in globvars.species_list:

        # Whole-atmosphere quanity [g] or [g d-1]
        varname = collection.strip() + "_" + spc
        q[spc + "_f"] = ds[varname].values * conv_factor[spc]

        # Shape of the data
        #q_shape = q[spc + "_f"].shape

        if "DryDep" in collection:
            # NOTE: DryDep is by nature trop-only
            # Therefore, data arrays don't have a lev dimension,
            # so special handling must be done.
            if globvars.is_gchp:
                sum_axes = (1, 2, 3)
            else:
                sum_axes = (1, 2)

        else:
            # Otherwise, expect a lev dimension
            if globvars.is_gchp:
                sum_axes = (1, 2, 3, 4)
            else:
                sum_axes = (1, 2, 3)

            # Trop-only quantities [g] or [g d-1]
            q[spc + "_t"] = np.ma.masked_array(q[spc + "_f"],
                                               globvars.tropmask)

        # Compute monthly averages, weighted by # of days in month
        # Special handling for DryDep, which lacks a "lev" dimension.
        if "DryDep" in collection:
            q_sum_f = np.sum(q[spc + "_f"], axis=sum_axes) \
                * globvars.d_per_mon
            q_sum_t = q_sum_f
            q_sum_s = 0.0
        else:
            q_sum_f = np.sum(q[spc + "_f"], axis=sum_axes) \
                * globvars.d_per_mon
            q_sum_t = np.sum(q[spc + "_t"], axis=sum_axes) \
                * globvars.d_per_mon
            q_sum_s = q_sum_f - q_sum_t

        # Take annual averages
        result[spc + "_f"] = np.sum(q_sum_f) / globvars.d_per_yr
        result[spc + "_t"] = np.sum(q_sum_t) / globvars.d_per_yr
        result[spc + "_s"] = np.sum(q_sum_s) / globvars.d_per_yr

    return result


def annual_average_sources(globvars):
    """
    Computes the annual average of radionuclide sources.

    Args:
        globvars: obj of type _GlobVars
            Global variables needed for budget computations.

     Returns:
        result: dict
            Source totals in strat, trop, and strat+trop regimes.
    """

    # Initialize
    q = {}
    q_sum_f = np.zeros(globvars.N_MONTHS)
    q_sum_t = np.zeros(globvars.N_MONTHS)
    q_sum_s = np.zeros(globvars.N_MONTHS)
    result = {}

    # Pb210 source is from Rn222 decay, in the RadioNuclide collection
    q["Pb210_f"] = globvars.ds_dcy["PbFromRnDecay"].values \
        * globvars.kg_s_to_g_d["Pb210"]

    # Determine the shape of the data
    q_shape = q["Pb210_f"].shape     # e.g. (time,lev,lat,lon)
    n_levs = q_shape[1]              # Number of levels

    # Determine which array axes to sum over for monthly sums
    if globvars.is_gchp:
        area_var = "Met_AREAM2"
        sum_axes = (1, 2, 3, 4)
    else:
        area_var = "AREA"
        sum_axes = (1, 2, 3)

    # Be7 and Be10 sources are in the HEMCO diagnostics collection
    q["Be7_f"] = np.zeros(q_shape)
    q["Be10_f"] = np.zeros(q_shape)

    # Convert Be7 and Be10 sources from kg/m2/s to g/day
    # NOTE: This is a kludgey way to do it but it works and
    # preserves the shape of the data as (time,lev,lat,lon).
    for t in range(globvars.N_MONTHS):
        for k in range(n_levs):
            if globvars.is_gchp:
                q["Be7_f"][t, k, :, :, :] = \
                    globvars.ds_hco["EmisBe7_Cosmic"].isel(time=t, lev=k) * \
                    globvars.ds_met[area_var].isel(time=t) * \
                    globvars.kg_s_to_g_d["Be7"]
                q["Be10_f"][t, k, :, :, :] = \
                    globvars.ds_hco["EmisBe10_Cosmic"].isel(time=t, lev=k) * \
                    globvars.ds_met[area_var].isel(time=t) * \
                    globvars.kg_s_to_g_d["Be10"]
            else:
                q["Be7_f"][t, k, :, :] = \
                    globvars.ds_hco["EmisBe7_Cosmic"].isel(time=t, lev=k) * \
                    globvars.ds_met[area_var].isel(time=t) * \
                    globvars.kg_s_to_g_d["Be7"]
                q["Be10_f"][t, k, :, :] = \
                    globvars.ds_hco["EmisBe10_Cosmic"].isel(time=t, lev=k) * \
                    globvars.ds_met[area_var].isel(time=t) * \
                    globvars.kg_s_to_g_d["Be10"]

    # Now proceed to computing the annual averages
    for spc in globvars.species_list:

        # Tropospheric-only quantities
        q[spc + "_t"] = np.ma.masked_array(q[spc + "_f"], globvars.tropmask)

        # Take monthly sums, weighted by the # of days in the month
        q_sum_f = np.sum(q[spc + "_f"], axis=sum_axes) \
            * globvars.d_per_mon
        q_sum_t = np.sum(q[spc + "_t"], axis=sum_axes) \
            * globvars.d_per_mon
        q_sum_s = q_sum_f - q_sum_t

        # Take annual averages
        result[spc + "_f"] = np.sum(q_sum_f) / globvars.d_per_yr
        result[spc + "_t"] = np.sum(q_sum_t) / globvars.d_per_yr
        result[spc + "_s"] = np.sum(q_sum_s) / globvars.d_per_yr

    return result


def trop_residence_time(globvars):
    """
    Computes the tropospheric residence time of radionuclides.

    Args:
        globvars: obj of type _GlobVars
            Global variables needed for budget computations.

    Returns:
        result: dict
            Tropopsheric residence time for all species.
    """

    # Initialize
    result = {}

    # Loop over species
    for spc in globvars.species_list:

        # Initialize
        result[spc + "_t"] = 0.0

        # Concentration [g]
        var = "SpeciesConc_" + spc
        q_cnc = globvars.ds_cnc[var].values * globvars.vv_to_g[spc]
        q_cnc = np.ma.masked_array(q_cnc, globvars.tropmask)

        # DryDep [g d-1]
        var = "DryDep_" + spc
        q_dry = globvars.ds_dry[var].values * globvars.mcm2s_to_g_d[spc]

        # Convective wet scavenging [g d-1]
        var = "WetLossConv_" + spc
        q_wcv = globvars.ds_wcv[var].values * globvars.kg_s_to_g_d[spc]
        q_wcv = np.ma.masked_array(q_wcv, globvars.tropmask)

        # Large-scale wet scavenging [g d-1]
        var = "WetLossLS_" + spc
        q_wls = globvars.ds_wls[var].values * globvars.kg_s_to_g_d[spc]
        q_wls = np.ma.masked_array(q_wls, globvars.tropmask)

        # Set a flag to denote if this is GCHP
        # and also denote which axes we should sum over
        if globvars.is_gchp:
            sum_axes = (1, 2, 3, 4)
            sum_dryd = (1, 2, 3)
        else:
            sum_axes = (1, 2, 3)
            sum_dryd = (1, 2)

        # Compute monthly averages [g]
        q_cnc_sum = np.sum(q_cnc, axis=sum_axes)
        q_dry_sum = np.sum(q_dry, axis=sum_dryd)
        q_wcv_sum = np.sum(q_wcv, axis=sum_axes)
        q_wls_sum = np.sum(q_wls, axis=sum_axes)

        # Compute lifetime
        for t in range(globvars.N_MONTHS):
            num = q_cnc_sum[t]
            denom = q_dry_sum[t] + q_wcv_sum[t] + q_wls_sum[t]
            result[spc + "_t"] += 1.0 / (num / denom)

        # Convert residence time [d]
        result[spc + "_t"] = 1.0 \
            / (result[spc + "_t"] / globvars.N_MONTHS_FLOAT)

    return result


def print_budgets(globvars, data, key):
    """
    Prints the trop+strat budget file.

    Args:
        globvars: object of type _GlobVars
            Global variables needed for budget computations.
        data: dict
            Nested dictionary containing budget info.
        key: list of str
            One of "_f", (full-atmosphere) "_t" (trop-only),
            or "_s" (strat-only).
    """
    # Create the plot directory hierarchy if it doesn't already exist
    if os.path.isdir(globvars.dst) and not globvars.overwrite:
        err_str = "Pass overwrite=True to overwrite files in that directory"
        print("Directory {} exists. {}".format(globvars.dst, err_str))
        return
    elif not os.path.isdir(globvars.dst):
        os.makedirs(globvars.dst)

    # Filename to print
    if "_f" in key:
        filename = "{}/Pb-Be_budget_trop_strat.txt".format(globvars.dst)
    elif "_t" in key:
        filename = "{}/Pb-Be_budget_troposphere.txt".format(globvars.dst)
    elif "_s" in key:
        filename = "{}/Pb-Be_budget_stratosphere.txt".format(globvars.dst)

    # Common title string
    title = "Annual Average Global Budgets of 210Pb, 7Be, and 10Be\n        "

    # Open file and print budgets
    with open(filename, "w+") as f:
        if "_f" in key:
            print(
                "Table 1. {} in the Troposphere + Stratosphere in {} for year {}\n".format(
                    title, globvars.devstr, globvars.y0_str), file=f)
        elif "_t" in key:
            print("Table 2. {} in the Troposphere in {} for year {}\n".format(
                title, globvars.devstr, globvars.y0_str), file=f)
        elif "_s" in key:
            print("Table 3. {} in the Stratosphere in {} for year {}\n".format(
                title, globvars.devstr, globvars.y0_str), file=f)
        print(
            "                                210Pb          7Be         10Be",
            file=f)
        print(
            "                          -----------  -----------  -----------",
            file=f)

        vals = [v[1] for v in list(data["burden"].items()) if key in v[0]]
        print("  Burden, g               {:11.4f}  {:11.4f}  {:11.4f}".format(
            *vals), file=f)
        print(file=f)

        if "_t" in key:
            vals = [v[1]
                    for v in list(data["res_time"].items()) if key in v[0]]

            print(
                "  Residence time, d       {:11.4f}  {:11.4f}  {:11.4f}".format(
                    *vals), file=f)
            print(file=f)

        vals = [v[1] for v in list(data["src_total"].items()) if key in v[0]]
        print("  Sources, g d-1          {:11.4f}  {:11.4f}  {:11.4f}".format(
            *vals), file=f)
        print(file=f)

        vals = [v[1] for v in list(data["snk_total"].items()) if key in v[0]]
        print("  Sinks, g d-1            {:11.4f}  {:11.4f}  {:11.4f}".format(
            *vals), file=f)

        vals = [v[1] for v in list(data["snk_drydep"].items()) if key in v[0]]
        print("    Dry deposition        {:11.4f}  {:11.4f}  {:11.4f}".format(
            *vals), file=f)
        print("    Wet deposition", file=f)

        vals = [v[1] for v in list(data["snk_wetls"].items()) if key in v[0]]
        print("      Stratiform          {:11.4f}  {:11.4f}  {:11.4f}".format(
            *vals), file=f)

        vals = [v[1] for v in list(data["snk_wetconv"].items()) if key in v[0]]
        print("      Convective          {:11.4f}  {:11.4f}  {:11.4f}".format(
            *vals), file=f)

        vals = [v[1] for v in list(data["snk_decay"].items()) if key in v[0]]
        print("    Radioactive decay     {:11.7f}  {:11.7f}  {:11.8f}".format(
            *vals), file=f)
        print(file=f)

        vals = [v[1] for v in list(data["src_minus_snk"].items())
                if key in v[0]]
        print("  Sources - Sinks, g d-1  {:11.6f}  {:11.6f}  {:11.6f}".format(
            *vals), file=f)

        print(file=f)
        print("  Accumulation Term", file=f)

        vals = [v[1] for v in list(data["accum_init"].items()) if key in v[0]]
        print("    Initial mass, g       {:11.4f}  {:11.4f}  {:11.4f}".format(
            *vals), file=f)

        vals = [v[1] for v in list(data["accum_final"].items()) if key in v[0]]
        print("    Final mass, g         {:11.4f}  {:11.4f}  {:11.4f}".format(
            *vals), file=f)

        vals = [v[1] for v in list(data["accum_diff"].items()) if key in v[0]]
        print("    Difference, g         {:11.7f}  {:11.7f}  {:11.7f}".format(
            *vals), file=f)
        f.close()


def transport_tracers_budgets(
        devstr,
        devdir,
        devrstdir,
        year,
        dst='./1yr_benchmark',
        is_gchp=False,
        overwrite=True,
        spcdb_dir=os.path.dirname(__file__)):
    """
    Main program to compute TransportTracersBenchmark budgets

    Args:
        maindir: str
            Top-level benchmark folder
        devstr: str
            Denotes the "Dev" benchmark version.
        year: int
            The year of the benchmark simulation (e.g. 2016).

    Keyword Args (optional):
        dst: str
            Directory where budget tables will be created.
            Default value: './1yr_benchmark'
        is_gchp: bool
            Denotes if data is from GCHP (True) or GCC (false).
            Default value: False
        overwrite: bool
            Denotes whether to ovewrite existing budget tables.
            Default value: True
        spcdb_dir: str
            Directory where species_database.yml is stored.
            Default value: GCPy directory
    """

    # Store global variables in a private class
    globvars = _GlobVars(devstr, devdir, devrstdir, year,
                         dst, is_gchp, overwrite, spcdb_dir)

    # Data structure for budgets
    data = {}
    # ==================================================================
    # Get the accumulation term (init & final masses) from the restarts
    # ==================================================================

    # Get initial mass from restart file
    ds = globvars.ds_ini
    tropmask = globvars.tropmask[0, :, :, :]
    data["accum_init"] = mass_from_rst(globvars, ds, tropmask)

    # Get initial mass from restart file
    ds = globvars.ds_end
    tropmask = globvars.tropmask[globvars.N_MONTHS - 1, :, :, :]
    data["accum_final"] = mass_from_rst(globvars, ds, tropmask)

    # Take the difference final - init
    data["accum_diff"] = dict_diff(data["accum_init"],
                                   data["accum_final"])

    # ==================================================================
    # Burdens [g]
    # ==================================================================
    data["burden"] = annual_average(globvars,
                                    ds=globvars.ds_cnc,
                                    collection='SpeciesConc',
                                    conv_factor=globvars.vv_to_g)

    # ==================================================================
    # Sources and sinks [g d-1]
    # ==================================================================
    data["src_total"] = annual_average_sources(globvars)

    # Radioactive decay [g d-1]
    data["snk_decay"] = annual_average(globvars,
                                       ds=globvars.ds_dcy,
                                       collection="RadDecay",
                                       conv_factor=globvars.kg_s_to_g_d)

    # Drydep fluxes
    data["snk_drydep"] = annual_average(globvars,
                                        ds=globvars.ds_dry,
                                        collection="DryDep",
                                        conv_factor=globvars.mcm2s_to_g_d)

    # Convective wet scavenging
    data["snk_wetconv"] = annual_average(globvars,
                                         ds=globvars.ds_wcv,
                                         collection="WetLossConv",
                                         conv_factor=globvars.kg_s_to_g_d)

    # Large-scale wet scavenging [g d-1]
    data["snk_wetls"] = annual_average(globvars,
                                       ds=globvars.ds_wls,
                                       collection="WetLossLS",
                                       conv_factor=globvars.kg_s_to_g_d)

    # Total sinks
    data["snk_total"] = total(globvars,
                              [data["snk_drydep"],
                               data["snk_wetls"],
                               data["snk_wetconv"],
                               data["snk_decay"]])

    # Sources - sinks
    data["src_minus_snk"] = dict_diff(data["snk_total"],
                                      data["src_total"])

    # Tropospheric residence time
    data["res_time"] = trop_residence_time(globvars)

    # ==================================================================
    # Print budgets
    # ==================================================================
    for key in ["_f", "_t", "_s"]:
        print_budgets(globvars, data, key)
