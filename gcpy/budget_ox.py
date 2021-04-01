#!/usr/bin/env python

"""
Computes the budget of Ox from 1-year GEOS-Chem Classic
or GCHP benchmark simulations.
"""

# ======================================================================
# Imports etc.
# ======================================================================

import os
import warnings
from calendar import monthrange
import numpy as np
import xarray as xr
from yaml import load as yaml_load_file
import gcpy.constants as constants
from gcpy.grid import get_troposphere_mask
from gcpy.util import rename_and_flip_gchp_rst_vars, reshape_MAPL_CS

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

    def __init__(
            self,
            devstr,
            devdir,
            devrstdir,
            year,
            dst,
            is_gchp,
            overwrite,
            spcdb_dir
    ):
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
        # --------------------------------------------------------------
        # Arguments from outside
        # --------------------------------------------------------------
        self.devstr = devstr
        self.devdir = devdir
        self.devrstdir = devrstdir
        self.dst = dst
        self.is_gchp = is_gchp
        self.overwrite = overwrite

        # ---------------------------------------------------------------
        # Benchmark year
        # ---------------------------------------------------------------
        self.y0 = int(year)
        self.y1 = self.y0 + 1
        self.y0_str = "{}".format(self.y0)
        self.y1_str = "{}".format(self.y1)

        # --------------------------------------------------------------
        # Collection file lists
        # --------------------------------------------------------------

        # Restarts
        if self.is_gchp:
            RstInit = os.path.join(
                self.devrstdir,
                "gcchem_internal_checkpoint.restart.{}{}".format(
                    self.y0_str,
                    "0101_000000.nc4"
                )
            )
            RstFinal = os.path.join(
                self.devrstdir,
                "gcchem_internal_checkpoint.restart.{}{}".format(
                    self.y1_str,
                    "0101_000000.nc4"
                )
            )
        else:
            RstInit = os.path.join(
                self.devrstdir,
                "GEOSChem.Restart.{}0101_0000z.nc4".format(self.y0_str)
            )
            RstFinal = os.path.join(
                self.devrstdir,
                "GEOSChem.Restart.{}0101_0000z.nc4".format(self.y1_str)
            )
        # Diagnostics
        DryDep = os.path.join(
            self.devdir,
            "*.DryDep.{}*.nc4".format(self.y0_str)
        )
        ProdLoss = os.path.join(
            self.devdir,
            "*.ProdLoss.{}*.nc4".format(self.y0_str)
        )
        StateMet = os.path.join(
            self.devdir,
            "*.StateMet.{}*.nc4".format(self.y0_str)
        )
        WetLossConv = os.path.join(
            self.devdir,
            "*.WetLossConv.{}*.nc4".format(self.y0_str)
        )
        WetLossLS = os.path.join(
            self.devdir,
            "*.WetLossLS.{}*.nc4".format(self.y0_str)
        )

        # --------------------------------------------------------------
        # Read data collections
        # --------------------------------------------------------------

        # Restarts
        skip_vars = constants.skip_these_vars
        extra_kwargs = {}

        # Restart files at start & end of benchmark
        self.ds_ini = xr.open_dataset(
            RstInit,
            drop_variables=skip_vars,
            **extra_kwargs
        )
        self.ds_end = xr.open_dataset(
            RstFinal,
            drop_variables=skip_vars,
            **extra_kwargs
        )

        # Change the restart datasets into format similar to GCC, and flip
        # vertical axis.  Also test if the restart files have the BXHEIGHT
        # variable contained within them.
        if is_gchp:
            self.ds_ini = rename_and_flip_gchp_rst_vars(self.ds_ini)
            self.ds_end = rename_and_flip_gchp_rst_vars(self.ds_end)

        # Diagnostics
        self.ds_pl = xr.open_mfdataset(
            ProdLoss,
            drop_variables=skip_vars,
            combine="nested",
            concat_dim="time",
            **extra_kwargs
        )
        self.ds_dry = xr.open_mfdataset(
            DryDep,
            drop_variables=skip_vars,
            combine="nested",
            concat_dim="time",
            **extra_kwargs
        )
        self.ds_wcv = xr.open_mfdataset(
            WetLossConv,
            drop_variables=skip_vars,
            combine="nested",
            concat_dim="time",
            **extra_kwargs
        )
        self.ds_wls = xr.open_mfdataset(
            WetLossLS,
            drop_variables=skip_vars,
            combine="nested",
            concat_dim="time",
            **extra_kwargs
        )

        # Met fields
        self.ds_met = xr.open_mfdataset(
            StateMet,
            drop_variables=skip_vars,
            combine="nested",
            concat_dim="time",
            **extra_kwargs
        )

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
        self.tropmask = get_troposphere_mask(
            self.ds_met,
        )

        # --------------------------------------------------------------
        # Months and days
        # --------------------------------------------------------------
        self.N_MONTHS = 12
        self.N_MONTHS_FLOAT = self.N_MONTHS * 1.0

        # Days per month in the benchmark year
        self.d_per_m = np.zeros(self.N_MONTHS)
        self.s_per_m = np.zeros(self.N_MONTHS)
        for t in range(self.N_MONTHS):
            self.d_per_m[t] = monthrange(self.y0, t + 1)[1] * 1.0
            self.s_per_m[t] = self.d_per_m[t] * 86400.0

        # Days and seconds in the benchmark year
        self.d_per_a = np.sum(self.d_per_m)
        self.s_per_a = np.sum(self.s_per_m)

        # Fraction of year occupied by each month
        self.frac_of_a = np.zeros(self.N_MONTHS)
        for t in range(self.N_MONTHS):
            self.frac_of_a[t] = self.d_per_m[t] / self.d_per_a

        # --------------------------------------------------------------
        # Species info
        # --------------------------------------------------------------

        # List of species (and subsets for the trop & strat)
        self.species_list = ["HNO3", "O3", "NO2", "NO3", "PAN", "PPN"]

        # Read the species database
        path = os.path.join(spcdb_dir, "species_database.yml")
        spcdb = yaml_load_file(open(path))

        # Molecular weights [kg mol-1], as taken from the species database
        self.mw = {}
        for v in self.species_list:
            self.mw[v] = spcdb[v]["MW_g"] * 1.0e-3
        self.mw["Air"] = constants.MW_AIR_g * 1.0e-3

        # kg/s --> Tg/d
        self.kg_s_to_tg_a_value = 86400.0 * self.d_per_a * 1e-9

        # Box volume [cm3]
        self.vol_cm3 = self.ds_met["Met_AIRVOL"].values * 1.0e6


# ======================================================================
# Functions
# ======================================================================


def init_and_final_mass(
        globvars
):
    """
    Computes global species mass from the initial & final restart files.

    Args:
        globvars: obj of type _GlobVars
            Global variables needed for budget computations.

    Returns:
        result: dict
            Contains initial & final tropospheric mass of O3.
    """

    # Meterological quantities (from restart file met)
    if 'time' in globvars.ds_ini["Met_DELPDRY"].dims:
        deltap_ini = globvars.ds_ini["Met_DELPDRY"].isel(time=0).values
    else:
        deltap_ini = globvars.ds_ini["Met_DELPDRY"].values
    if 'time' in globvars.ds_end["Met_DELPDRY"].dims:
        deltap_end = globvars.ds_end["Met_DELPDRY"].isel(time=0).values
    else:
        deltap_end = globvars.ds_end["Met_DELPDRY"].values
    g100 = 100.0 / constants.G
    airmass_ini = (deltap_ini * globvars.area_m2.values) * g100
    airmass_end = (deltap_end * globvars.area_m2.values) * g100

    # Conversion factors
    mw_ratio = globvars.mw["O3"] / globvars.mw["Air"]
    kg_to_tg = 1.0e-9
    vv_to_tg_ini = airmass_ini * (mw_ratio * kg_to_tg)
    vv_to_tg_end = airmass_end * (mw_ratio * kg_to_tg)

    # Determine if restart file variables have a time dimension
    # (if one does, they all do)
    v = "SpeciesRst_O3"
    if 'time' in globvars.ds_ini[v].dims:
        mass_ini = globvars.ds_ini[v].isel(time=0).values * vv_to_tg_ini
    else:
        mass_ini = globvars.ds_ini[v].values * vv_to_tg_ini
    if 'time' in globvars.ds_end[v].dims:
        mass_end = globvars.ds_end[v].isel(time=0).values * vv_to_tg_end
    else:
        mass_end = globvars.ds_end[v].values * vv_to_tg_end

    # Compute tropospheric masses [Tg Ox]
    tropmass_ini = np.ma.masked_array(
        mass_ini,
        globvars.tropmask[0, :, :, :]
    )
    tropmass_end = np.ma.masked_array(
        mass_end,
        globvars.tropmask[globvars.N_MONTHS-1, :, :, :]
    )

    # Create a dict to return values
    result = {}
    result["Ini"] = np.sum(tropmass_ini)
    result["End"] = np.sum(tropmass_end)
    result["Acc"] = result["End"] - result["Ini"]

    return result


def annual_average_prodloss(
        globvars
):
    """
    Args:
        globvars: obj of type _GlobVars
            Global variables needed for budget computations.

    Returns:
        result: dict
            Contains annual monthly-weighted average of Prod_Ox
            and Loss_Ox.
    """

    # Conversion factors
    mw_avo = globvars.mw["O3"] / constants.AVOGADRO
    kg_to_tg = 1.0e-9

    # Tropospheric P(Ox) and L(Ox) [molec/s]
    prod_trop = np.ma.masked_array(
        globvars.ds_pl["Prod_Ox"].values * globvars.vol_cm3,
        globvars.tropmask
    )
    loss_trop = np.ma.masked_array(
        globvars.ds_pl["Loss_Ox"].values * globvars.vol_cm3,
        globvars.tropmask
    )

    # Compute monthly-weighted averages [Tg Ox]
    prod_tot = 0.0
    loss_tot = 0.0
    for t in range(globvars.N_MONTHS):
        prod_tot += np.sum(prod_trop[t, :, :, :]) * globvars.frac_of_a[t]
        loss_tot += np.sum(loss_trop[t, :, :, :]) * globvars.frac_of_a[t]
    prod_tot *= globvars.s_per_a * kg_to_tg * mw_avo
    loss_tot *= globvars.s_per_a * kg_to_tg * mw_avo

    # Make a dict for returning results
    result = {}
    result["POx"] = np.sum(prod_tot)
    result["LOx"] = np.sum(loss_tot)
    result["POx-LOx"] = result["POx"] - result["LOx"]

    return result


def annual_average_drydep(
        globvars
):
    """
    Args:
        globvars: obj of type _GlobVars
            Global variables needed for budget computations.

    Returns:
        result : float
            Annual average dry deposition [Tg Ox]
    """

    # Conversion factors and area
    mw_avo = (globvars.mw["O3"] / constants.AVOGADRO)
    kg_to_tg = 1.0e-9
    area_cm2 = globvars.area_cm2.values

    # Get P(Ox) AND L(Ox) [kg/s]
    dry = globvars.ds_dry["DryDep_HNO3"].values \
        + globvars.ds_dry["DryDep_NO2"].values  \
        + globvars.ds_dry["DryDep_O3"].values   \
        + globvars.ds_dry["DryDep_PAN"].values  \
        + globvars.ds_dry["DryDep_PPN"].values

    # Monthly-weighted conv & LS wet losses [kg HNO3/s]
    dry_tot = 0.0
    for t in range(globvars.N_MONTHS):
        dry_tot += np.sum(dry[t, :, :] * area_cm2) * globvars.frac_of_a[t]
    result = dry_tot * globvars.s_per_a * kg_to_tg * mw_avo

    return result


def annual_average_wetdep(globvars):
    """
    Args:
        globvars: obj of type _GlobVars
            Global variables needed for budget computations.

    Returns:
        result : float
            Annual average wet deposition [Tg Ox]
    """

    # Conversion factors
    hno3_to_ox = globvars.mw["O3"] / globvars.mw["HNO3"]
    kg_to_tg = 1.0e-9

    # Tropospheric wet loss (convective & large scale)
    # Take only tropospheric values
    wetcv_trop = np.ma.masked_array(
        globvars.ds_wcv["WetLossConv_HNO3"].values,
        globvars.tropmask
    )
    wetls_trop = np.ma.masked_array(
        globvars.ds_wls["WetLossLS_HNO3"].values,
        globvars.tropmask
    )

    # Monthly-weighted conv & LS wet losses [kg HNO3/s]
    wetcv_tot = 0.0
    wetls_tot = 0.0
    for t in range(globvars.N_MONTHS):
        wetcv_tot += np.sum(wetcv_trop[t, :, :, :]) * globvars.frac_of_a[t]
        wetls_tot += np.sum(wetls_trop[t, :, :, :]) * globvars.frac_of_a[t]

    # Convert [kg HNO3/s] to [Tg Ox/a]
    wetcv_tot *= globvars.s_per_a * kg_to_tg  * hno3_to_ox
    wetls_tot *= globvars.s_per_a * kg_to_tg  * hno3_to_ox

    # Create a dict to return the results
    result = {}
    result["CV"] = wetcv_tot
    result["LS"] = wetls_tot
    result["Total"] = result["CV"] + result["LS"]

    return result


def print_budget(
        globvars,
        mass,
        prodloss,
        wetdep,
        drydep,
        dyn,
        net
):
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
    if not os.path.isdir(globvars.dst):
        os.makedirs(globvars.dst)

    # Filename
    filename = "{}/{}.O3_budget_troposphere_{}.txt".format(
        globvars.dst,
        globvars.devstr,
        globvars.y0_str
    )

    # Open file and print budgets
    with open(filename, "w+") as f:
        print("="*50, file=f)
        print("Annual Average Global Ox Budget", file=f)
        print("for GEOS-Chem {} 1-year benchmark\n".format(
            globvars.devstr), file=f)
        print("Start: {}-01-01 00:00 UTC".format(globvars.y0_str), file=f)
        print("End:   {}-01-01 00:00 UTC".format(globvars.y1_str), file=f)
        print("="*50, file=f)
        print("\n", file=f)
        print("  MASS ACCUMULATION       Tg Ox a-1", file=f)
        print("  -----------------      ----------\n", file=f)
        v = mass["Ini"]
        print("    Initial mass        {:11.4f}".format(v), file=f)
        v = mass["End"]
        print("    Final mass          {:11.4f}".format(v), file=f)
        print("                         ----------", file=f)
        v = mass["Acc"]
        print("    Difference          {:11.7f}".format(v), file=f)
        print("\n", file=f)
        print("  SOURCES AND SINKS       Tg Ox a-1", file=f)
        print("  -----------------      ----------\n", file=f)

        print("  * Chemistry", file=f)
        v = prodloss["POx"]
        print("      Total POx         {:11.4f}".format(v), file=f)
        v = prodloss["LOx"]
        print("      Total LOx         {:11.4f}".format(v), file=f)
        print("                         ----------", file=f)
        v = prodloss["POx-LOx"]
        print("      Net POx - LOx     {:11.6f}\n".format(v), file=f)
        v = drydep
        print("  * Dry deposition      {:11.4f}\n".format(v), file=f)
        print("  * Wet deposition", file=f)
        v = wetdep["CV"]
        print("      Convective        {:11.4f}".format(v), file=f)
        v = wetdep["LS"]
        print("      Large-Scale       {:11.4f}".format(v), file=f)
        print("                         ----------", file=f)
        v = wetdep["Total"]
        print("      Total Wetdep      {:11.6f}\n".format(v), file=f)
        v = dyn
        print("  * Dynamics", file=f)
        print("      Strat Flux(Ox)    {:11.4f}\n".format(v), file=f)
        v = net
        print("  * DELTA Ox", file=f)
        print("      Chem+Dyn-Dry-Wet  {:11.6f}".format(v), file=f)
        print(file=f)
        f.close()


def global_ox_budget(
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
    globvars = _GlobVars(
        devstr,
        devdir,
        devrstdir,
        year,
        dst,
        is_gchp,
        overwrite,
        spcdb_dir
    )

    # ==================================================================
    # Compute Ox budget [Tg a-1]
    # ==================================================================

    # Mass from initial & final restart file
    mass = init_and_final_mass(
        globvars
    )

    # Sources and sinks
    prodloss = annual_average_prodloss(
        globvars
    )
    wetdep = annual_average_wetdep(
        globvars
    )
    drydep = annual_average_drydep(
        globvars
    )

    # Stratospheric Ox is residual source
    dyn = mass["Acc"] - (prodloss["POx-LOx"] - wetdep["Total"] - drydep)

    # Net Ox: (Chem + Dyn) - (Wet + Dry)
    net = (prodloss["POx-LOx"] + dyn) - (wetdep["Total"] + drydep)

    # ==================================================================
    # Print budgets to file
    # ==================================================================
    print_budget(
        globvars,
        mass,
        prodloss,
        wetdep,
        drydep,
        dyn,
        net
    )

# For testing
if "__main__" in __name__:
    devstr = "12.9.0"
    devdir = "/n/holyscratch01/jacob_lab/ryantosca/B1yr/12.9.0/OutputDir"
    devrstdir = devdir
    year = 2016
    dst = "/n/holyscratch01/jacob_lab/ryantosca/B1yr/Benchmark_Results"
    ox_budget(devstr, devdir, devrstdir, year, dst)
