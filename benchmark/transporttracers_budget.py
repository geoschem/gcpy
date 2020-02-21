#!/usr/bin/env python

"""
Computes the budget of Pb and Be7 from the TransportTracers benchmarks.

NOTE: This works for GC-Classic, but may need modifications for GCHP.
 -- Bob Yantosca (21 Jan 2020)
"""

# ======================================================================
# Imports etc.
# ======================================================================

from calendar import monthrange
import numpy as np
import os
from os.path import join
import gcpy.constants as constants
from gcpy.benchmark import get_troposphere_mask
import warnings
import xarray as xr

# Suppress harmless run-time warnings (mostly about underflow in division)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# ======================================================================
# GLOBAL VARIABLES: Configurables (MUST EDIT)
# ======================================================================

# Main data directory
maindir     = "/path/to/main/data/dir"

# Version string
devstr      = "dev_version_str"

# Directory where budget tables will be created
plotsdir    = join(maindir, devstr, "Plots")

# Restart collections
rstdir      = join(maindir, devstr, "restarts",              )
RstInit     = join(rstdir,  "GEOSChem.Restart.2016*.nc4"     )
RstFinal    = join(rstdir,  "GEOSChem.Restart.2017*.nc4"     )

# Data collections
datadir     = join(maindir, devstr, "OutputDir"              )
HemcoDiag   = join(datadir, "HEMCO_diagnostics.2016*.nc"     )
DryDep      = join(datadir, "GEOSChem.DryDep.2016*.nc4"      )
RadioNucl   = join(datadir, "GEOSChem.RadioNuclide.2016*.nc4")
StateMet    = join(datadir, "GEOSChem.StateMet.2016*.nc4"    )
SpeciesConc = join(datadir, "GEOSChem.SpeciesConc.2016*.nc4" )
WetLossConv = join(datadir, "GEOSChem.WetLossConv.2016*.nc4" )
WetLossLS   = join(datadir, "GEOSChem.WetLossLS.2016*.nc4"   )

# ======================================================================
# GLOBAL VARIABLES: Months and days
# ======================================================================

# Number of months (nominally 12)
N_MONTHS = 12
N_MONTHS_FLOAT = N_MONTHS * 1.0

# Days per month in 2016
d_per_mon = np.zeros(N_MONTHS)
for t in range(N_MONTHS):
    d_per_mon[t] = monthrange(2016,t+1)[1] * 1.0

# Days in the year 2016
d_per_yr = np.sum(d_per_mon)

# Fraction of year occupied by each month
frac_of_yr = np.zeros(N_MONTHS)
for t in range(N_MONTHS):
    frac_of_yr[t] = d_per_mon[t] / d_per_yr
    
# ======================================================================
# GLOBAL VARIABLES: Met fields, constants, conversion factors
# ======================================================================

# Read the StateMet collection
ds = xr.open_mfdataset(StateMet)

# Read certain met data variables
area_m2 = ds["AREA"].isel(time=0)
area_cm2 = area_m2 * 1.0e4
tropmask = get_troposphere_mask(ds)

# List of species (and subsets for the trop & strat)
species_list = ["Pb210", "Be7", "Be10" ]
full_atm = [v + "_f" for v in species_list]
trop_only = [v + "_t" for v in species_list]
strat_only = [v + "_s" for v in species_list]

# Molecular weights
mw = { "Pb210": 210.0, "Be7": 7.0, "Be10": 10.0, "Air": 28.9644}

# kg/s --> g/day
kg_s_to_g_d_value= 86400.0 * 1000.0

# Conversion factors
kg_per_mol = {}
vv_to_g = {}
mcm2s_to_g_d  = {} 
kg_s_to_g_d = {}

for spc in species_list:

    # kg/s --> g/day (same for all species)
    kg_s_to_g_d[spc] = kg_s_to_g_d_value
    
    # kg/mole for each species
    kg_per_mol[spc] = constants.AVOGADRO / (mw[spc]  * 1e-3)

    # v/v dry --> g
    vv_to_g[spc] = ds["Met_AD"].values * (mw[spc] / mw["Air"]) * 1000.0

    # molec/cm2/s --> g/day
    mcm2s_to_g_d[spc] = area_cm2.values / kg_per_mol[spc] * kg_s_to_g_d[spc]

# ======================================================================
# Functions
# ======================================================================

def diff(dict0, dict1):
    """
    Function to take the difference of two dict objects.
    Assumes that both objects have the same keys.
    """
    result = {}
    for key, value in dict0.items():
        result[key] = dict1[key] - dict0[key]

    return result


def total(dict_list):
    """
    Function to take the difference of two dict objects.
    Assumes that all objects have the same keys.
    """
    # Initialize
    result = {}
    for spc in species_list:
        result[spc + "_f"] = 0.0   # full-atmosphere
        result[spc + "_t"] = 0.0   # trop-only
        result[spc + "_s"] = 0.0   # strat-only

    # Sum over all dictionaries
    for d in dict_list:
        for k, v in d.items():
            result[k] += v

    return result


def mass_from_rst(ds, tropmask):
    """
    Computes global species mass from a restart file.
    """
    # Initialize
    vv_to_g = {}
    rst_f = {}
    rst_t = {}
    result = {}
    
    # Conversion factors based on restart-file met fields
    g100    = 100.0 / constants.G
    airmass = ds["Met_DELPDRY"].isel(time=0) * ds["AREA"] * g100
    airmass = airmass.values

    # Loop over species
    for spc in species_list:

        # Conversion factor from mixing ratio to g, w/ met from rst file
        vv_to_g[spc] = airmass * (mw[spc] / mw["Air"]) * 1000.0

        # Whole-atmosphere mass
        rst_f[spc] = ds["SpeciesRst_" + spc].isel(time=0).values * vv_to_g[spc]

        # Troposphere-only mass
        rst_t[spc] = np.ma.masked_array(rst_f[spc], tropmask)

        # Sums
        result[spc + "_f"] = np.sum(rst_f[spc])
        result[spc + "_t"] = np.sum(rst_t[spc])
        result[spc + "_s"] = result[spc + "_f"] - result[spc + "_t"]

    return result


def annual_average(ds, collection, tropmask, conv_factor):
    """Take the annual average of a quantity"""

    # Initialize
    q = {}
    q_sum_f = np.zeros(N_MONTHS)
    q_sum_t = np.zeros(N_MONTHS)
    q_sum_s = np.zeros(N_MONTHS)
    result = {}
    
    for spc in species_list:
        
        # Whole-atmosphere and trop-only quantities [g]
        # NOTE: DryDep is by nature trop-only
        varname = collection.strip() + "_" + spc
        q[spc + "_f"] = ds[varname].values * conv_factor[spc]
        if "DryDep" not in collection:
            q[spc + "_t"] = np.ma.masked_array(q[spc + "_f"], tropmask)

        # Compute monthly averages, weighted by # of days in month
        # Special handling for Drydep, which is a 3-D array
        for t in range(N_MONTHS):
            if "DryDep" in collection:
                q_sum_f[t] = np.sum(q[spc + "_f"][t,:,:]) * d_per_mon[t]
                q_sum_t[t] = q_sum_f[t]
                q_sum_s[t] = 0.0
            else:
                q_sum_f[t] = np.sum(q[spc + "_f"][t,:,:,:]) * d_per_mon[t]
                q_sum_t[t] = np.sum(q[spc + "_t"][t,:,:,:]) * d_per_mon[t] 
                q_sum_s[t] = q_sum_f[t] - q_sum_t[t]

        # Take annual averages
        result[spc + "_f"] = np.sum(q_sum_f) / d_per_yr
        result[spc + "_t"] = np.sum(q_sum_t) / d_per_yr
        result[spc + "_s"] = np.sum(q_sum_s) / d_per_yr

    return result


def annual_average_sources(ds_hco, ds_dcy, tropmask):
    """
    Take the annual average of a sources -- special handling
    because Pb210 source is in one collection and Be7/Be10 sources
    are in another collection.
    """
    
    # Initialize
    q = {}
    q_sum_f  = np.zeros(N_MONTHS)
    q_sum_t = np.zeros(N_MONTHS)
    q_sum_s = np.zeros(N_MONTHS)
    result = {}
    
    # Pb210 source is from Rn222 decay, in the RadioNuclide collection
    q["Pb210_f"] = ds_dcy["PbFromRnDecay"].values * kg_s_to_g_d["Pb210"]
    q_shape = q["Pb210_f"].shape

    # Be7 and Be10 sources are in the HEMCO diagnostics collection
    q["Be7_f"] = np.zeros(q_shape)
    q["Be10_f"] = np.zeros(q_shape)

    # Convert Be7 and Be10 sources from kg/m2/s to g/day
    # NOTE: This is a kludgey way to do it but it works and
    # preserves the shape of the data as (time,lev,lat,lon).
    for t in range(N_MONTHS):
        for k in range(q_shape[1]):
            q["Be7_f"][t,k,:,:]  = \
                ds_hco["EmisBe7_Cosmic"].isel(time=t, lev=k) * \
                ds_hco["AREA"].isel(time=t) * \
                kg_s_to_g_d["Be7"]
            q["Be10_f"][t,k,:,:] = \
                ds_hco["EmisBe10_Cosmic"].isel(time=t, lev=k) * \
                ds_hco["AREA"].isel(time=t) * \
                kg_s_to_g_d["Be10"]

    # Now proceed to computing the annual averages
    for spc in species_list:

        # Tropospheric-only quantities
        q[spc + "_t"] = np.ma.masked_array(q[spc + "_f"], tropmask)

        # Compute monthly averages, weighted by # of days in month
        for t in range(N_MONTHS):
            q_sum_f[t] = np.sum(q[spc + "_f"][t,:,:,:]) * d_per_mon[t]
            q_sum_t[t] = np.sum(q[spc + "_t"][t,:,:,:]) * d_per_mon[t]
            q_sum_s[t] = q_sum_f[t] - q_sum_t[t]

        # Take annual averages
        result[spc + "_f"] = np.sum(q_sum_f) / d_per_yr
        result[spc + "_t"] = np.sum(q_sum_t) / d_per_yr
        result[spc + "_s"] = np.sum(q_sum_s) / d_per_yr

    return result


def trop_residence_time(ds_cnc, ds_dry, ds_wcv, ds_wls):
    """Take the annual average of a quantity"""

    # Initialize
    result = {}

    # Loop over species
    for spc in species_list:

        # Initialize
        result[spc + "_t"] = 0.0
        
        # Concentration [g]
        var = "SpeciesConc_" + spc
        q_cnc = ds_cnc[var].values * vv_to_g[spc]
        q_cnc = np.ma.masked_array(q_cnc, tropmask)

        # DryDep [g d-1]
        var = "DryDep_" + spc
        q_dry = ds_dry[var].values * mcm2s_to_g_d[spc]

        # Convective wet scavenging [g d-1]
        var = "WetLossConv_" + spc
        q_wcv = ds_wcv[var].values * kg_s_to_g_d[spc]
        q_wcv = np.ma.masked_array(q_wcv, tropmask)
        
        # Large-scale wet scavenging [g d-1]
        var = "WetLossLS_" + spc
        q_wls = ds_wls[var].values * kg_s_to_g_d[spc]
        q_wls = np.ma.masked_array(q_wls, tropmask)

        # Loop over months
        for t in range(N_MONTHS):

            # Compute monthly averages [g]
            q_cnc_sum = np.sum(q_cnc[t,:,:,:])
            q_dry_sum = np.sum(q_dry[t,:,:]  )
            q_wcv_sum = np.sum(q_wcv[t,:,:,:])
            q_wls_sum = np.sum(q_wls[t,:,:,:])

            # Compute lifetime
            num = q_cnc_sum
            denom = q_dry_sum + q_wcv_sum + q_wls_sum 
            result[spc + "_t"] += 1.0 / (num / denom)

        # Convert residence time [d]
        result[spc + "_t"] = 1.0 / ( result[spc + "_t"] / N_MONTHS_FLOAT)
            
    return result


def print_budgets(data, key):
    """Prints the trop+strat budget file"""
    
    # Filename to print
    if "_f" in key:
        filename = "{}/{}-TransportTracers.Pb-Be_budget_trop_strat.txt".format(
            plotsdir, devstr)
    elif "_t" in key:
        filename = "{}/{}-TransportTracers.Pb-Be_budget_troposphere.txt".format(
            plotsdir, devstr)
    elif "_s"in key:
        filename = \
            "{}/{}-TransportTracers.Pb-Be_budget_stratosphere.txt".format(
             plotsdir, devstr)


    # Open file and print budgets
    with open(filename, "w+") as f:
        if "_f" in key:
            print("Table 1. Annual Average Global Budgets of 210Pb, 7Be, and 10Be\n         in the Troposphere + Stratosphere for 2016\n", file=f)
        elif "_t" in key:
            print("Table 2. Annual Average Global Budgets of 210Pb, 7Be, and 10Be\n         in the Troposphere for 2016\n", file=f)
        elif "_s" in key:
            print("Table 3. Annual Average Global Budgets of 210Pb, 7Be, and 10Be\n         in the Stratosphere for 2016\n", file=f)
        print("                                210Pb          7Be         10Be",
              file=f)
        print("                          -----------  -----------  -----------",
              file=f)

        vals = [v[1] for v in list(data["burden"].items()) if key in v[0]]
        print("  Burden, g               {:11.4f}  {:11.4f}  {:11.4f}".format(
            *vals), file=f)
        print(file=f),

        if "_t" in key:
            vals = [v[1] for v in list(data["res_time"].items()) if key in v[0]]

            print("  Residence time, d       {:11.4f}  {:11.4f}  {:11.4f}".format(*vals), file=f)
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
        
        
def transport_tracers_budgets():
    """Driver program."""

    # Initialize
    data = {}
    
    # ==================================================================
    # Get the accumulation term (init & final masses) from the restarts
    # ==================================================================
    
    # Get initial mass from restart file
    ds = xr.open_mfdataset(RstInit)
    data["accum_init"] = mass_from_rst(ds, tropmask[0,:,:,:])

    # Get initial mass from restart file
    ds = xr.open_mfdataset(RstFinal)
    data["accum_final"] = mass_from_rst(ds, tropmask[11,:,:,:])

    # Take the difference final - init
    data["accum_diff"] = diff(data["accum_init"], data["accum_final"])

    # ==================================================================
    # Burdens [g]
    # ==================================================================
    ds_cnc = xr.open_mfdataset(SpeciesConc)
    data["burden"] = annual_average(ds_cnc, "SpeciesConc", tropmask, vv_to_g)

    # ==================================================================
    # Sources and sinks [g d-1]
    # ==================================================================

    # Sources
    ds_hco = xr.open_mfdataset(HemcoDiag)
    ds_dcy = xr.open_mfdataset(RadioNucl)
    data["src_total"] = annual_average_sources(ds_hco, ds_dcy, tropmask)
    
    # Radioactive decay [g d-1]
    data["snk_decay"] = annual_average(ds_dcy, "RadDecay",
                                       tropmask, kg_s_to_g_d)
    
    # Drydep fluxes
    ds_dry = xr.open_mfdataset(DryDep)
    data["snk_drydep"]= annual_average(ds_dry, "DryDep",
                                       tropmask, mcm2s_to_g_d)

    # Convective wet scavenging
    ds_wcv = xr.open_mfdataset(WetLossConv)
    data["snk_wetconv"] = annual_average(ds_wcv, "WetLossConv",
                                         tropmask, kg_s_to_g_d)

    # Large-scale wet scavenging [g d-1]
    ds_wls = xr.open_mfdataset(WetLossLS)
    data["snk_wetls"] = annual_average(ds_wls, "WetLossLS",
                                       tropmask, kg_s_to_g_d)

    # Total sinks and sources - sinks
    data["snk_total"] = total([data["snk_drydep"], data["snk_wetls"], 
                               data["snk_wetconv"], data["snk_decay"]])
    data["src_minus_snk"] = diff(data["snk_total"], data["src_total"])

    # Tropospheric residence time
    data["res_time"] = trop_residence_time(ds_cnc, ds_dry, ds_wcv, ds_wls)
    
    # ==================================================================
    # Print budgets
    # ==================================================================
    for key in ["_f", "_t", "_s"]:
        print_budgets(data, key)
        

if __name__ == "__main__":
    transport_tracers_budgets()

