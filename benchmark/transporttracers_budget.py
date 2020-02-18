#!/usr/bin/env python

"""
Computes the budget of Pb and Be7 from the TransportTracers benchmarks.

NOTE: In the future we can perhaps refactor this to use data structures
such as classes.  The initial development was intended to get this feature
working ASAP for the benchmarks. 
  -- Bob Yantosca (11 Feb 2020)
"""

# ======================================================================
# Imports etc.
# ======================================================================

from calendar import monthrange
import numpy as np
import os
from os.path import join
from gcpy.constants import AVOGADRO
from gcpy.benchmark import get_troposphere_mask
import warnings
import xarray as xr

# Suppress harmless run-time warnings (mostly about underflow in division)
warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', category=UserWarning)

# ======================================================================
# Configurables (MUST EDIT)
# ======================================================================

maindir     = '/path/to/data/dir'

# Version string
devstr      = 'version_string'

# Restart collections
rstdir      = join(maindir, 'restarts',                      )
RstInit     = join(rstdir,  'GEOSChem.Restart.2016*.nc4'     )
RstFinal    = join(rstdir,  'GEOSChem.Restart.2017*.nc4'     )

# Data collections
datadir     = join(maindir, 'OutputDir'                      )
HemcoDiag   = join(datadir, 'HEMCO_diagnostics.2016*.nc'     )
DryDep      = join(datadir, 'GEOSChem.DryDep.2016*.nc4'      )
RadioNucl   = join(datadir, 'GEOSChem.RadioNuclide.2016*.nc4')
StateMet    = join(datadir, 'GEOSChem.StateMet.2016*.nc4'    )
SpeciesConc = join(datadir, 'GEOSChem.SpeciesConc.2016*.nc4' )
WetLossConv = join(datadir, 'GEOSChem.WetLossConv.2016*.nc4' )
WetLossLS   = join(datadir, 'GEOSChem.WetLossLS.2016*.nc4'   )

# ======================================================================
# Meteorological variables, constants, conversion factors
# ======================================================================

# Number of months (normally=12, but can use 1 for testing)
N_MONTHS             = 12
N_MONTHS_FLOAT       = N_MONTHS * 1.0

# Read the StateMet collection
ds                   = xr.open_mfdataset(StateMet)

# Read certain met data variables
area_m2              = ds['AREA'] #.isel(time=0)
area_cm2             = area_m2 * 1.0e4
tropmask             = get_troposphere_mask(ds)

# Molecular weights
MW_AIR               = 28.9644       
MW_Pb210             = 210.0
MW_Be7               = 7.0
MW_Be10              = 10.0

# Conversion factors
KG_PER_MOL_Pb210     = AVOGADRO / (MW_Pb210 * 1e-3)
KG_PER_MOL_Be7       = AVOGADRO / (MW_Be7   * 1e-3)
KG_PER_MOL_Be10      = AVOGADRO / (MW_Be10  * 1e-3)
KG_S_TO_G_DAY        = 86400.0 * 1000.0      
Pb210_VV_TO_G        = ds['Met_AD'] * (MW_Pb210 / MW_AIR) * 1000.0
Be7_VV_TO_G          = ds['Met_AD'] * (MW_Be7   / MW_AIR) * 1000.0         
Be10_VV_TO_G         = ds['Met_AD'] * (MW_Be10  / MW_AIR) * 1000.0
Pb210_MCM2S_TO_G_DAY = area_cm2 / KG_PER_MOL_Pb210 * KG_S_TO_G_DAY
Be7_MCM2S_TO_G_DAY   = area_cm2 / KG_PER_MOL_Be7   * KG_S_TO_G_DAY
Be10_MCM2S_TO_G_DAY  = area_cm2 / KG_PER_MOL_Be10  * KG_S_TO_G_DAY

# Create arrays for monthly sums
Pb210_sum            = np.zeros(N_MONTHS)
Be7_sum              = np.zeros(N_MONTHS)
Be10_sum             = np.zeros(N_MONTHS)
Pb210_sum_tr         = np.zeros(N_MONTHS)
Be7_sum_tr           = np.zeros(N_MONTHS)
Be10_sum_tr          = np.zeros(N_MONTHS)

# ======================================================================
# Compute initial and final mass [g] from the restart files
# ======================================================================

# ------------------------------------
# Initial mass [mol/mol dry] -> [g]
# ------------------------------------

# Read initial restart data
ds             = xr.open_mfdataset(RstInit)

# Full-atmosphere initial mass arrays [g]
Pb210_init     = ds['SpeciesRst_Pb210'].isel(time=0) * Pb210_VV_TO_G
Be7_init       = ds['SpeciesRst_Be7'  ].isel(time=0) * Be7_VV_TO_G
Be10_init      = ds['SpeciesRst_Be10' ].isel(time=0) * Be10_VV_TO_G

# Trop-only initial mass arrays [g]
Pb210_init_tr  = np.ma.masked_array(Pb210_init.values, tropmask)
Be7_init_tr    = np.ma.masked_array(Be7_init.values,   tropmask)
Be10_init_tr   = np.ma.masked_array(Be10_init.values,  tropmask)

# Full-atmosphere annual average mass [g] (recycle variable names)
Pb210_init     = np.sum(Pb210_init.values) / N_MONTHS_FLOAT
Be7_init       = np.sum(Be7_init.values  ) / N_MONTHS_FLOAT
Be10_init      = np.sum(Be10_init.values ) / N_MONTHS_FLOAT

# Trop-only annual average mass [g] (recycle variable names)
Pb210_init_tr  = np.sum(Pb210_init_tr) / N_MONTHS_FLOAT
Be7_init_tr    = np.sum(Be7_init_tr  ) / N_MONTHS_FLOAT
Be10_init_tr   = np.sum(Be10_init_tr ) / N_MONTHS_FLOAT

# Strat-only annual average mass [g]
Pb210_init_st  = Pb210_init - Pb210_init_tr
Be7_init_st    = Be7_init   - Be7_init_tr
Be10_init_st   = Be10_init  - Be10_init_tr

# ------------------------------------
# Final mass [mol/mol dry] -> [g]
# ------------------------------------

# Read final restart file
ds             = xr.open_mfdataset(RstFinal)

# Full-atmosphere mass arrays [g]
Pb210_final    = ds['SpeciesRst_Pb210'].isel(time=0) * Pb210_VV_TO_G
Be7_final      = ds['SpeciesRst_Be7'  ].isel(time=0) * Be7_VV_TO_G
Be10_final     = ds['SpeciesRst_Be10' ].isel(time=0) * Be10_VV_TO_G

# Trop-only mass arrays [g]
Pb210_final_tr = np.ma.masked_array(Pb210_final.values, tropmask)
Be7_final_tr   = np.ma.masked_array(Be7_final.values,   tropmask)
Be10_final_tr  = np.ma.masked_array(Be10_final.values,  tropmask)

# Full-atmosphere annual average mass sums [g] (recycle variable names)
Pb210_final    = np.sum(Pb210_final.values) / N_MONTHS_FLOAT
Be7_final      = np.sum(Be7_final.values  ) / N_MONTHS_FLOAT
Be10_final     = np.sum(Be10_final.values ) / N_MONTHS_FLOAT

# Trop-only annual average mass [g] (recycle variable names)
Pb210_final_tr = np.sum(Pb210_final_tr) / N_MONTHS_FLOAT
Be7_final_tr   = np.sum(Be7_final_tr  ) / N_MONTHS_FLOAT
Be10_final_tr  = np.sum(Be10_final_tr ) / N_MONTHS_FLOAT

# Strat-only annual average mass [g]
Pb210_final_st = Pb210_final - Pb210_final_tr
Be7_final_st   = Be7_final   - Be7_final_tr
Be10_final_st  = Be10_final  - Be10_final_tr

# ------------------------------------
# Final mass - Initial mass [g]
# ------------------------------------
Pb210_diff     = Pb210_final    - Pb210_init
Pb210_diff_tr  = Pb210_final_tr - Pb210_init_tr
Pb210_diff_st  = Pb210_final_st - Pb210_init_st
Be7_diff       = Be7_final      - Be7_init     
Be7_diff_tr    = Be7_final_tr   - Be7_init_tr  
Be7_diff_st    = Be7_final_st   - Be7_init_st  
Be10_diff      = Be10_final     - Be10_init     
Be10_diff_tr   = Be10_final_tr  - Be10_init_tr
Be10_diff_st   = Be10_final_st  - Be10_init_st

# ======================================================================
# Burdens [mol/mol dry] -> [g]
# ======================================================================

# Read initial restart data
ds                  = xr.open_mfdataset(SpeciesConc)

# Full-atmosphere atmospheric burdens [g]
Pb210_burden        = ds['SpeciesConc_Pb210'] * Pb210_VV_TO_G
Be7_burden          = ds['SpeciesConc_Be7'  ] * Be7_VV_TO_G
Be10_burden         = ds['SpeciesConc_Be10' ] * Be10_VV_TO_G

# Trop-only atmospheric burdens [g]
Pb210_burden_tr     = np.ma.masked_array(Pb210_burden.values, tropmask)
Be7_burden_tr       = np.ma.masked_array(Be7_burden.values,   tropmask)
Be10_burden_tr      = np.ma.masked_array(Be10_burden.values,  tropmask)

# Compute monthly sums [g]
for t in range(N_MONTHS):
    Pb210_sum[t]    = np.sum(Pb210_burden.isel(time=t).values)
    Be7_sum[t]      = np.sum(Be7_burden.isel(time=t).values)
    Be10_sum[t]     = np.sum(Be10_burden.isel(time=t).values)
    Pb210_sum_tr[t] = np.sum(Pb210_burden_tr[t,:,:,:])
    Be7_sum_tr[t]   = np.sum(Be7_burden_tr[t,:,:,:]  )
    Be10_sum_tr[t]  = np.sum(Be10_burden_tr[t,:,:,:] )

# Save trop-only sums [g] for residence time computation
Pb210_sum_burden    = Pb210_sum_tr
Be7_sum_burden      = Be7_sum_tr
Be10_sum_burden     = Be10_sum_tr

# Compute annual average burdens [g] (recycle variable names)
Pb210_burden        = np.mean(Pb210_sum)
Pb210_burden_tr     = np.mean(Pb210_sum_tr)
Pb210_burden_st     = Pb210_burden - Pb210_burden_tr
Be7_burden          = np.mean(Be7_sum)
Be7_burden_tr       = np.mean(Be7_sum_tr)
Be7_burden_st       = Be7_burden - Be7_burden_tr
Be10_burden         = np.mean(Be10_sum)
Be10_burden_tr      = np.mean(Be10_sum_tr)
Be10_burden_st      = Be10_burden - Be10_burden_tr

# ======================================================================
# Drydep [molec/cm2/s] -> [g/day]
# NOTE: By nature, drydep is only in the troposphere
# ======================================================================

# Read initial restart data
ds               = xr.open_mfdataset(DryDep)

# Drydep flux [g/day]
Pb210_dry        = ds['DryDep_Pb210'] * Pb210_MCM2S_TO_G_DAY
Be7_dry          = ds['DryDep_Be7']   * Be7_MCM2S_TO_G_DAY
Be10_dry         = ds['DryDep_Be10']  * Be7_MCM2S_TO_G_DAY

# Compute monthly sums [g]
for t in range(N_MONTHS):
    Pb210_sum[t] = np.sum(Pb210_dry.isel(time=t).values)
    Be7_sum[t]   = np.sum(Be7_dry.isel(time=t).values)
    Be10_sum[t]  = np.sum(Be10_dry.isel(time=t).values)

# Save trop-only sums [g] for residence time computation
Pb210_sum_dry    = Pb210_sum
Be7_sum_dry      = Be7_sum
Be10_sum_dry     = Be10_sum

# Compute annual average drydep flux [g/day] (recycle variable names)
Pb210_dry        = np.mean(Pb210_sum)
Be7_dry          = np.mean(Be7_sum)
Be10_dry         = np.mean(Be10_sum)

# ======================================================================
# Large scale wetdep [kg/s] -> [g day]
# ======================================================================

# Read wetdep data
ds                  = xr.open_mfdataset(WetLossLS)

# Large-scale wet scavenging [g/day]
Pb210_wls           = ds['WetLossLS_Pb210'] * KG_S_TO_G_DAY
Be7_wls             = ds['WetLossLS_Be7']   * KG_S_TO_G_DAY
Be10_wls            = ds['WetLossLS_Be10']  * KG_S_TO_G_DAY

# Large-scale wet scavenging [g/day], troposphere only
Pb210_wls_tr        = np.ma.masked_array(Pb210_wls.values, tropmask)
Be7_wls_tr          = np.ma.masked_array(Be7_wls.values,   tropmask)
Be10_wls_tr         = np.ma.masked_array(Be10_wls.values,  tropmask)

# Create monthly sums [g/day]
for t in range(N_MONTHS):
    Pb210_sum[t]    = np.sum(Pb210_wls.isel(time=t).values)
    Be7_sum[t]      = np.sum(Be7_wls.isel(time=t).values)
    Be10_sum[t]     = np.sum(Be10_wls.isel(time=t).values)
    Pb210_sum_tr[t] = np.sum(Pb210_wls_tr[t,:,:,:])
    Be7_sum_tr[t]   = np.sum(Be7_wls_tr[t,:,:,:])
    Be10_sum_tr[t]  = np.sum(Be10_wls_tr[t,:,:,:])

# Save trop-only sums [g] for residence time computation
Pb210_sum_wls       = Pb210_sum_tr
Be7_sum_wls         = Be7_sum_tr
Be10_sum_wls        = Be10_sum_tr

# Annual average large-scale wet-scavenging [g/day] (recycle variable names)
Pb210_wls           = np.mean(Pb210_sum)
Pb210_wls_tr        = np.mean(Pb210_sum_tr)
Be7_wls             = np.mean(Be7_sum)
Be7_wls_tr          = np.mean(Be7_sum_tr)
Be10_wls            = np.mean(Be10_sum)
Be10_wls_tr         = np.mean(Be10_sum_tr)

# ======================================================================
# Convective wet scavenging [kg/s] -> [g day]
# ======================================================================

# Read convective wet scavenging data
ds                  = xr.open_mfdataset(WetLossConv)

# Convective wet scavenging [g/day], full-atmosphere
Pb210_wcv           = ds['WetLossConv_Pb210'] * KG_S_TO_G_DAY
Be7_wcv             = ds['WetLossConv_Be7']   * KG_S_TO_G_DAY
Be10_wcv            = ds['WetLossConv_Be10']  * KG_S_TO_G_DAY

# Convective wet scavenging [g/day], troposphere only
Pb210_wcv_tr        = np.ma.masked_array(Pb210_wcv.values, tropmask)
Be7_wcv_tr          = np.ma.masked_array(Be7_wcv.values,   tropmask)
Be10_wcv_tr         = np.ma.masked_array(Be10_wcv.values,  tropmask)

# Compute monthly sums [g]
for t in range(N_MONTHS):
    Pb210_sum[t]    = np.sum(Pb210_wcv.isel(time=t).values)
    Be7_sum[t]      = np.sum(Be7_wcv.isel(time=t).values)
    Be10_sum[t]     = np.sum(Be10_wcv.isel(time=t).values)
    Pb210_sum_tr[t] = np.sum(Pb210_wcv_tr[t,:,:,:])
    Be7_sum_tr[t]   = np.sum(Be7_wcv_tr[t,:,:,:])
    Be10_sum_tr[t]  = np.sum(Be10_wcv_tr[t,:,:,:])

# Save trop-only sums [g] for residence time computation
Pb210_sum_wcv       = Pb210_sum_tr
Be7_sum_wcv         = Be7_sum_tr
Be10_sum_wcv        = Be10_sum_tr

# Annual average convective wet scavenging [g/day] (recycle variable names)
Pb210_wcv           = np.mean(Pb210_sum)
Pb210_wcv_tr        = np.mean(Pb210_sum_tr)
Be7_wcv             = np.mean(Be7_sum)
Be7_wcv_tr          = np.mean(Be7_sum_tr)
Be10_wcv            = np.mean(Be10_sum)
Be10_wcv_tr         = np.mean(Be10_sum_tr)

# ======================================================================
# Emissions [kg/m2/s] -> [g day]
# ======================================================================

# Read sources and sinks data
ds                  = xr.open_mfdataset(HemcoDiag)
ds_dcy              = xr.open_mfdataset(RadioNucl)

# Sources [g/day], full-atmosphere
Pb210_src           = ( ds_dcy['PbFromRnDecay']         ) * KG_S_TO_G_DAY
Be7_src             = ( ds['EmisBe7_Cosmic']  * area_m2 ) * KG_S_TO_G_DAY
Be10_src            = ( ds['EmisBe10_Cosmic'] * area_m2 ) * KG_S_TO_G_DAY

# Sources [g/day], troposphere-only
Pb210_src_tr        = np.ma.masked_array(Pb210_src.values, tropmask)
Be7_src_tr          = np.ma.masked_array(Be7_src.values,   tropmask)
Be10_src_tr         = np.ma.masked_array(Be10_src.values,  tropmask)

# Compute monthly sums
for t in range(N_MONTHS):
    Pb210_sum[t]    = np.sum(Pb210_src.isel(time=t).values)
    Be7_sum[t]      = np.sum(Be7_src.isel(time=t).values)
    Be10_sum[t]     = np.sum(Be10_src.isel(time=t).values)
    Pb210_sum_tr[t] = np.sum(Pb210_src_tr[t,:,:,:])
    Be7_sum_tr[t]   = np.sum(Be7_src_tr[t,:,:,:])
    Be10_sum_tr[t]  = np.sum(Be10_src_tr[t,:,:,:])

# Annual average sources [g/day] (recycle variable names)
Pb210_src           = np.mean(Pb210_sum)
Pb210_src_tr        = np.mean(Pb210_sum_tr)
Pb210_src_st        = Pb210_src - Pb210_src_tr
Be7_src             = np.mean(Be7_sum)
Be7_src_tr          = np.mean(Be7_sum_tr)
Be7_src_st          = Be7_src - Be7_src_tr
Be10_src            = np.mean(Be10_sum)
Be10_src_tr         = np.mean(Be10_sum_tr)
Be10_src_st         = Be10_src - Be10_src_tr

# ======================================================================
# Radioactive decay [kg/s] -> [g day]
# ======================================================================

# Radioactive decay [g/day], full-atmosphere
Pb210_dcy           = ds_dcy['RadDecay_Pb210'] * KG_S_TO_G_DAY
Be7_dcy             = ds_dcy['RadDecay_Be7']   * KG_S_TO_G_DAY
Be10_dcy            = ds_dcy['RadDecay_Be10']  * KG_S_TO_G_DAY

# Radioactive decay [g/day], troposphere only
Pb210_dcy_tr        = np.ma.masked_array(Pb210_dcy.values, tropmask)
Be7_dcy_tr          = np.ma.masked_array(Be7_dcy.values,   tropmask)
Be10_dcy_tr         = np.ma.masked_array(Be10_dcy.values,  tropmask)

# Compute monthly sums
for t in range(N_MONTHS):
    Pb210_sum[t]    = np.sum(Pb210_dcy.isel(time=t).values)
    Be7_sum[t]      = np.sum(Be7_dcy.isel(time=t).values)
    Be10_sum[t]     = np.sum(Be10_dcy.isel(time=t).values)
    Pb210_sum_tr[t] = np.sum(Pb210_dcy_tr[t,:,:,:])
    Be7_sum_tr[t]   = np.sum(Be7_dcy_tr[t,:,:,:])
    Be10_sum_tr[t]  = np.sum(Be10_dcy_tr[t,:,:,:])

# Annual average sinks [g/day] (recycle variable names)
Pb210_dcy           = np.mean(Pb210_sum)
Pb210_dcy_tr        = np.mean(Pb210_sum_tr)
Pb210_dcy_st        = Pb210_dcy - Pb210_dcy_tr
Be7_dcy             = np.mean(Be7_sum)
Be7_dcy_tr          = np.mean(Be7_sum_tr)
Be7_dcy_st          = Be7_dcy - Be7_dcy_tr
Be10_dcy            = np.mean(Be10_sum)
Be10_dcy_tr         = np.mean(Be10_sum_tr)
Be10_dcy_st         = Be10_dcy - Be10_dcy_tr

# ======================================================================
# Compute total sinks [g/day]
# ======================================================================
Pb210_snk           = Pb210_dry + Pb210_wls + Pb210_wcv + Pb210_dcy
Be7_snk             = Be7_dry   + Be7_wls   + Be7_wcv   + Be7_dcy
Be10_snk            = Be10_dry  + Be10_wls  + Be10_wcv  + Be10_dcy
Pb210_snk_tr        = Pb210_dry + Pb210_wls + Pb210_wcv + Pb210_dcy_tr
Be7_snk_tr          = Be7_dry   + Be7_wls   + Be7_wcv   + Be7_dcy_tr
Be10_snk_tr         = Be10_dry  + Be10_wls  + Be10_wcv  + Be10_dcy_tr

# ======================================================================
# Compute tropospheric residence time [days]
# ======================================================================
tau_Pb210       = 0.0
tau_Be7         = 0.0
tau_Be10        = 0.0

for t in range(N_MONTHS):
    days_in_mon = monthrange(2016,t+1)[1] * 1.0

    # Lifetime of Pb210 [days]
    num         = Pb210_sum_burden[t]
    denom       = Pb210_sum_dry[t] + Pb210_sum_wcv[t] + Pb210_sum_wls[t] 
    tau_Pb210  += 1.0 / (num / denom * days_in_mon)

    # Lifetime of Be7 [days]
    num         = Be7_sum_burden[t]
    denom       = Be7_sum_dry[t] + Be7_sum_wcv[t] + Be7_sum_wls[t] 
    tau_Be7    += 1.0 / (num / denom * days_in_mon)

    # Lifetime of Be210 [days]
    num         = Be10_sum_burden[t]
    denom       = Be10_sum_dry[t] + Be10_sum_wcv[t] + Be10_sum_wls[t] 
    tau_Be10   += 1.0 / (num / denom * days_in_mon)

# Residence time [days]
res_Pb210       = 1.0 / ( tau_Pb210 / N_MONTHS_FLOAT )
res_Be7         = 1.0 / ( tau_Be7   / N_MONTHS_FLOAT )
res_Be10        = 1.0 / ( tau_Be10  / N_MONTHS_FLOAT )

# ======================================================================
# Table 1: Trop + Strat budgets
# ======================================================================
filename = "{}-TransportTracers.Pb-Be_budget_trop_strat.txt".format(devstr)
with open(filename, "w+") as f:
    print("Table 1. Annual Average Global Budgets of 210Pb, 7Be, and 10Be",
          file=f)
    print("         in the Troposphere + Stratosphere for 2016\n", file=f)
    print("                                210Pb          7Be         10Be",
          file=f)
    print("                          -----------  -----------  -----------",
          file=f)
    print("  Burden, g               {:11.4f}  {:11.4f}  {:11.4f}".format(
        Pb210_burden,  Be7_burden, Be10_burden), file=f)
    print("    Stratosphere          {:11.4f}  {:11.4f}  {:11.4f}".format(
        Pb210_burden_st, Be7_burden_st, Be10_burden_st), file=f)
    print("    Troposphere           {:11.4f}  {:11.4f}  {:11.4f}".format(
        Pb210_burden_tr, Be7_burden_tr, Be10_burden_tr), file=f)    
    print(file=f),
    print("  Sources, g d-1          {:11.4f}  {:11.4f}  {:11.4f}".format(
        Pb210_src, Be7_src, Be10_src), file=f)
    print("    Stratosphere          {:11.6f}  {:11.6f}  {:11.6f}".format(
        Pb210_src_st, Be7_src_st, Be10_src_st), file=f)
    print("    Troposphere           {:11.4f}  {:11.4f}  {:11.4f}".format(
        Pb210_src_tr, Be7_src_tr, Be10_src_tr), file=f)
    print(file=f)
    print("  Sinks, g d-1            {:11.4f}  {:11.4f}  {:11.4f}".format(
        Pb210_snk, Be7_snk, Be10_snk), file=f)
    print("    Dry deposition        {:11.4f}  {:11.4f}  {:11.4f}".format(
        Pb210_dry, Be7_dry, Be10_dry), file=f)
    print("    Wet deposition", file=f)
    print("      Stratiform          {:11.4f}  {:11.4f}  {:11.4f}".format(
        Pb210_wls, Be7_wls, Be10_wls), file=f)
    print("      Convective          {:11.4f}  {:11.4f}  {:11.4f}".format(
        Pb210_wcv, Be7_wcv, Be10_wcv), file=f)
    print("    Radioactive decay     {:11.7f}  {:11.7f}  {:11.8f}".format(
        Pb210_dcy, Be7_dcy, Be10_dcy), file=f)
    print("      Stratosphere        {:11.7f}  {:11.7f}  {:11.8f}".format(
        Pb210_dcy_st, Be7_dcy_st, Be10_dcy_st), file=f)
    print("      Tropoosphere        {:11.7f}  {:11.7f}  {:11.8f}".format(
        Pb210_dcy_tr, Be7_dcy_tr, Be10_dcy_tr), file=f)

    print(file=f)
    print("  Accumulation Term", file=f)
    print("    Initial mass, g       {:11.4f}  {:11.4f}  {:11.4f}".format(
        Pb210_init, Be7_init, Be10_init), file=f)
    print("      Stratosphere        {:11.4f}  {:11.4f}  {:11.4f}".format(
        Pb210_init_st, Be7_init_st, Be10_init_st), file=f)
    print("      Troposphere         {:11.4f}  {:11.4f}  {:11.4f}".format(
        Pb210_init_tr, Be7_init_tr, Be10_init_tr), file=f)
    print(file=f)
    print("    Final mass, g         {:11.4f}  {:11.4f}  {:11.4f}".format(
        Pb210_final, Be7_final, Be10_final), file=f)
    print("      Stratosphere        {:11.4f}  {:11.4f}  {:11.4f}".format(
        Pb210_final_st, Be7_final_st, Be10_final_st), file=f)
    print("      Troposphere         {:11.4f}  {:11.4f}  {:11.4f}".format(
        Pb210_final_tr, Be7_final_tr, Be10_final_tr), file=f)
    print(file=f)
    print("    Difference, g         {:11.7f}  {:11.7f}  {:11.7f}".format(
        Pb210_diff, Be7_diff, Be10_diff), file=f)
    print("      Stratosphere        {:11.7f}  {:11.7f}  {:11.7f}".format(
        Pb210_diff_st, Be7_diff_st, Be10_diff_st), file=f)
    print("      Troposphere         {:11.7f}  {:11.7f}  {:11.7f}".format(
        Pb210_diff_tr, Be7_diff_tr, Be10_diff_tr), file=f)
    f.close()

# ======================================================================
# Table 2: Trop only budgets
# ======================================================================
filename = "{}-TransportTracers.Pb-Be_budget_troposphere.txt".format(devstr)
with open(filename, "w+") as f:
    print("Table 1. Annual Average Global Budgets of 210Pb, 7Be, and 10Be",
          file=f)
    print("         in the Troposphere for 2016\n", file=f)
    print("                                210Pb          7Be         10Be",
          file=f)
    print("                          -----------  -----------  -----------",
          file=f)
    print("  Burden, g               {:11.4f}  {:11.4f}  {:11.4f}".format(
        Pb210_burden_tr,  Be7_burden_tr, Be10_burden_tr), file=f)
    print(file=f),
    print("  Residence time, d       {:11.4f}  {:11.4f}  {:11.4f}".format(
        res_Pb210, res_Be7, res_Be10), file=f)
    print(file=f)
    print("  Sources, g d-1", file=f)
    print("    From stratosphere     {:11.6f}  {:11.6f}  {:11.6f}".format(
        Pb210_src_st, Be7_src_st, Be10_src_st), file=f)
    print("    Within troposphere    {:11.4f}  {:11.4f}  {:11.4f}".format(
        Pb210_src_tr, Be7_src_tr, Be10_src_tr), file=f)
    print(file=f)
    print("  Sinks, g d-1            {:11.4f}  {:11.4f}  {:11.4f}".format(
        Pb210_snk, Be7_snk, Be10_snk), file=f)
    print("    Dry deposition        {:11.4f}  {:11.4f}  {:11.4f}".format(
        Pb210_dry, Be7_dry, Be10_dry), file=f)
    print("    Wet deposition", file=f)
    print("      Stratiform          {:11.4f}  {:11.4f}  {:11.4f}".format(
        Pb210_wls, Be7_wls, Be10_wls), file=f)
    print("      Convective          {:11.4f}  {:11.4f}  {:11.4f}".format(
        Pb210_wcv, Be7_wcv, Be10_wcv), file=f)
    print("    Radioactive decay     {:11.7f}  {:11.7f}  {:11.8f}".format(
        Pb210_dcy_tr, Be7_dcy_tr, Be10_dcy_tr), file=f)
    print(file=f)
    print("  Accumulation Term", file=f)
    print("    Initial mass, g       {:11.4f}  {:11.4f}  {:11.4f}".format(
        Pb210_init_tr, Be7_init_tr, Be10_init_tr), file=f)
    print("    Final mass, g         {:11.4f}  {:11.4f}  {:11.4f}".format(
        Pb210_final_tr, Be7_final_tr, Be10_final_tr), file=f)
    print("    Difference, g         {:11.7f}  {:11.7f}  {:11.7f}".format(
        Pb210_diff_tr, Be7_diff_tr, Be10_diff_tr), file=f)
    f.close()
