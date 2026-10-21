#!/usr/bin/env python
"""
Plots GEOS-Chem species timeseries at a single (lat, lon) point,
overlaying a Ref and (optionally) a Dev simulation.

This script reads GEOS-Chem Classic History output on a lat/lon grid
(GEOSChem.<collection>.*.nc4; GCHP cubed-sphere output is not
supported), selects the grid box nearest to the requested point at the
requested model level, and plots one column per species:

   top row    : Ref and Dev concentrations
   bottom row : (Dev - Ref) / Ref [%]  (only when --dev is given)

Mixing ratios (mol/mol) are displayed in ppbv, or in pptv for small
values.  The extracted timeseries can optionally be written to a CSV
file.

This is useful for checking how a Ref vs. Dev change behaves at a
single site.  For example, comparing level 1 with levels higher in the
boundary layer helps separate changes in surface emissions or dry
deposition from changes in boundary-layer mixing.

Examples
--------

Run the following commands:

.. code-block:: console

   $ conda activate gcpy_env
   (gcpy_env) $ python -m gcpy.examples.timeseries.plot_point_timeseries \\
                  --ref /path/to/gc_ref_rundir --dev /path/to/gc_dev_rundir \\
                  --lat 35.0 --lon -85.0 --lev 1 \\
                  --species O3 HNO3 ISOP --out point_35N_85W.png

Add --figsize W H (inches) to set the figure size, e.g. --figsize 11 8.5
for a landscape letter page.  By default the width scales with the
number of species.  Use --help to list all options.

References
----------

Author: Bob Yantosca (GitHub: @yantosca), Oct 2026, with assistance
from Claude Code.
"""
import argparse
import glob
import importlib.util
import os
import sys

import matplotlib
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

# Series colors (Ref = blue, Dev = orange)
REF_COLOR = "#2a78d6"
DEV_COLOR = "#eb6834"
INK = "#3a3a38"
MUTED = "#8a8a84"
GRID = "#e4e3dc"

# Gas constant for dry air / gravity [m/K], U.S. Standard Atmosphere 1976
R_OVER_G = 287.053 / 9.80665


def find_files(rundir, collection):
    """
    Returns the sorted list of History files for a collection.  Accepts
    either a run directory (looks in OutputDir/) or the output dir itself.

    Parameters
    ----------
    rundir : str
        Path to a GEOS-Chem Classic run directory or its OutputDir.
    collection : str
        History collection name (e.g. "SpeciesConc").

    Returns
    -------
    files : list of str
        Sorted paths of the GEOSChem.<collection>.*.nc4 files.

    Raises
    ------
    FileNotFoundError
        If no matching files exist in rundir or rundir/OutputDir.
    """
    pattern = f"GEOSChem.{collection}.*.nc4"
    for path in (os.path.join(rundir, "OutputDir"), rundir):
        files = sorted(glob.glob(os.path.join(path, pattern)))
        if files:
            return files
    raise FileNotFoundError(
        f"No files matching {pattern} in {rundir} or {rundir}/OutputDir"
    )


def resolve_varname(ds, spc):
    """
    Returns the History variable name for a species.  Tries
    SpeciesConcVV_<spc> first, then SpeciesConc_<spc>, then <spc> as
    given (so a full variable name can also be passed).

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset read from a History file.
    spc : str
        Species name (e.g. "O3") or full variable name.

    Returns
    -------
    name : str
        Name of the matching variable in ds.

    Raises
    ------
    KeyError
        If none of the candidate names is in ds.
    """
    for name in (f"SpeciesConcVV_{spc}", f"SpeciesConc_{spc}", spc):
        if name in ds.data_vars:
            return name
    raise KeyError(
        f"Could not find a variable for species '{spc}' "
        f"(tried SpeciesConcVV_{spc}, SpeciesConc_{spc}, {spc})"
    )


def wrap180(dlon):
    """
    Wraps a longitude difference into -180..180.

    Parameters
    ----------
    dlon : float or numpy.ndarray
        Longitude difference(s) [degrees].

    Returns
    -------
    dlon_wrapped : float or numpy.ndarray
        dlon wrapped into the range -180..180.
    """
    return ((dlon + 180.0) % 360.0) - 180.0


def nearest_indices(ds0, lat, lon):
    """
    Returns the (lat, lon) indices of the grid box nearest the requested
    point.  Longitude distances are taken modulo 360, so this works for
    grids in either -180..180 or 0..360 and across the dateline.

    Parameters
    ----------
    ds0 : xarray.Dataset
        Dataset with "lat" and "lon" coordinates.
    lat, lon : float
        Requested point [degrees].

    Returns
    -------
    ilat, ilon : int
        Indices of the nearest grid box.
    """
    ilat = int(np.argmin(np.abs(ds0["lat"].values - lat)))
    ilon = int(np.argmin(np.abs(wrap180(ds0["lon"].values - lon))))
    return ilat, ilon


def box_edges(ds0, ilat, ilon):
    """
    Returns the edges of grid box (ilat, ilon).  Uses the CF bounds
    written by GEOS-Chem History when present, otherwise half the grid
    spacing on either side of the box center.

    Parameters
    ----------
    ds0 : xarray.Dataset
        Dataset with "lat" and "lon" coordinates, and optionally
        "lat_bnds" and "lon_bnds".
    ilat, ilon : int
        Indices of the grid box.

    Returns
    -------
    edges : tuple of float or None
        (lat_lo, lat_hi, lon_lo, lon_hi), or None if the edges cannot
        be determined.
    """
    blat = float(ds0["lat"][ilat])
    blon = float(ds0["lon"][ilon])

    if "lat_bnds" in ds0 and "lon_bnds" in ds0:
        lat_lo, lat_hi = sorted(float(v) for v in ds0["lat_bnds"][ilat])
        lon_lo, lon_hi = (float(v) for v in ds0["lon_bnds"][ilon])
        return lat_lo, lat_hi, lon_lo, lon_hi

    if ds0.sizes["lat"] < 2 or ds0.sizes["lon"] < 2:
        return None
    dlat = float(np.median(np.abs(np.diff(ds0["lat"].values))))
    dlon = float(np.median(np.abs(np.diff(ds0["lon"].values))))
    return (blat - dlat / 2.0, blat + dlat / 2.0,
            blon - dlon / 2.0, blon + dlon / 2.0)


def check_point_in_box(lat, lon, edges, label):
    """
    Warns (to stderr) if the requested point lies outside the selected
    grid box, e.g. when the output is a regional subset that does not
    contain the point.  Longitudes are compared modulo 360.

    Parameters
    ----------
    lat, lon : float
        Requested point [degrees].
    edges : tuple of float or None
        (lat_lo, lat_hi, lon_lo, lon_hi) as returned by box_edges.
        If None, no check is done.
    label : str
        Run label (e.g. "Ref") used in the warning message.
    """
    if edges is None:
        return
    lat_lo, lat_hi, lon_lo, lon_hi = edges
    lon_ctr = 0.5 * (lon_lo + lon_hi)
    half_w = 0.5 * abs(lon_hi - lon_lo)
    tol = 1.0e-6
    inside = (lat_lo - tol <= lat <= lat_hi + tol
              and abs(wrap180(lon - lon_ctr)) <= half_w + tol)
    if not inside:
        print(f"WARNING: {label}: requested point ({lat:.2f}, {lon:.2f}) "
              f"is outside the selected grid box (lat {lat_lo:.2f}.."
              f"{lat_hi:.2f}, lon {lon_lo:.2f}..{lon_hi:.2f}).  The output "
              "may be a regional subset that does not contain this point; "
              "check LON_RANGE/LAT_RANGE in HISTORY.rc.", file=sys.stderr)


def std_atm_height_km(p_hpa):
    """
    Returns the approximate geometric height for a pressure, using the
    U.S. Standard Atmosphere (1976) layers up to 32 km and extending the
    32 km layer above that.  Surface pressure is 1013.25 hPa.

    Parameters
    ----------
    p_hpa : float
        Pressure [hPa].

    Returns
    -------
    height : float
        Approximate height [km].
    """
    if p_hpa >= 226.321:                     # Troposphere (0-11 km)
        return 44.3308 * (1.0 - (p_hpa / 1013.25) ** (R_OVER_G * 0.0065))
    if p_hpa >= 54.7489:                     # Isothermal layer (11-20 km)
        return 11.0 + R_OVER_G * 216.65 * np.log(226.321 / p_hpa) / 1000.0
    # Warming layer (20-32 km, lapse rate -1 K/km), extended upward
    return 20.0 + 216.65 * ((p_hpa / 54.7489) ** (-R_OVER_G * 0.001)
                            - 1.0)


def load_level_table():
    """
    Returns the GEOS-Chem 72-level reference table shipped with GCPy
    (gcpy/benchmark/modules/GC_72_vertical_levels.csv).  The package is
    located without importing it, since importing gcpy loads many heavy
    modules.

    Returns
    -------
    table : pandas.DataFrame or None
        The level table, or None if GCPy or the file cannot be found.
    """
    try:
        spec = importlib.util.find_spec("gcpy")
    except (ImportError, ValueError):
        return None
    if spec is None or not spec.submodule_search_locations:
        return None
    path = os.path.join(spec.submodule_search_locations[0], "benchmark",
                        "modules", "GC_72_vertical_levels.csv")
    if not os.path.isfile(path):
        return None
    return pd.read_csv(path, encoding="utf-8-sig")


def level_height_km(ds0, lev, table):
    """
    Returns the approximate height of the midpoint of a model level.

    The midpoint pressure is computed from hyam and hybm with a surface
    pressure of 1013.25 hPa.  If it matches a level of the GCPy 72-level
    table to within 0.5%, that level's number and altitude are used, so
    files saved with a "levels" subset still report the true model
    level.  Otherwise the height comes from the standard atmosphere.

    Parameters
    ----------
    ds0 : xarray.Dataset
        Dataset read from a History file.
    lev : int
        Level as a 1-based index into the file's lev dimension.
    table : pandas.DataFrame or None
        The level table returned by load_level_table.

    Returns
    -------
    height_km : float or None
        Approximate height of the level midpoint [km].
    model_level : int or None
        True model level, if it was found in the level table.
    source : str or None
        "GC 72-level table" or "standard atmosphere".

    All three are None if the file has no hybrid coefficients
    (e.g. 2-D collections) or lev is out of range.
    """
    if "hyam" not in ds0 or "hybm" not in ds0 or "lev" not in ds0.dims:
        return None, None, None
    if not 1 <= lev <= ds0.sizes["lev"]:
        return None, None, None
    hyam = float(ds0["hyam"][lev - 1])
    if ds0["hyam"].attrs.get("units", "hPa").strip() == "Pa":
        hyam /= 100.0
    p_hpa = hyam + float(ds0["hybm"][lev - 1]) * 1013.25

    if table is not None:
        p_tab = table["Pressure (hPa)"].to_numpy(dtype=float)
        idx = int(np.argmin(np.abs(p_tab - p_hpa)))
        if abs(p_tab[idx] - p_hpa) <= 0.005 * p_hpa:
            return (float(table["Altitude (km)"].iloc[idx]),
                    int(table["Lev"].iloc[idx]), "GC 72-level table")

    return float(std_atm_height_km(p_hpa)), None, "standard atmosphere"


def read_point(rundir, collection, species, lat, lon, lev, label="Ref"):
    """
    Reads timeseries for each species at the grid box nearest (lat, lon)
    and model level lev.

    Parameters
    ----------
    rundir : str
        Path to a GEOS-Chem Classic run directory or its OutputDir.
    collection : str
        History collection name (e.g. "SpeciesConc").
    species : list of str
        Species names or full variable names.
    lat, lon : float
        Requested point [degrees].
    lev : int
        Model level, 1-based (1 = surface).  Ignored for 2-D data.
    label : str, optional
        Run label used in warning messages (default: "Ref").

    Returns
    -------
    data : dict of xarray.DataArray
        Timeseries at the selected grid box, keyed by species.
    sel_lat, sel_lon : float
        Center of the selected grid box [degrees].
    height : tuple
        (height_km, model_level, source) as returned by level_height_km.

    Raises
    ------
    ValueError
        If lev is outside the range of the file's lev dimension.
    """
    files = find_files(rundir, collection)
    with xr.open_dataset(files[0]) as ds0:
        varnames = [resolve_varname(ds0, spc) for spc in species]
        ilat, ilon = nearest_indices(ds0, lat, lon)
        edges = box_edges(ds0, ilat, ilon)
        height = level_height_km(ds0, lev, load_level_table())
    check_point_in_box(lat, lon, edges, label)

    def keep_vars(ds):
        return ds[varnames]

    ds = xr.open_mfdataset(
        files,
        combine="by_coords",
        preprocess=keep_vars,
        data_vars="minimal",
        coords="minimal",
        compat="override",
        parallel=False,
    )

    # GEOS-Chem "lev" is a hybrid coordinate value, not an index,
    # so select the level by position (lev=1 is the surface).
    sel = ds.isel(lat=ilat, lon=ilon)
    if "lev" in sel.dims:
        nlev = sel.sizes["lev"]
        if not 1 <= lev <= nlev:
            raise ValueError(f"--lev {lev} outside 1..{nlev}")
        sel = sel.isel(lev=lev - 1)
    sel = sel.load()

    out = {spc: sel[vn] for spc, vn in zip(species, varnames)}
    return out, float(sel["lat"]), float(sel["lon"]), height


def display_scale(da_list):
    """
    Picks a display unit.  Mixing ratios (mol/mol) are shown in ppbv,
    or in pptv when the values are small; other units are left alone.

    Parameters
    ----------
    da_list : list of xarray.DataArray
        Data to be plotted together.  The units are taken from the
        first element.

    Returns
    -------
    scale : float
        Factor to multiply the data by.
    unit : str
        Display unit.
    """
    units = da_list[0].attrs.get("units", "")
    if "mol mol-1" in units or units.strip() in ("v/v", "mol/mol"):
        vmax = max(float(np.nanmax(np.abs(da.values))) for da in da_list)
        if vmax < 1.0e-9:
            return 1.0e12, "pptv"
        return 1.0e9, "ppbv"
    return 1.0, units


def style_axis(ax):
    """
    Keeps the grid and spines recessive and sets concise date ticks.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to style.
    """
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=INK, labelsize=9)
    loc = mdates.AutoDateLocator()
    ax.xaxis.set_major_locator(loc)
    ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(loc))


def fmt_lat(lat):
    """Formats a latitude for display, e.g. 35.00°N."""
    return f"{abs(lat):.2f}°{'N' if lat >= 0 else 'S'}"


def fmt_lon(lon):
    """Formats a longitude wrapped to -180..180 for display, e.g. 85.00°W."""
    lon = ((lon + 180.0) % 360.0) - 180.0
    return f"{abs(lon):.2f}°{'E' if lon >= 0 else 'W'}"


def main():
    """
    Parses command-line arguments, reads the Ref (and Dev) point
    timeseries, and writes the figure (and optional CSV file).
    """
    # Plot to a file, not the screen
    matplotlib.use("Agg")

    parser = argparse.ArgumentParser(
        description="Plot Ref/Dev species timeseries at a (lat, lon) point."
    )
    parser.add_argument("--ref", required=True,
                        help="Ref run directory (or its OutputDir)")
    parser.add_argument("--dev", default=None,
                        help="Dev run directory (or its OutputDir); optional")
    parser.add_argument("--refname", default="Ref",
                        help="Label for the Ref run (default: Ref)")
    parser.add_argument("--devname", default="Dev",
                        help="Label for the Dev run (default: Dev)")
    parser.add_argument("--lat", type=float, required=True,
                        help="Latitude in degrees north")
    parser.add_argument("--lon", type=float, required=True,
                        help="Longitude, either -180..180 or 0..360")
    parser.add_argument("--lev", type=int, default=1,
                        help="Model level, 1 = surface (default: 1)")
    parser.add_argument("--species", nargs="+",
                        default=["O3", "HNO3", "ISOP"],
                        help="Species names or full variable names")
    parser.add_argument("--collection", default="SpeciesConc",
                        help="History collection name (default: SpeciesConc)")
    parser.add_argument("--out", default="point_timeseries.png",
                        help="Output figure (.png or .pdf)")
    parser.add_argument("--csv", default=None,
                        help="Optionally write the extracted series to CSV")
    parser.add_argument("--figsize", nargs=2, type=float, default=None,
                        metavar=("W", "H"),
                        help="Figure size in inches, e.g. 11 8.5 for "
                             "landscape letter (default: scales with the "
                             "number of species)")
    args = parser.parse_args()

    # ------------------------------------------------------------------
    # Read data
    # ------------------------------------------------------------------
    ref, rlat, rlon, height = read_point(
        args.ref, args.collection, args.species,
        args.lat, args.lon, args.lev, label=args.refname)
    print(f"{args.refname}: nearest grid box lat={rlat:.3f}, "
          f"lon={rlon:.3f}, lev={args.lev}")
    height_km, model_lev, height_src = height
    if height_km is not None:
        print(f"Level height: ~{height_km:.3f} km ({height_src}"
              + (f", model level {model_lev}" if model_lev else "") + ")")

    dev = None
    if args.dev:
        dev, dlat, dlon, _ = read_point(args.dev, args.collection, args.species,
                                     args.lat, args.lon, args.lev,
                                     label=args.devname)
        print(f"{args.devname}: nearest grid box lat={dlat:.3f}, "
              f"lon={dlon:.3f}, lev={args.lev}")
        if not (np.isclose(rlat, dlat) and np.isclose(rlon, dlon)):
            print("WARNING: Ref and Dev selected different grid boxes!",
                  file=sys.stderr)

    # ------------------------------------------------------------------
    # Plot
    # ------------------------------------------------------------------
    nspc = len(args.species)
    nrow = 2 if dev is not None else 1
    fig, axes = plt.subplots(
        nrow, nspc,
        figsize=(tuple(args.figsize) if args.figsize
                 else (5.2 * nspc, 3.4 * nrow + 0.6)),
        squeeze=False,
        sharex="col",
        gridspec_kw={"height_ratios": [2, 1]} if nrow == 2 else None,
    )

    table = {}
    for col, spc in enumerate(args.species):
        r = ref[spc]
        d = dev[spc] if dev is not None else None
        if d is not None:
            r, d = xr.align(r, d, join="inner")
            if r.sizes["time"] == 0:
                sys.exit(f"No common times between Ref and Dev for {spc}")

        scale, unit = display_scale([r] if d is None else [r, d])
        times = pd.to_datetime(r["time"].values)

        # Concentrations
        ax = axes[0, col]
        ax.plot(times, r.values * scale, color=REF_COLOR, lw=1,
                label=args.refname)
        if d is not None:
            ax.plot(times, d.values * scale, color=DEV_COLOR, lw=1,
                    label=args.devname)
            ax.legend(frameon=False, fontsize=9, labelcolor=INK)
        ax.set_title(spc, color=INK, fontsize=12, loc="left",
                     fontweight="bold")
        ax.set_ylabel(unit, color=INK, fontsize=10)
        style_axis(ax)

        table[f"{args.refname}_{spc}_{unit}"] = \
            r.values.astype(np.float64) * scale

        # Percent difference
        if d is not None:
            rv = r.values.astype(np.float64)
            dv = d.values.astype(np.float64)
            with np.errstate(divide="ignore", invalid="ignore"):
                pct = np.where(rv != 0.0, (dv - rv) / rv * 100.0, np.nan)

            ax = axes[1, col]
            ax.axhline(0.0, color=MUTED, lw=1)
            ax.plot(times, pct, color=INK, lw=1)
            ax.set_ylabel(f"({args.devname}−{args.refname})/"
                          f"{args.refname} [%]", color=INK, fontsize=9)
            if np.isfinite(pct).any():
                ax.text(0.01, 0.95,
                        f"mean {np.nanmean(pct):+.2f}%   "
                        f"max |Δ| {np.nanmax(np.abs(pct)):.2f}%",
                        transform=ax.transAxes, va="top", ha="left",
                        fontsize=8.5, color=INK)
            style_axis(ax)

            table[f"{args.devname}_{spc}_{unit}"] = dv * scale
            table[f"pctdiff_{spc}"] = pct

    t0 = pd.to_datetime(ref[args.species[0]]["time"].values[0])
    t1 = pd.to_datetime(ref[args.species[0]]["time"].values[-1])
    level_str = f"level {model_lev or args.lev}"
    if height_km is not None:
        level_str += f", ≈{height_km:.3f} km"
    fig.suptitle(
        f"{fmt_lat(rlat)}, {fmt_lon(rlon)}, {level_str}   "
        f"({t0:%Y-%m-%d %H:%M} to {t1:%Y-%m-%d %H:%M})",
        color=INK, fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"Wrote {args.out}")

    # ------------------------------------------------------------------
    # Optional CSV output
    # ------------------------------------------------------------------
    if args.csv:
        r0 = ref[args.species[0]]
        if dev is not None:
            r0, _ = xr.align(r0, dev[args.species[0]], join="inner")
        df = pd.DataFrame(table, index=pd.to_datetime(r0["time"].values))
        df.index.name = "time"
        df.to_csv(args.csv)
        print(f"Wrote {args.csv}")


if __name__ == "__main__":
    main()
