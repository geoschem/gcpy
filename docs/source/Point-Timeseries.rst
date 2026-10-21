.. _point-timeseries:

##########################
Plot timeseries at a point
##########################

This example plots GEOS-Chem Classic species timeseries at the grid box nearest to a requested (latitude, longitude) point and model level, overlaying a Ref and (optionally) a Dev simulation. When a Dev simulation is given, the percent difference :literal:`(Dev - Ref) / Ref` is plotted underneath. This is useful for checking how a model change behaves at a single site. For example, comparing level 1 with levels higher in the boundary layer can help separate changes in surface emissions or dry deposition from changes in boundary-layer mixing.

.. _point-timeseries-code:

===========
Source code
===========

.. list-table::
   :header-rows: 1

   * - Description
     - Script location
   * - :mod:`gcpy.examples.timeseries.plot_point_timeseries`
     - `gcpy/examples/timeseries/plot_point_timeseries.py <https://github.com/geoschem/gcpy/blob/main/gcpy/examples/timeseries/plot_point_timeseries.py>`_

.. _point-timeseries-usage:

=====
Usage
=====

Activate your GCPy environment with conda:

.. code-block:: console

   $ conda activate gcpy_env

Then run the script as a module. For example, to plot O3, HNO3, and ISOP in the surface level at 35°N, 85°W:

.. code-block:: console

   $ python -m gcpy.examples.timeseries.plot_point_timeseries \
       --ref /path/to/gc_ref_rundir \
       --dev /path/to/gc_dev_rundir \
       --lat 35.0 --lon -85.0 --lev 1 \
       --species O3 HNO3 ISOP \
       --out point_35N_85W.png

The :literal:`--ref` and :literal:`--dev` arguments can point either to a GEOS-Chem Classic run directory (the script then looks in its :file:`OutputDir/` subfolder) or to the folder containing the History files themselves. Omit :literal:`--dev` to plot the Ref simulation alone.

.. _point-timeseries-args:

======================
Command-line arguments
======================

.. list-table::
   :header-rows: 1
   :widths: 25 50 25

   * - Argument
     - Description
     - Default
   * - :literal:`--ref`
     - Ref run directory (or its :file:`OutputDir`). Required.
     - none
   * - :literal:`--dev`
     - Dev run directory (or its :file:`OutputDir`).
     - none (plot Ref only)
   * - :literal:`--refname`
     - Label for the Ref run, used in the legend, axis labels, and CSV column names.
     - :literal:`Ref`
   * - :literal:`--devname`
     - Label for the Dev run.
     - :literal:`Dev`
   * - :literal:`--lat`
     - Latitude in degrees north. Required.
     - none
   * - :literal:`--lon`
     - Longitude in degrees, either -180..180 or 0..360. Required.
     - none
   * - :literal:`--lev`
     - Model level as a 1-based index into the file's :literal:`lev` dimension (1 = surface). Ignored for 2-D collections.
     - :literal:`1`
   * - :literal:`--species`
     - One or more species names (e.g. :literal:`O3`) or full variable names.
     - :literal:`O3 HNO3 ISOP`
   * - :literal:`--collection`
     - History collection name, i.e. files named :file:`GEOSChem.<collection>.*.nc4`.
     - :literal:`SpeciesConc`
   * - :literal:`--out`
     - Output figure (:file:`.png` or :file:`.pdf`).
     - :file:`point_timeseries.png`
   * - :literal:`--csv`
     - Also write the extracted timeseries to this CSV file.
     - none
   * - :literal:`--figsize W H`
     - Figure width and height in inches, e.g. :literal:`11 8.5` for a landscape letter page.
     - scales with the number of species

.. _point-timeseries-output:

======
Output
======

The figure has one column per species:

- **Top row:** Ref (blue) and Dev (orange) concentrations. Mixing ratios (mol/mol) are displayed in ppbv, or in pptv when the values are below 1 ppbv. Other units are shown as they are.
- **Bottom row** (only with :literal:`--dev`): :literal:`(Dev - Ref) / Ref` in percent, annotated with the mean and maximum absolute percent difference.

The figure title gives the center of the selected grid box, the model level and its approximate height, and the time range plotted.

If :literal:`--csv` is given, the CSV file has a :literal:`time` column plus one column per species and run, named :literal:`<refname>_<species>_<unit>` and :literal:`<devname>_<species>_<unit>`, and a :literal:`pctdiff_<species>` column for each species when a Dev run is given.

.. _point-timeseries-notes:

=====
Notes
=====

- The variable for each species is found by trying :literal:`SpeciesConcVV_<species>`, then :literal:`SpeciesConc_<species>`, then the name exactly as given.
- The nearest grid box is chosen with longitudes compared modulo 360, so it works for grids in either -180..180 or 0..360 and across the dateline.
- A warning is printed if the requested point lies outside the selected grid box. This usually means the output is a regional subset (see :literal:`LON_RANGE` and :literal:`LAT_RANGE` in :file:`HISTORY.rc`) that does not contain the point. A warning is also printed if Ref and Dev select different grid boxes.
- The approximate height of the level is computed from the midpoint pressure (from the hybrid coefficients :literal:`hyam` and :literal:`hybm` with a surface pressure of 1013.25 hPa). If this pressure matches a level of the GEOS-Chem 72-level table (:file:`gcpy/benchmark/modules/GC_72_vertical_levels.csv`) to within 0.5%, that level's number and altitude are reported. This means the true model level is shown even for output saved with a :literal:`levels` subset. Otherwise the height is estimated from the 1976 U.S. Standard Atmosphere.
- When a Dev run is given, only the times present in both Ref and Dev are plotted.
- Only GEOS-Chem Classic output on a latitude/longitude grid is supported. GCHP cubed-sphere output is not.
