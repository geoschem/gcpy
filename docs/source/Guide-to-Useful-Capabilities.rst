Overview of Capabilities
============================

This page outlines the capabilities of GCPy with links to detailed
function documentation.



Spatial Plotting
----------------

One hallmark of GCPy is easy-to-use spatial plotting of GEOS-Chem data.
Available plotting falls into two layouts: single panel (one map of one
variable from a dataset) and six panel (six maps comparing a variable
between two datasets). The maps in these plots can display data at a
single vertical level of your input dataset or in a zonal mean for all
layers of the atmosphere.

Single Panel Plots
~~~~~~~~~~~~~~~~~~

Single panel plots are generated through the ``plot.single_panel()``
function. ``plot.single_panel()`` uses Matplotlib and Cartopy plotting
capabilities while handling certain behind the scenes operations that
are necessary for plotting GEOS-Chem data, particularly for cubed-sphere
and/or zonal mean data.

.. code:: python

    import xarray as xr
    import gcpy.plot as gcplot
    import matplotlib.pyplot as plt
    ds = xr.open_dataset('GEOSChem.Restart.20160701_0000z.nc4')
    #plot surface Ozone over the North Pacific
    gcplot.single_panel(ds['SpeciesRst_O3'].isel(lev=0), title='Surface Ozone over the North Pacific', extent=[80, -90, -10, 60])
    plt.show()


.. image:: images/single\_panel\_single\_level.png

.. code:: python

    #plot global zonal mean of Ozone
    gcplot.single_panel(ds['SpeciesRst_O3'], plot_type='zonal_mean', title='Global Zonal Mean of Ozone')
    plt.show()

.. image:: images/single\_panel\_zonal\_mean.png

`Click here <Single_panel.html>`__ for an example single panel plotting script.
`Click here <Plotting.html#single-panel>`__ for detailed documentation for ``single_panel()``.

Six Panel Comparison Plots
~~~~~~~~~~~~~~~~~~~~~~~~~~

Six panel plots are used to compare results across two different model
runs. Single level and zonal mean plotting options are both available.
The two model runs do not need to be the same resolution or even the
same grid type (GEOS-Chem Classic and GCHP output can be mixed at will).

.. code:: python

    import xarray as xr
    import gcpy.plot as gcplot
    import matplotlib.pyplot as plt
    gcc_ds = xr.open_dataset('GEOSChem.SpeciesConc.20160701_0000z.nc4')
    gchp_ds = xr.open_dataset('GCHP.SpeciesConc.20160716_1200z.nc4')
    #Plot comparison of surface ozone over the North Pacific
    gcplot.compare_single_level(gcc_ds, 'GEOS-Chem Classic', gchp_ds, 'GCHP', varlist=['SpeciesConc_O3'], extra_title_txt='Surface')
    plt.show()


.. image:: images/six\_panel\_single\_level.png

.. code:: python

    #Plot comparison of global zonal mean ozone
    gcplot.compare_zonal_mean(gcc_ds, 'GEOS-Chem Classic', gchp_ds, 'GCHP', varlist=['SpeciesConc_O3'])
    plt.show()

.. image:: images/six\_panel\_zonal\_mean.png

`Click here <Six_panel.html>`__ for an example six panel plotting script.
`Click here <Plotting.html#compare-single-level-and-compare-zonal-mean>`__
for complete documentation for ``compare_single_level()`` and ``compare_zonal_mean()``.

Comprehensive Benchmark Plotting
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The GEOS-Chem Support Team uses comprehensive plotting functions from
``benchmark.py`` to generate full plots of benchmark diagnostics.
Functions like ``benchmark.make_benchmark_conc_plots`` by default create
plots for every variable in a given collection (e.g. ``SpeciesConc``) at
multiple vertical levels (surface, 500hPa, zonal mean) and divide plots
into separate folders based on category (e.g. Chlorine, Aerosols). The
GCST uses full benchmark plotting / table scripts similar to `this example <benchmark_plotting.html>`__ 
to produce plots and tables for official model benchmarks. Full documentation for the
benchmark plotting functions can be found
`here <Plotting.html#benchmark-plotting-functions>`__.

Table Creation
--------------

GCPy has several dedicated functions for tabling GEOS-Chem output data
in text file format. These functions and their outputs are primarily
used for model benchmarking purposes.

Budget Tables
~~~~~~~~~~~~~

Currently, budget tables can be created for "operations" (table shows
change in mass after each category of model operation, as contained in
the GEOS-Chem ``Budget`` diagnostics) or in overall averages for
different aerosols or the Transport Tracers simulation.

Operations budget tables are created using the
``benchmark.make_benchmark_operations_budget`` function and appear as
follows:

.. image:: images/budget\_table.png

Full documentation for operations budget table creation can be found
`here <Tabling.html#operations-budget-tables>`__.

Mass Tables
~~~~~~~~~~~

The ``benchmark.make_benchmark_mass_tables`` function uses species
concentrations and info from meteorology files to generate the total
mass of species in certain segments of the atmosphere (currently global or only the
troposphere). An example table is shown below:

.. image:: images/mass\_table.png

Full documentation for mass table creation can be found
`here <Tabling.html#mass-tables>`__.

Emissions Tables
~~~~~~~~~~~~~~~~

The ``benchmark.make_benchmark_emis_tables`` function creates tables of
total emissions categorized by species or by inventory. Examples of both
emissions table types are shown below:


.. image:: images/emissions\_totals.png 

.. image:: images/inventory\_totals.png

Full documentation for emissions table creation can be found `here <Tabling.html#emissions-tables>`__.

Regridding
----------

General Regridding Rules
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GCPy supports regridding between all horizontal GEOS-Chem grid types, including
latitude/longitude grids (the grid format of GEOS-Chem
Classic), standard cubed-sphere (the standard grid format of GCHP), and stretched-grid
(an optional grid format in GCHP). GCPy contains several horizontal regridding functions
built off of xESMF. GCPy automatically handles most regridding needs when plotting GEOS-Chem data.

``gcpy.file_regrid`` allows you to regrid NetCDF files between different grid types / resolutions and can be called from the command line or as a function.

The 72-level and 47-level vertical grids are pre-defined in GCPy. Other vertical grids can also
be defined if you provide `the A and B coefficients of the hybrid vertical grid <wiki.seas.harvard.edu/geos-chem/index.php/GEOS-Chem_vertical_grids>`__.

When plotting data of differing grid types or horizontal resolutions using
``compare_single_level`` or ``compare_zonal_mean``, you can specify a comparison resolution using the ``cmpres`` argument.
This resolution will be used for the difference panels in each plot (the bottom four panels rather than the top two raw data panels).
If you do not specify a comparison resolution, GCPy will automatically choose one.

For more extensive regridding information, visit the `detailed regridding documentation <Regridding.html>`__.
