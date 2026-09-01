.. _plot:

########
Plotting
########

This page describes in depth the general plotting capabilities of GCPy,
including possible argument values for every plotting function.

For information about GCPy functions that are specific to the
GEOS-Chem benchmark workflow, please see our :ref:`Benchmarking <bmk>`
chapter.

.. _plot-six-panel:

==========================
Six-panel comparison plots
==========================

Functions :func:`gcpy.plot.compare_single_level` and
:func:`gcpy.plot.compare_zonal_mean` generate six-panel plots
comparing each variable between two datasets. These plots can either
be saved to PDFs or generated sequentially for visualization in the
Matplotlib GUI using :func:`matplotlib.pyplot.show`. Each plot uses
data passed from a reference (:literal:`Ref`) dataset and a
development (:literal:`Dev`) dataset.  Both functions share
significant structural overlap both in output  appearance and code
implementation.

You can import these routines into your code with these statements:

.. code-block:: python

   from gcpy.plot.compare_single_level import compare_single_level
   from gcpy.plot.compare_zonal_mean import compare_zonal_mean

Each panel has a title describing the type of panel, a colorbar for
the values plotted in that panel, and the units of the data plotted in
that panel. The upper two panels of each plot show actual values from
the :literal:`Ref` (left) and :literal:`Dev` (right) datasets for a
given variable. The middle two panels show the difference
(:literal:`Dev - Ref`) between the values in the :literal:`Dev`
dataset and the values in the :literal:`Ref` dataset. The left middle
panel uses a full dynamic color map, while the right middle panel caps
the color map at the 5th and 95th percentiles.  The bottom two panels
show the ratio (:literal:`Dev/Ref`) between the values in the Dev
dataset and the values in the Ref Dataset. The left bottom panel uses
a full dynamic color map, while the right bottom panel caps the color
map at 0.5 and 2.0.

.. _plot-flat-colorbars:

Panels with nothing to show
---------------------------

When there is no meaningful structure to plot in a panel, GCPy
replaces its color scale with a flat one and its colorbar with a
single label explaining why.  This is deliberate: without it, a
colorbar would be stretched across nothing but numerical noise, which
renders as spurious color "striping" (see `GitHub issue #330
<https://github.com/geoschem/gcpy/issues/330>`_).  The labels are:

.. list-table::
   :header-rows: 1

   * - Colorbar label
     - Meaning
   * - :literal:`Zero throughout domain`
     - Every value in the panel is exactly zero.
   * - :literal:`Undefined throughout domain`
     - Every value in the panel is :literal:`NaN`.
   * - :literal:`Differences negligible throughout domain`
     - The :literal:`Dev - Ref` differences in this panel are
       negligible compared to the magnitude of the :literal:`Ref` and
       :literal:`Dev` data themselves, so they are numerical noise
       (e.g. floating-point round-off or regridding error) rather
       than a real signal.
   * - :literal:`Ref and Dev equal throughout domain`
     - The :literal:`Dev/Ref` ratio is 1 everywhere.
   * - :literal:`Ref is zero throughout domain`
     - Shown on a ratio panel when :literal:`Ref` is zero everywhere
       but :literal:`Dev` is not.  :literal:`Dev/Ref` is then a
       division by zero at every point, so no ratio can be drawn even
       though the difference panels above may show a real change.
   * - :literal:`Dev is zero throughout domain`
     - Shown on a ratio panel when :literal:`Dev` is zero everywhere
       but :literal:`Ref` is not, so the ratio is zero everywhere.
   * - :literal:`Zero within the 5th-95th percentile range`
     - Shown on a restricted-range difference panel when the field is
       zero over most of the domain (e.g. aircraft emissions), so that
       its 5th and 95th percentiles are both zero.  The dynamic-range
       panel beside it may still show real differences.

Gray cells in a ratio panel mark places where no meaningful
:literal:`Dev/Ref` ratio exists: either :literal:`Ref` is zero there,
or both :literal:`Ref` and :literal:`Dev` are negligible compared to
the magnitude of the field.  The latter matters for fields that are
zero over much of the domain, such as aircraft emissions: regridding
leaves a tiny residue rather than an exact zero, and dividing one
residue by another yields an arbitrary ratio that would otherwise
saturate the color scale.  A cell in which :literal:`Ref` is
negligible but :literal:`Dev` is not represents a real change, so its
ratio is kept.

The negligible-difference threshold is relative, not absolute: a
difference is suppressed only if it is smaller than
:code:`gcpy.plot.core.NOISE_REL_TOL` (currently
:math:`1 \times 10^{-5}`) times the magnitude of the data being
differenced, i.e. unless Ref and Dev agree to better than about
5 parts per million.  A field with very small absolute values, such as an
aircraft emissions flux of order :math:`10^{-13}` kg m\ :sup:`-2`
s\ :sup:`-1`, will therefore still have its real differences plotted.

.. _plot-csl:

Function :code:`compare_single_level`
-------------------------------------

This function generates a comparison plot such as:

.. image:: _static/images/six\_panel\_single\_level.png
   :align: center

For a list of input parameters, click on this link:
:func:`gcpy.plot.compare_single_level`.

.. _plot-czm:

Function :code:`compare_zonal_mean`
-----------------------------------

This function generates a comparison plot such as:

.. image:: _static/images/six\_panel\_zonal\_mean.png
   :align: center

For a list of input parameters, click on this link:
:func:`gcpy.plot.compare_single_level`.

.. _plot-shared:

Shared structure
----------------

Both :func:`gcpy.plot.compare_single_level` and
:func:`gcpy.plot.compare_zonal_mean` have four positional (required)
arguments.

.. option:: refdata <xarray.Dataset>

   Dataset used as reference in comparison

.. option:: refstr <str> | <list of str>

   String description for reference data to be used in plots OR list
   containing [ref1str, ref2str] for diff-of-diffs plots

.. option:: devdata : xarray.Dataset

   Dataset used as development in comparison

.. option:: devstr <str> | <list of str>

   String description for development data to be used in plots
   OR list containing [dev1str, dev2str] for diff-of-diffs plots

:option:`refstr` and :option:`devstr` title the top two panels of
each six panel plot.

Functions :func:`gcpy.plot.compare_single_level` and
:func:`gcpy.plot.compare_zonal_mean` share many arguments. Some of
these arguments are plotting options that change the format of the plots:

For example, you may wish to convert units to :math:`\mu`\ g/m\ :sup:`3` when
generating comparison plots of aerosol species.  Activate this option
by setting the keyword argument :literal:`convert_to_ugm3=True`.

Other arguments are necessary to achieve a correct plot depending on
the format of :literal:`refdata` and :literal:`devdata` and require
you to know certain traits of your input data. For example, you must
specify if one of the datasets should be flipped vertically if Z
coordinates in that dataset do not denote decreasing pressure as Z
index increases, otherwise the vertical coordinates between your two
datasets may be misaligned and result in an undesired plotting
outcome.  This may be done with by setting the boolean options
:literal:`flip_ref=True` and/or :literal:`flip_dev=True`.

The :literal:`n_job` argument governs the parallel plotting settings
of :func:`gcpy.plot.compare_single_level` and
:func:`gcpy.plot.compare_zonal_mean`. GCPy uses the JobLib library to
create plots in parallel. Due to limitations with matplotlib, this
parallelization creates plots (pages) in parallel rather than
individual panels on a single page. Parallel plot creation is not
enabled when you do not save to a PDF. The default value of
:literal:`n_job=-1` allows the function call to automatically scale up
to, at most, the number of cores available on your system.

.. note::

   On systems with higher (12+) core counts, the maximum number of
   cores is not typically reached because of the process handling
   mechanics of JobLib. However, on lower-end systems with lower core
   counts or less available memory, it is advantageous to use
   :literal:`n_job` to limit the max number of processes.

   Due to how Python handles memory management on Linux systems, using
   more cores may result in memory not returned to the system after
   the plots are created.  Requesting fewer cores with
   :literal:`n_job` may help to avoid this situation.

.. _plot-six-panel-example:

Example script
--------------

Here is a basic script that calls both :func:`gcpy.plot.compare_zonal_mean` and
:func:`gcpy.plot.compare_single_level`:

.. code-block:: python

   #!/usr/bin/env python

   import xarray as xr
   import matplotlib.pyplot as plt
   from gcpy.plot.compare_single_level import compare_single_level
   from gcpy.plot.compare_zonal_mean import compare_zonal_mean

   file1 = '/path/to/ref'
   file2 = '/path/to/dev'
   ds1 = xr.open_dataset(file1)
   ds2 = xr.open_dataset(file2)
   compare_zonal_mean(ds1, 'Ref run', ds2, 'Dev run')
   plt.show()
   compare_single_level(ds1, 'Ref run', ds2, 'Dev run')
   plt.show()

.. _plot-single-panel:

==================
Single panel plots
==================

Function :func:`gcpy.plot.single_panel` is used to create plots
containing only one panel of GEOS-Chem data.  This function is used
within :func:`gcpy.plot.compare_single_level` and
:func:`gcpy.plot.compare_zonal_mean` to generate each panel plot. It
can also be called directly on its own to quickly plot GEOS-Chem data
in zonal mean or single level format.

The :func:`gcpy.plot.single_panel` function expects data with a 1-length
(or non-existent) :literal:`T` (time) dimension, as well as a 1-length
or non-existent :literal:`Z` (vertical level) dimension.  It also
contains a few amenities to help with plotting GEOS-Chem data,
including automatic grid detection for lat/lon or standard
cubed-sphere :func:`xarray.DataArray`-s. You can also pass NumPy
arrays to plot, though you'll need to manually pass grid info in this
case (with the :literal:`gridtype`, :literal:`pedge`, and
:literal:`pedge_ind` keyword arguments).

The sample script shown below shows how you can data at a single level
and timestep from an :func:`xarray.DataArray` object.

.. code-block:: python

   #!/usr/bin/env python

   import xarray as xr
   import matplotlib.pyplot as plt
   from gcpy.plot.single_panel import single_panel

   # Read data from a file into an xr.Dataset object
   dset = xr.open_dataset('GEOSChem.SpeciesConc.20160701_0000z.nc4')

   # Extract ozone (v/v) from the xr.Dataset object,
   # for time=0 (aka first timestep) and lev=0 (aka surface)
   sfc_o3 = dset['SpeciesConcVV_O3'].isel(time=0).isel(lev=0)

   # Plot the data!
   single_panel(sfc_o3)
   plt.show()
