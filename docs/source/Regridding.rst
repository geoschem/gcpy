.. |br| raw:: html

   <br/>

.. _regrid:

##########
Regridding
##########

This page describes the regridding capabilities of GCPy. GCPy
currently supports regridding of data from GEOS-Chem restarts and
output NetCDF files. Regridding is supported across any horizontal
resolution and any grid type available in GEOS-Chem, including lat/lon
(global or non-global), global standard cubed-sphere, and global
stretched-grid. GCPy also supports arbitrary vertical regridding
across different vertical resolutions.

Regridding with GCPy is currently undergoing an overhaul. As of the current
release, regridding is split into two different categories - regridding
GEOS-Chem Classic format files (lat/lon), and regridding GCHP format files
(standard cubed-sphere, stretched cubed-sphere).

.. _regrid-classic:

====================================
Regridding Files - GEOS-Chem Classic
====================================

You can regrid existing GEOS-Chem Classic restart or output diagnostic files
between lat/lon resolutions using :code:`gcpy.file_regrid`.
:code:`gcpy.file_regrid` can either be called directly from the command line
using :code:`python -m gcpy.file_regrid` or as a function
(:code:`gcpy.file_regrid.file_regrid()`) from a Python script or interpreter.
The syntax of :code:`file_regrid` is as follows:

.. code-block:: python

   def file_regrid(fin, fout, dim_format_in, dim_format_out, ll_res_out='0x0'):
   """
   Regrids an input file to a new horizontal grid specification and saves it
   as a new file.
   """

Required Arguments:
-------------------

.. option:: fin : str

      The input filename

.. option:: fout : str

      The output filename (file will be overwritten if it already exists)

.. option:: dim_format_in : str

      Format of the input file's dimensions (set this to 'classic' - denoting
      a GEOS-Chem Classic file with a lat/lon grid)

.. option:: dim_format_out : str

      Format of the output file's dimensions (set this to 'classic' - denoting
      a GEOS-Chem Classic file with a lat/lon grid)

Optional arguments:
-------------------

.. option:: ll_res_out : str

      The lat/lon resolution of the output dataset.

      Default value: '0x0'

There is now only one grid format supported for regridding files using the
:code:`gcpy.file_regrid` method: :literal:`classic`. You must specify
:literal:`classic` as the value of both :code:`dim_format_in` and
:code:`dim_format_out`, as well as specifying a resolution as the value of
:code:`ll_res_out`.

As stated previously, you can either call
:code:`file_regrid.file_regrid()` directly or call it from the command
line using :code:`python -m gcpy.file_regrid ARGS`. An example command
line call (separated by line for readability) for regridding a 2x2.5 lat/lon
restart file to a 4x5 lat/lon grid looks like:

.. code-block::

   python -m gcpy.file_regrid                     \
         --filein initial_GEOSChem_rst.2x2.5.nc   \
         --dim_format_in classic                  \
         --fileout GEOSChem_rst.4x5.nc            \
         --ll_res_out 4x5                         \
         --dim_format_out classic

.. _regrid-gchp:

=======================
Regridding Files - GCHP
=======================

GCHP regridding is where the first steps of the overhaul in GCPy regridding have
happened. We are moving towards an integrated approach for all GEOS-Chem grid
types using `gridspec <https://github.com/liambindle/gridspec>`_ and
`sparselt <https://github.com/liambindle/sparselt>`_. For now, this is only
supported for GCHP grid formats, but in a later GCPy this will be the single
method for regridding all GEOS-Chem grid formats.

Currently, this method is only available from the command line. The syntax of
:code:`regrid_restart_file` is as follows:

Required Arguments:
-------------------

.. option:: file_to_regrid : str

      The GCHP restart file to be regridded

.. option:: regridding_weights_file : str

      Regridding weights to be used in the regridding transformation, generated
      by :literal:`ESMF_RegridWeightGen`

.. option:: template_file : str

      The GCHP restart file to use as a template for the regridded restart
      file - attributes, dimensions, and variables for the output file will be
      taken from this template. Typically this will be the same file as the file
      you are regridding!

Optional arguments:
-------------------

.. option:: --stretched-grid : switch

      A switch to indicate that the target grid is a stretched cubed-sphere grid

.. option:: --stretch-factor : float

      The grid stretching factor for the target stretched grid. Only takes
      effect when :code:`--stretched-grid` is set. See the
      `GCHP documentation <https://gchp.readthedocs.io/en/latest/supplement/stretched-grid.html#choose-stretching-parameters>`_
      for more information

.. option:: --target-latitude : float

      The latitude of the centre point for stretching the target grid. Only
      takes effect when :code:`--stretched-grid` is set. See the
      `GCHP documentation <https://gchp.readthedocs.io/en/latest/supplement/stretched-grid.html#choose-stretching-parameters>`_
      for more information

.. option:: --target-longitude : float

      The longitude of the centre point for stretching the target grid. Only
      takes effect when :code:`--stretched-grid` is set. See the
      `GCHP documentation <https://gchp.readthedocs.io/en/latest/supplement/stretched-grid.html#choose-stretching-parameters>`_
      for more information

.. _regrid-gchp-firsttime:

First Time Setup
-----------------

Until GCPy contains a complete regridding implementation that works for all
GEOS-Chem grid formats, we recommend that you create a small
`conda <https://docs.conda.io/en/latest/>`_ environment in which to carry out
your GCHP regridding.

The following conda `environment file <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file>`_
will get you set up with an environment for regridding with
:literal:`gridspec` and :literal:`sparselt`:

.. code-block:: yaml

   name: gchp_regridding
   channels:
     - conda-forge
   dependencies:
     - python=3.9
     - esmf
     - gridspec
     - numpy
     - requests
     - sparselt
     - xarray
     - xesmf

.. tip::

   For your convenience, we have placed a copy of the above
   environment file at the path
   :file:`docs/environment/gchp_regridding.yml`.

After installing and switching to this new conda environment, you should have
the :literal:`gridspec` commands available to you at the command line.

.. _regrid-gchp-procedure:

Regridding
----------

Regridding with :literal:`gridspec` and :literal:`sparselt` is a three stage
process:

#. Create grid specifications for the source and target grids using
   :literal:`gridspec`

#. Create regridding weights for the transformation using
   :literal:`ESMF_RegridWeightGen`

#. Run the regridding operation using the new :code:`regrid_restart_file`
   submodule of GCPy


Standard Cubed-Sphere Regridding
--------------------------------

We will use the example of regridding the out-of-the-box
:literal:`GEOSChem.Restart.20190701_0000z.c48.nc4` restart file from C48 to
C60 to demonstrate the standard cubed-sphere regridding process:

#. Create a source grid specification using :code:`gridspec-create`.

   .. code-block:: console

      $ gridspec-create gcs 48

   This will produce 7 files - :literal:`c48_gridspec.nc` and
   :literal:`c48.tile[1-6].nc`

#. Create a target grid specification using :code:`gridspec-create`.

   .. code-block:: console

      $ gridspec-create gcs 60

   Again, this will produce 7 files - :literal:`c60_gridspec` and
   :literal:`c60.tile[1-6].nc`

#. Create the regridding weights for the regridding transformation using
   :code:`ESMF_RegridWeightGen`.

   .. code-block:: console

      $ ESMF_RegridWeightGen            \
          --source c48_gridspec.nc      \
          --destination c60_gridspec.nc \
          --method conserve             \
          --weight c48_to_c60_weights.nc

   This will produce a log file, :literal:`PET0.RegridWeightGen.Log`, and our
   regridding weights, :literal:`c48_to_c60_weights.nc`

#. Finally, use the grid weights produced in step 3 to complete the regridding. You will need to activate your GCPy python environment for this step.

   .. code-block:: console

      $ python -m gcpy.regrid_restart_file        \
          GEOSChem.Restart.20190701_0000z.c48.nc4 \
          c48_to_c60_weights.nc                   \
          GEOSChem.Restart.20190701_0000z.c48.nc4

   This will produce a single file, :literal:`new_restart_file.nc`, regridded
   from C48 to C60, that you can rename and use as you please.

Stretched Cubed-Sphere Regridding
---------------------------------

We will use the example of regridding the out-of-the-box
:literal:`GEOSChem.Restart.20190701_0000z.c48.nc4` restart file from C48 to
a C120 base resolution stretched grid with a stretch factor of 4.0 over Bermuda
to demonstrate the stretched cubed-sphere regridding process:

#. Create a source grid specification using :code:`gridspec-create`.

   .. code-block:: console

      $ gridspec-create gcs 48

   This will produce 7 files - :literal:`c48_gridspec.nc` and
   :literal:`c48.tile[1-6].nc`

#. Create a target grid specification using :code:`gridspec-create`.

   .. code-block:: console

      $ gridspec-create sgcs 120 -s 4.0 -t 32.0 -64.0

   Here, the :code:`-s` option denotes the stretch factor and the :code:`-t`
   option denotes the latitude / longitude of the centre point of the grid
   stretch.

   Again, this will produce 7 files - :literal:`c120_..._gridspec.nc` and
   :literal:`c120_..._tile[1-6].nc`, where :literal:`...` denotes randomly
   generated characters.

#. Create the regridding weights for the regridding transformation using
   :code:`ESMF_RegridWeightGen`, replacing :literal:`c120_..._gridspec.nc`
   with the actual name of the file created in the previous step.

   .. code-block:: console

      $ ESMF_RegridWeightGen                 \
          --source c48_gridspec.nc           \
          --destination c120_..._gridspec.nc \
          --method conserve                  \
          --weight c48_to_c120_stretched_weights.nc

   This will produce a log file, :literal:`PET0.RegridWeightGen.Log`, and our
   regridding weights, :literal:`c48_to_c120_stretched_weights.nc`

#. Finally, use the grid weights produced in step 3 to complete the regridding.
   You will need to switch to your GCPy python environment for this step.

   .. code-block:: console

      $ python -m gcpy.regrid_restart_file        \
          --stretched-grid                        \
          --stretch-factor 4.0                    \
          --target-latitude 32.0                  \
          --target-longitude -64.0                \
          GEOSChem.Restart.20190701_0000z.c48.nc4 \
          c48_to_c120_stretched_weights.nc        \
          GEOSChem.Restart.20190701_0000z.c48.nc4

   This will produce a single file, :literal:`new_restart_file.nc`, regridded
   from C48 to C120, with a stretch factor of 4.0 over 32.0N, -64.0E, that you
   can rename and use as you please. It is generally a good idea to rename the
   file to include the grid resolution, stretch factor, and target lat/lon for
   easy reference.

   .. code-block:: console

      $ mv new_restart_file.nc GEOSChem.Restart.20190701_0000z.c120.s4_32N_64E.nc

.. _regrid-plot:

===============================
Regridding for Plotting in GCPy
===============================

When plotting in GCPy (e.g. through :code:`compare_single_level()` or
:code:`compare_zonal_mean()`), the vast majority of regridding is
handled internally. You can optionally request a specific
horizontal comparison resolution in :code:`compare_single_level()``
and :code:`compare_zonal_mean()`.  Note that all regridding in these
plotting functions only applies to the comparison panels (not the top
two panels which show data directly from each dataset). There are only
two scenarios where you will need to pass extra information to GCPy to
help it determine grids and to regrid when plotting.

Pass stretched-grid file paths
------------------------------

Stretched-grid parameters cannot currently be automatically determined
from grid coordinates. If you are plotting stretched-grid data in
:code:`compare_single_level()` or :code:`compare_zonal_mean()` (even
if regridding to another format), you need to use the
:code:`sg_ref_path` or :code:`sg_dev_path` arguments to pass the path
of your original stretched-grid restart file to GCPy.
If using :code:`single_panel()`, pass the file path using
:code:`sg_path`. Stretched-grid restart files created using GCPy
contain the specified stretch factor, target longitude, and
target latitude in their metadata.  Currently, output files from
stretched-grid runs of GCHP do not contain any metadata that specifies
the stretched-grid used.

Pass vertical grid parameters for non-72/47-level grids
-------------------------------------------------------

GCPy automatically handles regridding between different vertical grids
when plotting except when you pass a dataset that is not on the
typical 72-level or 47-level vertical grids. If using a different
vertical grid, you will need to pass the corresponding `grid
parameters
<http://wiki.seas.harvard.edu/geos-chem/index.php/GEOS-Chem_vertical_grids#Reference_section_for_vertical_grids>`_
using the :code:`ref_vert_params` or :code:`dev_vert_params` keyword
arguments.

Automatic regridding decision process
-------------------------------------

When you do not specify a horizontal comparison resolution using the
:code:`cmpres` argument in :code:`compare_single_level()` and
:code:`compare_zonal_mean()`, GCPy follows several steps to determine
what comparison resolution it should use:

- If both input grids are lat/lon, use the highest resolution between
  them (don't regrid if they are the same resolution).
- Else if one grid is lat/lon and the other is cubed-sphere (standard
  or stretched-grid), use a 1x1.25 lat/lon grid.
- Else if both grids are cubed-sphere and you are plotting zonal
  means, use a 1x1.25 lat/lon grid.
- Else if both grids are standard cubed-sphere, use the highest
  resolution between them (don't regrid if they are the same
  resolution).
- Else if one or more grids is a stretched-grid, use the grid of the
  ref dataset.

For differing vertical grids, the smaller vertical grid is currently
used for comparisons.
