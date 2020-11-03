Regridding
==========

This page describes the regridding capabilities of GCPy. GCPy currently supports regridding of data from GEOS-Chem restarts and output NetCDF files.
Regridding is supported across any horizontal resolution and any grid type available in GEOS-Chem, including lat/lon (global or non-global), global 
standard cubed-sphere, and global stretched-grid. GCPy also supports arbitrary vertical regridding across different vertical resolutions.
   

Regridding for Plotting in GCPy
-------------------------------

When plotting in GCPy (e.g. through ``compare_single_level()`` or ``compare_zonal_mean()``), the vast majority of regridding is handled
internally. You can optionally request a specific horizontal comparison resolution in ``compare_single_level()`` and ``compare_zonal_mean()``.
Note that all regridding in these plotting functions only applies to the comparison panels (not the top two panels which show data directly from
each dataset).
There are only two scenarios where you will need to pass extra information to GCPy to help it determine grids and to regrid when plotting.

Pass stretched-grid file paths
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Stretched-grid parameters cannot currently be automatically determined from grid coordinates. If you are plotting stretched-grid data in 
``compare_single_level()`` or ``compare_zonal_mean()`` (even if regridding to another format), 
you need to use the ``sg_ref_path`` or ``sg_dev_path`` arguments to pass the path of your original stretched-grid restart file to GCPy. 
If using ``single_panel()``, pass the file path using ``sg_path``.
Stretched-grid restart files created using GCPy contain the specified stretch factor, target longitude, and target latitude in their metadata.
Currently, output files from stretched-grid runs of GCHP do not contain any metadata that specifies the stretched-grid used.

Pass vertical grid parameters for non-72/47-level grids
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GCPy automatically handles regridding between different vertical grids when plotting except when you pass a dataset that is not on 
the typical 72-level or 47-level vertical grids. If using a different vertical grid, you will need to pass the corresponding
`grid parameters <http://wiki.seas.harvard.edu/geos-chem/index.php/GEOS-Chem_vertical_grids#Reference_section_for_vertical_grids>`_ 
using the ``ref_vert_params`` or ``dev_vert_params`` keyword arguments. 



Automatic regridding decision process
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When you do not specify a horizontal comparison resolution using the ``cmpres`` argument in ``compare_single_level()`` and ``compare_zonal_mean()``,
GCPy follows several steps to determine what comparison resolution it should use:

- If both input grids are lat/lon, use the highest resolution between them (don't regrid if they are the same resolution).
- Else if one grid is lat/lon and the other is cubed-sphere (standard or stretched-grid), use a 1x1.25 lat/lon grid.
- Else if both grids are cubed-sphere and you are plotting zonal means, use a 1x1.25 lat/lon grid.
- Else if both grids are standard cubed-sphere, use the highest resolution between them (don't regrid if they are the same resolution).
- Else if one or more grids is a stretched-grid, use the grid of the ref dataset.


For differing vertical grids, the smaller vertical grid is currently used for comparisons.


Regridding Files
----------------

You can regrid existing GEOS-Chem restart or output diagnostic files between lat/lon and cubed-sphere formats using ``gcpy.file_regrid``. 
``gcpy.file_regrid`` can either be called directly from the command line using ``python -m gcpy.file_regrid`` 
or as a function (``gcpy.file_regrid.file_regrid()``) from a Python script or interpreter. The syntax of ``file_regrid`` is as follows: 

.. code-block:: python

   def file_regrid(fin, fout, dim_format_in, dim_format_out, cs_res_out=0, ll_res_out='0x0', 
   sg_params_in=[1.0, 170.0, -90.0], sg_params_out=[1.0, 170.0, -90.0]
   ):
   """
   Regrids an input file to a new horizontal grid specification and saves it
   as a new file.
   """


Required Arguments:
~~~~~~~~~~~~~~~~~~~

.. option:: fin : str

      The input filename

.. option:: fout : str

      The output filename (file will be overwritten if it already exists)

.. option:: dim_format_in : str

      Format of the input file's dimensions (choose from: classic, checkpoint, diagnostic),
      where classic denotes lat/lon and checkpoint / diagnostic are cubed-sphere formats

.. option:: dim_format_out : str

      Format of the output file's dimensions (choose from: classic, checkpoint, diagnostic),
      where classic denotes lat/lon and checkpoint / diagnostic are cubed-sphere formats

Optional arguments:
~~~~~~~~~~~~~~~~~~~

.. option:: cs_res_out : int

      The cubed-sphere resolution of the output dataset. Not used if dim_format_out is classic

      Default value: 0

.. option:: ll_res_out : str

      The lat/lon resolution of the output dataset. Not used if dim_format_out is not classic

      Default value: '0x0'

.. option:: sg_params_in : list[float, float, float]

      Input grid stretching parameters [stretch-factor, target longitude, target latitude].
      Not used if dim_format_in is classic

      Default value: [1.0, 170.0, -90.0] (No stretching)

.. option:: sg_params_out : list[float, float, float]

      Output grid stretching parameters [stretch-factor, target longitude, target latitude].
      Not used if dim_format_out is classic

      Default value: [1.0, 170.0, -90.0] (No stretching)

There are three dimension formats available for regridding: `classic` (GEOS-Chem Classic lat/lon format), `checkpoint` (GCHP restart file format),
and `diagnostic` (GCHP output file format). You can regrid between any of these formats using ``file_regrid``, as well as between different resolutions
and/or grid-types within each dimension format (e.g. standard cubed-sphere checkpoint to stretched-grid checkpoint). Note that although the ``cs_res_out``
and ``ll_res_out`` parameters are technically optional in the function, you must specify at least one of these in your call to ``file_regrid``.

As stated previously, you can either call ``file_regrid.file_regrid()`` directly or call it from the command line using ``python -m gcpy.file_regrid ARGS``.
An example command line call (separated by line for readability) for regridding a C90 cubed-sphere restart file to a C48 stretched-grid 
with a stretch factor of 3, a target longitude of 260.0, and a target latitude of 40.0 looks like:

.. code-block::

   python -m gcpy.file_regrid             \
         -i initial_GEOSChem_rst.c90_standard.nc   \
         --dim_format_in checkpoint      \
         -o sg_restart_c48_3_260_40.nc       \
         --cs_res_out 48            \
         --sg_params_out 3.0 260.0 40.0      \
         --dim_format_out checkpoint
