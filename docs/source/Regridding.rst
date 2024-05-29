.. |br| raw:: html

   <br/>

.. _regrid:

##########
Regridding
##########

:program:`GCPy` currently supports regridding of data from:

#. GEOS-Chem Classic restart files
#. GEOS-Chem Classic diagnostic files
#. GCHP restart files
#. GCHP diagnostic files
#. HEMCO restart files
#. HEMCO diagnostic files
#. As well as any netCDF file adhering to `COARDS
   <https://ferret.pmel.noaa.gov/Ferret/documentation/coards-netcdf-conventions>`_
   or `CF <https://cfconventions.org/>`_  conventions.

Regridding is supported across any horizontal resolution and any grid
type available in GEOS-Chem, including lat/lon (global or non-global),
global standard cubed-sphere, and global stretched-grid. GCPy also
supports arbitrary vertical regridding across different vertical
resolutions.

Regridding with GCPy is currently undergoing an overhaul. As of the
current release, regridding is split into two different
categories:

#. Regridding between lat-lon grids using regridding weights computed
   on the fly by GCPy, and
#. Regridding either lat-lon or cubed-sphere using regridding weights
   computed as a preprocessing step.

The latter method may be used for creating GCHP standard grid
and stretched grid restart files from either GCHP or GEOS-Chem Classic
restart files.

.. _regrid-classic:

===============================
Using Online Regridding Weights
===============================

You can regrid existing GEOS-Chem restart or diagnostic files using
GCPy function :code:`gcpy.file_regrid`. This function can called
directly from the command line (:ref:`see the examples below
<regrid-classic-example>`) or from a Python script or
interpreter (:code:`gcpy.file_regrid.file_regrid()`)

.. note::

   For regridding to or from GCHP stretched-grid restart files, we
   recommend using the :ref:`offline regridding weights method
   <regrid-gchp>`.

The syntax of :code:`file_regrid` is as follows:

.. code-block:: python

   def file_regrid(
           fin,
           fout,
           dim_format_in,
           dim_format_out,
           cs_res_out=0,
           ll_res_out='0x0',
           sg_params_in=None,
           sg_params_out=None,
           vert_params_out=None,
   ):
       """
       Regrids an input file to a new horizontal grid specification
       and saves it as a new file.
       """

gcpy.file_regrid required arguments:
------------------------------------

.. option:: fin : str

      The input filename

.. option:: fout : str

      The output filename (file will be overwritten if it already exists)

.. option:: dim_format_in : str

      Format of the input file's dimensions.  Accepted values are:

      - :literal:`classic`: For GEOS-Chem Classic restart & diagnostic files
      - :literal:`checkpoint` : For GCHP checkpoint & restart files
      - :literal:`diagnostic`: For GCHP diagnostic files

.. option:: dim_format_out : str

      Format of the output file's dimensions.  Accepted values are:

      - :literal:`classic`: For GEOS-Chem Classic restart & diagnostic files
      - :literal:`checkpoint` : For GCHP checkpoint & restart files
      - :literal:`diagnostic`: For GCHP diagnostic files

gcpy.file_regrid optional arguments:
------------------------------------

.. option:: sg_params_in : list of float

      Stretching parameters (:literal:`stretch-factor`,
      :literal:`target-longitude`, :literal:`target-latitude`) for the
      input grid.  Only needed when the data contained in file
      :option:`fin` is on a GCHP stretched grid.

      Default value: :literal:`[1.0, 170.0, -90.0]`

.. option:: sg_params_out : list of float

      Stretching parameters (:literal:`stretch-factor`,
      :literal:`target-longitude`, :literal:`target-latitude`) for the
      output  grid.  Only needed when the data to be contained in file
      :option:`fout` is to be placed on a GCHP stretched grid.

      Default value: :literal:`[1.0, 170.0, -90.0]`

.. option:: cs_res_out : int

      Cubed-sphere resolution of the output dataset.  Only needed when
      the data in file :option:`fin` is on a GCHP cubed-sphere grid.

      Default value: :code:`0`

.. option:: ll_res_out : str

      The lat/lon resolution of the output dataset.  Only needed when
      the data to be contained in file :option:`fout` is to be placed
      on a GEOS-Chem Classic lat-lon grid.

      Default value: :code:`"0x0"`.

.. option:: vert_params_out : list of float

      Hybrid grid parameter :math:`A` (in :literal:`hPa` and :math:`B`
      (:literal:`unitless`), returned in list format: :code:`[A, B]`

      Default value: :code:`None`

.. _regrid-classic-example:

Examples
--------

As stated previously, you can call
:code:`gcpy.file_regrid.file_regrid()` from a Python script, or from
the command line.  Here we shall focus on command-line examples.

#. Regrid a 4x5 GEOS-Chem Classic restart or diagnostic file to a
   GEOS-Chem Classic 2x2.5 file:

   .. code-block::

      $ python -m gcpy.file_regrid              \
        --filein         /path/to/file_4x5.nc4  \
        --dim_format_in  classic                \
        --fileout        /path/to/file_2x25.nc4 \
        --ll_res_out     2x2.5                  \
        --dim_format_out classic

   |br|

#. Regrid a 4x5 GEOS-Chem Classic restart or diagnostic file to a
   GCHP C24 restart file:

   .. code-block::

      $ python -m gcpy.file_regrid              \
        --filein         /path/to/file_4x5.nc4  \
        --dim_format_in  classic                \
        --fileout        /path/to/file_c24.nc4  \
        --cs_res_out     24                     \
        --dim_format_out checkpoint

   |br|

#. Regrid a GCHP C48 restart file to a GCHP stretched grid C48 restart
   file.  The stretch parameters are:

   - stretch-factor: 5
   - target-longitude: -72
   - target-latitude: 41

   .. code-block::

      $ python -m gcpy.file_regrid                 \
        --filein         /path/to/file_c48.nc4     \
        --dim_format_in  checkpoint                \
        --fileout        /path/to/file_c48_sg.nc4  \
        --cs_res_out     48                        \
        --dim_format_out checkpoint                \
	--sg_params_out  5 -72 41

   |br|

#. Regrid the GCHP stretched grid C48 restart file from Example 3
   above to a GCHP C24 diagnostic file.

   .. code-block::

      $ python -m gcpy.file_regrid                 \
        --filein         /path/to/file_c48_sg.nc4  \
	--sg_params_in   5 -72 41                  \
        --dim_format_in  checkpoint                \
        --fileout        /path/to/file_c24.nc4     \
        --cs_res_out     24                        \
        --dim_format_out diagnostic

.. _regrid-gchp:

================================
Using Offline Regridding Weights
================================

This approach requires generating regridding weights using python
packages `gridspec <https://github.com/liambindle/gridspec>`_ and
`sparselt <https://github.com/liambindle/sparselt>`_. Regridding with
:literal:`GCPy`, :literal:`gridspec` and :literal:`sparselt` is a
three stage process:

#. Create grid specifications for the source and target grids using
   :literal:`gridspec`.
#. Create regridding weights for the transformation using
   :literal:`ESMF_RegridWeightGen`.
#. Run the regridding operation using the :code:`regrid_restart_file`
   submodule of GCPy.

.. note::

   As of GCPy 1.4.0, the :ref:`default GCPy environment
   <install>` (aka :literal:`gcpy_env`) now contains
   :literal:`gridspec` and :literal:`sparselt` packages.  You no
   longer need to use the separate :literal:`gchp_regridding`
   environment as in prior versions.

.. _regrid-gchp-args:

gcpy.regrid_restart_file required arguments:
--------------------------------------------

There are three arguments required by the GCPy function
:literal:`regrid_restart_file`:

.. option:: file_to_regrid : str

      The GEOS-Chem Classic or GCHP data file to be regridded.

.. option:: regridding_weights_file : str

      Regridding weights to be used in the regridding transformation,
      generated by :literal:`ESMF_RegridWeightGen`

.. option:: template_file : str

      The GC-Classic or GCHP restart file to use as a template for the
      regridded restart file. Attributes, dimensions, and variables
      for the output file will be taken from this template.

gcpy.regrid_restart_file optional arguments:
--------------------------------------------

There are four optional arguments, all of which are for regridded to a
stretched cubed-sphere grid.

.. option:: --stretched-grid : switch

      A switch to indicate that the target grid is a stretched
      cubed-sphere grid.

.. option:: --stretch-factor : float

      The grid stretching factor for the target stretched grid. Only
      takes  effect when :code:`--stretched-grid` is set. See the
      `GCHP documentation
      <https://gchp.readthedocs.io/en/latest/supplement/stretched-grid.html#choose-stretching-parameters>`_
      for more information. Make sure this value exactly matches the
      value you plan to use in GCHP configuration file
      :file:`setCommonRunSettings.sh`.

.. option:: --target-latitude : float

      The latitude of the centre point for stretching the target
      grid. Only takes effect when :code:`--stretched-grid` is
      set. See the `GCHP documentation
      <https://gchp.readthedocs.io/en/latest/supplement/stretched-grid.html#choose-stretching-parameters>`_
      for more information. Make sure this value exactly matches the
      value you plan to use in GCHP configuration file
      :file:`setCommonRunSettings.sh`.

.. option:: --target-longitude : float

      The longitude of the centre point for stretching the target
      grid. Only takes effect when :code:`--stretched-grid` is
      set. See the `GCHP documentation <https://gchp.readthedocs.io/en/latest/supplement/stretched-grid.html#choose-stretching-parameters>`_
      for more information. Make sure this value exactly matches the
      value you plan to use in GCHP configuration file
      :file:`setCommonRunSettings.sh`.

.. _regrid-gchp-example-1:

Example 1: Standard Lat-Lon to Cubed-Sphere Regridding
------------------------------------------------------

This example will show regridding a GC-Classic 4x5 restart file to a
GCHP c24 restart file.

#. Activate your GCPy environment.

   .. code-block:: console

      $ mamba activate gcpy_env  # Or whatever your environment's name is

   |br|

#. Create a lat-lon source grid specification using
   :code:`gridspec-create`.

   .. code-block:: console

      $ gridspec-create latlon --pole-centered --half-polar 46 72

   This will produce 1 file: :file:`regular_lat_lon_46x72.nc`. |br|
   |br|

#. Create a target grid specification using :code:`gridspec-create`.

   .. code-block:: console

      $ gridspec-create gcs 24

   This will produce 7 files: :file:`c24_gridspec.nc` and
   :file:`c24.tile[1-6].nc` |br|
   |br|

#. Create the regridding weights for the regridding transformation
   (46x72 to C24) using :code:`ESMF_RegridWeightGen`.

   .. code-block:: console

      $ ESMF_RegridWeightGen                   \
        --source      regular_lat_lon_46x72.nc \
        --destination c24_gridspec.nc          \
        --method      conserve                 \
        --weight      46x72_to_c24_weights.nc

   This will produce a log file, :file:`PET0.RegridWeightGen.Log`, and our
   regridding weights, :file:`46x72_to_c24_weights.nc` |br|
   |br|

#. Use the grid weights produced in previous steps to complete the
   regridding.

   .. code-block:: console

      $ python -m gcpy.regrid_restart_file          \
        GEOSChem.Restart.20190701_0000z.nc4         \
        46x72_to_c24_weights.nc                     \
        GEOSChem.Restart.20190701_0000z.c24_old.nc4

   The arguments to :code:`gcpy.regrid_restart_file`
   :ref:`are described above <regrid-gchp-args>`.  In this example
   (lat-lon to cubed-sphere) we need to use a GEOS-Chem Classic
   restart file as the file to be regridded and a GCHP restart file as
   the template file.

   .. note::

      The resolution of the template file does not matter as long as it
      contains all of the variables and attributes that you wish to
      include in the regridded restart file.

   After running :code:`gcpy.regrid_restart_file`, a single restart file
   named :file:`new_restart_file.nc` will be created.  You can rename
   this file and use it to initialize your GCHP C24 simulation. |br|
   |br|

#. Deactivate your GCPy environment when finished.

   .. code-block:: console

      $ mamba deactivate

.. _regrid-gchp-example-2:

Example 2: Standard Cubed-Sphere to Cubed-Sphere Regridding
-----------------------------------------------------------

We will use the example of regridding the out-of-the-box
:file:`GEOSChem.Restart.20190701_0000z.c48.nc4` restart file from
C48 to C60 to demonstrate the standard cubed-sphere regridding process:

#. Activate your GCPy environment.

   .. code-block:: console

      $ mamba activate gcpy_env  # Or whatever your environment's name is

   |br|

#. Create a source grid specification using :code:`gridspec-create`.

   .. code-block:: console

      $ gridspec-create gcs 48

   This will produce 7 files: :literal:`c48_gridspec.nc` and
   :literal:`c48.tile[1-6].nc` |br|
   |br|

#. Create a target grid specification using :code:`gridspec-create`.

   .. code-block:: console

      $ gridspec-create gcs 60

   Again, this will produce 7 files: :literal:`c60_gridspec.nc` and
   :literal:`c60.tile[1-6].nc` |br|
   |br|

#. Create the regridding weights for the regridding transformation
   (C48 to C60) using :code:`ESMF_RegridWeightGen`.

   .. code-block:: console

      $ ESMF_RegridWeightGen               \
        --source      c48_gridspec.nc      \
        --destination c60_gridspec.nc      \
        --method      conserve             \
        --weight      c48_to_c60_weights.nc

   This will produce a log file, :file:`PET0.RegridWeightGen.Log`,
   and our regridding weights, :file:`c48_to_c60_weights.nc` |br|
   |br|

#. Use the grid weights produced in earlier steps to complete the regridding.

   .. code-block:: console

      $ python -m gcpy.regrid_restart_file      \
        GEOSChem.Restart.20190701_0000z.c48.nc4 \
        c48_to_c60_weights.nc                   \
        GEOSChem.Restart.20190701_0000z.c48.nc4

   The arguments to :code:`gcpy.regrid_restart_file`
   :ref:`are described above <regrid-gchp-args>`.  Because we are
   regridding from one cubed-sphere grid to another cubed-sphere grid,
   we can use the file to be regridded as the template file.

   After running :code:`gcpy.regrid_restart_file`, a single restart
   file named :file:`new_restart_file.nc` will be created.  You can
   rename this file as you wish and use it for your GCHP C60
   simulation. |br|
   |br|

#. Deactivate your GCPy environment when you have finished.

   .. code-block:: console

      $ mamba deactivate

Example 3: Standard to Stretched Cubed-Sphere Regridding
--------------------------------------------------------

This example regrids the out-of-the-box c48 restart file
(:file:`GEOSChem.Restart.20190701_0000z.c48.nc4`) from a standard
cubed-sphere grid to a stretched grid. The base resolution will remain
the same at c48. The regridded file will have a stretch factor of 4.0
over Bermuda which means a regional grid resolution of c196 (4
times 48) in that area.

#. Activate your GCPy environment:

   .. code-block:: console

      $ mamba activate gcpy_env  # Or whatever your environment's name is

   |br|

#. Create a source grid specification using :code:`gridspec-create`.

   .. code-block:: console

      $ gridspec-create gcs 48

   This will produce 7 files: :file:`c48_gridspec.nc` and
   :file:`c48.tile[1-6].nc` |br|
   |br|

#. Create a target grid specification using :code:`gridspec-create`.
   This will be for the stretched grid.

   .. code-block:: console

      $ gridspec-create sgcs 48 -s 4.0 -t 32.0 -64.0

   Here, the :code:`-s` option denotes the stretch factor and the
   :code:`-t` option denotes the latitude / longitude of the centre
   point of the grid stretch.

   Again, this will produce 7 files: :file:`c48_..._gridspec.nc` and
   :file:`c48_..._tile[1-6].nc`, where :file:`...` denotes randomly
   generated characters. Be sure to look for these since you will need
   them in the next step. |br|
   |br|

#. Create the regridding weights for the regridding transformation
   (C48 to C48-stretched) using :code:`ESMF_RegridWeightGen`,
   replacing  :file:`c48_..._gridspec.nc` with the actual name of the
   file created in the previous step. An example is shown below.

   .. code-block:: console

      $ ESMF_RegridWeightGen                              \
        --source      c48_gridspec.nc                     \
        --destination c48_s4d00_tdtdqp9ktebm5_gridspec.nc \
        --method      conserve                            \
        --weight      c48_to_c48_stretched_weights.nc

   This will produce a log file, :file:`PET0.RegridWeightGen.Log`, and our
   regridding weights, :file:`c48_to_c48_stretched_weights.nc` |br|
   |br|

#. Use the grid weights produced in earlier steps to complete the
   regridding.

   .. code-block:: console

      $ python -m gcpy.regrid_restart_file       \
         --stretched-grid                        \
         --stretch-factor 4.0                    \
         --target-latitude 32.0                  \
         --target-longitude -64.0                \
         GEOSChem.Restart.20190701_0000z.c48.nc4 \
         c48_to_c48_stretched_weights.nc         \
         GEOSChem.Restart.20190701_0000z.c48.nc4

   The arguments to :code:`gcpy.regrid_restart_file`
   :ref:`are described above <regrid-gchp-args>`.  Because we are
   regridding from one cubed-sphere grid to another cubed-sphere grid,
   we can use the file to be regridded as the template file.

   This will produce a single file, :literal:`new_restart_file.nc`,
   regridded from C48 standard to C48 stretched with a stretch factor
   of 4.0 over 32.0N, -64.0E, that you can rename and use as you
   please.

   .. tip::

      It is generally a good idea to rename the file to include
      the grid  resolution, stretch factor, and target lat/lon for easy
      reference. You can copy it somewhere to keep long-term and link to
      it from the GCHP Restarts subdirectory in the run directory.

      .. code-block:: console

         $ mv new_restart_file.nc GEOSChem.Restart.20190701_0000z.c120.s4_32N_64E.nc

   You can also easily reference the file's stretch parameters by
   looking at the global attributes in the file. When using the
   file as a restart file in GCHP make sure that you use the exact
   same parameters in both  the file's global attributes and GCHP
   configuration file :file:`setCommonRunSettings.sh`.

#. Deactivate your GCPy environment when you have finished.

   .. code-block:: console

      $ mamba deactivate

.. _regrid-plot:

===============================
Regridding for Plotting in GCPy
===============================

When plotting in GCPy (e.g. through
:code:`gcpy.compare_single_level()` or
:code:`gcpy.compare_zonal_mean()`), the vast majority of regridding is
handled internally. You can optionally request a specific horizontal
comparison resolution in :code:`compare_single_level()` and
:code:`compare_zonal_mean()`.  Note that all regridding in these
plotting functions only applies to the comparison panels (not
the top two panels which show data directly from each dataset). There
are only two scenarios where you will need to pass extra information
to GCPy to help it determine grids and to regrid when plotting.

Pass stretched-grid file paths
------------------------------

Stretched-grid parameters cannot currently be automatically determined
from grid coordinates. If you are plotting stretched-grid data in
:code:`gcpy.compare_single_level()` or
:code:`gcpy.compare_zonal_mean()` (even if regridding to another
format), you need to use the :code:`sg_ref_path` or
:code:`sg_dev_path` arguments to pass the path of your original
stretched-grid restart file to GCPy. If using :code:`single_panel()`,
pass the file path using :code:`sg_path`. Stretched-grid restart files
created using GCPy contain the specified stretch factor, target
longitude, and target latitude in their metadata.  Currently, output
files from stretched-grid runs of GCHP do not contain any metadata
that specifies the stretched-grid used.

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
:code:`cmpres` argument in :code:`gcpy.compare_single_level()` and
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
