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

================
Regridding Files
================

You can regrid existing GEOS-Chem restart or output diagnostic files
between lat/lon and cubed-sphere formats using
:code:`gcpy.file_regrid`. :code:`gcpy.file_regrid` can either be
called directly from the command line using :code:`python -m
gcpy.file_regrid` or as a function
(:code:`gcpy.file_regrid.file_regrid()`) from a Python script or
interpreter. The syntax of :code:`file_regrid` is as follows:

.. code-block:: python

   def file_regrid(fin, fout, dim_format_in, dim_format_out,
   cs_res_out=0, ll_res_out='0x0',
   sg_params_in=[1.0, 170.0, -90.0], sg_params_out=[1.0, 170.0, -90.0]
   ):
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

      Format of the input file's dimensions (choose from: classic,
      checkpoint, diagnostic), where classic denotes lat/lon and
      checkpoint / diagnostic are cubed-sphere formats

.. option:: dim_format_out : str

      Format of the output file's dimensions (choose from: classic,
      checkpoint, diagnostic), where classic denotes lat/lon and
      checkpoint / diagnostic are cubed-sphere formats

Optional arguments:
-------------------

.. option:: cs_res_out : int

      The cubed-sphere resolution of the output dataset. Not used if
      dim_format_out is classic.

      Default value: 0

.. option:: ll_res_out : str

      The lat/lon resolution of the output dataset. Not used if
      dim_format_out is not classic/

      Default value: '0x0'

.. option:: sg_params_in : list[float, float, float]

      Input grid stretching parameters [stretch-factor, target
      longitude, target latitude]. Not used if dim_format_in is classic

      Default value: [1.0, 170.0, -90.0] (No stretching)

.. option:: sg_params_out : list[float, float, float]

      Output grid stretching parameters [stretch-factor, target
      longitude, target latitude].  Not used if dim_format_out is classic.

      Default value: [1.0, 170.0, -90.0] (No stretching)

There are three dimension formats available for regridding: :literal:`classic`
(GEOS-Chem Classic lat/lon format), :literal:`checkpoint` (GCHP restart file
format), and :literal:`diagnostic` (GCHP output file format). You can
regrid between any of these formats using :code:`file_regrid`, as well as
between different resolutions  and/or grid-types within each dimension
format (e.g. standard cubed-sphere checkpoint to stretched-grid
checkpoint). Note that although the :code:`cs_res_out` and
:code:`ll_res_out` parameters are technically optional in the
function, you must specify at least one of these in your call to
:code:`file_regrid`.

As stated previously, you can either call
:code:`file_regrid.file_regrid()` directly or call it from the command
line using :code:`python -m gcpy.file_regrid ARGS`. An example command
line call (separated by line for readability) for regridding a C90
cubed-sphere restart file to a C48 stretched-grid with a stretch
factor of 3, a target longitude of 260.0, and a target latitude of
40.0 looks like:

.. code-block::

   python -m gcpy.file_regrid             \
         -i initial_GEOSChem_rst.c90_standard.nc   \
         --dim_format_in checkpoint      \
         -o sg_restart_c48_3_260_40.nc       \
         --cs_res_out 48            \
         --sg_params_out 3.0 260.0 40.0      \
         --dim_format_out checkpoint

.. _regrid-sparselt:

=====================================
Regridding with gridspec and sparselt
=====================================

GCPy 1.3.0 and later supports regridding with the `gridspec <https://github.com/liambindle/gridspec>`_ and `sparselt <https://github.com/liambindle/sparselt>`_
utilities.

.. _regrid-sparselt-firsttime:

First-time setup
-----------------

#. Install command line tool gridspec in your bin directory

   .. code-block:: console

      $ pip install git+https://github.com/LiamBindle/gridspec.git

#. Make sure location of installation is added to path in your bashrc
   (or equivalent)

   .. code-block:: bash

      $ export PATH=/path/to/home/.local/bin:$PATH

#. Install sparselt as a python package.

   .. code-block:: console

      $ conda install -c conda-forge sparselt==0.1.3

.. _regrid-sparselt-gridcombo:

One-time setup per grid resolution combination
----------------------------------------------

#. Create a directory structure to store files that you will use in
   regridding. Ideally this would be in a shared location where all of
   the GCPy users at your institution coud access it.

   Navigate to this directory.

   .. code-block:: console

      $ mkdir /path/to/RegridInfo

#. Within this top level directory, create two directories that will
   store grid information and regridding weights.  Navigate to the
   grid information folder.

   .. code-block:: console

      $ mkdir Grids
      $ mkdir Weights
      $ cd Grids

#. Create tilefiles (if cubed-sphere) and grid spec file for each
   input and output grid resolution (see also gridspec README):

   For uniform cubed-sphere global grid, specify face side length.

   #. For simplicity, keep all cubed-sphere data in subdirectories
      of the Grids folder.

      .. code-block:: console

         $ mkdir c24
         $ gridspec-create gcs 24
         $ mv c24*.nc c24

         $ mkdir c48
         $ gridspec-create gcs 48
         $ mv c48*.nc c48

          ... etc for other grids ...

   #. For cubed-sphere stretched grid, specify face side length,
      stretch factor, and target latitude and longitude:

      .. code-block:: console

         $ mkdir sc24
         $ gridspec-create sgcs 24 -s 2 -t 40 -100
         $ mv *c24*.nc sc24

   #. For uniform global lat-lon grid, specify the number of latitude and
      longitude grid boxes. For a list of optional settings, run the
      command :command:`gridspec-create latlon --help`.

      Create a subdirectory named latlon and move all of your latlon grid
      specification files there.

      .. code-block:: console

         $ gridspec-create latlon 90 180                # Generic 1 x 1 grid
         $ gridspec-create latlon 46 72 -dc -pc -hp     # GEOS-Chem Classic 4 x 5
         $ gridspec-create latlon 91 144 -dc -pc -hp    # GEOS-Chem Classic 2 x 2.5
         $ gridspec-create latlon 361 576 -dc -pc -hp   # MERRA-2 0.5 x 0.625
         $ gridspec-create latlon 721 1172 -dc -pc -hp  # GEOS-FP 0.25 x  0.3125

         $ mkdir latlon
         $ mv regular_lat_lon*.nc latlon

#. (Optional) View contents of grid spec file:

   .. code-block:: console

      $ gridspec-dump c24/c24_gridspec.nc

      ... etc. for other grids ...

#. Initialize your GCPy conda environmnt (which includes ESMF as a
   dependency):

   .. code-block:: console

      $ conda activate gcpy_env

#. Navigate to the directory that will store the regridding
   weights. (Recall that we created this in created this in step #2.

   .. code-block:: console

      $ cd /path/to/RegridInfo/Weights

#. Generate regridding weights (see also sparselt sample data files
   README), specifying the following:

   - Path to input file horizontal resolution grid spec netcdf file
   - Path to output file horizontal resolution grid spec netcdf file
   - Regridding type, e.g. conserve for conservative (string)
   - Name of output regridding weights file (include input and output
     resolutions)
   - Name of directory containing grid spec tilefiles

   .. code-block:: console

      (gcpy_env) $ /ESMF_RegridWeightGen                                  \
                   -s /path/to/RegridInfo/Grids/c48/c48_gridspec.nc       \
                   -d /path/to/RegridInfo/Grids/regular_lat_lon_90x180.nc \
                   -m conserve                                            \
                   -w ./regrid_weights_c48_to_latlon90x180.nc             \
                   --tilefile_path /path/to/RegridInfo/Grids/c48

      ... etc. for other grid combinations ...

#. (Optional) Consider using a bash script such as the one shown below
   if you need to create regridding weights to/from several grids.

   .. code-block:: bash

      #!/bin/bash

      # Generates regridding weights with ESMF_RegridWeightGen

      # The top-level directory containing Grids and Weights subdirectories
      # (EDIT AS NEEDED)
      main_dir="/path/to/RegridInfo"

      # Subdirectories for grid specifications and regridding weights
      grids_dir="${main_dir}/Grids"
      weights_dir="${main_dir}/Weights"

      # GCHP cubed-sphere grids (EDIT AS NEEDED)
      cs_list=(c24 c48 c90 c180 c360)

      # GCClassic lat-lon grids (EDIT AS NEEDED)
      ll_list=(46x72 91x144 361x576 721x1172)

      # Loop over cubed-sphere grids
      for cs in ${cs_list[@]}; do

          # Cubed-sphere gridspec file
          cs_grid_info="${grids_dir}/${cs}/${cs}_gridspec.nc"
          if [[ ! -f ${cs_grid_info} ]]; then
              echo "Could not find ${cs_grid_info}!"
              exit 1
          fi

          # Loop over latlon grids
          for ll in ${ll_list[@]}; do

              # Latlon gridspec file
              ll_grid_info="${grids_dir}/latlon/regular_lat_lon_${ll}.nc"
              if [[ ! -f ${ll_grid_info} ]]; then
                  echo "Could not find ${ll_grid_info}!"
                  exit 1
              fi

              # Cubed-sphere -> latlon regridding
              echo "----"
              echo "Regridding from ${cs} to ${ll}"
              weightfile="${weights_dir}/regrid_weights_${cs}_to_latlon${ll}.nc"
              ESMF_RegridWeightGen                  \
                  -s ${cs_grid_info}                \
                  -d ${ll_grid_info}                \
                  -m conserve                       \
                  -w ${weightfile}                  \
                  --tilefile_path ${grids_dir}/${cs}
              unset weightfile

              # Latlon -> cubed-sphere regridding
              echo "----"
              echo "Regridding from ${ll} to ${cs}"
              weightfile="${weights_dir}/regrid_weights_latlon${ll}_to_${cs}.nc"
              ESMF_RegridWeightGen                  \
                  -s ${ll_grid_info}                \
                  -d ${cs_grid_info}                \
                  -m conserve                       \
                  -w ${weightfile}                  \
                  --tilefile_path ${grids_dir}/${cs}
              unset weightfile

          done
      done

.. _regrid-sparselt-regrid:

Sample regridding script
------------------------

Once you have created the tilefiles and regridding weights, you can
use them to regrid data files.  Shown below is a sample Python script
that you can modify.

.. code-block:: python

   #!/usr/bin/env python

   # Imports
   import xarray as xr
   import sparselt.esmf
   import sparselt.xr

   # Create a linear transform object from the regridding weights file
   # for the combination of source and target horizontal resolutions.
   transform = sparselt.esmf.load_weights(
       'path/to/RegridInfo/Weights/regrid_weights_c48_to_latlon90x180.nc',
        input_dims=[('nf', 'Ydim', 'Xdim'), (6, 48, 48)]
        output_dims=[('lat', 'lon'), (90, 180)],
   )

   # Open file to regrid as xarray DataSet.
   ds = xr.open_dataset('my_data_c48.nc')

   # Regrid the DataSet using the transform object.
   ds = sparselt.xr.apply(transform, ds)

   # Write xarray DataSet contents to netcdf file.
   ds.to_netcdf("my_data_latlon90x180.nc")
