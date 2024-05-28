.. _obspack:

##############################################################
Create a coordinates file for the GEOS-Chem ObsPack diagnostic
##############################################################

This community contribution script generates a netCDF file containing
longitude, latitude, and altitude values at which
`GEOS-Chem ObsPack diagnostic output
<https://geos-chem.readthedocs.io/en/latest/gcclassic-user-guide/obspack.html>`_
will be archived.

Please reach out to the author directly if you have questions about
using this script:

   +-----------------+---------------+
   | Author          | GitHub Handle |
   +=================+===============+
   | Alli Moon       | @alli-moon    |
   +-----------------+---------------+
   | Yuk Yun Chan    | @yc-chan      |
   +-----------------+---------------+

.. _obspack-code:

===========
Source code
===========

**Script location:** `gcpy/community/create_obspack_coords_file.py <https://github.com/geoschem/gcpy/blob/main/gcpy/community/create_obspack_coords_file.py>`_

.. _obspack-usage:

=====
Usage
=====

Navigate to a local directory where you will run the script.  Then
copy :file:`gcpy/community/create_obspack_coords_file.py` there:

.. code-block:: console

   $ cd /path/to/my/local/directory

   $ cp /path/to/gcpy/gcpy/community/create_obspack_coords_file.py .

Open :file:`create_obspack_coords_file.py` in your favorite editor and
edit the configurable settings for your particular use case:

.. code-block:: python

   # ==================================================================
   # Configurable settings -- edit for your own use case!
   DATASET_ID = 'GC-MODEL'
   OUTPUT_DIR = './'

   SAMPLING_DATES_SUMMER = pd.date_range(
       start='2022-05-31',
       end='2022-07-01',
       freq='d'
   )
   SAMPLING_DATES_WINTER = pd.date_range(
       start='2023-01-13',
       end='2023-2-13',
       freq='d'
   )

   SITE_LAT = 32.2644
   SITE_LON = -64.8767
   SITE_ALT = 0
   ASSUMED_INLET_HEIGHT_ABOVE_GROUND = 30 # Unit:m
   # ==================================================================

Once you have done this, you may run the script:

.. code-block:: console

   $ python -m create_obspack_coords_file
