.. _hco-fmt:

######################################
Format netCDF files for input to HEMCO
######################################

This community contribution script formats a netCDF file so that its
global and variable attributes are COARDS-compliant.  This is a
requirement for input to `GEOS-Chem
<https://geos-chem.readthedocs.io>`_ via the `Harmonized Emission
Component (HEMCO) <https://hemco.readthedocs.io>`_.

Please reach out to the author directly if you have questions about
using this script:

   +-----------------+-----------------+
   | Author          | GitHub Handle   |
   +=================+=================+
   | Hannah Nesser   | @hannahnesser   |
   +-----------------+-----------------+

.. _hco-fmt-code:

===========
Source code
===========

**Script location:** `gcpy/community/format_hemco_data.py <https://github.com/geoschem/gcpy/blob/main/gcpy/community/format_hemco_data.py>`_

**Related example script:**  `gcpy/examples/hemco/format_hemco_demo.py
<https://github.com/geoschem/gcpy/blob/main/gcpy/community/format_hemco_demo.py>`_

.. _hco-fmt-usage:

=====
Usage
=====

Navigate to a local directory where you will run the script.  Then
download this data file from the AWS :file:`s3://gcgrid` bucket:

.. code-block:: console

   $ wget https://gcgrid.s3.amazonaws.com/HEMCO/GCClassic_Output/14.0.0/2019/GEOSChem.ProdLoss.20190101_0000z.nc4

and rename it to :file:`HEMCO_demonstration_file.nc`:

.. code-block:: console

   $ mv GEOSChem.ProdLoss.20190101_0000z.nc4 HEMCO_demonstration_file.nc

Then you may run the demonstration script:

.. code-block:: console

   $ python -m gcpy.examples.hemco.format_hemco_demo

You can view the :file:`gcpy/examples/hemco/format_hemco_demo.py`
script in your favorite editor to see the individual commands that
are executed.
