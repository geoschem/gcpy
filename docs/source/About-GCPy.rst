.. |br| raw:: html

   <br/>

.. _about:

##########
About GCPy
##########

:program:`GCPy` is a Python-based toolkit containing useful functions for
working specifically with the :program:`GEOS-Chem` model of
atmospheric chemistry and composition.

GCPy aims to build on the well-established scientific
Python technical stack, leveraging tools like :program:`cartopy`,
:program:`numpy`, and :program:`xarray` to simplify the task of
working with GEOS-Chem model output and performing atmospheric
chemistry analyses.

.. _about-what-gcpy-does:

============================
What GCPy was intended to do
============================

#. Produce plots and tables from `GEOS-Chem
   <https://geos-chem.readthedocs.io>`_ output using simple function
   calls.
#. Generate the standard evaluation plots and tables from GEOS-Chem
   benchmark simulations.
#. Obtain GEOS-Chem's horizontal and vertical grid information.
#. Implement `GCHP <https://gchp.readthedocs.io>`_-specific regridding
   functionalities (e.g. cubed-sphere to lat-lon regridding).
#. Provide example scripts for creating specific types of plots or
   analysis from GEOS-Chem output.
#. Provide user-submitted scripts for specific applications related to
   GEOS-Chem and `HEMCO <https://hemco.readthedocs.io>`_.

.. _about-what-gcpy-doesnt-do:

================================
What GCPy was not intended to do
================================

#. General NetCDF file modification: (crop a domain, extract some variables):

   -  Instead, use netCDF tools such as:

      - `xarray <http://xarray.pydata.org>`_
      - `netCDF Operators (NCO) <https://nco.sourceforge.net/>`_
      - `Climate Data Operators (CDO) <https://mpimet.mpg.de/cdo>`_

   -  Also see our `Work with netCDF data
      <https://geos-chem.readthedocs.io/en/latest/geos-chem-shared-docs/supplemental-guides/netcdf-guide.html>`_
      guide at `geos-chem.readthedocs.io
      <https://geos-chem.readthedocs.io>`_.

#. Statistical analysis:

   -  Instead, use statistical tools such as:

      - `scipy <http://www.scipy.org>`_
      - `scikit-learn <https://scikit-learn.org>`_
      - `R <https://r-project.org>`_
      - etc.

#. Machine Learning:

   -  Instead, use machine learning tools such as:

      - `pytorch <https://pytorch.org>`_
      - `tensorflow <https://www.tensorflow.org>`_
      - `julia <https://julialang.org>`_
      - etc.

=======
License
=======

GCPy is distributed under the `MIT license
<https://opensource.org/license/mit/>`_.  Please see the `GCPy license
agreement  <https://github.com/geoschem/gcpy/blob/dev/LICENSE.txt>`_
and `List of GCPy developers
<https://github.com/geoschem/gcpy/blob/dev/AUTHORS.txt>`_ for more
information.

==================
Requesting support
==================

To report a bug or suggest a new feature, please see our `Support
Guidelines <https://github.com/geoschem/gcpy/blob/dev/SUPPORT.md>`_.

=======================
Submitting new features
=======================

If you are interested in submitting code to GCPy, please see our
`Contributing Guidelines <https://github.com/geoschem/gcpy/blob/dev/CONTRIBUTING.md>`_.
