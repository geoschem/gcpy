.. _about:

##########
About GCPy
##########

**GCPy** is a Python-based toolkit containing useful functions for
working specifically with the GEOS-Chem model of atmospheric chemistry
and composition.

GCPy aims to build on the well-established scientific Python technical
stack, leveraging tools like cartopy and xarray to simplify the task of
working with model output and performing atmospheric chemistry analyses.

.. _about-what-gcpy-does:

============================
What GCPy was intended to do
============================

#. Produce plots and tables from GEOS-Chem output using simple function
   calls.
#. Generate the standard evaluation plots and tables from GEOS-Chem
   benchmark output.
#. Obtain GEOS-Chem's horizontal/vertical grid information.
#. Implement GCHP-specific regridding functionalities (e.g. cubed-sphere
   to lat-lon regridding).
#. Provide example scripts for creating specific types of plots or
   analysis from GEOS-Chem output.

.. _about-what-gcpy-doesnt-do:

================================
What GCPY was not intended to do
================================

#. General NetCDF file modification: (crop a domain, extract some variables):

   -  Use `xarray <http://xarray.pydata.org>`__ instead.
   -  Also see `Work with netCDF data files
      <https://geos-chem.readthedocs.io/en/latest/geos-chem-shared-docs/supplemental-guides/netcdf-guide.html>`_
      at the GEOS-Chem ReadTheDocs site.

#. Statistical analysis:

   -  Use `scipy <http://www.scipy.org>`__ and `scikit-learn
      <https://scikit-learn.org>`__ tools instead.

#. Machine Learning:

   -  Use the standard machine learning utilities
      (`pytorch <https://pytorch.org>`__,
      `tensorflow <https://www.tensorflow.org>`__,
      `julia <https://julialang.org>`__, etc.).
