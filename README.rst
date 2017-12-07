GCPy: Python toolkit for GEOS-Chem
==================================

.. image:: https://readthedocs.org/projects/gcpy/badge/?version=latest
    :target: http://gcpy.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

**GCPy** is a Python-based toolkit containing useful functions and routines for
working with GEOS-Chem_, meant to update and replace the IDL-based
GAMAP_ utility. While not a complete re-write of GAMAP_, **GCPy** aims to
build on the well-established scientific python technical stack, leveraging
tools like cartopy_ and xarray_ to simplify the task of working with model
output and performing atmospheric chemistry analyses.

This package is comprised of two major components:

1. A library of functions, implementing some of the core chemistry and
   thermodynamic calculations in GEOS-Chem_ for use outside of the model
2. Documentation including long-form articles, interactive notebooks, and short
   example snippets illustrating workflows using Python and **GCPy**

For more details, please visit `our documentation on ReadTheDocs <http://gcpy.readthedocs.io/en/latest/?badge=latest>`_.


Requirements
------------

**GCPy** is built on top of Python 3 and the scientific Python / NumPy
stack, including

- `Python >= 3.6 <https://www.python.org/>`_
- cartopy_
- `matplotlib <https://matplotlib.org/>`_
- `numpy <http://www.numpy.org/>`_
- `scipy <http://www.scipy.org/>`_
- xarray_

To create an environment for working with **GCPy**, we recommend using
the `Anaconda Python distribution <https://www.continuum.io/downloads>`_
or curating your own *virtualenv* or *conda* environment. Please
refer to `our documentation <http://gcpy.readthedocs.io/en/latest/getting_started.html>`_
for more details.


Installation
------------

At the moment, the easiest way to install **GCPy** is directly from
our bitbucket repository::

    $ pip install git+https://bitbucket.org/gcst/gcpy.git

For development purposes, we recommend cloning the repository and then
manually installing the package in editable mode via *setuptools*::

    $ git clone https://bitbucket.org/gcst/gcpy.git gcpy
    $ cd gcpy
    $ pip install -e .

Currently, **GCPy** is not available via conda-forge or PyPI, but we
anticipate posting early versions of the package to those resources
in conjunction with the official release of GEOS-Chem v11-02.


License
-------

To be determined.


Contact
-------

For more information, please contact the `GEOS-Chem Support Team <geos-chem-support@as.harvard.edu>`_
or reach out to us on the `GEOS-Chem Wiki <http://wiki.seas.harvard.edu/geos-chem/index.php/Main_Page>`_.

.. _cartopy: http://scitools.org.uk/cartopy/
.. _GAMAP: http://acmg.seas.harvard.edu/gamap/
.. _GEOS-Chem: http://acmg.seas.harvard.edu/geos/
.. _xarray: http://xarray.pydata.org/
