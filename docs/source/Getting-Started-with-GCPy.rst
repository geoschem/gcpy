.. _install:

###############
Installing GCPy
###############

GCPy and its dependencies can be installed using `Conda
<https://github.com/conda/conda>`__ in either standard user mode
(allow Conda to handle installation without Git support) or
development mode using Conda and `Conda-build
<https://github.com/conda/conda-build>`__ (install from a Git
clone). You can also manually install GCPy using a clone of the source
code, but this option requires you to add the package to your
:envvar:`PYTHONPATH` manually and to install properly versioned
dependencies on your own.

.. _install-reqs:

=====================
Requirements for GCPy
=====================

.. _install-reqs-software:

Software Prerequisites
----------------------

GCPy is currently supported for Linux and MacOS operating systems. Due
to a reliance on several packages without Windows support, **GCPy is
not currently supported for Windows**. You will receive an error
message if you attempt to install GCPy through Conda on Windows.

The only essential software package you need before installing GCPy is a
distribution of the Conda package manager, which is used to install GCPy
and its dependencies. It is also highly recommended that you create an
environment in Conda for working with GCPy. Steps to setup Conda are
described below:

#. Install `Miniconda or Anaconda <https://github.com/conda/conda>`__.
   Miniconda is much more lightweight and functions perfectly well for
   GCPy purposes, while Anaconda automatically includes many extra
   packages that are not directly relevant to GCPy.
#. After installing Miniconda or Anaconda, create a Conda environment
   for using GCPy. The basic usage (also found on the `Conda Github
   hompeage <https://github.com/conda/conda>`__) is:

   .. code-block:: bash

      # Navigate to the top-level GCPy folder
      cd /path/to/gcpy

      # Create a Conda environment for working with GCPy
      conda env create -n gcpy_env --file=environment.yml

      # Activate (enter) your new Conda environment
      $ conda activate gcpy_env

      # Deactivate (exit) your Conda environment
      $ conda deactivate

   From within your Conda environment, you can follow the instructions
   on `Installing normally through Conda (if you don't plan on
   modifying GCPy source code) <#installing-gcpy-for-non-developers-using-conda>`__ or `Installing in development
   mode through Conda-build (for developers) <#install_dev>`__.

.. _install-reqs-pydeps:

Python dependencies
-------------------

Conda handles the installation of all dependencies for GCPy
automatically. Most dependencies have minimum version
requirements. We recommend using GCPy with Python 3.9. The list of
dependencies (not including sub-dependencies) that are installed by
Conda includes:

-  `Python 3.9 <https://www.python.org/>`_
-  `cartopy <https://scitools.org.uk/cartopy/docs/latest/>`_
-  `matplotlib <https://matplotlib.org/>`_
-  `numpy <http://www.numpy.org/>`_
-  `scipy <http://www.scipy.org/>`_
-  `xarray <http://xarray.pydata.org>`_
-  `xesmf <https://xesmf.readthedocs.io>`_
-  `esmpy <https://www.earthsystemcog.org/projects/esmpy/>`_
-  `pypdf2 <https://pythonhosted.org/PyPDF2/>`_
-  `joblib <https://joblib.readthedocs.io/en/latest/>`_
-  `xbpch <https://github.com/darothen/xbpch>`_
-  `pandas <https://pandas.pydata.org/docs/>`_
-  `sparselt >= 0.1.3 <https://github.com/liambindle/sparselt>`_

A full list of package version requirements may be found in
:file:`docs/source/environment.yml`. There is also a symbolic link to
this file from the top-level gcpy folder.

.. _install-non-devs:

==============================================
Installing GCPy for non-developers using Conda
==============================================

GCPy is available through the :code:`conda-forge` channel under the
name :code:`geoschem-gcpy`. Installing GCPy in your Conda environment
requires two commands:

.. code:: console

   $ conda config --add channels conda-forge
   $ conda install geoschem-gcpy

Conda will handle the installation of all dependencies and
sub-dependencies for GCPy, which includes many Python packages and
several non-Python libraries.

.. _install-devs:

==============================
Installing GCPy for developers
==============================

If you wish to make changes to the GCPy source code with the goal of
contributing to GCPy development, you will need to install GCPy from a
clone of the GCPy Git repository:

.. code:: console

   $ git clone https://github.com/geoschem/gcpy.git
   $ cd gcpy
   $ conda config --add channels conda-forge
   $ conda install geoschem-gcpy --only-deps
   $ pip install -e .

Conda will handle the installation of dependencies when you install
from this clone, and pip will point all GCPy links to this directory.
