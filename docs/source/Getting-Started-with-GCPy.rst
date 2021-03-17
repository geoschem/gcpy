Getting Started with GCPy
=========================

This page describes the installation process for GCPy and how to test to
make sure GCPY is installed correctly.


Installing GCPy
---------------

GCPy and its dependencies can be installed using
`Conda <https://github.com/conda/conda>`__ in either standard user mode
(allow Conda to handle installation without Git support) or development
mode using Conda and
`Conda-build <https://github.com/conda/conda-build>`__ (install from a
Git clone). You can also manually install GCPy using a clone of the
source code, but this option requires you to add the package to your
``PYTHONPATH`` manually and to install properly versioned dependencies
on your own.

Software Prerequisites
~~~~~~~~~~~~~~~~~~~~~~

GCPy is currently supported for Linux and MacOS operating systems. Due to a reliance on several packages without Windows support,
**GCPy is not currently supported for Windows**. You will receive an error message if you attempt to install GCPy through Conda
on Windows.

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

   .. code:: bash

       # Create a Conda environment for working with GCPy
       conda create -n gcpy_env

       # Activate (enter) your new Conda environment
       $ conda activate gcpy_env #Linux / MacOS
       > activate gcpy_env        #Windows

       # Deactivate (exit) your Conda environment
       $ conda deactivate #Linux / MacOS
       > deactivate        #Windows

   From within your Conda environment, you can follow the instructions
   on `Installing normally through Conda (if you don't plan on
   modifying GCPy source code) <#installing-gcpy-for-non-developers-using-conda>`__ or `Installing in development
   mode through Conda-build (for developers) <#install_dev>`__.

Python dependencies
^^^^^^^^^^^^^^^^^^^

Conda handles the installation of all dependencies for GCPy
automatically. Most dependencies have minimum version requirements. GCPy has been tested with Python 3.6,
3.7, and 3.8. The list of dependencies (not including
sub-dependencies) that are installed by Conda includes:

-  `Python >= 3.6 and < 3.9 <https://www.python.org/>`__
-  `cartopy <https://scitools.org.uk/cartopy/docs/latest/>`__
-  `matplotlib <https://matplotlib.org/>`__
-  `numpy <http://www.numpy.org/>`__
-  `scipy <http://www.scipy.org/>`__
-  `xarray <http://xarray.pydata.org>`__
-  `xesmf <https://xesmf.readthedocs.io>`__
-  `esmpy <https://www.earthsystemcog.org/projects/esmpy/>`__
-  `pypdf2 <https://pythonhosted.org/PyPDF2/>`__
-  `joblib <https://joblib.readthedocs.io/en/latest/>`__
-  `xbpch <https://github.com/darothen/xbpch>`__
-  `pandas <https://pandas.pydata.org/docs/>`__

A full list of package version requirements can be found in
``docs/environment.yml``.

Installing GCPy for non-developers using Conda
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GCPy is available through the ``conda-forge`` channel under the name
``geoschem-gcpy``. Installing GCPy in your Conda environment requires two commands:

.. code:: bash

   conda config --add channels conda-forge
   conda install geoschem-gcpy


Conda will handle the installation of all dependencies and
sub-dependencies for GCPy, which includes many Python packages and
several non-Python libraries.

Installing GCPy for developers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you wish to make changes to the GCPy source code with the goal of
contributing to GCPy development, you will need to install GCPy from a
clone of the GCPy Git repository:

.. code:: bash

    git clone https://github.com/geoschem/gcpy.git
    cd gcpy
    conda config --add channels conda-forge
    conda install geoschem-gcpy --only-deps
    pip install -e .
	
	
Conda will handle the installation of dependencies when you install
from this clone, and pip will point all GCPy links to this directory.


Manual install using source code (pre-1.0.0)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Versions of GCPy prior to 1.0.0 do not support installation through
Conda. However, you can still use Conda to install requisite
dependencies by `creating a Conda environment from the sample
environment
file <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file>`__
at ``docs/environment_files/gcpy_min/environment.yml``. Then clone the GCPy repository using
``git clone https://github.com/geoschem/gcpy.git``. You will also need
to add the GCPy directory to the Python path using
``export PYTHONPATH=/path/to/gcpy:$PYTHONPATH``, where
``/path/to/gcpy/`` is the top-level directory of the GCPy repository.

Optional extra Python libraries
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The GCPy repository contains a few different ``environment.yml`` files for creating
new Conda environments. ``docs/environment_files`` features three different options:
``gcpy_min``, ``gcpy_full``, and ``gcpy_extra``. 

-  ``gcpy_min`` contains only the libraries necessary for executing all GCPy functions, and is equivalent to the environment generated by running ``conda install geoschem-gcpy``.
-  ``gcpy_full`` contains everything in ``gcpy_min`` as well as Jupyter (for working with / developing Jupyter notebook examples) and IPython.
-  ``gcpy_extras`` contains everything in ``gcpy_full`` as well as extra libraries for scientific analysis in Python outside of GCPy, such as scikit-learn.


Testing your GCPy installation
------------------------------

Once you've installed GCPy using one of the methods installed above, you
should make sure the package functions correctly. From within your Conda
environment, type:

::

    $    python
    >>>  import gcpy

If no errors appear, congratulations! GCPy and its dependencies are probably properly
installed. If you run into any problems, feel free to open an issue at
`the GCPy Issues page on
Github <https://github.com/geoschem/gcpy/issues>`__.


