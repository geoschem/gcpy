.. |br| raw:: html

   <br/>

.. _install:

###############
Installing GCPy
###############

In this chapter we will walk you through the steps of installing GCPy
on your computer system

.. _install-reqs:

============
Requirements
============

.. _install-reqs-platforms:

Platforms
---------

:program:`GCPy` is currently supported on the following platforms:

#. Linux (x86_64)
#. `Windows Subsystem for Linux
   <https://learn.microsoft.com/en-us/windows/wsl/about>`_ under Microsoft Windows 11
#. MacOS

.. _install-reqs-pydeps:

Software dependencies
---------------------

GCPy requires several other Python packages, which are listed below.

.. list-table:: GCPy dependencies
   :header-rows: 1
   :align: center

   * - Package
     - Version (Python 3.12)
     - Version (Python 3.13)
   * - `cartopy <https://scitools.org.uk/cartopy/docs/latest/>`_
     - 0.23.0
     - 0.24.0
   * - cf_xarray
     - 0.9.1
     - 0.10.0
   * - dask
     - 2025.3.0
     - 2025.3.0
   * - esmf [#A]_
     - 8.6.1
     - 8.8.1
   * - `esmpy <https://www.earthsystemcog.org/projects/esmpy/>`_ [#A]_
     - 8.6.1
     - 8.8.1
   * - gridspec
     - 0.1.0
     - 0.1.0
   * - ipython
     - 8.25.0
     - 9.0.0
   * - `joblib <https://joblib.readthedocs.io/en/latest/>`_
     - 1.4.2
     - 1.4.2
   * - jupyter
     - 1.0.0
     - 1.1.1
   * - `matplotlib <https://matplotlib.org/>`_
     - 3.8.4
     - 3.10.1
   * - netcdf4
     - 1.6.5
     - 1.7.2
   * - netcdf-fortran
     - 4.6.1
     - 4.6.1
   * - `numpy <http://www.numpy.org/>`_
     - 1.26.4
     - 2.1.3
   * - `pandas <https://pandas.pydata.org/docs/>`_
     - 2.2.2
     - 2.2.3
   * - pip
     - 24.0
     - 25.0.1
   * - pylint
     - 3.2.2
     - 3.3.4
   * - pyproj
     - 3.6.1
     - 3.7.1
   * - `python <https://www.python.org/>`_
     - 3.12.0
     - 3.13.0
   * - pypdf
     - 4.2.0
     - 5.3.1
   * - requests
     - 2.32.3
     - 2.32.3
   * - `scipy <http://www.scipy.org/>`_
     - 1.13.1
     - 1.15.2
   * - `sparselt <https://github.com/liambindle/sparselt>`_
     - 0.1.3
     - 0.1.3
   * - tabulate
     - 0.9.0
     - 0.9.0
   * - tk
     - 8.6.13
     - 8.6.13
   * - `xarray <http://xarray.pydata.org>`_
     - 2024.5.0
     - 2025.1.2
   * - `xesmf <https://xesmf.readthedocs.io>`_
     - 0.8.5
     - 0.8.5

.. rubric:: Notes

.. [#A] GCPy requires the :program:`esmf` and :program:`esmpy`
	packages, which are not supported on Microsoft Windows.  Thus,
	the only way to use GCPy on a Windows PC is from within a
	Windows Subsystem for Linux (WSL) instance.
       
The default GCPy environment uses Python 3.12.  In the root folder you
will find a symbolic link :file:`environment.yml`, which points to the
file :file:`docs/environment_files/gcpy_environment_py312.yml` file.
Ths file contains the package specifications listed above under the
**Version (Python 3.12)** column.

We also maintain GCPy environment files for newer Python versions
(e.g. :file:`docs/environment_files/gcpy_environment_py313.yml`, which
is based on Python 3.13). This will allow you to build a GCPy
environment based on a newer Python version, which is often necessary
for testing and development.  For most GCPy users, it should be
sufficient to use the default environment based on Python 3.12.

.. _install-methods:

===========================
Methods for installing GCPy
===========================

You can choose among the following installation methods:

.. list-table:: GCPy installation methods
   :header-rows: 1
   :align: center

   * - Method
     - Complexity
     - Who should use it
   * - :ref:`install-pip`
     - Simple
     - Most GCPy users
   * - :ref:`install-conda-forge`
     - Medium
     - GCPy users who have experience building mamba/conda environments
   * - :ref:`install-dev`
     - Complex
     - GCPy developers

Unless you are going to be actively developing GCPy, you should
install from conda-forge.

.. _install-pip:

======================
Install GCPy from PyPI
======================

If you only plan on using GCPy for visualization of GEOS-Chem
simulation results, you can install GCPy from the :program:`Python
Package Index (PyPi)` using the `Pip installer
<https://pypi.org/project/pip/>`_.

If your system does not already have Pip installed, you may install it
with the `get-pip.py
<https://pip.pypa.io/en/stable/installation/#get-pip-py>`_ script.

.. _install-pip-first:

First-time installation with Pip
--------------------------------

Once you are sure that Pip is installed, you may proceed to download
GCPy with this command:

.. code-block:: console

   $ pip install geoschem-gcpy

To validate the installation, we recommend running the
:ref:`test-plot` example script.

.. _install-pip-update:

Updating to a newer version with Pip
------------------------------------

Use this command to update an existing GCPy installation to a newer version:

.. code-block:: console

   $ pip install -U geoschem-gcpy

You may now skip ahead to the :ref:`mpl-backend` chapter.

.. _install-conda-forge:

=============================
Install GCPy from conda-forge
=============================

GCPy is available through the :code:`conda-forge` channel under the
name :code:`geoschem-gcpy`. :program:`Mamba` or :program:`Conda`
will handle the installation of all dependencies and sub-dependencies
for GCPy, which includes many Python packages and several non-Python
libraries.  If you do not already have a version of :program:`Mamba`
or :program:`Conda` on your system, please see our
:ref:`install-mamba-conda` Supplemental Guide.

.. _install-conda-forge-mamba:

Installing GCPy with Mamba
--------------------------

Use these :program:`Mamba` commands to create a Python environment
named :literal:`gcpy_env` and to install GCPy into this environment.

.. code-block:: console

   $ mamba create -n gcpy_env
   $ mamba activate gcpy_env
   $ mamba install geoschem-gcpy

After you have installed GCPy, check if the installation was
successful by running a test program:

.. code-block:: console

   $ export MPLBACKEND=tkagg   # Sets the matplotlib backend to Tk/Tcl
   $ python -m gcpy.examples.plotting.create_test_plot

If a plot appears on your screen, you have installed GCPy
successfully.  Close the plot window (click the close button or type
:command:`q`) and then deactivate the environment:

.. code-block:: console

   $ mamba deactivate

.. _install-conda-forge-conda:

Installing GCPy with Conda
--------------------------

Use these :program:`Mamba` commands to create a Python environment
named :literal:`gcpy_env` and to install GCPy into this environment.

.. code-block:: console

   $ conda create -n gcpy_env
   $ conda activate gcpy_env
   $ conda install geoschem-gcpy

After you have installed GCPy, check if the installation was
successful by running a test program:

.. code-block:: console

   $ export MPLBACKEND=tkagg   # Sets the matplotlib backend to Tk/Tcl
   $ python -m gcpy.examples.plotting.create_test_plot

If a plot appears on your screen, you have installed GCPy
successfully.  Close the plot window (click the close button or type
:command:`q`) and then deactivate the environment:

.. code-block:: console

   $ conda deactivate

.. _install-dev:

=============================================================
Download GCPy with Git and build a Python virtual environment
=============================================================

If you plan on actively developing GCPy, we recommend that you install
GCPy from Git and create a :program:`Mamba` or :program:`Conda`
environment. If you do not already have a version of :program:`Mamba`
or :program:`Conda` on your system, please see our
:ref:`install-mamba-conda` Supplemental Guide.


Install GCPy and its dependencies
---------------------------------

Once you have made sure that :ref:`a Mamba or Conda installation
exists on your system <install-mamba-conda-check>`, you may create a
Python environment for GCPy. Follow these steps:

#. **Download the GCPy source code.**

   Create and go to the directory in which you would like to store GCPy. In
   this example we will store GCPy in your :file:`$HOME/python/`
   path, but you can store it wherever you wish.  You can also name
   the GCPy download whatever you want. In this example the GCPy
   directory is called :file:`GCPy`.

   .. code-block:: console

      $ cd $HOME/python
      $ git clone https://github.com/geoschem/gcpy.git GCPy
      $ cd GCPy

   |br|

#. **Create a new Python virtual environment for GCPy.**

   A Python virtual environment is a named set of Python installs,
   e.g. packages, that are independent of other virtual
   environments. Using an environment dedicated to GCPy is useful to
   maintain a set of package dependencies compatible with GCPy without
   interfering with Python packages you use for other work. You can
   create a Python virtual environment from anywhere on your
   system. It will be stored in your :program:`Mamba` (or
   :program:`Conda` installation rather than the directory from which
   you create it).

   You can create a Python virtual environment using a file that lists
   all packages and their versions to be included in the environment.
   GCPy includes such as file, :file:`environment.yml`, located in the
   top-level directory of the package.
   
   Run one of the following commands at the command prompt to create a virtual
   environment for use with GCPy. You can name environment whatever you
   wish. This example names it :file:`gcpy_env`.

   .. code-block:: console

      $ mamba env create -n gcpy_env --file=environment.yml   # If using Mamba

      $ conda env create -n gcpy_env --file=environment.yml   # If using Conda

   A list of packages to be downloaded will be displayed.  A
   confirmation message will ask you if you really wish to install all
   of the listed packages.  Type :command:`Y` to proceed or
   :command:`n` to abort.

   Once successfully created you can activate the environment with
   one of these commands:

   .. code-block:: console

      $ mamba activate gcpy_env   # If using Mamba

      $ conda activate gcpy_env   # If using Conda

   To exit the environment, use one of these commands:

   .. code-block:: console

      $ mamba deactivate   # If using Mamba

      $ conda deactivate   # If using Conda

   |br|

#. **Add GCPy to** :envvar:`PYTHONPATH`

   The environment variable :envvar:`PYTHONPATH` specifies the
   locations of Python libraries on your system that were not
   installed by :program:`Mamba`.

   Add the path to your GCPy source code folder :file:`~/.bashrc` file:

   .. code-block:: bash

      export PYTHONPATH=$PYTHONPATH:$HOME/python/GCPy

   and then use

   .. code-block:: console

      $ source ~/.bashrc

   to apply the change. |br|
   |br|

#. **Set the** :envvar:`MPLBACKEND` **environment variable**

   The environment variable :envvar:`MPLBACKEND` specifies the X11
   backend that the Matplotlib package will use to render plots to the
   screen.

   Add this line to your :file:`~/.bashrc` file on your local PC/Mac
   and on any remote computer systems where you will use GCPy:

   .. code-block:: bash

      export MPLBACKEND=tkagg

   And then use:

   .. code-block:: console

      $ source ~/.bashrc

   to apply the change. |br|
   |br|

#. **Perform a simple test:**

   Make sure that you have specified the proper :ref:`mpl-backend` for
   your system.  Then run the following commands in your terminal:

   .. code-block:: console

      $ source $HOME/.bashrc                      # Alternatively close and reopen your terminal
      $ echo $PYTHONPATH                          # Check it contains path to your GCPy clone
      $ mamba activate gcpy_env
      $ mamba list                                # Check it contains contents of gcpy env file
      $ python -m gcpy.examples.create_test_plot  # Create a test plot

If the plot appears on your screen, then the GCPy installation was successful.

If no error messages are displayed, you have successfully installed
GCPy and its dependencies.

.. _install-dev-upgrade:

Upgrading GCPy versions
-----------------------

Sometimes the GCPy dependency list changes with a new GCPy version,
either through the addition of new packages or a change in the minimum
version. You can always update to the latest GCPy version from within
you GCPy clone, and then update your virtual environment using the
environment.yml file included in the package.

Run the following commands to update both your GCPy version to the
latest available.

.. code-block:: console

   $ cd $HOME/python/GCPy
   $ git fetch -p
   $ git checkout main
   $ git pull

You can also checkout an older version by doing the following:

.. code-block:: console

   $ cd $HOME/python/GCPy
   $ git fetch -p
   $ git tag
   $ git checkout tags/version_you_want

Once you have the version you wish you use you can do the following
commands to then update your virtual environment:

.. code-block:: console

   $ mamba activate gcpy_env
   $ cd $HOME/python/GCPy
   $ mamba env update --file environment.yml --prune
