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
     - 6.4.0
     - 6.4.0
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

The default GCPy environment uses Python 3.13.  In the root folder you
will find a symbolic link :file:`environment.yml`, which points to the
file :file:`docs/environment_files/gcpy_environment_py313.yml` file.
Ths file contains the package specifications listed above under the
**Version (Python 3.13)** column.

We also maintain a GCPy environment file for backwards compatibility
with Python 3.12
(:file:`docs/environment_files/gcpy_environment_py312.yml`).  However,
we have not yet been able to construct a GCPy environment using Python
3.14, due to conflicts in installing the :program:`esmf` and
:program:`xesmf` packages.  We hope to be able to resolve this in the
near future.

.. _install-methods:

===========================
Methods for installing GCPy
===========================

You can choose among the following installation methods:

.. list-table:: GCPy installation methods
   :header-rows: 1
   :align: center
   :widths: 40 10 50

   * - Method
     - Complexity
     - Who should use it
   * - :ref:`install-conda-forge` (Recommended)
     - Simple
     - Most GCPy users
   * - :ref:`install-pip`
     - Simple
     - Those who do not wish to install a Conda environment
   * - :ref:`install-dev`
     - Complex
     - GCPy developers

Unless you are going to be actively developing GCPy, you should
:ref:`install from conda-forge <install-conda-forge>`.

.. _install-conda-forge:

=============================
Install GCPy from conda-forge
=============================

GCPy is available through the :code:`conda-forge` channel under the
name :code:`geoschem-gcpy`. :program:`Conda` will handle the
installation of all dependencies and sub-dependencies for GCPy, which
includes many Python packages and several non-Python libraries.  If
you do not already have a version of :program:`Conda` on your system,
please see our :ref:`install-conda` Supplemental Guide.

1. Create and activate a Python virtual environment
----------------------------------------------------

You will use :program:`Conda` to create a Python virtual environment
named :literal:`gcpy_env`, into which you will install GCPy and its
dependencies.  A :program:`Python virtual environment` is a named set
of Python packages that are independent of other virtual
environments. Using an environment dedicated to GCPy is useful to
maintain a set of package dependencies compatible with GCPy without
interfering with Python packages you use for other work.

To create and activate the GCPy Python environment, use these
commands:

.. code-block:: console

   $ conda create -n gcpy_env
   $ conda activate gcpy_env
   (gcpy_env) $

Activating the environment adds the :literal:`(gcpy_env)` prefix added
to the command prompt.  This is a visual cue to remind you that the
environment is active.

2. Install GCPy with Conda
--------------------------

Next, install GCPy from conda-forge:

.. code-block:: console

   (gcpy_env) $ conda install geoschem-gcpy

Conda will ask you to confirm that you want to proceed with the
installation.  Type :program:`y` to proceed.

3. Verify the installation
--------------------------

Run the :ref:`test-plot` example script to validate the installation.
If you see a plot created on your screen, GCPy has been installed
successfully.

4. Deactivate the Python environment
------------------------------------

Next, deactivate the :literal:`gcpy_env` environment:

.. code-block:: console

   (gcpy_env) $ conda deactivate
   $

This also removes the :literal:`(gcpy_env)` prefix from the command
prompt.

You may now skip ahead to the :ref:`next chapter <mpl-backend>`.

.. _install-pip:

======================
Install GCPy from PyPI
======================

If you do not wish to create a Python virtual environment, you may
install GCPy from the :program:`Python Package Index (PyPi)` using the
`Pip installer <https://pypi.org/project/pip/>`_.  If your system does
not already have Pip installed, you may install it with the
`get-pip.py
<https://pip.pypa.io/en/stable/installation/#get-pip-py>`_ script.

.. attention::

   If you install GCPy from PyPI, you will lose the ability to keep
   the Python packages needed for GCPy separate from other Python
   packages. This can potentially lead to package conflicts.  For this
   reason we recommend :ref:`installing GCPy from conda-forge
   <install-conda-forge>`.

.. _install-pip-first:

1. Install GCPy with Pip
------------------------

Once you are sure that Pip is installed, you may proceed to download
GCPy with this command:

.. code-block:: console

   $ pip install geoschem-gcpy

2. Validate the installation
----------------------------

Run the :ref:`test-plot` example script to validate the installation.
If you see a plot created on your screen, GCPy has been installed
successfully.

.. _install-pip-update:

(Optional) Update to a newer GCPy version
-----------------------------------------

Use this command to update an existing GCPy installation to a newer version:

.. code-block:: console

   $ pip install -U geoschem-gcpy

You may now skip ahead to the :ref:`next chapter <mpl-backend>`.

.. _install-dev:

=============================================================
Download GCPy with Git and build a Python virtual environment
=============================================================

If you plan on actively developing GCPy, we recommend that you install
GCPy from Git and create a :program:`Python virtual environment` with
:program:`Conda`.  If :program:`Conda` is not already installed on
your system, you may follow the instructions in our
:ref:`install-conda` Supplemental Guide.

.. _install-dev-gcpy-install:

Install GCPy and its dependencies
---------------------------------

Once you have made sure that :ref:`a Conda installation exists on your
system <install-conda-check>`, you may create a Python environment for
GCPy. Follow these steps:

1. Download the GCPy source code
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create and go to the directory in which you would like to store GCPy. In
this example we will store GCPy in your :file:`$HOME/python/`
path, but you can store it wherever you wish.  You can also name
the GCPy download whatever you want. In this example the GCPy
directory is called :file:`GCPy`.

.. code-block:: console

   $ cd $HOME/python
   $ git clone https://github.com/geoschem/gcpy.git GCPy
   $ cd GCPy

2. Create a new Python virtual environment for GCPy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A :program:`Python virtual environment` is a named set of Python
packages that are independent of other virtual environments. Using an
environment dedicated to GCPy is useful to maintain a set of package
dependencies compatible with GCPy without interfering with Python
packages you use for other work.

GCPy ships with YAML files that specify Python environments for Python
3.12 and 3.13.  These are located in the
:file:`docs/environment_files` folder.  The symbolic link
:file:`environment.yml`, located in the top-level directory of the
package, points to the default environment (which is based on Python
3.13).

Use this command to create a Python environment named :file:`gcpy_env`
for use with GCPy:

.. code-block:: console

   $ conda env create -n gcpy_env --file=environment.yml

A list of packages to be downloaded will be displayed.  A
confirmation message will ask you if you really wish to install all
of the listed packages.  Type :command:`Y` to proceed or
:command:`n` to abort.

Once successfully created you can activate the environment with
one of these commands:

.. code-block:: console

   $ conda activate gcpy_env
   (gcpy_env) $

Activating the GCPy Python environment adds the
:literal:`(gcpy_env)` added to the command prompt.  This is a visual
cue to remind you that the environment is active.

To exit the environment, use this command:

.. code-block:: console

   (gcpy_env) $ conda deactivate
   $

Deactivating the environment removes :literal:`(gcpy_env)` from the
command prompt.

3. Add GCPy to :envvar:`PYTHONPATH`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The environment variable :envvar:`PYTHONPATH` specifies the
locations of Python libraries on your system that were not
installed by :program:`Conda`.

Add the path to your GCPy source code folder :file:`~/.bashrc` file:

.. code-block:: bash

   export PYTHONPATH=$PYTHONPATH:$HOME/python/GCPy

and then use

.. code-block:: console

   $ source ~/.bashrc

 to apply the change.

4. Set the :envvar:`MPLBACKEND` environment variable
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The environment variable :envvar:`MPLBACKEND` specifies the X11
backend that the Matplotlib package will use to render plots to the
screen.

Add this line to your :file:`~/.bashrc` file on your local PC/Mac
and on any remote computer systems where you will use GCPy:

.. code-block:: bash

   export MPLBACKEND=tkagg   # or MacOSX if you are on a Mac

And then use:

.. code-block:: console

   $ source ~/.bashrc

to apply the change.

5. Perform a simple test
~~~~~~~~~~~~~~~~~~~~~~~~

Make sure that you have specified the proper :ref:`mpl-backend` for
your system.  Then run the following commands in your terminal:

.. code-block:: console

   $ source $HOME/.bashrc                                 # Or, close and reopen your terminal

   $ echo $PYTHONPATH                                     # Check it contains path to your GCPy clone

   $ conda activate gcpy_env                              # Activate the environment

   (gcpy_env) $ conda list                                # Check it contains contents of gcpy env file

You may now run the :ref:`test-plot` example script.  If you see a
plot created on your screen, the installation was successful.

Then deactivate the environment with:

.. code-block:: console

   (gcpy_env) conda deactivate
   $

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

   $ conda activate gcpy_env
   $ cd $HOME/python/GCPy
   $ conda env update --file environment.yml --prune
