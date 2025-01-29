.. |br| raw:: html

   <br/>

###############
Installing GCPy
###############

.. _install:

In this chapter we will walk you through the steps of installing GCPy
on your computer system.

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
     - Required version
     - Package
     - Required version
   * - `cartopy <https://scitools.org.uk/cartopy/docs/latest/>`_
     - 0.23.0
     - cf_xarray
     - 0.9.1
   * - dask
     - 2024.5.2
     - esmf
     - 8.6.1
   * - `esmpy <https://www.earthsystemcog.org/projects/esmpy/>`_
     - 8.6.1
     - gridspec
     - 0.1.0
   * - ipython
     - 8.25.0
     - `joblib <https://joblib.readthedocs.io/en/latest/>`_
     - 1.4.2
   * - jupyter
     - 1.0.0
     - `matplotlib <https://matplotlib.org/>`_
     - 3.8.4
   * - netcdf4
     - 1.6.5
     - netcdf-fortran
     - 4.6.1
   * - `numpy <http://www.numpy.org/>`_
     - 1.26.4
     - `pandas <https://pandas.pydata.org/docs/>`_
     - 2.2.2
   * - pip
     - 24.0
     - pylint
     - 3.2.2
   * - pyproj
     - 3.6.1
     - `python <https://www.python.org/>`_
     - 3.12.0
   * - pypdf
     - 4.2.0
     - requests
     - 2.32.3
   * - `scipy <http://www.scipy.org/>`_
     - 1.13.1
     - `sparselt <https://github.com/liambindle/sparselt>`_
     - 0.1.3
   * - tabulate
     - 0.9.0
     - tk
     - 8.6.13
   * - `xarray <http://xarray.pydata.org>`_
     - 2024.5.0
     - `xesmf <https://xesmf.readthedocs.io>`_
     - 0.8.5

These packages are also listed the environment file
:file:`docs/environment_files/gcpy_environment.yml`, which is
symlinked to the GCPy root directory as :file:`environment.yml`.

.. note::

   GCPy requires the :program:`esmf` and :program:`esmpy` packages,
   which are not supported on Microsoft Windows.  Thus, the only way
   to use GCPy on a Windows PC is from within a Windows Subsystem for
   Linux (WSL) instance.

.. _install-reqs-mamba-conda:

Mamba or Conda
--------------

In addition to the :ref:`Python packages listed above
<install-reqs-pydeps>` , you will need either the `Mamba
<https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html>`_
or `Conda <https://anaconda.org/anaconda/conda>`_ package manager to
install GCPy.  :program:`Mamba` is a fast drop-in replacement for the
widely-used :program:`Conda` package manager.  We recommend using
:program:`Mamba` if possible, but if your system already has
:program:`Conda` installed, feel free to use it instead.

.. _install-check-mamba-conda:

============================================
Check if Mamba or Conda is already installed
============================================

Follow these instructions to check if you already have a version of
Mamba or Conda installed on your computer system.

First check if :program:`Mamba` has been installed:

.. code-block:: console

   $ mamba --version

If a :program:`Mamba` version exists, you will see output such as:

.. code-block:: console

   mamba version X.Y.Z
   conda version A.B.C

where :literal:`X.Y.Z` and :literal:`A.B.C` are the version numbers.
If you see this output, you may skip ahead to the :ref:`install-methods`
section.

Next, check if :program:`Conda` has been installed:

.. code-block:: console

   $ conda --version

If a :program:`Conda` version exists, you will see its version number
printed to the screen:

.. code-block:: console

   conda version A.B.C

.. note::

   If your :program:`Conda` version is earlier than 23.7, you will
   need to do the following additional steps.

   .. code-block:: console

      $ conda install -n base conda-libmamba-solver
      $ conda config --set solver libmamba

   This will install the fast :program:`Mamba` environment solver into
   your :program:`Conda` base environment. Using the :program:`Mamba`
   solver within :program:`Conda` will considerably speed up the
   Python environment creation.

If a :program:`Conda` version exists, you may skip ahead to the
:ref:`install-methods` section.

If neither :program:`Conda` or :program:`Mamba` are installed, we
recommend installing the :program:`Mamba` package manager yourself, as
described below.

.. _install-mamba:

=============
Install Mamba
=============

This section will walk you through installation of the
:program:`Mamba` package manager.

.. _install-mamba-install:

Install the MambaForge distribution
-----------------------------------

We recommend installing the :program:`MambaForge`, distribution, which
is a full implementation of :program:`Mamba` (as opposed to the
minimal :program:`MicroMamba` distribution).

Follow the instructions below to install :program:`MambaForge`:

MacOS
~~~~~

#. Install :program:`MambaForge` with `Homebrew <https://brew.sh/>`_:

   .. code-block:: console

      $ brew install mambaforge

   |br|

#. Initialize :program:`Mamba` for your shell.  Type one of the
   following commands:

   .. code-block:: console

      $ mamba init bash    # If you use the bash shell (recommended!)
      $ mamba init zsh     # If you use the zsh shell
      $ mamba init fish    # If you use the fish shell

   :program:`Mamba` will add some code to your :file:`~/.bash_profile`
   startup script that will tell your shell where to look for
   Python environments.

   |br|

#. Exit your current terminal session and open a new terminal
   session.  This will apply the changes.

You may now skip ahead to the :ref:`install-methods` section.

Linux and Windows Subsystem for Linux
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#. Download the :program:`MambaForge` installer script from the
   `conda-forge GitHub releases page
   <https://github.com/conda-forge/miniforge/releases>`_:

   .. code-block:: console

      $ wget https://github.com/conda-forge/miniforge/releases/download/24.11.3-0/Miniforge3-24.11.3-0-Linux-x86_64.sh

   This will download the :program:`MambaForge` installer script
   :file:`Mambaforge-24.11.3-0-Linux-x86_64.sh` to your computer.

   .. note::

      As of this writing (January 2025), the latest
      :program:`MambaForge` version is :literal:`24.11.3-0`.  If you
      find that the version has since been updated, simply replace the
      version number :literal:`24.11.3-0` in the above command with the
      most recent version number.

   |br|

#. Change the permission of the :program:`MambaForge` installer script
   so that it is executable.

   .. code-block:: console

      $ chmod 755 Mambaforge-24.11.3-0-Linux-x86_64.sh

   |br|

#. Execute the :program:`Mambaforge` installer script.

   .. code-block::

      $ ./Mambaforge-24.11.3-0-Linux-x86_64.sh

   To update an older version of :program:`Mamba`,  add the
   :literal:`-u` option to the above command.  |br|
   |br|

#. Review and accept the license agreement.

   .. code-block:: console

      In order to continue the installation process, please review the license
      agreement.
      Please, press ENTER to continue
      >>>

   Press :literal:`ENTER` and then :literal:`SPACE` until you reach
   the end of the license agreement.  Then you will be asked:

   .. code-block:: console

      Do you accept the license terms? [yes|no]
      [no] >>>

   Type :literal:`yes` and hit :literal:`ENTER`. |br|
   |br|


#. Specify the root installation path for :program:`MambaForge`.

   .. code-block::

      Mambaforge will now be installed into this location:
     /home/YOUR-USER-NAME/mambaforge

     - Press ENTER to confirm the location
     - Press CTRL-C to abort the installation
     - Or specify a different location below
     [/home/YOUR-USER-NAME/mambaforge] >>>

   In most cases, it should be OK to accept the default installation
   location.  But on some systems, users may be encouraged to install
   software into a different location (e.g. if there is a faster
   filesystem available than the home directory filesystem).
   Consult your sysadmin or IT staff if you are unsure where to
   install :program:`MambaForge`.

   Press the :literal:`ENTER` key to accept the default installation
   path or type a new path and then press :literal:`ENTER`.

   .. code-block:: console

      :program:`MambaForge` will downlad and install Python software
      packages into the  :file:`pkgs` subfolder of the root
      installation path.  Similarly, when you :ref:`create Python
      environments <install-dev-gcpy-install>`, these will be
      installed to the :file:`envs` subfolder of the root installation
      path.

   |br|

#. You may see this warning:

   .. code-block:: console

      WARNING:
       You currently have a PYTHONPATH environment variable set. This may cause
       unexpected behavior when running the Python interpreter in Mambaforge.
       For best results, please verify that your PYTHONPATH only points to
       directories of packages that are compatible with the Python interpreter
       in Mambaforge: /home/YOUR-USER-NAMEb/mambaforge

   As long as your :envvar:`PYTHONPATH` environment variable only
   contains the path to the root-level GCPy folder, you may safely
   ignore this.  (More on :envvar:`PYTHONPATH` :ref:`later
   <install-dev>`.) |br|
   |br|

#. Tell the installer to initialize :program:`MambaForge`.

   .. code-block:: console

      Do you wish the installer to initialize Mambaforge
      by running conda init? [yes|no]
      [no] >>>

   Type :literal:`yes` and then :literal:`ENTER`.  The installer
   script will add some code to your :file:`~/.bashrc` system startup
   file that will tell your shell where to find Python
   environments. |br|
   |br|


#. Exit your current terminal session.  Start a new terminal session
   to apply the updates.  You are now ready to install GCPy.

.. _install-methods:

===========================
Methods for installing GCPy
===========================

Now that you have ensured that a version of :program:`Mamba` or
:program:`Conda` has been installed on your system, you can proceed to
installing GCPy.  There are two different installation methods that
you can use.

.. list-table:: GCPy installation methods
   :header-rows: 1
   :align: center

   * - Method
     - Complexity
     - Who should use it
   * - :ref:`install-conda-forge`
     - Simple
     - GCPy users
   * - :ref:`install-dev`
     - Medium
     - GCPy developers

Unless you are going to be actively developing GCPy, you should
install from conda-forge.

.. _install-conda-forge:

================================
Installing GCPy from conda-forge
================================

GCPy is available through the :code:`conda-forge` channel under the
name :code:`geoschem-gcpy`. :program:`Mamba` or :program:`Conda`
will handle the installation of all dependencies and sub-dependencies
for GCPy, which includes many Python packages and several non-Python
libraries.

Installing GCPy with Mamba
--------------------------

Use these :program:`Mamba` commands to create a Python environment
named :literal:`gcpy_env` and to install GCPy into this environment.

.. code-block:: console

   $ mamba env create -n gcpy_env
   $ mamba activate gcpy_env
   $ mamba config --add channels conda-forge
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

Installing GCPy with Conda
--------------------------

Use these :program:`Mamba` commands to create a Python environment
named :literal:`gcpy_env` and to install GCPy into this environment.

.. code-block:: console

   $ conda env create -n gcpy_env
   $ conda activate gcpy_env
   $ conda config --add channels conda-forge
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
environment.

Install GCPy and its dependencies
---------------------------------

Once you have made sure that :ref:`a Mamba or Conda installation
exists on your system <install-check-mamba-conda>`, you may create a
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
