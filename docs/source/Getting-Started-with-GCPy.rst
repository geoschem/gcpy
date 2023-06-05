.. |br| raw:: html

   <br/>

.. _install:

###############
Installing GCPy
###############

.. _requirements:

============
Requirements
============

:program:`GCPy` is currently supported for Linux and MacOS operating
systems. Due to a reliance on several packages without Windows
support, **GCPy is not currently supported for Windows**. You will
receive an error message if you attempt to use GCPy on Windows.

.. tip::

   Windows 11 (and some later builds of Windows 10) support the
   `Windows Subsystem for Linux (WSL)
   <https://learn.microsoft.com/en-us/windows/wsl/install>`_. If your
   Windows version is WSL-compatible, you can install GCPy into a
   Linux instance (such as Ubuntu 22.04) running under Windows.  At
   present, this is the only way to use GCPy locally on a Windows
   computer.

The only essential software you need before installing GCPy is a
distribution of the :program:`Conda` package manager. This is used to
create a Python environment for GCPy containing all of its software
dependences, including what version of Python you use. You must
using GCPy with Python version 3.9.

You can check if you already have Conda installed by running the
following command:

.. code-block:: console

   $ conda --version

.. attention::

   You must use Conda 4.12.0 or earlier to install GCPy and its
   dependencies.  Newer versions of Conda than this will install
   Python package versions that are incompatible with GCPy. See
   :ref:`Installing Conda 4.12.0 with Miniconda <conda412_install>`
   below.

   In the future we hope to be able to resolve this installation issue
   so that you can use the latest Conda version.

If Conda is not already installed, you must use :program:`Miniconda`
to install Conda 4.12.0.  Miniconda is a minimal installer for Conda
that generally includes many fewer packages in the base environment
than are available for download. This provides a lightweight Conda
installation from which you can create custom Python environments with
whatever Python packages you wish to use, including an environment
with GCPy dependencies.

.. _conda412_install:

============================================
Steps to install Conda 4.12.0 with Miniconda
============================================

If you already have a Conda version prior to 4.12.0 installed on your
system, you may skip this step and proceed to the section entitled
:ref:`gcpy_install`.

If you need to install Conda 4.12.0, follow these steps:

#. Download the Miniconda installer script for your operating system
   as shown below. The script will install Conda version 4.12.0 using
   Python 3.9.

   **Linux (x86_64 CPUs)**

   .. code-block:: console

      $ wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh

   **MacOS (M1 CPUs)**

   .. code-block:: console

      $ wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-MacOSX-arm64.sh

   **MacOS (x86_64 CPUs)**

   .. code-block:: console

      $ wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-MacOSX-x86_64.sh

   .. tip::

      If you do not have :program:`wget` installed on MacOS, you can
      download it with the :program:`Homebrew` package manager:

      .. code-block::

	 $ brew install wget

   In the steps that follow, we will walk through installation using
   the Linux installer script.  The steps are the same for MacOS; just
   substitute the appropriate MacOS script name for the Linux script
   name in steps 2 and 3 below. |br|
   |br|


#. Change the permission of the Miniconda installer script so that it
   is executable:

   .. code-block:: console

      $ chmod 755 Miniconda3-py39_4.12.0-Linux-x86_64.sh

   |br|

#. Run the Miniconda installer script.

   .. code-block:: console

      $ ./Miniconda3-py39_4.12.0-Linux-x86_64.sh

   |br|

#. Accept the license agreement.

   When the installer script starts, you will be prompted to accept
   the Miniconda license agreement:

   .. code-block:: console

     Welcome to Miniconda3 py39_4.12.0

     In order to continue the installation process, please review the license
     agreement.
     Please, press ENTER to continue
     >>>

   When you press :literal:`ENTER`, you will see the license agreement
   in all of its gory legalese detail.  Press the space bar repeatedly
   to scroll down ot the end. You will then see this prompt:

   .. code-block:: console

      Do you accept the license terms? [yes|no]
      [no] >>>

   Type :literal:`yes` and hit :literal:`ENTER` to accept. |br|
   |br|


#. Specify the installation path.

   You will then be prompted to provide a directory path for the
   installation:

   .. code-block:: console

      Miniconda3 will now be installed into this location:
      /home/YOUR-USERNAME/miniconda3

      - Press ENTER to confirm the location
      - Press CTRL-C to abort the installation
      - Or specify a different location below

      [/home/YOUR-USERNAME/miniconda3] >>>

   Press :literal:`ENTER` to continue, or specify a new path and then
   press :literal:`ENTER`.

   .. tip::

      If a previous Conda installation is already installed to the
      default path, you may choose to delete the previous installation
      folder, or install Conda 4.12.0 to a different path.

   The script will then start installing the Conda 4.12.0 package
   manager. |br|
   |br|


#. Specify post-installation options.

   You will see this text at the bottom of the screen printout upon
   successful installation:

   .. code-block:: console

      Preparing transaction: done
      Executing transaction: done
      installation finished.
      Do you wish the installer to initialize Miniconda3
      by running conda init? [yes|no]
      [no] >>>

   Type :literal:`yes` and press :literal:`ENTER`.  You will see
   output similar to this:

   .. code-block:: console

      no change     /home/bob/miniconda3/condabin/conda
      no change     /home/bob/miniconda3/bin/conda
      no change     /home/bob/miniconda3/bin/conda-env
      no change     /home/bob/miniconda3/bin/activate
      no change     /home/bob/miniconda3/bin/deactivate
      no change     /home/bob/miniconda3/etc/profile.d/conda.sh
      no change     /home/bob/miniconda3/etc/fish/conf.d/conda.fish
      no change     /home/bob/miniconda3/shell/condabin/Conda.psm1
      no change     /home/bob/miniconda3/shell/condabin/conda-hook.ps1
      no change     /home/bob/miniconda3/lib/python3.9/site-packages/xontrib/conda.xsh
      no change     /home/bob/miniconda3/etc/profile.d/conda.csh
      no change     /home/bob/.bashrc
      No action taken.
      If you'd prefer that conda's base environment not be activated on startup,
         set the auto_activate_base parameter to false:

      conda config --set auto_activate_base false

      Thank you for installing Miniconda3!

   |br|

#. Disable the base Conda environment from being activated at startup

   Close the terminal window that you used to install Conda 4.12.0 and
   open a new terminal window.  You will see this prompt:

   .. code-block:: console

      (base) $

   By default, Conda will open the :literal:`base` environment each
   time that you open a new terminal window.  to disable this
   behavior, type:

   .. code-block:: console

      (base) $ conda config --set auto_activate_base false

   The next time you open a terminal window, you will just see the
   regular prompt, such as;

   .. code-block:: console

      $

   (or whatever you have defined your prompt to be in your startup scripts).

Now that you have installed Conda 4.12.0, you may proceed to creating
a new Conda environment for GCPy, as shown below.

.. _gcpy_install:

==========================================
Steps to install GCPy and its dependencies
==========================================

#. Install Conda if it is not already installed.

   If Conda 4.12.0 or prior is already installed on your system, you
   may skip this step.  Otherwise, please follow the instructions
   listed in :ref:`conda412_install`. |br|
   |br|

#. Download the GCPy source code.

   Create and go to the directory in which you would like to store GCPy. In
   this example we will store GCPy in a :file:`python/packages`
   subdirectory in your home directory, but you can store it wherever
   you wish. You can also name the GCPy download whatever you want. In
   this example the GCPy directory is called :file:`GCPy`.

   .. code-block:: console

      $ cd $HOME/python/packages
      $ git clone https://github.com/geoschem/gcpy.git GCPy
      $ cd GCPy

   |br|

#. Create a new Python virtual environment for GCPy.

   A Python virtual environment is a named set of Python installs,
   e.g. packages, that are independent of other virtual
   environments. Using an environment dedicated to GCPy is useful to
   maintain a set of package dependencies compatible with GCPy without
   interfering with Python packages you use for other work. You can
   create a Python virtual environment from anywhere on your
   system. It will be stored in your Conda installation rather than
   the directory from which you create it.

   You can create a Python virtual environment using a file that lists
   all packages and their versions to be included in the environment.
   GCPy includes such as file, environment.yml, located in the
   top-level directory of the package.

   Run the following command at the command prompt to create a virtual
   environment for use with GCPy. You can name environment whatever you
   wish. This example names it :file:`gcpy_env`.

   .. code-block:: console

      $ conda env create -n gcpy_env --file=environment.yml

   Once successfully created you can load the environment by running the
   following command, specifying the name of your environment.

   .. code-block:: console

      $ conda activate gcpy_env

   To exit the environment do the following:

   .. code-block:: console

      $ conda deactivate

   |br|

#. Add GCPy to Python path.

   The environment variable :envvar:`PYTHONPATH` specifies the
   locations of Python libraries on your system that are not included
   in your conda environment. If GCPy is included in
   :envvar:`PYTHONPATH` then Python will recognize its existence
   when you try to use. Add the following line to your startup script,
   e.g. :file:`.bashrc`, and edit the path to where you are storing
   GCPy.

   .. code-block:: bash

      export PYTHONPATH=$PYTHONPATH:$HOME/python/packages/GCPy

   |br|

#. Perform a simple test.

   Run the following commands in your terminal to check if the
   installation was succcesful.

   .. code-block:: console

      $ source $HOME/.bashrc     # Alternatively close and reopen your terminal
      $ echo $PYTHONPATH         # Check it contains path to your GCPy clone
      $ conda activate gcpy_env
      $ conda list               # Check it contains contents of gcpy env file
      $ python
      >>> import gcpy

If no error messages are displayed, you have successfully installed
GCPy and its dependencies.

=======================
Upgrading GCPy versions
=======================

Sometimes the GCPy dependency list changes with a new GCPy version,
either through the addition of new packages or a change in the minimum
version. You can always update to the latest GCPy version from within
you GCPy clone, and then update your virtual environment using the
environment.yml file included in the package.

Run the following commands to update both your GCPy version to the
latest available.

.. code-block:: console

   $ cd $HOME/python/packages/GCPy
   $ git fetch -p
   $ git checkout main
   $ git pull

You can also checkout an older version by doing the following:

.. code-block:: console

   $ cd $HOME/python/packages/GCPy
   $ git fetch -p
   $ git tag
   $ git checkout tags/version_you_want

Once you have the version you wish you use you can do the following
commands to then update your virtual environment:

.. code-block:: console

   $ source activate gcpy_env
   $ cd $HOME/python/packages/GCPy
   $ conda env update --file environment.yml --prune
