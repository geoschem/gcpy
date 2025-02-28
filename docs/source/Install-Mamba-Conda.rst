.. |br| raw:: html

   <br/>

.. _install-mamba-conda:

#############################################
Install Mamba or Conda Python package manager
#############################################

As we learned in the :ref:`install` chapter, there are multiple
installation methods for GCPy.  Some of these installation use either
the `Mamba
<https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html>`_
or `Conda <https://anaconda.org/anaconda/conda>`_ package manager.
:program:`Mamba` is a fast drop-in replacement for the
widely-used :program:`Conda` package manager.  We recommend using
:program:`Mamba` if possible, but if your system already has
:program:`Conda` installed, feel free to use it instead.

.. _install-mamba-conda-check:

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

.. _install-mamba-conda-mambaforge:

===================================
Install the MambaForge distribution
===================================

We recommend installing the :program:`MambaForge`, distribution, which
is a full implementation of :program:`Mamba` (as opposed to the
minimal :program:`MicroMamba` distribution).

Follow the instructions below to install :program:`MambaForge`:

MacOS
-----

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
-------------------------------------

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
   to apply the updates.  You are now ready to :ref:`install GCPy
   <install>`.
