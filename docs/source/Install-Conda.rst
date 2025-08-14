.. |br| raw:: html

   <br/>

.. _install-mamba-conda:

########################################
Install the Conda Python package manager
########################################

As we learned in the :ref:`install` chapter, there are multiple
installation methods for GCPy.  Some of these installation methods use
the `Conda <https://anaconda.org/anaconda/conda>`_ package manager.
If :program:`Conda` has not been already installed on your system, you
can install it following the instructions in this chapter.
not already present on your system.

.. important::

   Previous versions of this documentation encouraged users to install
   the :program:`Mamba` package mamager.  However, :program:`Mamba` has
   been deprecated as of August 2024 and is slated for removal as a
   stand-alone package. :program:`Mamba` functionality has since been
   incorporated into :program:`Conda` version 2.24 and later.

   We have updated these installation instructions accordingly, and
   now direct users to install :program:`Conda` via the
   :program:`Miniforge` distribution.


.. _install-mamba-conda-check:

===================================
Check if Conda is already installed
===================================

To check if a version of :program:`Conda` has already been installed
on your system, type:

.. code-block:: console

   $ conda --version

If a :program:`Conda` version exists, you will see its version number
printed to the screen:

.. code-block:: console

   conda version A.B.C

If a :program:`Conda` version exists, you may skip ahead to the
:ref:`install-methods` section.  If not, then proceed as described
below.

.. _install-mamba-conda-mambaforge:

==================================
Install the Miniforge distribution
==================================

We recommend installing :program:`Conda` from the :program:`Miniforge` distribution.


#. Download the :program:`Miniforge` installer script with the
   :program:`wget` download utility:

   .. code-block:: console

      $ wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"


   Your shell will run the :program:`uname` command to fill in your OS
   type and processor type automatically.

   |br|


#. Execute the :program:`Miniforge` installer script.

   .. code-block::

      $ bash Miniforge3-$(uname)-$(uname -m).sh

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


#. Specify the root installation path:

   .. code-block:: console

      Miniforge3 will now be installed into this location:
      /home/YOUR-USER-NAME/miniforge3

      - Press ENTER to confirm the location
      - Press CTRL-C to abort the installation
      - Or specify a different location below

      [/home/YOUR-USER-NAME/miniforge3] >>>

   In most cases, it should be OK to accept the default installation
   location.  But on some systems, users may be encouraged to install
   software into a different location (e.g. if there is a faster
   filesystem available than the home directory filesystem).
   Consult your sysadmin or IT staff if you are unsure where to
   install :program:`Miniforge`.

   Press the :literal:`ENTER` key to accept the default installation
   path or type a new path and then press :literal:`ENTER`.


   :program:`Miniforge` will downlad and install Python software
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
       unexpected behavior when running the Python interpreter in Miniforge3.
       For best results, please verify that your PYTHONPATH only points to
       directories of packages that are compatible with the Python interpreter
       in Miniforge3: /path/to/miniforge3
    

   As long as your :envvar:`PYTHONPATH` environment variable only
   contains the path to the root-level GCPy folder, you may safely
   ignore this.  (More on :envvar:`PYTHONPATH` :ref:`later
   <install-dev>`.) |br|
   |br|

#. Tell the installer to initialize :program:`Miniforge`:

   .. code-block:: console

      Do you wish the installer to initialize Miniforge
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
