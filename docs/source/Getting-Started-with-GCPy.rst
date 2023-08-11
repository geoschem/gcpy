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
distribution of the :program:`Mamba` package manager. :program:`Mamba`
is a drop-in replacement for the widely-used :program:`Conda`
package manager.  You will use :program:`Mamba` to create a Python
environment for GCPy, which contains a Python version (currently
3.9.6) plus specific versions of various dependent packages.

Check if you already have :program:`Mamba` by using this command:

.. code-block:: console

   $ mamba --version

If :program:`Mamba` has already been installed on your system, **skip
ahead to the** :ref:`gcpy_install` **section**.

.. _mamba_install:

==============================
Installing Mamba as MicroMamba
==============================

:program:`MicroMamba` is a minimal :program:`Mamba` version that is
easier to install than the full :program:`MambaForge` distribution.
Follow the installation steps listed below:

#. Download and initialize :program:`MicroMamba`.

   Execute this command:

   .. code-block:: console

      $ "${SHELL}" <(curl -L micro.mamba.pm/install.sh)

   This will download the :program:`MicroMamba` via the
   :program:`curl` utility.  You will be then asked several questions:

   .. code-block:: console

      Micromamba binary folder? [~/.local/bin]

   This prompt is asking where you prefer to install the
   :program:`MicroMamba` executable. Press :command:`ENTER` to accept
   the default or type a new location, then press :command:`ENTER`.    

   .. code-block::   
      
      Prefix location? [~/micromamba]

   This prompt is asking where the Python packages and environments
   created with Micromamba should be installed.  The default location
   is your :file:`$HOME/micromamba` folder.  Press :command:`ENTER` to
   accept the default or specify a new location and then press
   :command:`ENTER`.

   .. code-block::
   
      Init shell? [Y/n]

   This prompt is asking if you would like :program:`MicroMamba` to
   add some code into your :file:`$HOME/.bashrc` startup script.
   Press :command:`Y` and then :command:`ENTER`.

   .. code-block:: 

      Configure conda-forge? [Y/n]

   This prompt is asking if you would like :program:`MicroMamba` to
   have access to the :literal:`conda-forge` repository. Press
   :command:`Y` then :command:`ENTER`. 

   Then :program:`MicroMamba` installer will print out the following
   information to the screen:
   
   .. code-block:: console

      Modifying RC file "/path/to/.bashrc"
      Generating config for root prefix "/path/to/root/prefix"
      Setting mamba executable to: "/path/to/mamba/executable-dir/micromamba
      Adding (or replacing) the following in your "/path/to/.bashrc" file

   |br|

#. Tell your shell where it can find the :program:`MicroMamba` executable.

   If you have not done so already, add the following
   line to your :file:`$HOME/.bashrc` startup script:

   .. code-block:: bash

      export PATH="/path/to/mamba/executable-dir:$PATH"

   where :file:`/path/to/mamba/executable-dir` is the same text as
   displayed in the previous step.
      
   .. note::

      Some shared computer systems prefer that users place
      modifications not into the :file:`$HOME/.bashrc` file, but
      instead to a different script (e.g. :file:`$HOME/.bash_aliases`)
      that is executed by :file:`$HOME/.bashrc`.  Ask your system
      administrator for more information pertaining to your particular
      setup.

   Apply the change with this command:

   .. code-block:: console

      $ source $HOME/.bashrc

   This will tell your shell to look for executable files in your
   :file:`$HOME/bin` folder before it looks through the rest of your
   search path.  |br|
   |br|

#. Define the :literal:`mamba` convenience alias.

   Add the following lines to your :file:`$HOME/.bashrc` file

   .. code-block:: bash

      # Invoke micromamba as "mamba"
      alias mamba="micromamba"

   Apply the change with this command:

   .. code-block:: console

      $ source ~/.bashrc

   This will allow you to invoke :program:`MicroMamba` by typing
   :literal:`mamba`. |br|
   |br|

   You are now ready to use :program:`Mamba` (installed as
   :program:`MicroMamba`)!

.. _gcpy_install:

=================================
Install GCPy and its dependencies
=================================

Once :program:`Mamba` has been installed, you may proceed use it to
create a Python environment for GCPy.  (Please return to
:ref:`mamba_install` if you have not yet installed :program:`Mamba`.)

#. Download the GCPy source code.

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

#. Create a new Python virtual environment for GCPy.

   A Python virtual environment is a named set of Python installs,
   e.g. packages, that are independent of other virtual
   environments. Using an environment dedicated to GCPy is useful to
   maintain a set of package dependencies compatible with GCPy without
   interfering with Python packages you use for other work. You can
   create a Python virtual environment from anywhere on your
   system. It will be stored in your :program:`Mamba` installation
   rather than the directory from which you create it.

   You can create a Python virtual environment using a file that lists
   all packages and their versions to be included in the environment.
   GCPy includes such as file, :file:`environment.yml`, located in the
   top-level directory of the package.

   Run the following command at the command prompt to create a virtual
   environment for use with GCPy. You can name environment whatever you
   wish. This example names it :file:`gcpy_env`.

   .. code-block:: console

      $ mamba env create -n gcpy_env --file=environment.yml

   A list of packages to be downloaded will be displayed.  A
   confirmation message will ask you if you really wish to install all
   of the listed packages.  Type :command:`Y` to proceed or
   :command:`n` to abort.

   Once successfully created you can activate the environment with
   this command:

   .. code-block:: console

      $ mamba activate gcpy_env

   To exit the environment, use this command:

   .. code-block:: console

      $ mamba deactivate

   |br|

#. Add GCPy to Python path.

   The environment variable :envvar:`PYTHONPATH` specifies the
   locations of Python libraries on your system that are not included
   in your conda environment. If GCPy is included in
   :envvar:`PYTHONPATH` then Python will recognize its existence.

   Add the path to your GCPy source code folder  :file:`~/.bashrc` file:

   .. code-block:: bash

      export PYTHONPATH=$PYTHONPATH:$HOME/python/GCPy

   and then use

   .. code-block:: console

      $ source ~/.bashrc

   to apply the change. |br|
   |br|

#. Perform a simple test:

   Run the following commands in your terminal to check if the
   installation was succcesful.

   .. code-block:: console

      $ source $HOME/.bashrc     # Alternatively close and reopen your terminal
      $ echo $PYTHONPATH         # Check it contains path to your GCPy clone
      $ mamba activate gcpy_env
      $ mamba list               # Check it contains contents of gcpy env file
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
