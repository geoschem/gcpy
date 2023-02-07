.. _install:

###############
Installing GCPy
###############

============
Requirements
============

GCPy is currently supported for Linux and MacOS operating systems. Due
to a reliance on several packages without Windows support, **GCPy is
not currently supported for Windows**. You will receive an error
message if you attempt to use GCPy on Windows.

The only essential software you need before installing GCPy is a
distribution of the Conda package manager. This is used to create a
python environment for GCPy containing all of its software dependences,
including what version of python you use. We recommend using GCPy with
python version 3.9.

You can check if you already have conda installed by running the
following command:

.. code:: console

   $ conda --version

If conda is not already installed then we recommend using Miniconda to
install it. Miniconda is a minimal installer for conda that generally
includes many fewer packages in the base environment than are available
for download. This provides a lightweight conda install from which you
can create custom python environments with whatever python packages you
wish to use, including an environment with GCPy dependencies. To install
Miniconda follow instructions in the  `Miniconda docs <https://docs.conda.io/en/latest/miniconda.html>`__. We recommend using Python 3.9.

==========================================
Steps to install GCPy and its dependencies
==========================================

#. Step 0: Install conda if not already installed.

See the Requirements section above.

#. Step 1: Download GCPy

Create and go to the directory in which you would like to store GCPy. In
this example we will store GCPy in a python/packages subdirectory in the
home directory, but you can store it wherever you wish. You can also name
the GCPy download whatever you want. In this example the GCPy directory
is called GCPy.

.. code:: console

   $ cd $HOME/python/packages
   $ git clone https://github.com/geoschem/gcpy.git GCPy
   $ cd GCPy

#. Step 2: Create new python virtual environment for GCPy

A python virtual environment is a named set of python installs,
e.g. packages, that are independent of other virtual environments.
Using an environment dedicated to GCPy is useful to maintain a set
of package dependencies compatible with GCPy without interfering with
python packages you use for other work. You can create a python virtual
environment from anywhere on your system. It will be stored in your
conda install rather than the directory from which you create it.

You can create a python virtual environment using a file that lists
all packages and their versions to be included in the environment.
GCPy includes such as file, environment.yml, located in the top-level
directory of the package.

Run the following command at the command prompt to create a virtual
environment for use with GCPy. You can name environment whatever you
wish. This example names it gcpy_env.

.. code:: console

   $ conda env create -n gcpy_env --file=environment.yml

Once successfully created you can load the environment by running the
following command, specifying the name of your environment.

.. code:: console

   $ conda activate gcpy_env

To exit the environment do the following:

.. code:: console

   $ conda deactivate

#. Step 3: Add GCPy to python path

The environment variable PYTHONPATH specifies the locations of python
libraries on your system that are not included in your conda environment.
If GCPy is included in PYTHONPATH then python will recognize its
existence when you try to use. Add the following line to your startup
script, e.g. .bashrc, and edit the path to where you are storing GCPy.

.. code:: console

   PYTHONPATH=$PYTHONPATH:$HOME/python/packages/GCPy

#. Step 4: Perform a simple test

Run the following commands in your terminal to check if the 
installation was succcesful.

.. code:: console

   $ source $HOME/.bashrc     # Alternatively close and reopen your terminal
   $ echo $PYTHONPATH         # Check it contains path to your GCPy clone
   $ conda activate gcpy_env    
   $ conda list               # Check it contains contents of gcpy env file
   $ python
   \>>> import gcpy

If no errors were encountered then you successfully installed GCPy and
its dependencies.

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

.. code:: console

   $ cd $HOME/python/packages/GCPy
   $ git fetch -p
   $ git checkout main
   $ git pull

You can also checkout an older version by doing the following:

.. code:: console

   $ cd $HOME/python/packages/GCPy
   $ git fetch -p
   $ git tag
   $ git checkout tags/version_you_want

Once you have the version you wish you use you can do the following
commands to then update your virtual environment:

.. code:: console

   $ source activate gcpy_env
   $ cd $HOME/python/packages/GCPy
   $ conda env update --file environment.yml --prune
