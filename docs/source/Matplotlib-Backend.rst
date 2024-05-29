.. _mpl-backend:

##############################
Specify backend for Matplotlib
##############################

GCPy uses the `matplotlib <https://matplotlib.org/>`_ for creating the
various types of plots as described in this manual.  You will need to
specify the proper backend for the Python matplotlib library.  This
will make sure that plots appear on your screen properly.

You can set the :envvar:`MPLBACKEND` environment variable to specify
the `backend
<https://matplotlib.org/stable/users/explain/figure/backends.html#backends>`_
that you wish to use.

You can place one of these commands in your :file:`~/.bashrc` startup
script so the

==========
For MacOSX
==========

Use this command if you are running GCPy on a Mac:

.. code-block:: console

   $ export MPLBACKEND=MacOSX

=========
For Linux
=========

Use this command if you are running GCPy on Linux (or on a Windows PC
using the `Windows Subsystem for Linux
<https://learn.microsoft.com/en-us/windows/wsl/>`_): 

.. code-block:: console

   $ export MPLBACKEND=tkagg
