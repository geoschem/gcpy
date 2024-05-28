.. |br| raw:: html

   <br/>

.. _single-panel:

#####################
Single Panel Plotting
#####################

This example demonstrates GCPy's single-panel plotting capabilities.
Following the example below will generate a series of plots such as
this surface-level ozone concentrations plot:

.. image:: _static/images/single\_panel\_single\_level.png
   :align: center

|br|

and this zonal mean plot:

.. image:: _static/images/single\_panel\_zonal\_mean.png
   :align: center

|br|

.. _single-panel-code:

===========
Source code
===========

**Script location:** `gcpy/examples/plotting/plot_single_panel.py <https://github.com/geoschem/gcpy/blob/main/gcpy/examples/plotting/plot_single_panel.py>`_

.. _single-panel-usage:

=====
Usage
=====

Make sure that you have :ref:`specified the proper Matplotlib backend
<mpl-backend>` for your system. Then run the example with:

.. code-block:: console

   $ python -m gcpy.examples.plotting.plot_single_panel \
     --infile   /path/to/ref/diagnostic/or/restart/file \
     --varname  variable-name-to-plot \
     --level    level-to-plot

Note that the :code:`level-to-plot` starts from 0.

For example, to plot July 2019 ozone concentrations, you would use
the command below.  Data from the surface level is used for the
single-level plot.

.. code-block:: console

   $ python -m gcpy.examples.plotting.plot_single_panel \
     --infile  /path/to/GEOSChem.SpeciesConc.20190701_0000z.nc4 \
     --varname SpeciesConcVV_O3 \
     --level 0
