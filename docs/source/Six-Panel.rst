.. _six-panel:

##################
Six Panel Plotting
##################

This example demonstrates GCPy's comparison plotting capabilities.
Following the example below will generate a plot similar to this:

.. image:: _static/images/six\_panel\_single\_level.png
   :align: center

.. _six-panel-code:

===========
Source code
===========

**Script location:** `gcpy/examples/plotting/plot_comparisons.py <https://github.com/geoschem/gcpy/blob/main/gcpy/examples/plotting/plot_comparisons.py>`_

.. _six-panel-call:

================
Calling sequence
================

Make sure that you have :ref:`specified the proper Matplotlib backend
<mpl-backend>` for  your system. Then run the example with:

.. code-block:: console

   $ python -m gcpy.examples.plotting.plot_comparisons \
     --ref   /path/to/ref/diagnostic/or/restart/file \
     --dev   /path/to/dev/diagnostic/or/restart/file \
     --var   variable-name-to-plot \
     --level level-to-plot

Note that the :code:`level-to-plot` starts from 0.

For example, to plot July 2019 ozone concentrations at the surface level
from two different GEOS-Chem simulations, you would use this command:

.. code-block:: console

   $ python -m gcpy.examples.plottings.plot_comparisons \
     --ref   /path/to/ref/GEOSChem.SpeciesConc.20190701_0000z.nc4 \
     --dev   /path/to/dev/GEOSChem.SpeciesConc.20190701_0000z.nc4 \
     --var   SpeciesConcVV_O3 \
     --level 0
