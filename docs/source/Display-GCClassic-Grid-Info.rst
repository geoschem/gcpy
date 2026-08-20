.. _display-gcc-grid-info:

####################################
Display GEOS-Chem Classic grid info
####################################

GEOS-Chem Classic supports several horizontal grid resolutions,
including both global and regional ("nested") domains. This example
demonstrates how you can print out the metadata (resolution,
longitude/latitude ranges, and grid-cell centers/edges) for any of
these supported grids, which is useful when configuring a run
directory or writing a custom regridding script.

.. _display-gcc-grid-info-code:

===========
Source code
===========

.. list-table::
   :header-rows: 1

   * - Description
     - Script location
   * - :mod:`gcpy.examples.grids.display_gcclassic_grid_info`
     - `gcpy/examples/grids/display_gcclassic_grid_info.py <https://github.com/geoschem/gcpy/blob/main/gcpy/examples/grids/display_gcclassic_grid_info.py>`_

.. _display-gcc-grid-info-usage:

=====
Usage
=====

Activate your GCPy environment with mamba or conda:

.. code-block:: console

   $ mamba activate gcpy_env    # If using mamba, or

   $ conda activate gcpy_env    # If using conda

Then run the script as a module:

.. code-block:: console

   $ python -m gcpy.examples.grids.display_gcclassic_grid_info

You will be prompted to select a grid from a numbered menu of the
global and nested-domain resolutions that GEOS-Chem Classic supports
(4.0° x 5.0°, 2.0° x 2.5°, 0.5° x 0.625°, 0.25° x 0.3125°, and
0.125° x 0.15625°, each with their available nested regions). After
you make a selection, the script will print the grid's name,
resolution, longitude/latitude ranges, and the longitude/latitude
values at the grid-cell centers and edges.
