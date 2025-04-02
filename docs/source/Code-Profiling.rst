.. |br| raw:: html

   <br/>

.. _code-profiling:

###################################
Plotting output from code profilers
###################################

GCPy contains functions to plot output from the :program:`gprofng` and
:program:`Intel VTune` code profilers, which can help you to identify
computational bottlenecks (aka "hotspots") in GEOS-Chem
Classic, GCHP, and HEMCO.

.. _code-profiling-gprofng:

=======
gprofng
=======

This example demonstrates how you can plot function profiles generated
by the :program:`gprofng` performance profiler.

.. image:: _static/images/gprofng\_functions.png
   :align: center


.. _code-profiling-gprofng-source:

Source code
-----------

**Script location:** `gcpy/examples/gprofng/plot_functions.py <https://github.com/geoschem/gcpy/blob/main/gcpy/examples/gprofng/plot_functions.py>`_

.. _code-profiling-gprofng-usage:

Usage
-----

First, generate a function profile with :program:`gprofng`:

.. code-block:: console

   $ gprofng collect app /path/to/executable/file

where :code:`/path/to/executable/file` is the path to the program that
you wish to profile.  For example, to profile GEOS-Chem Classic, you
would use this command:

.. code-block:: console

   $ gprofng collect app ./gcclassic

:program:`Gprofng` will send profiling output to a folder named
:file:`test.N.er`, where :code:`N` is an integer index and :code:`er`
stands for  "experiment record".

Next, send function profiling information to a file:

.. code-block:: console

   $ echo functions | gprofng display text test.N.er > functions_profile.txt

Here is a sample :file:`functions_profile.txt` for GEOS-Chem Classic.

.. code-block:: text

   (gp-display-text) Functions sorted by metric: Exclusive Total CPU Time

   Excl. Total     Incl. Total      Name
   CPU             CPU
      sec.      %     sec.      %
   508.406 100.00  508.406 100.00   <Total>
    62.594  12.31   62.594  12.31   __unitconv_mod_MOD_convertbox_kgm2_to_kg
    61.653  12.13   61.653  12.13   __unitconv_mod_MOD_convertbox_kg_to_kgm2
    58.401  11.49   58.401  11.49   <static>@0x7696c (<libm-2.28.so>)
    46.833   9.21   64.075  12.60   __tomas_mod_MOD_mnfix
    41.609   8.18   41.609   8.18   __gckpp_linearalgebra_MOD_kppdecomp
    29.941   5.89   29.941   5.89   __gckpp_linearalgebra_MOD_kppsolve
    26.609   5.23   26.609   5.23   __gckpp_function_MOD_fun_split
    ... etc ...

The :program:`Excl. Total` (total Exclusive Time) metric is useful for
identifying computational bottlenecks.  This represents the amount of
time spent in a subroutine, excluding time spent in any subroutines
called by the subroutine.

Make sure that you have :ref:`specified the proper Matplotlib backend
<mpl-backend>` for your system. Then run the example script with the
following command:

.. code-block:: console

   $ python -m gcpy.examples.profile.gprofng_functions functions_profile.txt 1 40

This will create a plot similar to that shown above, which shows the
top 40 functions sorted by exclusive time. To see the next 40
functions sorted by exclusive time, use this command:

.. code-block:: console

   $ python -m gcpy.examples.profile.gprofng_functions functions_profile.txt 41 80

etc.  You should display fewer than 50 functions in order to prevent
the time and percent labels from overlapping.

.. _code-profiling-vtune:

===========
Intel VTune
===========

To be added
