.. |br| raw:: html

   <br/>

.. _test-plot:

##################
Create a Test Plot
##################

This example demonstrates how you can create a simple test plot with
GCPy.  You can use this test plot to verify that GCPy was installed
properly onto your system.

.. image:: _static/images/create\_test\_plot.png
   :align: center

|br|

.. _test-plot-code:

===========
Source code
===========

**Script location:** `gcpy/examples/plotting/create_test_plot.py <https://github.com/geoschem/gcpy/blob/main/gcpy/examples/plotting/create_test_plot.py>`_

.. _test-plot-usage:

=====
Usage
=====

.. code-block:: console

   $ conda activate gcpy_env

   (gcpy_env) $ export MPLBACKEND=tkagg   # Or MacOSX if you are on a Mac

   (gcpy_env) $ python -m gcpy.examples.plotting.create_test_plot

At this point you should see the plot above on your screen.  To close
the plot window you may either :program:`X` button or type
:command:`q`.  You may then deactivate the Python environment:

.. code-block:: console

   (gcpy_env) $ conda deactivate
