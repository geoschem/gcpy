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

This example script may be found at `gcpy/examples/plotting/create_test_plot.py <https://github.com/geoschem/gcpy/blob/main/gcpy/examples/plotting/create_test_plot.py>`_.

.. _test-plot-call:

================
Calling sequence
================

Make sure that you :ref:`specified the proper Matplotlib backend
<mpl-backend>` for  your system. Then run the example with:

.. code-block:: console

   $ python -m gcpy.examples.plotting.create_test_plot
