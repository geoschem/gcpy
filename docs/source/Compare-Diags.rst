.. _comp-diags:

##########################
Compare diagnostic outputs
##########################

This example demonstrates GCPy's diagnostic comparison capabilities. 
Following the example below will generate a table commparing the sums
of individual variables from two GEOS-Chem diagnostic or restart
files.  

.. code-block:: console

   Using configuration file compare_diags.yml
   ... Printing totals and differences
   Variable               Ref=GCC_ref              Dev=GCC_dev              Dev - Ref
   AREA                 : 510065600000000.0      | 510065600000000.0      | 0.0 
   SpeciesConcVV_A3O2   : 9.399016e-10           | 9.399016e-10           | 0.0 
   SpeciesConcVV_ACET   : 6.726078e-05           | 6.726078e-05           | 0.0 
   SpeciesConcVV_ACTA   : 5.329012e-06           | 5.329012e-06           | 0.0 
   SpeciesConcVV_AERI   : 7.7059624e-08          | 7.7059624e-08          | 0.0 
   SpeciesConcVV_ALD2   : 5.2878436e-06          | 5.2878436e-06          | 0.0 
   SpeciesConcVV_ALK4   : 5.894393e-06           | 5.894393e-06           | 0.0 
   SpeciesConcVV_AONITA : 2.7138583e-07          | 2.7138583e-07          | 0.0 
   SpeciesConcVV_AROMP4 : 1.7938361e-09          | 1.7938361e-09          | 0.0 
   SpeciesConcVV_AROMP5 : 1.0746459e-09          | 1.0746459e-09          | 0.0 
   SpeciesConcVV_AROMRO2 : 4.7994303e-10          | 4.7994303e-10          | 0.0 
   SpeciesConcVV_ASOA1  : 3.987758e-08           | 3.987758e-08           | 0.0 
   SpeciesConcVV_ASOA2  : 1.0983177e-08          | 1.0983177e-08          | 0.0 
   SpeciesConcVV_ASOA3  : 3.7467963e-08          | 3.7467963e-08          | 0.0 
   SpeciesConcVV_ASOAN  : 2.9784314e-07          | 2.9784314e-07          | 0.0 
   SpeciesConcVV_ASOG1  : 8.251855e-08           | 8.251855e-08           | 0.0 
   . . .
		
as well as optional :ref:`six-panel plots <six-panel>`.  This
allows you to determine if two GEOS-Chem simulations have yielded
identical results or not.

.. _comp-diags-code:

===========
Source code
===========

**Script location:** `gcpy/examples/diagnostics/compare_diags.py
<https://github.com/geoschem/gcpy/blob/main/gcpy/examples/plotting/plot_comparisons.py>`_

**Related configuration file:**
`gcpy/examples/diagnostics/compare_diags.yml
<https://github.com/geoschem/gcpy/blob/main/gcpy/examples/diagnostics/compare_diags.yml>`_ 

.. _comp-diags-usage:

=====
Usage
=====

Make sure that you :ref:`specified the proper Matplotlib backend
<mpl-backend>` for  your system.

First, copy the :file:`compare_diags.yml` file to your local folder.

.. code-block:: console

   $ cp /path/to/GCPy/gcpy/examples/diagnostics/compare_diags.yml .

.. tip::

   You can rename your copy of the file if you wish.  This may be
   useful if you intend to do multiple comparisons.
   
Next, customize the :file:`compare_diags.yml` file so that it contains
the proper directory paths to your GEOS-Chem output files.  You can
also decide whether or not to create the optional single-level and
zonal mean plots.

.. code-block:: yaml

   ---
   paths:
     main_dir: /path/to/your/data   # Add the path to your output here
     plots_dir: ./Results
     weights_dir: /path/to/regridding/weights/folder
   
   data:
     ref:
       label: "GCC_ref"
       dir: GCC_ref
       subdir: OutputDir
       file: GEOSChem.SpeciesConc.20190701_0000z.nc4
     dev:
       label: "GCC_dev"
       dir: GCC_dev
       subdir: OutputDir
       file: GEOSChem.SpeciesConc.20190701_0000z.nc4
   
   options:
     verbose: False
     restrict_vars: []
     level_plot:
       create_plot: True
       pdfname: single_level_comparison.pdf
       level_to_plot: 0
     zonal_mean:
       create_plot: True
       pdfname: zonal_mean_comparison.pdf
     totals_and_diffs:
       create_table: True
       diff_type: absdiff             # Values: percent, pctdiff, %, abs, absdiff
       print_to_screen: True
       filename: ''
       skip_small_diffs: True
       small_diff_threshold: 0.0000
     n_cores: -1
		
Then, run the script with:

.. code-block:: console

   $ python -m gcpy.examples.diagnostics.compare_diags compare_diags.yml
