.. |br| raw:: html

   <br/>

.. _kppsa:

###########################################
Visualizing KPP-Standalone box model output
###########################################

This page describes how you can plot results from the `KPP-Standalone
box model
<https://geos-chem.readthedocs.io/en/stable/geos-chem-shared-docs/supplemental-guides/using-kpp-standalone.html>`_
with GCPy.

.. _kppsa-usage:

=========================================
Creating plots from KPP-Standalone output
=========================================

.. _kppsa-usage-quick-look:

Quick-look plots
----------------

You can use GCPy :func:`gcpy.kpp.kppsa_quick_look` to
generate a vertical profile for a given species at a given observation
site.  This can be useful as a quick-look sanity check.

Type this at the command line:

.. code-block:: console

    $ conda activate gcpy_env
    (gcpy_env) $ python -m gcpy.kpp.kppsa_quick_look \
     --dirname /path/to/KPP-Standalone/output        \
     --label   Rosenbrock                            \
     --pattern Beijing*20190701_0040.log             \
     --species O3

which will generate this plot:

.. image:: _static/images/kppsa\_quick\_look.png
   :align: center

.. _kppsa-usage-vert-prof:

Vertical profiles at multiple observation sites
-----------------------------------------------

You can use GCPy module :func:`gcpy.kpp.kppsa_plot_sites` to
generate vertical profile plots for a given species at multiple
observation locations.

Type this at the command line:

.. code-block:: console

    conda activate gcpy_env
    (gcpy_env) $ python -m gcpy.kpp.kppsa_plot_sites  \
     --refdir   /path/to/KPP-Standalone/Ref/log/files \
     --reflabel Rosenbrock                            \
     --devdir   /path/to/KPP-Standalone/Dev/log/files \
     --devlabel Backwards Euler                       \
     --pattern  20190701_0040.log                     \
     --species  O3                                    \
     --pdfname  KPP-Standalone-O3-20190701-0040.pdf

This will create a PDF file with several pages.  Each page will look
similar to this:

.. image:: _static/images/kppsa\_plot\_sites.png
   :align: center

In some instances (such as in the plot above), the Dev simulation (blue line)
will overlap the Ref simulation (red line).   This indicates that Dev
has produced identical or nearly-identical results to Ref.

.. _kppsa-ref:

=====================
Source code reference
=====================

.. _kppsa-ref-list:

List of functions
-----------------

.. list-table:: Functions for reading and plotting KPP-Standalone box
		model output
   :header-rows: 1
   :align: center

   * - Function
     - Description
   * - :func:`gcpy.kpp.kppsa_quick_look`
     - Creates a "quick-look" vertical profile plot.  Useful as a
       sanity check.
   * - :func:`gcpy.kpp.kppsa_plot_sites`
     - Creates vertical profile plots of a given species at various
       locations.
   * - :func:`gcpy.kpp.kpp_utils.kppsa_get_file_list`
     - Returns a list of KPP-Standalone log files matching a search
       criteria.
   * - :func:`gcpy.kpp.kppsa_utils.kppsa_read_one_csv_file`
     - Reads a single log file (in CSV format) from the KPP
       standalone box model into a pandas.DataFrame object.
   * - :func:`gcpy.kpp.kppsa_utils.kppsa_read_csv_files`
     - Reads all KPP standalone log files for a given site
       in a given directory.
   * - :func:`gcpy.kpp.kppsa_utils.kppsa_prepare_site_data`
     - Returns a pd.DataFrame object containing data for a given
       species, and observation site, as well as the corresponding
       top-of-plot title.
   * - :func:`gcpy.kpp.kppsa_utils.kppsa_plot_single_site`
     - Plots observation data vs. model data at a single station
       site.
   * - :func:`gcpy.kpp.kppsa_utils.kppsa_plot_one_page`
     - Plots a single page of models vs. observations.
   * - :func:`gcpy.kpp.kppsa_utils.kppsa_get_unique_site_names`
     - Returns a list of unique sites where KPP-Standalone box model
       output has been archived.
