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

You can use GCPy module :file:`gcpy/kpp/kppsa_quick_look.py` to
generate a vertical profile for a given species at a given observation
site.  This can be useful as a quick-look sanity check.

First, activate your Python virtual environment

.. code-block:: console

    $ mamba activate gcpy_env   # If using Mamba, or

    $ conda activate gcpy_env   # If using Conda

Then run the :file:`kppsa_quick_look.py` module from the command line:

.. code-block:: console

    (gcpy_env) $ python -m gcpy.kpp.kppsa_quick_look \
     --dirname /path/to/KPP-Standalone/output        \
     --label   Rosenbrock                            \
     --pattern Beijing*20190701_0040.log             \
     --species O3

This will create the following plot:

.. image:: _static/images/kppsa\_quick\_look.png
   :align: center

.. _kppsa-usage-vert-prof:

Vertical profiles at multiple observation sites
-----------------------------------------------

You can use GCPy module :file:`gcpy/kpp/kppsa_plot_sites.py` to
generate vertical profile plots for a given species at multiple
observation locations.

First, activate your Python virtual environment

.. code-block:: console

    $ mamba activate gcpy_env   # If using Mamba, or

    $ conda activate gcpy_env   # If using Conda

Then run the :file:`kppsa_plot_sites.py` module from the command line:

.. code-block:: console

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
   :widths: 180 150 220

   * - Function
     - In file
     - Description
   * - :ref:`kppsa_make_quick_look_plot <kppsa-ref-quick-look>`
     - :file:`kppsa_quick_look.py`
     - Creates a "quick-look" vertical profile plot.  Useful as a
       sanity check.
   * - :ref:`kppsa_plot_species_at_sites <kppsa-ref-sites>`
     - :file:`kppsa_plot_sites.py`
     - Creates vertical profile plots of a given species at various
       locations.
   * - :ref:`kppsa_get_file_list <kppsa-ref-file-list>`
     - :file:`kppsa_utils.py`
     - Returns a list of KPP-Standalone log files matching a search
       criteria.
   * - :ref:`kppsa_read_one_csv_file <kppsa-ref-one-csv>`
     - :file:`kppsa_utils.py`
     - Reads a single log file (in CSV format) from the KPP
       standalone box model into a pandas.DataFrame object.
   * - :ref:`kppsa_read_csv_files <kppsa-ref-csv>`
     - :file:`kppsa_utils.py`
     - Reads all KPP standalone log files for a given site
       in a given directory.
   * - :ref:`kppsa_prepare_site_data <kppsa-ref-prepare>`
     - :file:`kppsa_utils.py`
     - Returns a pd.DataFrame object containing data for a given
       species, and observation site, as well as the corresponding
       top-of-plot title.
   * - :ref:`kppsa_plot_single_site <kppsa-ref-one-site>`
     - :file:`kppsa_utils.py`
     - Plots observation data vs. model data at a single station
       site.
   * - :ref:`kppsa_plot_one_page <kppsa-ref-one-page>`
     - :file:`kppsa_utils.py`
     - Plots a single page of models vs. observations.
   * - :ref:`kppsa_get_unique_site_names <kppsa-ref-unique>`
     - :file:`kppsa_utils.py`
     - Returns a list of unique sites where KPP-Standalone box model
       output has been archived.

.. _kppsa-ref-quick-look:


Function :file:`kppsa_make_quick_look_plot`
-------------------------------------------

.. code-block:: python

   def kppsa_make_quick_look_plot(file_list, label, species):
       """
       Creates a quick-look plot from KPP-Standalone box model output.

       Args
       file_list : list : List of KPP-Standalone log files
       site_name : str  : Name of the site that you wish to plot
       label     : str  : Descriptive label for the data
       species   : str  : Name of the species that you wish to plot
       """

.. _kppsa-ref-sites:

Function :code:`kppsa_plot_species_at_sites`
--------------------------------------------

.. code-block:: python

   def kppsa_plot_species_at_sites(
           ref_file_list,
           ref_label,
           dev_file_list,
           dev_label,
           species,
           pdfname,
   ):
       """
       Creates vertical profile plots of a given species
       from KPP-Standalone box model output.

       Args
       ref_file_list : list : KPP-Standalone log files for "Ref" version
       ref_label     : str  : Label for the "Ref" version
       dev_file_list : list : KPP-Standalone log files for "Dev" version
       dev_label     : str  : Label for the "Dev" version
       species       : str  : Name of the species to plot
       pdfname       : str  : Name of the output PDF file
       """

.. _kppsa-ref-file-list:

Function :code:`kppsa_get_file_list`
------------------------------------

.. code-block:: python

   def kppsa_get_file_list(
           input_dir,
           pattern=""
   ):
       """
       Returns a list of KPP-Standalone log files matching
       a search criteria.

       Args
       input_dir : str  : Directory with KPP-Standalone log files
       pattern   : str  : Read files matching this pattern (Default = "")

       Returns
       file_list : list : List of files matching the criteria
       """

.. _kppsa-ref-one-csv:

Function :code:`kppsa_read_one_csv_file`
----------------------------------------

.. code-block:: python

   def kppsa_read_one_csv_file(file_name):
       """
       Reads a single log file (in CSV format) from the KPP
       standalone box model into a pandas.DataFrame object.
   
       Args
       file_name : str          : File to be read
   
       Returns
       dframe    : pd.DataFrame : DataFrame with the results
       """

.. _kppsa-ref-csv:

Function :code:`kppsa_read_csv_files`
-------------------------------------

.. code-block:: python

   def kppsa_read_csv_files(file_list):
       """
       Reads all KPP standalone log files for a given site
       in a given directory.

       Args
       input_dir  : str          : Directory to search
       site       : str          : KPP standalone site name

       Returns
       dframe_all : pd.DataFrame : Observations at all levels
       """

.. _kppsa-ref-prepare:

Function :code:`kppsa_prepare_site_data`
----------------------------------------

.. code-block:: python

   def kppsa_prepare_site_data(
           dframe,
           site_name,
           species,
   ):
       """
       Returns a pd.DataFrame object containing data for a given species,
       and observation site, as well as the corresponding top-of-plot
       title.  Species data is limited from the surface to 500 hPa.

       Args
       dframe     : pd.DataFrame : KPP-Standalone output data
       site_name  : str          : Name of site to plot
       species    : species      : Name of species to plot

       Returns
       site_data  : pd.DataFrame : Data for the given site & species
       site_title : str          : Corresponding plot title string
       """

.. _kppsa-ref-one-site:

Function :code:`kppsa_plot_single_site`
---------------------------------------

.. code-block:: python

   def kppsa_plot_single_site(
           fig,
           rows_per_page,
           cols_per_page,
           subplot_index,
           subplot_title,
           ref_data,
           ref_label,
           dev_data,
           dev_label,
           species,
           font_scale,
   ):
       """
       Plots observation data vs. model data at a single station site.

       Args:
       fig            : mpl.figure.Figure : Figure object for the plot
       rows_per_page  : int               : # of rows to plot on a page
       cols_per_page  : int               : # of columns to plot on a page
       subplot_index  : int               : Index of each subplot
       subplot_title  : str               : Title for each subplot
       ref_data       : pd.DataFrame      : Observations at each station site
       ref_label      : str               : Label for the Ref model data
       dev_data       : pd.DataFrame      :
       dev_label      : str               : Label for the Dev model data
       site_name      : str               : Name of the station site
       species        : pd.Series         : Data from the Ref model version
       font_scale     : float             : Scale fac to increase font size
       """

.. _kppsa-ref-one-page:

Function :code:`kppsa_plot_one_page`
------------------------------------

.. code-block:: python

   def kppsa_plot_one_page(
           pdf,
           site_names,
           ref_dframe,
           ref_label,
           dev_dframe,
           dev_label,
           species="O3",
           rows_per_page=3,
           cols_per_page=3,
           font_scale=1.0,
   ):
       """
       Plots a single page of models vs. observations.

       Args:
       pdf             : pdf          : PDF object
       ref_dframe      : pd.DataFrame : Observations at each station site.
       ref_label       : str          : Label for the observational data
       dev_dframe      : pd.DataFrame : Data from the Ref model version
       dev_label       : str          : Label for the Ref model data
       species         : str          : Name of the species to plot
       dev_dataarray   : xr.DataArray : Data from the Dev model version
       dev_label       : str          : Label for the Dev model data
       dev_cs_grid     : str|None     : Metadata for Dev cubed-sphere grid
       gc_levels       : pd.DataFrame : Metadata for model vertical levels
       rows_per_page   : int          : Number of rows to plot on a page
       cols_per_page   : int          : Number of cols to plot on a page
       font_scale      : float        : PDF output file name
       """

.. _kppsa-ref-unique:

Function :code:`kppsa_get_unique_site_names`
--------------------------------------------

.. code-block:: python

   def kppsa_get_unique_site_names(dframe):
       """
       Returns a list of unique sites where KPP-Standalone box model
       output has been archived.
   
       Args
       dframe     : pd.DataFrame : Object containing KPP-Standalone output
   
       Returns
       site_names : list of str  : List of unique site names
       """
