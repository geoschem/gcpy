Tabling
========

This page describes the tabling capabilities of GCPy, including possible argument values for every tabling function.
These functions are primarily used for model benchmarking purposes. All tables are printed to text files.


Emissions tables
----------------

.. code-block:: python


    def make_benchmark_emis_tables(reflist, refstr, devlist,
        devstr, dst="./benchmark", refmet=None, devmet=None,
        overwrite=False, ref_interval=[2678400.0], dev_interval=[2678400.0],
        spcdb_dir=os.path.dirname(__file__)
    ):
        """
        Creates a text file containing emission totals by species and
        category for benchmarking purposes.
        """

Arguments:
~~~~~~~~~~

.. option:: reflist: list of str
  
   List with the path names of the emissions file or files
   (multiple months) that will constitute the "Ref"
   (aka "Reference") data set.

.. option:: refstr : str

   A string to describe ref (e.g. version number)

.. option:: devlist : list of str

   List with the path names of the emissions file or files
   (multiple months) that will constitute the "Dev"
   (aka "Development") data set

.. option:: devstr : str

   A string to describe dev (e.g. version number)


Keyword arguments:
~~~~~~~~~~~~~~~~~~

.. option:: dst : str

   A string denoting the destination folder where the file
   containing emissions totals will be written.

   Default value: ./benchmark

.. option:: refmet : str

   Path name for ref meteorology

   Default value: None

.. option:: devmet : str

   Path name for dev meteorology  

   Default value: None

.. option:: overwrite : bool

   Set this flag to True to overwrite files in the
   destination folder (specified by the dst argument).

   Default value: False

.. option:: ref_interval : list of float

   The length of the ref data interval in seconds. By default, interval
   is set to [2678400.0], which is the number of seconds in July
   (our 1-month benchmarking month).

   Default value: [2678400.0]

.. option:: dev_interval : list of float

   The length of the dev data interval in seconds. By default, interval
   is set to [2678400.0], which is the number of seconds in July
   (our 1-month benchmarking month).

   Default value: [2678400.0]

.. option:: spcdb_dir : str

   Directory of species_datbase.yml file

   Default value: Directory of GCPy code repository


``gcpy.benchmark.make_benchmark_emis_tables()`` generates tables of total emissions categorized by species or by inventory.
These tables contain total global emissions over the lengths of the Ref and Dev datasets, as well as the differences between
totals across the two datasets. Passing a list of datasets as Ref or Dev (e.g. multiple months of emissions files) will result
in printing totals emissions summed across all files in the list. Make sure to update the ``ref_interval`` and/or ``dev_interval``
arguments if you pass input that does not correspond with 1 31 day month.   


Mass Tables
-----------

.. code-block:: python

   def make_benchmark_mass_tables(ref, refstr, dev, devstr,
      varlist=None, dst="./benchmark", subdst=None, overwrite=False,
      verbose=False, label="at end of simulation", spcdb_dir=os.path.dirname(__file__),
      ref_met_extra='', dev_met_extra=''
   ):
      """
      Creates a text file containing global mass totals by species and
      category for benchmarking purposes.
      """

Arguments:
~~~~~~~~~~

.. option:: reflist : str

   Pathname that will constitute
   the "Ref" (aka "Reference") data set.

.. option:: refstr : str

   A string to describe ref (e.g. version number)

.. option:: dev : list of str

   Pathname that will constitute
   the "Dev" (aka "Development") data set.  The "Dev"
   data set will be compared against the "Ref" data set.

.. option:: devstr : str

   A string to describe dev (e.g. version number)


Keyword arguments:
~~~~~~~~~~~~~~~~~~

.. option:: varlist : list of str

   List of variables to include in the list of totals.
   If omitted, then all variables that are found in either
   "Ref" or "Dev" will be included.  The varlist argument
   can be a useful way of reducing the number of
   variables during debugging and testing.

   Default value: None

.. option:: dst : str

   A string denoting the destination folder where the file
   containing emissions totals will be written.

   Default value: ./benchmark

.. option:: subdst : str

   A string denoting the sub-directory of dst where PDF
   files containing plots will be written.  In practice,
   subdst is only needed for the 1-year benchmark output,
   and denotes a date string (such as "Jan2016") that
   corresponds to the month that is being plotted.

   Default value: None

.. option:: overwrite : bool

   Set this flag to True to overwrite files in the
   destination folder (specified by the dst argument).

   Default value: False

.. option:: verbose : bool

   Set this flag to True to print extra informational output.

   Default value: False.

.. option:: spcdb_dir : str

   Directory of species_datbase.yml file

   Default value: Directory of GCPy code repository

.. option:: ref_met_extra : str

   Path to ref Met file containing area data for use with restart files
   which do not contain the Area variable.
   Default value : ''

.. option:: dev_met_extra : str

   Path to dev Met file containing area data for use with restart files
   which do not contain the Area variable.

   Default value: ''         


``gcpy.benchmark.make_benchmark_mass_tables`` is used to create global mass tables of GEOS-Chem species from a ``Restart`` file.
This function will create one table of total mass by species from the earth's surface to the top of the stratosphere and one table for only the troposphere.
The tables contain total mass for each of the ref and dev datasets in Gg, as well as absolute and percentage difference between the two datasets.
If your restart files do not contain an Area variable (``"AREA`` for GEOS-Chem Classic or ``"Met_AREAM2`` for GCHP) then you will need to use the 
``ref_met_extra`` and/or ``dev_met_extra`` arguments to pass the paths of NetCDF files containing the corresponding area variables 
(usually contained in meteorology diagnostic output).


Operations Budget Tables
------------------------

.. code-block:: python

   def make_benchmark_operations_budget(refstr, reffiles, devstr,
      devfiles, ref_interval, dev_interval, benchmark_type=None,
      label=None, col_sections=["Full", "Trop", "PBL", "Strat"],
      operations=["Chemistry","Convection","EmisDryDep","Mixing",
   "Transport","WetDep"], compute_accum=True,
      require_overlap=False, dst='.', species=None, overwrite=True
   ):
      """
      Prints the "operations budget" (i.e. change in mass after
      each operation) from a GEOS-Chem benchmark simulation.
      """


Arguments:
~~~~~~~~~~

.. option:: refstr : str

   Labels denoting the "Ref" versions

.. option:: reffiles : list of str

   Lists of files to read from the "Ref" version.

.. option:: devstr : str

   Labels denoting the "Dev" versions

.. option:: devfiles : list of str

   Lists of files to read from "Dev" version.

.. option:: interval : float

   Number of seconds in the diagnostic interval.


Keyword arguments:

.. option:: benchmark_type : str

   "TransportTracersBenchmark" or "FullChemBenchmark".

   Default value: None

.. option:: label : str

   Contains the date or date range for each dataframe title.

   Default value: None

.. option:: col_sections : list of str

   List of column sections to calculate global budgets for. May
   include Strat eventhough not calculated in GEOS-Chem, but Full
   and Trop must also be present to calculate Strat.

   Default value: ["Full", "Trop", "PBL", "Strat"]

.. option:: operations : list of str

   List of operations to calculate global budgets for. Accumulation
   should not be included. It will automatically be calculated if
   all GEOS-Chem budget operations are passed and optional arg
   compute_accum is True.

   Default value: ["Chemistry","Convection","EmisDryDep",
          "Mixing","Transport","WetDep"]

.. option:: compute_accum : bool

   Optionally turn on/off accumulation calculation. If True, will 
   only compute accumulation if all six GEOS-Chem operations budgets
   are computed. Otherwise a message will be printed warning that
   accumulation will not be calculated.

   Default value: True

.. option:: require_overlap : bool

   Whether to calculate budgets for only species that are present in
   both Ref or Dev.

   Default value: False

.. option:: dst : str

   Directory where plots & tables will be created.

   Default value: '.' (directory in which function is called)

.. option:: species : list of str

   List of species for which budgets will be created.

   Default value: None (all species)

.. option:: overwrite : bool

   Denotes whether to overwrite existing budget file.

   Default value: True



``gcpy.benchmark.make_benchmark_operations_budget()`` creates tables of budgets for species separated by model operation.
The tables show budgets for each of the ref and dev datasets in Gg, as well as absolute and percentage difference between the two datasets.
Note that total accumulation across all operations will only be printed if you set ``compute_accum==True`` and
all operations are included in ``operations``. Note also that when using the non-local mixing scheme (default), ``'Mixing'`` 
includes emissions and dry deposition applied below the PBL. ``'EmisDryDep'`` therefore only captures fluxes above the PBL.
When using full mixing, ``'Mixing'`` and ``'EmisDryDep'`` are fully separated.

Aerosol Budgets and Burdens
---------------------------

.. code-block:: python

   def make_benchmark_aerosol_tables(devdir, devlist_aero, devlist_spc,
      devlist_met, devstr, year, days_per_mon, dst='./benchmark',
         overwrite=False, is_gchp=False, spcdb_dir=os.path.dirname(__file__)
   ):
      """
      Compute FullChemBenchmark aerosol budgets & burdens
      """


Arguments:
~~~~~~~~~~
        
.. option:: devdir: str

   Path to development ("Dev") data directory

.. option:: devlist_aero : list of str

   List of Aerosols collection files (different months)

.. option:: devlist_spc : list of str

   List of SpeciesConc collection files (different months)

.. option:: devlist_met : list of str

   List of meteorology collection files (different months)

.. option:: devstr : str

   Descriptive string for datasets (e.g. version number)

.. option:: year : str

   The year of the benchmark simulation (e.g. '2016'). 

.. option:: days_per_mon : list of int

   List of number of days per month for all months


Keyword arguments:
~~~~~~~~~~~~~~~~~~

.. option:: dst : str 

   Directory where budget tables will be created.

   Default value: './benchmark'         

.. option:: overwrite : bool

   Overwrite burden & budget tables? (default=True)

   Default value: False

.. option:: is_gchp : bool

   Whether datasets are for GCHP

   Default value: False

.. option:: spcdb_dir : str

   Directory of species_datbase.yml file

   Default value: Directory of GCPy code repository


``gcpy.benchmark.make_benchmark_aerosol_tables()`` generates two different tables using output from a single dataset. One contains annual mean aerosol burdens in Tg in the stratosphere, 
troposphere, and combined stratosphere and troposphere. The other table shows annual global mean AOD in the stratosphere, troposphere, and combined
stratosphere and troposphere. Aerosol species used are pre-defined in ``aod_species.yml``: BCPI, OCPI, SO4, DST1, SALA, and SALC.

