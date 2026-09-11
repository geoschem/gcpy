.. _gcst-examples:

##############################################
GEOS-Chem Support Team (GCST) utility scripts
##############################################

The :file:`gcpy/examples/gcst` folder contains standalone utility
scripts that the `GEOS-Chem Support Team
<https://geoschem.github.io/support-team.html>`_ (GCST) uses for
routine tasks such as building `GCHP <https://gchp.readthedocs.io>`_
configuration file entries and summarizing integration test results.
This page describes how to use each of these scripts.

.. _gcst-speciesconcvv:

==================================================
Generate SpeciesConcVV entries for GCHP HISTORY.rc
==================================================

All diagnostic quantities must be individually specified in the `GCHP
<https://gchp.readthedocs.io>`_' :file:`HISTORY.rc` run directory
configuration file.  Typing these manually is both tedious and error
prone.  Instead, you can generate this automatically from the list of
species indices that is written out to the GEOS-Chem log file.

.. _gcst-speciesconcvv-source:

Source code
-----------

.. list-table::
   :header-rows: 1

   * - Description
     - Script location
   * - :mod:`gcpy.examples.gcst.generate_gchp_speciesconcvv_list`
     - `gcpy/examples/gcst/generate_gchp_speciesconcvv_list.py <https://github.com/geoschem/gcpy/blob/main/gcpy/examples/gcst/generate_gchp_speciesconcvv_list.py>`_

.. _gcst-speciesconcvv-usage:

Usage
-----

First, run GCHP and send output to a log file (such as
:file:`gchp.YYYYMMDD_hhmmss.log`, where :file:`YYYYMMDD` and
:file:`hhmmss` are the date and time at the start of your
simulation. This log file will contains a section that lists each
GEOS-Chem species and its correspoinding indices.  For example:

.. code-block:: text

   ===============================================================================
   SPECIES NAMES AND INDICES

   Name               ModelId  DryDepId  WetDepId  PhotolId HygGrthId  KppSpcId
   -------------------------------------------------------------------------------
   ACET                     1         1         -         1         -       299
   ACTA                     2         2         1         -         -       286
   ACR                      3         3         -         -         -       236
   ...
   O2                     394         -         -       150         -       360
   ===============================================================================

Next, activate your GCPy environment with:

.. code-block:: console

   $ conda activate gcpy_env

Then run script :mod:`gcpy.examples.gcst.generate_gchp_speciesconcvv_list`
and pass it the log file and the name of a file where output will be sent:

.. code-block:: console

   (gcpy_env) $ python -m gcpy.examples.gcst.generate_gchp_speciesconcvv_list \
                gchp.YYYYMMDD_hhmmss.log \
                speciesconcvv_entries.txt

The output file will contain one comma-terminated entry per species,
in the order expected in the GCHP :file:`HISTORY.rc` file (reversed
alphabetical order):

.. code-block:: text

                                'SpeciesConcVV_O2             ', "GCHPchem",
                                'SpeciesConcVV_N2             ', "GCHPchem",
                                'SpeciesConcVV_H2             ', "GCHPchem",
                                'SpeciesConcVV_OH             ', "GCHPchem",
                                'SpeciesConcVV_HO2            ', "GCHPchem",
                                ...

You can then paste these lines directly into your GCHP
:file:`HISTORY.rc` file.

Deactivate the python environment when you are finished.

.. code-block:: console

   (gcpy_env) $ conda deactivate

.. _gcst-inttest-report:

===================================
Generate an integration test report
===================================

GEOS-Chem Classic and GCHP integration tests produce an execution log
(listing which tests passed, failed, or are still pending) and a diff
log (listing which tests produced output that differs from a
previous integration test).  Reading through both logs by hand to
compile a summary is time-consuming.

Script :mod:`gcpy.examples.gcst.generate_inttest_report` parses both
logs and produces a concise, plain-text (GitHub-Markdown-friendly)
report that summarizes the execution and diff results, suitable for
pastingdirectly into a pull request or issue comment. 

.. _gcst-inttest-report-source:

Source code
-----------

.. list-table::
   :header-rows: 1

   * - Description
     - Script location
   * - :mod:`gcpy.examples.gcst.generate_inttest_report`
     - `gcpy/examples/gcst/generate_inttest_report.py <https://github.com/geoschem/gcpy/blob/main/gcpy/examples/gcst/generate_inttest_report.py>`_

.. _gcst-inttest-report-usage:

Usage
-----

Navigate to the integration test folder.  Then run the difference test
script and pipe its output to a log.

.. code-block:: console

   $ cd /path/to/integration/test/folder
   $ CodeDir/test/difference/diffTest.sh ../path/to/previous/int/test . | tee diff.log

Next, activate your GCPy environment:

.. code-block:: console

   $ conda activate gcpy_env

Then run :mod:`gcpy.examples.gcst.generate_inttest_report`, supplying the
execution-test log, the diff-test log, and a label identifying the
reference version that the tests were compared against:

.. code-block:: console

   (gcpy_env) $ python -m gcpy.examples.gcst.generate_inttest_report \
       --exec-log  logs/results.execute.log                          \
       --diff-log  diff.log                                          \
       --ref-label "Name of previous integration test"               \
       --output    report.txt

:literal:`--exec-log` and :literal:`--diff-log` are required and
accept either GEOS-Chem Classic or GCHP log formats (the model type is
auto-detected from the execution log).

:literal:`--ref-label` is required and should match the
reference-version label as it appears in the diff-log file
paths (e.g. :literal:`gcc.X.Y.Z-alpha.N` or
:literal:`gchp.X.Y.Z-alpha.N`).

:literal:`--output` is optional; if
omitted, the report is printed to :literal:`stdout` instead of being
written to a file.

If all execution tests passed and no differences were found, the
report is a compact confirmation:

.. code-block:: text

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%  All execution tests passed!  %%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   All GEOS-Chem Classic tests were zero-diff w/r/t <ref-label>

Otherwise, the report lists the failed/pending tests and the tests
with differences, along with the verbatim header block and failure
lines from the execution log.

You may now deactivate the GCPy Python environment:

.. code-block:: console

   $ conda deactivate

Run the script once per model type (GC-Classic and GCHP) and paste the
output from each model into a GitHub issue or pull request.
