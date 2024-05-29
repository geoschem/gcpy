.. |br| raw:: html

   <br/>

.. _release-guide:

######################
Releasing new versions
######################

This page describes some of the steps required for releasing new
versions of GCPy on Github, PyPi, and conda-forge. 

#. Update :file:`CHANGELOG.md` as necessary. |br|
   |br|

#. For clarity, update version numbers to the new release
   (:literal:`X.Y.Z`) in the following locations:

   - :file:`CHANGELOG.md`
   - :file:`setup.py`
   - :file:`gcpy/_version.py`
   - :file:`docs/source/conf.py`
   - :file:`gcpy/benchmark/run_benchmark.py`
   - :file:`gcpy/benchmark/modules/run_1yr_fullchem_benchmark.py`
   - :file:`gcpy/benchmark/modules/run_1yr_tt_benchmark.py`

   This can be done by running the utility script:

   .. code-block:: console

      $ cd /path/to/GCPy/.release
      $ ./changeVersionNumbers.sh X.Y.Z
   
#. Merge **dev** into **main** |br|
   |br|
   
#. Publish the release on Github. |br|
   |br|

#. A GitHub Action will push the :literal:`geoschem-gcpy` package to
   the :program:`Python Package Index (PyPi)`.
