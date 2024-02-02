.. |br| raw:: html

   <br/>

.. _release-guide:

######################
Releasing new versions
######################

This page describes some of the steps required for releasing new
versions of GCPy on Github, PyPi, and conda-forge. 

#. For clarity, update version numbers to the new release in the
   following locations:

   - :file:`setup.py`
   - :file:`gcpy/_version.py`
   - :file:`docs/source/conf.py`
   - :file:`gcpy/benchmark/run_benchmark.py`
   - :file:`gcpy/benchmark/modules/run_1yr_fullchem_benchmark.py`
   - :file:`gcpy/benchmark/modules/run_1yr_tt_benchmark.py`
     |br|

#. Update :file:`CHANGELOG.md` |br|
   |br|

   
#. Merge **dev** into **main** |br|
   |br|
   
#. Publish the release on Github. |br|
   |br|
   
#. Install :code:`twine` using :code:`pip install twine` (if you
   haven't done this before). |br|
   |br|
   
#. To package GCPy for publication to PyPi, run the following from the
   root of your local GCPy repository:

   .. code-block:: console
     
      $ conda activate gcpy_env   # or whatever your conda env is named
      $ python setup.py sdist bdist_wheel
      $ twine check dist/*
      $ run twine upload --repository-url https://test.pypi.org/legacy/ dist/*

   Enter your login credentials for :file:`test.pypi.org` as
   requested. Publishing to test.pypi ensures there are no issues with
   packaging the new release before publication to the primary
   PyPi database. |br|
   |br|

#. Publish to PyPi by running :code:`run twine upload dist/*`, and enter
   your login information for pypi.org as requested. |br|
   |br|

#. Verify the new release is visible at
   https://pypi.org/project/geoschem-gcpy/ (may take a few
   minutes). |br|
   |br|

#. After a period of time (around an hour), you will be notified of a
   new PR at https://github.com/conda-forge/geoschem-gcpy-feedstock
   indicating conda-forge has detected a new release on PyPi. You
   should be able to merge this PR without any additinal interference
   once all checks have passed. |br|
   |br|

#. Once the feedstock PR has been merged and after another period of
   waiting, you should see builds for the new release when running
   :code:`conda search -f geoschem-gcpy`.  This indicates the new
   version is publicly available for installation through
   conda-forge. 
