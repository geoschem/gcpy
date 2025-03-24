.. |br| raw:: html

   <br/>

.. _release-guide:

##################################
Releasing new versions (GCST only)
##################################

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

#. Update branches accordingly:

   .. code-block:: console

      $ git checkout main
      $ git merge dev                             # Merge dev into main
      $ git checkout dev
      $ git merge main                            # Update dev with main
      $ git checkout docs/dev
      $ git merge main
      $ git checkout main                         # Update docs/dev with main
      $ git tag X.Y.Z                             # Create tag
      $ git push origin main dev docs/dev X.Y.Z   # Push everything to origin

#. Publish the release on Github. |br|
   |br|

#. A GitHub Action will push the :literal:`geoschem-gcpy` package to
   the :program:`Python Package Index (PyPi)`. |br|
   |br|

#. Verify the new release is visible at
   https://pypi.org/project/geoschem-gcpy/.  This may take a few
   minutes). |br|
   |br|

#. After a period of time (around an hour), you will be notified of a
   new PR at https://github.com/conda-forge/geoschem-gcpy-feedstock
   indicating conda-forge has detected a new release on PyPi. You
   should be able to merge this PR without any additinal interference
   once all checks have passed. |br|
   |br|

#. If any checks on the feedstock should fail, this typically
   indicates Python package incompatibilities.  Follow these steps to
   resolve these:

   #. Clone a local copy of
      https://github.com/geoschem/geoschem-gcpy-feedstock and check
      out the feedstock PR's branch. |br|
      |br|

   #. Look at the errors on the PR's GitHub Actions tab and the Azure
      Pipelines page.  You should be able to figure out which package
      is the one causing the issue. |br|
      |br|

   #. Edit the Python packages in the file :file:`recipes/meta.yaml`
      and also increment the build number setting by one.  If packages
      have been pegged to a specific version, try removing the
      peg. |br|
      |br|

   #. Push your updates to the feedstock PR branch and wait for all
      checks to pass. |br|
      |br|

   #. If checks are still failing, repeat steps 2-4 above and then
      wait to see if the checks pass. |br|
      |br|

   #. Merge the feedstock PR manually after all checks have passed.
      |br|

#. Once the feedstock PR has been merged and after another period of
   waiting, you should see builds for the new release when running
   :code:`conda search -f geoschem-gcpy` or by checking
   https://anaconda.org/conda-forge/geoschem-gcpy.  This indicates the
   new version is publicly available for installation through
   conda-forge.
