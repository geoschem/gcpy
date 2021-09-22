Releasing new versions
======================

This page describes some of the steps required for releasing new versions of GCPy on Github, PyPi, and conda-forge.



   1. For clarity, update version numbers to the new release in ``setup.py``, _version.py, and the benchmark scripts.
   2. Update CHANGELOG.md
   3. Merge main with dev
   4. Publish the release on Github.
   5. Install ``twine`` using ``pip install twine``.
   6. To package GCPy for publication to PyPi, run the following from the root of your local GCPy repository::

        run python setup.py sdist bdist_wheel
        run twine check dist/*
        run twine upload --repository-url https://test.pypi.org/legacy/ dist/*

      Enter your login credentials for test.pypi.org as requested. Publishing to test.pypi ensures there are no issues with packaging the new release
      before publication to the primary PyPi database.


	
   7. Publish to PyPi by running ``run twine upload dist/*``, and enter your login information for pypi.org as requested.
   8. Verify the new release is visible at https://pypi.org/project/geoschem-gcpy/ (may take a few minutes).
   9. After a period of time (around an hour), you will be notified of a new PR at https://github.com/conda-forge/geoschem-gcpy-feedstock indicating conda-forge has 
      detected a new release on PyPi. You should be able to merge this PR without any additinal interference once all checks have passed.
   10. Once the feedstock PR has been merged and after another period of waiting, you should see builds for the new release when running ``conda search -f geoschem-gcpy``.
      This indicates the new version is publicly available for installation through conda-forge.