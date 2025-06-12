#!/usr/bin/env python

import os
import warnings

from setuptools import setup, find_packages

from textwrap import dedent

DESCRIPTION = ""
LONG_DESCRIPTION = """\
Python toolkit for working with GEOS-Chem output.
"""

DISTNAME = "geoschem-gcpy"
AUTHOR = "GEOS-Chem Support Team"
AUTHOR_EMAIL = "geos-chem-support@g.harvard.edu"
URL = "https://github.com/geoschem/gcpy"
LICENSE = "MIT"

CLASSIFIERS = [
    'Development Status :: 5 - Production/Stable',
    'License :: OSI Approved :: MIT License',
    'Operating System :: OS Independent',
    'Intended Audience :: Science/Research',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3.12',
    'Topic :: Scientific/Engineering',
]

MAJOR = 1
MINOR = 6
MICRO = 2
EXTRA = '' # for alpha (aN), beta (bN), rc (rcN) versions

VERSION = f"{MAJOR}.{MINOR}.{MICRO}{EXTRA}"
"""
#DEV format (using git hash) is intriguing but incompatible with PEP 440
#No hashes can be used in version field
DEV = True


# Correct versioning with git info if DEV

if DEV:
    import subprocess

    pipe = subprocess.Popen(
        ['git', "describe", "--always", "--match", "v[0-9]*"],
        stdout=subprocess.PIPE)
    so, err = pipe.communicate()

    if pipe.returncode != 0:
        # no git or something wrong with git (not in dir?)
        warnings.warn("WARNING: Couldn't identify git revision, using generic version string")
        VERSION += ".dev"
    else:
        git_rev = so.strip()
        git_rev = git_rev.decode('ascii') # necessary for Python >= 3

        VERSION += ".dev-{}".format(git_rev)
"""

def _write_version_file():

    fn = os.path.join(os.path.dirname(__file__), 'gcpy', '_version.py')

    version_str = dedent(
        """
        __version__ = '{}'
        """)

    # Write version file
    with open(fn, 'w') as version_file:
        version_file.write(version_str.format(VERSION))

# Write version and install
_write_version_file()

setup(
    name = DISTNAME,
    author = AUTHOR,
    author_email = AUTHOR_EMAIL,
    maintainer = AUTHOR,
    maintainer_email = AUTHOR_EMAIL,
    description = DESCRIPTION,
    long_description = LONG_DESCRIPTION,
    license = LICENSE,
    url = URL,
    version = VERSION,
    packages = find_packages(),
    include_package_data=True,
    install_requires=[
        "cartopy==0.23.0",
        "cf_xarray==0.9.1",
        "dask==2025.3.0",
        "esmf==8.6.1",
        "esmpy==8.6.1",
        "gridspec==0.1.0",
        "ipython==8.25.0",
        "joblib==1.4.2",
        "jupyter==1.0.0",
        "matplotlib==3.8.4",
        "netcdf4==1.6.5",
        "netcdf-fortran==4.6.1",
        "numpy==1.26.4",
        "pandas==2.2.2",
        "pip==24.0",
        "pylint==3.2.2",
        "pyproj==3.6.1",
        "python==3.12.0",
        "pypdf==4.2.0",
        "requests==2.32.4",
        "scipy==1.13.1",
        "sparselt==0.1.3",
        "tabulate==0.9.0",
        "tk==8.6.13",
        "xarray==2024.5.0",
        "xesmf==0.8.5",
    ],
    classifiers = CLASSIFIERS
)


'''
Instructions for publishing new version to TestPyPi (mimics PyPi) and PyPi:

1. Update version numbers above and in _version.py and commit before releasing
2. Install twine if not available (pip install twine)
3. run python setup.py sdist bdist_wheel
4. run twine check dist/*
5. run twine upload --repository-url https://test.pypi.org/legacy/ dist/*
6. Test installation in a virtual environment using pip install --extra-index-url https://test.pypi.org/simple/ geoschem-gcpy
7. To publish to PyPi, run twine upload dist/*

'''
