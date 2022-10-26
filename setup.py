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
    'Programming Language :: Python :: 3.9',
    'Topic :: Scientific/Engineering',
]

MAJOR = 1
MINOR = 3
MICRO = 2
EXTRA = '' # for alpha (aN), beta (bN), rc (rcN) versions

VERSION = "{}.{}.{}{}".format(MAJOR, MINOR, MICRO, EXTRA)
'''
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
'''

def _write_version_file():

    fn = os.path.join(os.path.dirname(__file__), 'gcpy', '_version.py')

    version_str = dedent("""
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
        "cartopy==0.19.0.post1",
        "dask==2021.7.1",          
        "docutils==0.16",      
        "esmf==8.1.1",
        "esmpy==8.1.1",
        "h5netcdf==0.11.0",
        "h5py==3.3.0",
        "jinja2==3.0.3",
        "joblib==1.0.1",
        "matplotlib==3.4.2",
        "pandas==1.3.1",
        "netcdf4==1.6.0",
        "netcdf-fortran==4.5.4",
        "numpy==1.21.1",
        "pypdf2==1.26.0",
        "recommonmark==0.7.1",
        "requests==2.26.0",
        "scipy==1.7.0",
        "sparselt>=0.1.3",
        "sphinx==3.5.4",
        "sphinx-autoapi==1.9.0",
        "sphinx-autobuild==2021.3.14",
        "sphinxcontrib-bibtex==2.2.0",
        "sphinx_rtd_theme==0.5.2",
        "tabulate==0.8.9",
        "xarray==0.17.0",
        "xesmf==0.5.1"
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
