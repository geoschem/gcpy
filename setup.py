#!/usr/bin/env python

import os
import warnings

from setuptools import setup, find_packages

from textwrap import dedent

DESCRIPTION = ""
LONG_DESCRIPTION = """\
Python toolkit for working with GEOS-Chem output
"""

DISTNAME = "gcpy"
AUTHOR = "GEOS-Chem Support Team"
AUTHOR_EMAIL = "geos-chem-support@as.harvard.edu"
URL = "http://gcpy.readthedocs.io/en/latest/"
LICENSE = "MIT"

CLASSIFIERS = [
    'Development Status :: 4 - Beta',
    'License :: OSI Approved :: MIT License',
    'Operating System :: OS Independent',
    'Intended Audience :: Science/Research',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
    'Topic :: Scientific/Engineering',
]

MAJOR = 0
MINOR = 0
MICRO = 2
VERSION = "{}.{}.{}".format(MAJOR, MINOR, MICRO)
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


def _write_version_file():

    fn = os.path.join(os.path.dirname(__file__), DISTNAME, 'version.py')

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
    package_data = {},

    classifiers = CLASSIFIERS
)