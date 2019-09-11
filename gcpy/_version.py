'''
Module for obtaining version information.
'''
import os
import platform
import sys
import json
import importlib


def show_versions(as_json=False):
    '''
    Prints a list of GCPy dependencies and their version numbers.
    Can also return the list of dependencies as a JSON object.

    Keyword Args (optional):
    ------------------------
    as_json : bool
        Set this switch to True to return the list of versions
        as a JSON object instead of printing them to stdout.
        Default value: False

    Example:
    --------
    >>> import gcpy
    >>> gcpy.show_versions()
    '''
    # List of dependent packages
    deps = ['bottleneck', 'cartopy', 'cython', 'dask',
            'esmf', 'esmpy', 'graphviz', 'future', 'gcpy',
            'h5netcdf', 'h5py', 'h5pyd', 'IPython', 'jupyter',
            'matplotlib', 'netCDF4', 'notebook', 'numpy', 'pandas',
            'pip', 'pycodestyle', 'pyresample', 'pytest', 'python',
            'scipy', 'seaborn', 'setuptools', 'six', 'sphinx',
            'xbpch', 'xarray', 'xesmf']

    # Make a directory of version numbers for each module
    versions = {}
    for modname in deps:
        versions[modname] = get_version_number(modname)

    # If as_json is True, return a JSON of the version numbers
    # of each module.  Otherwise print this info to stdout.
    if as_json:
        return json.dumps(versions)
    else:
        print('\nINSTALLED VERSIONS')
        print('------------------')
        print('')
        for mod, ver in versions.items():
            print('%s: %s' % (mod.ljust(16), ver))


def get_version_number(modname):
    '''
    Returns the version number for a given module.

    Args:
    ----
    modname : str
       Name of a Python module

    Returns:
    --------
    version : str
       Version number corresponding to module.

    Example:
    --------
    >>> import json
    >>> print(get_version_number('json')
    >>> 2.0.9
    '''

    # Search in this order: Python, system modules, imported modules
    if modname.strip() == 'python':
        return platform.python_version()
    elif modname in sys.modules:
        try:
            mod = sys.modules[modname]
            return mod.__version__
        except (AttributeError, ModuleNotFoundError):
            return None
    else:
        try:
            mod = importlib.import_module(modname)
            return mod.__version__
        except (AttributeError, ModuleNotFoundError):
            return None


def main():
    '''
    Main program to call the show_versions() function.
    '''
    if len(sys.argv) == 2:
        return show_versions(sys.argv[1])
    else:
        return show_versions()


if __name__ == '__main__':
    sys.exit(main())
