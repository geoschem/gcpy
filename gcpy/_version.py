'''
Module for obtaining version information.
'''
import platform
import sys
import yaml
import importlib


def show_versions(as_yaml=False):
    '''
    Prints a list of GCPy dependencies and their version numbers.
    Can also return the list of dependencies as a YAML object.

    Keyword Args (optional):
    ------------------------
    as_yaml : bool
        Set this switch to True to return the list of versions
        as a YAML object instead of printing them to stdout.
        Default value: False
    '''
    
    # List of dependent packages
    deps = ['platform', 'python', 'bottleneck', 'cartopy', 'cython',
            'dask', 'esmf', 'esmpy', 'graphviz', 'future', 'gcpy',
            'h5netcdf', 'h5py', 'h5pyd', 'IPython', 'jupyter',
            'matplotlib', 'netCDF4', 'notebook', 'numpy', 'pandas',
            'pip', 'pycodestyle', 'pyresample', 'pytest',
            'scipy', 'seaborn', 'setuptools', 'six', 'sphinx',
            'xbpch', 'xarray', 'xesmf']
    
    # Make a directory of version numbers for each module
    # Skip packages for which we cannot get 
    versions = {}
    for modname in deps:
        versions[modname] = get_version_number(modname)
        if versions[modname] is None:
            del versions[modname]

    # If as_yaml is True, return a YAML of the version numbers
    # of each module.  Otherwise print this info to stdout.
    if as_yaml:
        return yaml.dump(versions)
    else:
        print('\nSYSTEM INFORMATION')
        print('------------------')
        for mod in ['platform', 'python']:
            print('{}: {}'.format(mod.ljust(16), versions[mod]))

        del versions['python']
        del versions['platform']
            
        print('\nVERSION NUMBERS FOR GCPy DEPENDENCIES')
        print('-------------------------------------')
        print('')
        for mod, ver in versions.items():
            print('{}: {}'.format(mod.ljust(16), ver))


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
    '''

    # Search in this order: Python, system modules, imported modules
    if modname.strip() == 'python':
        return platform.python_version()
    elif modname.strip() == 'platform':
        return platform.platform()
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
