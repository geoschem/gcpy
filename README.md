[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://github.com/geoschem/gcpy/blob/master/LICENSE.txt)

# GCPy: Python toolkit for GEOS-Chem

**GCPy** is a Python-based toolkit containing useful functions for working specifically with the GEOS-Chem model of atmospheric chemistry and composition.

GCPy aims to build on the well-established scientific Python technical stack, leveraging tools like cartopy and xarray to simplify the task of working with model output and performing atmospheric chemistry analyses.

## What GCPy was intended to do:

1. Generate the standard evaluation plots from GEOS-Chem benchmark output.
2. Obtain GEOS-Chem's horizontal/vertical grid information.
3. Implement GCHP-specific regridding functionalities (e.g. cubed-sphere to lat-lon regridding)
4. Provide example scripts for creating specific types of plots or analysis from GEOS-Chem output.

## What GCPY was not intended to do:

1. NetCDF file modification: (crop a domain, extract some variables):
    * Use [xarray](http://xarray.pydata.org) instead.
    * Also see [our *Working with netCDF data files* wiki page](http://wiki.geos-chem.org/Working_with_netCDF_data_files).
2. Simple plotting on lat-lon grids:
    * Can be done directly with [cartopy](https://scitools.org.uk/cartopy/docs/latest/), [matplotlib](https://matplotlib.org/), etc.
    * See our [GEOS-Chem Python tutorial](https://github.com/geoschem/GEOSChem-python-tutorial) for more examples!
3. Statistical analysis:
    * Use [scipy](http://www.scipy.org)/[scikit-learn](https://scikit-learn.org) tools instead
4. Machine Learning:
    * Use the standard machine learning utilities ([pytorch](https://pytorch.org), [tensorflow](https://www.tensorflow.org), [julia](https://julialang.org), etc.)

## Requirements:
**GCPy** is built on top of Python 3 and the scientific Python / NumPy
stack, including

- [Python >= 3.6](https://www.python.org/)
- [cartopy](https://scitools.org.uk/cartopy/docs/latest/)
- [matplotlib](https://matplotlib.org/)
- [numpy](http://www.numpy.org/)
- [scipy](http://www.scipy.org/)
- [xarray](http://xarray.pydata.org)
- [xesmf](https://xesmf.readthedocs.io)
- [esmpy](https://www.earthsystemcog.org/projects/esmpy/)
- [pypdf2](https://pythonhosted.org/PyPDF2/)

To create an environment for working with **GCPy**, we recommend using the [Anaconda Python distribution](https://www.continuum.io/downloads) or curating your own *virtualenv* or *conda* environment. Please see
```gcpy/docs/environment.yml``` for an example.


## Installation
At the moment, the easiest way to install **GCPy** is directly from our GitHub repository.

    $ git clone https://github.com/geoschem/gcpy.git gcpy

You will need to manually set your PYTHONPATH to include the top-level GCPy directory:

    $ export PYTHONPATH=/path/to/gcpy:$PYTHONPATH

Currently, **GCPy** is not available via conda-forge or PyPI, but will be starting with the release of version 1.0.0.


## License
GCPy is distributed under the MIT license.  Please read the license documents LICENSE.txt and AUTHORS.txt, which are located in the root folder.


## Contact
To contact us, please [open a new issue on the issue tracker connected to this repository](https://github.com/geoschem/gcpy/issues/new/choose). You can ask a question, report a bug, or request a new feature.