[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://github.com/geoschem/gcpy/blob/master/LICENSE.txt)

# GCPy: Python toolkit for GEOS-Chem

**GCPy** is a Python-based toolkit containing useful functions for working specifically with the GEOS-Chem model of atmospheric chemistry and composition.

GCPy aims to build on the well-established scientific Python technical stack, leveraging tools like cartopy and xarray to simplify the task of working with model output and performing atmospheric chemistry analyses.



## What GCPy was intended to do:

1. Produce plots and tables from GEOS-Chem output using simple function calls.
2. Generate the standard evaluation plots and tables from GEOS-Chem benchmark output.
3. Obtain GEOS-Chem's horizontal/vertical grid information.
4. Implement GCHP-specific regridding functionalities (e.g. cubed-sphere to lat-lon regridding)
5. Provide example scripts for creating specific types of plots or analysis from GEOS-Chem output.

## What GCPY was not intended to do:

1. General NetCDF file modification: (crop a domain, extract some variables):
    * Use [xarray](http://xarray.pydata.org) instead.
    * Also see [our *Working with netCDF data files* wiki page](http://wiki.geos-chem.org/Working_with_netCDF_data_files).
2. Statistical analysis:
    * Use [scipy](http://www.scipy.org)/[scikit-learn](https://scikit-learn.org) tools instead
3. Machine Learning:
    * Use the standard machine learning utilities ([pytorch](https://pytorch.org), [tensorflow](https://www.tensorflow.org), [julia](https://julialang.org), etc.)


## Documentation:

For more information on installing and using GCPy, visit the official documentation at [gcpy.readthedocs.io](https://gcpy.readthedocs.io/).


## License

GCPy is distributed under the MIT license.  Please read the license documents LICENSE.txt and AUTHORS.txt, which are located in the root folder.


## Contact

To contact us, please [open a new issue on the issue tracker connected to this repository](https://github.com/geoschem/gcpy/issues/new/choose). You can ask a question, report a bug, or request a new feature.
