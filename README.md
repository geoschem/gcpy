# GCPy: Python toolkit for GEOS-Chem

<p>
   <a href="https://github.com/geoschem/gcpy/releases"><img src="https://img.shields.io/github/v/release/geoschem/gcpy?label=Latest%20Stable%20Release"></a>
   <a href="https://anaconda.org/conda-forge/geoschem-gcpy"> <img src="https://anaconda.org/conda-forge/geoschem-gcpy/badges/version.svg" /> </a>
   <a href="https://img.shields.io/pypi/v/geoschem-gcpy"><img alt="PyPI - Version" src="https://img.shields.io/pypi/v/geoschem-gcpy"></a>
   <a href="https://github.com/geoschem/gcpy/releases/"><img src="https://img.shields.io/github/release-date/geoschem/gcpy"></a>
   <br />
   <a href="https://anaconda.org/conda-forge/geoschem-gcpy"> <img src="https://anaconda.org/conda-forge/geoschem-gcpy/badges/platforms.svg" /> </a>
   <a href="https://doi.org/10.5281/zenodo.3689589"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.3689589.svg" alt="DOI"></a>
   <a href="https://github.com/geoschem/gcpy/blob/main/LICENSE.txt"><img src="https://img.shields.io/badge/License-MIT-blue.svg"></a>
   <br />
   <a href="https://gcpy.readthedocs.io/en/latest/"><img src="https://img.shields.io/readthedocs/gcpy?label=ReadTheDocs"></a>
   <a href="https://github.com/geoschem/gcpy/actions/workflows/build-gcpy-environment.yml"><img src="https://github.com/geoschem/gcpy/actions/workflows/build-gcpy-environment.yml/badge.svg"></a>
   <a href="https://anaconda.org/conda-forge/geoschem-gcpy"> <img src="https://anaconda.org/conda-forge/geoschem-gcpy/badges/downloads.svg" /> </a>  
</p>

**GCPy** is a Python-based toolkit containing useful functions for working specifically with the GEOS-Chem model of atmospheric chemistry and composition.

**GCPy** aims to build on the well-established scientific Python technical stack, leveraging tools like **cartopy**, **numpy**, and **xarray** to simplify the task of working with GEOS-Chem model output and performing atmospheric chemistry analyses.


## What GCPy was intended to do:

1. Produce plots and tables from [GEOS-Chem](https://geos-chem.readthedocs.io) output using simple function calls.
2. Generate the standard evaluation plots and tables for GEOS-Chem benchmark simulations.
3. Obtain GEOS-Chem's horizontal and vertical grid information.
4. Implement [GCHP](https://gchp.readthedocs.io)-specific regridding functionalities (e.g. cubed-sphere to lat-lon regridding)
5. Provide example scripts for creating specific types of plots or  analysis from GEOS-Chem output.
6. Provide user-submitted scripts for specific applications related to GEOS-Chem and [HEMCO](https://hemco.readthedocs.io).

## What GCPy was not intended to do:

1. General NetCDF file modification: (crop a domain, extract some variables):
    * Instead, use netCDF tools such as:
	  * [xarray](http://xarray.pydata.org)
	  * [netCDF operators (NCO)](https://nco.sourceforge.net)
	  * [Climate Data Operators](https://mpimet.mpg.de/cdo) instead.
    * Also see our [*Work with netCDF files* guide](https://geos-chem.readthedocs.io/en/latest/geos-chem-shared-docs/supplemental-guides/netcdf-guide.html) at [geos-chem.readthedocs.io](https://geos-chem.readthedocs.io)

2. Statistical analysis:
    * Instead, use statistical tools such as:
	  * Use [scipy](http://www.scipy.org)
	  * [scikit-learn](https://scikit-learn.org)
	  * [R](https://r-project.org)
	  * etc

3. Machine Learning:
    * Instead, use machine learning tools such as:
	  * [pytorch](https://pytorch.org),
	  * [tensorflow](https://www.tensorflow.org)
	  * [julia](https://julialang.org)
	  * etc.

## Documentation:

For more information on installing and using GCPy, visit the official documentation at [gcpy.readthedocs.io](https://gcpy.readthedocs.io/).

## License

GCPy is distributed under the MIT license.  Please see the [GCPy license agreement](https://github.com/geoschem/gcpy/blob/dev/LICENSE.txt) and [List of GCPy developers](https://github.com/geoschem/gcpy/blob/dev/AUTHORS.txt) for more information.

## Requesting support

To report a bug or suggest a new feature, please see our [Support
Guidelines](https://github.com/geoschem/gcpy/blob/dev/SUPPORT.md).

## Submitting new features

If you are interested in submitting code to GCPy, please see our
[Contributing Guidelines](https://github.com/geoschem/gcpy/blob/dev/CONTRIBUTING.md).
