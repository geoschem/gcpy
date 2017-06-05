.. _python-overview:

Scientific Python Overview
--------------------------

Python_ is a widely-used programming language that has merged as a popular
environment to perform data analysis and visualization. This popularity is based
around a large user community that develops and maintains a large number of free,
open-source packages for many different applications in the worlds of science,
finance, business, system programming, and more.

Although Python_ is provided "with batteries included"
(it has a `comprehensive standard library <https://docs.python.org/3/>`_ with
built-in tools for many tasks including file I/O and processing, datetime
handling, command line tool integration, and much more), its use for scientific
visualization and analysis is built on a suite of community-built tools:

.. figure:: _static/state_of_the_stack_2015.png
    :alt: State of the Stack, Jake Vanderplas, 2015
    :align: center
    :width: 75%

    State of the "Scientific Python Stack" circa 2015 (courtesy `Jake
    VanderPlas`_)

*n*-dimensional Arrays
======================

Out of the box, Python does not have a high-performance multi-dimensional array
object. This sort of data structure is **critical** for numerical analysis. In
the Python ecosystem, multi-dimensional arrays are provided by the NumPy_
library, which contains a high-performance array object and related linear algebra
routines. The array "interface" provided by NumPy supports vector programming and
is very general, underpinning many of the numerical tools widely used in the Python world. For instance, the companion SciPy_ package implements many useful
functions foor optimization, statistics, interpolation, image processing, and
spatial mathematics. Many of these functions work with NumPy_ arrays having
arbitrary size, dimensions, and data types.

.. note::

    The version of Python you will most likely use has been implemented in C, and
    designed in such a way that it is very easy to work with other compiled
    codebases, especially those written in C/C++ or Fortran. In fact, libraries
    like NumPy_ heavily rely on optimized, compiled codes in these languages,
    which greatly improves their performance, reliability, and speed.

Another library which builds on the NumPy_ array interface is pandas_, which
extends the array with labeling and a powerful engine for the transformation and
analysis of structured datasets. The data structures provided by pandas -
particularly the ``DataFrame`` - greatly simplify timeseries analysis,
split-apply-combine workflows, and other common research processing tasks. The
xarray_ package extends pandas by providing addtional data structures to handle
*n*-dimensional labeled arrays (such as those contained in a NetCDF file). Most
importantly, the array and labeling semantics in NumPy, pandas, and xarray are
similar if not identical in the majority of cases, and because each package
sequentially builds on the others, they can all be used within the same
analysis context and data can easily be shuttled back-and-forth to whatever format
works best for a given task.

Visualization
=============

The core visualization library in the Python world is matplotlib_. Although it
originated as an emulation of the graphics capabilities of MATLAB, matplotlib has
grown into the defacto base layer for 2D graphics in Python. Matplotlib provides
fine-grained control of graphics, and works natively with both base Python objects
and NumPy_ array derived-types (including pandas ``Series`` and ``DataFrame``\s and xarray ``DataArray``\s). In fact, both pandas and xarray provide shim layers to help automate plotting numerical data through matplotlib.

Many libraries extend the core features of matplotlib. For instance, seaborn_
implements many useful statistical visualizations and leans particularly heavily
on pandas to help organize data for plotting. To plot geographical
information, one can use the cartopy_ library, which itself wraps several
open-source cartographic libraries and has support for geographical projections,
shapefiles, and more. Users coming from R who love ggplot2_ should be aware of
several upcoming grammar of graphics implementations based on matplotlib,
including `plotnine <https://plotnine.readthedocs.io/en/stable/>`_,
`ggplot <http://ggplot.yhathq.com/>`_, and `altair
<https://altair-viz.github.io/>`_

Domain-Specific Toolkits
========================

Python users can be found throughout the ranks of researchers in the natural
sciences, and many contribute specialized toolkits to help with their own niche
applications. Here is a short summary of tools that can be useful for different
research tasks in Python:

.. todo:: fix link bolding below

**`scikit-learn <http://scikit-learn.org/stable/documentation.html>`_**
    A toolkit implement a wide variety of algorithms for un/supervised machine
    learning tasks, including regressions, clustering, manifold learning,
    principal components, density estimation, and much more. It also provides
    many useful tools to help build "pipelines" for managing modeling tasks such
    as data processing/normalization, feature engineering, cross-validation,
    fitting, and prediction.

**`statsmodels <http://www.statsmodels.org/>`_**
    A module for fitting and estimating many different types of statistical models
    as well as performing hypothesis testing and exploratory data analysis. It
    features tools for fitting generalized linear models, survival analyses,
    and multi-variate statistics. Furthermore, it implements formula-based
    regression specification similar to R which natively works with pandas_
    data structures.

**`scikit-image <http://scikit-image.org/>`_**
    An image processing library featuring many common operations including
    convolutional mapping, filtering, edge detection, and image segmentation.

**`pyresample <https://pyresample.readthedocs.io/en/latest/>`_**
    A Python package for reprojecting earth observing satellite data, capable
    of handling both swath data from polar-orbitting satellites and gridded data
    from geostationary satellites.

**`shapely <http://toblerity.org/shapely/manual.html>`_**
    A spatial analysis library which extends Python to work as a fully-featured
    GIS environmental comparable to commercial software such as ArcGIS.

**`SymPy <http://www.sympy.org/en/index.html>`_**
    A full-featured computer algebra system (CAS) similar to Mathematica or Maple.
    SymPy powers an additional ecosystem of domain-specific tools used in
    pure mathematics research and which have many applications in physics.

**`pymc3 <https://pymc-devs.github.io/pymc3/>`_**
    A toolkit for Bayesian statistical modeling and probabilistic programming,
    including a suite of Markov chain Monte Carlo fitting algorithms.

Compilation/Optimization
========================

Python is slower than statically-typed, compiled languages such as C/C++ and
Fortran. However, it doesn't *have* to be slow. Vectorized programming through
NumPy and pandas can dramatically increase the performance of Python in executing
numerical analyses and calculations. However, in applications where vectorization
is non-trivial or inappropriate, Python's performance can be dramatically improved
by using one of several different approaches.

First, an optimising compiler called cython_ is available to compile your code
into high-performance, efficient C kernels. Cython will work on your normal Python
code with few modifications, and can often times increase its performance by 1-2x.
However, by incorporating a special set of static typing directives into your code
(similar to what you'd do in Fortran by declaring variable types), Cython can go
a step further and yield much more significant performance improvements. It also
trivializes the task of wrapping legacy code from C/C++ or Fortran applications.

Alternatively, one can use a `just-in-time (JIT)
<https://en.wikipedia.org/wiki/Just-in-time_compilation>`_ compiler to compile
code on-the-fly. One approach in the Python world implementing a JIT is the
`PyPy <https://pypy.org/>`_ project, which is an alternative implementation of
Python itself. A drawback to PyPy is that it does not totally support all of the
numerical libraries in the scientific Python stack. Instead, one can use
`Numba <http://numba.pydata.org/>`_ to target specific, expensive functions or
subroutines in a codebase. Numba-compiled functions can target multi-core
or GPU architectures when available.

Niche optimization tools also exist in the Python world. For instance, the
`PyCUDA <https://mathema.tician.de/software/pycuda/>`_ package helps to glue
together Python with high-performance GPGPU routines written in C/CUDA.
Meta-programming libraries for working with tensor mathematics are also widely
used, including `theano <http://deeplearning.net/software/theano/>`_ and
`TensorFlow <https://www.tensorflow.org/>`_.


Parallel Computing
==================

The most recent versions of Python include `modules and infrastructure for
asynchronous and coroutine programming
<https://docs.python.org/3/library/asyncio.html?highlight=asyncio#module-asyncio>`_
in the standard library. The programming model used in this paradigm (using
"futures" or other "delayed" objects representing a contract for future results
from calculations) is extended by several libraries in the Python ecosystem.

In particular, the `joblib <http://pythonhosted.org/joblib/>`_ library implements
a very lightweight, easy-to-use set of tools in this programming model. Joblib
strives to let you make simple modifications to your single-threaded code to
achieve parallelism when and only where it becomes necessary. An alternative
library is `dask <http://dask.pydata.org/en/latest/>`_, which provides a similar
API but works natively with array-like and DataFrame-like structures from NumPy
and pandas. Dask abstracts the parallel programming model one step further,
tracking a graph representing your computations; it optimizes this graph before
any calculations are actually performed, which allows it to optimize the amount
of data held in memory at any given time, or scale to arbitrary resources as
they become available.

Traditional `MPI <https://mpi4py.readthedocs.io/en/stable/>`_ tools also exist
in the Python ecosystem, although these tend to be very low-level and
"un-Pythonic."


.. _Anaconda: https://www.continuum.io/anaconda-overview
.. _cartopy: http://scitools.org.uk/cartopy/docs/latest/
.. _cython: http://cython.org/
.. _Enthought Python Distribution: https://www.enthought.com/products/epd/
.. _ggplot2: http://ggplot2.org/
.. _Jake Vanderplas: https://staff.washington.edu/jakevdp/
.. _matplotlib: http://matplotlib.org/
.. _NumPy: http://www.numpy.org/
.. _pandas: http://pandas.pydata.org/
.. _Python: https://www.python.org/
.. _SciPy: https://www.scipy.org/
.. _seaborn: http://seaborn.pydata.org/index.html
.. _xarray: https://xarray.pydata.org/
