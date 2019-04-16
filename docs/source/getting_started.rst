.. _getting_started:

Getting Started
---------------

Prelude
=======

The only hard-and-fast rule about scientific computing in Python is this:

    **Do not use your system Python installation!**

The version of Python that ships with operating systems such as Red Hat Linux and
macOS is usually outdated, but configured to support system functions. Although
it is entirely possible to install and use the packages mentioned in
:ref:`python-overview` using the system Python, it's much more practical on both
your local machine and any cluster you work with to curate a specialized Python
installation.

.. note::

    **Should I use Python 2 or Python 3**?

    You should use Python 3. The majority of the scientific Python packages are
    `moving to only support Python 3 <http://www.python3statement.org/>`_ in the
    near future without any backwards compatibility. The differences between
    Python 2 and Python 3 are mostly superficial, but large enough that it is
    cumbersome to mantain large codebases that are compatible with both. With the
    exception of a handful of packages you may encounter which do not support
    Python 3, there is no compelling reason to use Python 2 today.

.. _packages:

Packages
~~~~~~~~

Once you have a Python installation, setting up all the packages you want to use
is very easy thanks to the package management ecosystems available. The most
common way to install packages is to search for them on the official
`PyPI <https://pypi.python.org/pypi>`_ index. Once you've found the package you
want to install (you may have also just found it on github or elsewhere), you
simply execute from a terminal::

    $ pip install <package-name>

and it will fetch the source code, build it, and install it to wherever your
``$PYTHONPATH`` is set. This works in the vast majority of cases, particularly
when the code you're installing doesn't have any compiled dependencies. However,
because in the scientific Python world we care about performance and building
tools which interface with vetted third-party libraries, we sometimes have much
more complex dependencies. To deal with this situation, an open source package
management system called `conda <https://conda.io/>`_ was created. To install
conda, you can grab it from PyPI::

    $ pip install conda

Then, you can install packages from an official, curated set of packages which are
built and tested for a number of different system configurations on Linux,
Windows, and macOS::

    $ conda install <package-name>

Additionally, there is a `community-maintained collection of packages/recipes
<https://conda-forge.github.io/>`_ which are accessible through conda as a
channel::

    $ conda install -c conda-forge <package-name>

You can usually find bleeding-edge versions of packages on conda-forge. If you
can't find a package on either PyPI or conda-forge, you can always install it
directly from the source code. If the package is on github, ``pip`` already has
an alias to do this for you::

    $ pip install git+https://github.com/<user>/<package-name>.git

If all else fails, you can always download the source code and install it manually
like::

    $ wget https:/path/to/my/pkg/source.tar.gz
    $ tar -xvzf source.tar.gz
    $ cd source/
    $ python setup.py install

.. note::

    You can also use ``pip`` to install code you've downloaded::

        $ cd source/
        $ pip install -e .

    This will automatically call **setup.py** for you. The "**-e**" flag will
    install the package in "editable" mode, which means that any change you make
    to the source code will automatically be recognized when you load the package
    in Python; this is *especially* useful when you're developing code.

Finally, you don't *have* to go through this process of installing packages.
If you have code sitting on your disk somewhere, you can always modify the
environmental variable ``$PYTHONPATH`` to include a path to that code, and
Python will find it for you.
However, you *should not do this* if it can be avoided, because it is
extremely difficult (if not impossible) to be sure that any compiled code will
link against the correct libraries it needs, and it is very hard to debug errors
associated with mis-matched libraries/headers if you go this route.
Besides, using packages greatly improves transparency and reproducibility, so
you're already developing all your code as packages, right?

jupyter / IPython
~~~~~~~~~~~~~~~~~

Python is an interpreted language, which means you'll spend most of your time
inside an interactive shell/environment typing in commands or running scripts.
However, the default Python interpeter is quite barebones. For serious usage, you
should use a tool like `IPython <https://ipython.org/>`_, which extends the
default interpeter with all sorts of useful features, including:

1. As-you-go syntax highlighting
2. Documentation access
3. Multi-line command entry
4. Code-completition and history
5. Interoperability with system shell

The IPython project is part of a larger scientific computing project called
`Jupyter <https://jupyter.org/>`_, which provides even more sophisticated tools
and environments for working with your code. An extremely popular environment is
the **Notebook**, which provides a browser-based interface for working with
your code in a document which mixes code, markup-language/documentation, and
much more. Jupyter Notebooks can be converted to stand-alone documents, reports,
slideshows, or other multi-media.

.. note::

    As an example, the work supporting the LIGO discovery of gravitational
    waves is `fully documented using Jupyter Notebooks
    <https://losc.ligo.org/s/events/GW150914/GW150914_tutorial.html>`_. This is a
    major milestone in terms of reproducibility and open science.

    Let's be frank for a moment: these researchers will *undoubtedly* win a Nobel
    Prize in Physics for this work sometime in the next decade. If these tools are
    good enough for work leading to a *Nobel Prize*, then they're good enough for
    you to consider trying out, right?

To install IPython, Jupyter, and the Notebook environment, simply install their
packages, e.g. through conda::

    $ conda install jupyter notebook ipython

Once installed, you can open an IPython prompt by executing from your command
line::

    $ ipython

which will open up a prompt that looks something like this

.. parsed-literal::

    Python 3.5.2 \|Continuum Analytics, Inc.\| (default, Jul  2 2016, 17:52:12)
    Type "copyright", "credits" or "license" for more information.

    IPython 5.1.0 -- An enhanced Interactive Python.
    ?         -> Introduction and overview of IPython's features.
    %quickref -> Quick reference.
    help      -> Python's own help system.
    object?   -> Details about 'object', use 'object??' for extra details.

    In [1]:

You can then enter Python commands as if you were in a normal Python interpreter.


One Step to Scientific Python
=============================

The easiest way to set up a full-stack scientific Python deployment is to use a
**Python distribution**. This is an installation of Python with a set of curated
libraries. Two examples of such a distribution are the Anaconda_ distribution from
Continuum IO and the `Enthought Python Distribution`_ from Enthought. Both of
these distributions include one-click installers, and provide some graphical
utilities to help manage any packages you may want to install which are not
already included in the curated inclusion list.

Three Steps to Scientific Python
================================

Alternatively, the way I'd recommend to start up a scientific Python environment
is to follow these steps:

1. **Obtain a minimal Python installer**
    I like to use the `Miniconda <https://conda.io/miniconda.html>`_ installer;
    this provides a Python install for your operating system, plus the **conda**
    package manager. This way, you can install just the packages you want and
    need.

2. **Run the installer**
    You'll probably need to do this from the command line, e.g.::

        $ sh Miniconda3-latest-MacOSX-x86_64.sh

    Follow the instructions; you can choose where to place the installation
    (preferably somewhere you have write access without super-user/root access,
    like your home directory). At the end of this process, add this path to your
    \*rc configuration::

        $ echo "PATH=$PATH:/path/to/miniconda/bin" > ~/.bashrc

    If you do this, your ``$PYTHONPATH`` will be implicitly configured correctly
    and you will never have to touch it.

3. **Install any packages you want**
    As shown in :ref:`packages`, install whatever packages you want
    using **conda**.

That's all there is to it! In general, this is a better way to go because you can
quickly curate your own scientific Python installation on any external computer
resources you may wish to use (e.g. university cluster).

Environments
============

Python coupled with a package manager provides a way to make isolated,
reproducible *environments* where you have fine-tuned control over all packages
and configuration. One environment solution that works well with PyPI is
`virtualenv <https://virtualenv.pypa.io/en/stable/>`_; you can find many resources
on using virtualenv on the internet as it's widely used in web application
deployments.

For scientific Python, you can alternatively use **conda**\'s built in
environment management system. To create a conda environment, you simply
execute the following command::

    $ conda create --name my_environment python=3.6 numpy

This will create a special environment in ``$MINICONDA_HOME/envs/my_environment``
with only Python and numpy to begin with. Here, we've also told conda to install
Python version 3.6; you can specify exact versions or minima, and conda will
take care of figuring out all the compatibilties between versions for you. To use
this environment, simply "activate" it by executing::

    (my_environment) $ source activate my_environment

Regardless of your shell, you should now see the string ``(my_environment)``
prepended to your prompt. Now, if you execute any Python-related tool from the
command line, it will first search in ``$MINICONDA_HOME/envs/my_environment/bin``
to find them. You can deactivate your environment by typing::

    $ source deactivate

For extensive documentation on using environments, please see
`the conda documentation <https://conda.io/docs/using/envs.html#>`_. The most
important feature to review here is the ability to *share and export* your
environment; this is the basis for reproducibility in the scientific Python stack.
At any time from the shell, you can execute::

    $ conda list

to get a complete summary of all the packages installed in your environment, the
channel they were installed from, and their full version info. Using this info,
you can create an **environment file** in YAML syntax which documents the exact
contents of your environment. With that file, a new environment with the exact
configuration can be installed by executing::

    $ conda env create -f my_environment.yml


Geosciences Python Environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Combining all of the previous sections, we can very easily spin-up a
full-featured scientific Python environment with a set of packages curated for the
geosciences. Copy the ``environment.yml`` file located in gcpy/docs somewhere
on your local hard drive.

.. note::

    Installing this environment will also install many dependencies, including
    compiled libraries. This is totally fine; even if you have these libraries
    already installed through your system package manager, **conda** will install
    and link for use in the environment a configuration which should be guaranteed
    to play nicely and work with all of its components.

Create this environment through **conda**::

    $ conda env create -f /path/to/environment.yml

Activate this environment::

    $ source activate geo_scipy

You're now ready to reproduce any example analysis in this documentation.
