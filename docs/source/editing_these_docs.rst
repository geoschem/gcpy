.. |br| raw:: html

   <br />

.. _editing_this_user_guide:

#######################
Editing this User Guide
#######################

This user guide is generated with `Sphinx
<https://www.sphinx-doc.org/>`_.  Sphinx is an open-source Python
project designed to make writing software documentation easier.  The
documentation is written in :ref:`reStructuredText (reST)
<editing_this_user_guide_rest>`, a plaintext markup language that
Sphinx extends for software documentation. The source for the
documentation is the :file:`docs/source` directory in top-level of the
source code (and its subdirectories).

.. _editing_this_user_guide_quickstart:

===========
Quick start
===========

First-time setup: Install Sphinx
--------------------------------

To build this user guide on your local machine, you need to install
Sphinx and its dependencies, which are listed in the table below.

.. list-table::
   :header-rows: 1
   :widths: 30 50 20

   * - Package
     - Description
     - Version
   * - sphinx
     - Creates online user manual documentation from markup text files
     - 7.2.6
   * - `sphinx-autobuild <https://github.com/sphinx-doc/sphinx-autobuild>`_
     - Dynamically builds Sphinx documentation and displays it in a
       browser
     - 2021.3.14
   * - `sphinx_rtd_theme <https://github.com/readthedocs/sphinx_rtd_theme>`_
     - Sphinx theme for ReadTheDocs
     - 2.0.0
   * - `sphinxcontrib-bibtex <https://pypi.org/project/sphinxcontrib-bibtex/>`_
     - Inserts LaTeX-style bibliography citations into ReadTheDocs
       documentation
     - 2.6.2
   * - `docutils <https://docutils.sourceforge.io/>`_
     - Processes plaintext documentation into HTML and other formats
     - 0.20.1
   * - `recommonmark  <https://github.com/readthedocs/recommonmark>`_
     - Parses text for docutils
     - 0.7.1
   * - `jinja2 <https://jinja.palletsprojects.com/en/stable/>`_
     - Replaces tokenized strings with text
     - 3.1.6

We recommend that you create a standalone :program:`Conda` environment
to install Sphinx and its dependencies.  The YAML file
:file:`docs/environment_files/read_the_docs_environment.yaml` contains the proper package
specifications.  Use these commands:

.. code-block:: console

   $ cd docs
   $ conda env create -n rtd_env --file=environment_files/read_the_docs_environment.yml

This step only needs to be done once.

Build the documentation
-----------------------

#. Activate the :program:`Conda` environment containing
   :program:`Sphinx` and its dependencies:

   .. code-block:: console

      $ conda activate rtd_env

#. Navigate to the :file:`docs/` folder:

   .. code-block:: console

      (rtd_env) $ cd docs     # Skip if you are already in the docs folder

#. Check out the :file:`docs/dev` branch of this repository, as this
   is the branch from which the :program:`latest` ReadTheDocs version
   will be built:

   .. code-block:: console

      (rtd_env) $ git checkout docs/dev   # Skip if you are already on the docs/dev branch

#. Start the :command:`sphinx-autobuild` server:

   .. code-block:: console

      (rtd_env) $ sphinx-autobuild source build/html

#. Remove any HTML files (in :file:`docs/build/html`) that might be
   left behind from a previous build:

   .. code-block:: console

      (rtd_env) $ make clean

   This will parse the reST-format files in the :file:`docs/source/`
   directory tree and generate new HTML files in
   :file:`docs/build/html`. |br|
   |br|

#. Open a web browser and navigate to :file:`localhost:8000`. |br|
   |br|

#. Open your favorite text editor and start making changes to the
   reST-format documentation files in the :file:`docs/source`
   directory tree.  While :program:`sphinx-autobuild` is running, you
   will see your updates rendered in the web browser as soon as you
   soon as you save your changes to disk. |br|
   |br|

#. Once you are satisfied with your edits, commit your changes to Git
   and push the documentation to the :file:`docs/dev` remote branch of
   this repository, |br|
   |br|

#. Remove the generated HTML documentation files:

   .. code-block:: console

      (rtd_env) $ make clean

#. Halt the :program:`sphinx-autobuild` server by typing
   :program:`CTRL-C`. |br|
   |br|

#. Deactivate the :program:`Conda` environment:

   .. code-block:: console

      (rtd_env) $ conda deactivate

.. _editing_this_user_guide_rest:

==========================
Learning reStructured Text
==========================

ReadTheDocs documentation is generated from text files in **reStructured
Text (reST)**, which is an easy-to-read, what-you-see-is-what-you-get
plaintext markup language. It is the default markup language used by
Sphinx.

Writing reST can be tricky at first. Whitespace matters, and some
directives can be easily miswritten. Two important things you should
know right away are:

- Indents are 3-spaces
- "Things" are separated by 1 blank line. For example, a list or
  code-block following a paragraph should be separated from the
  paragraph by 1 blank line.

You should keep these in mind when you're first getting
started. Dedicating an hour to learning reST will save you time in the
long-run. Below are some good resources for learning reST.

- `reStructuredText primer
  <https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`_:
  (single best resource; however, it's better read than skimmed) |br|
  |br|

- Official `reStructuredText reference
  <https://docutils.sourceforge.io/docs/user/rst/quickref.html>`_
  (there is *a lot* of information here) |br|
  |br|

- `Presentation by Eric Holscher
  <https://www.youtube.com/watch?v=eWNiwMwMcr4>`_ (co-founder of Read
  The Docs) at DjangoCon US 2015 (the entire presentation is good, but
  reST is described from 9:03 to 21:04) |br|
  |br|

- `YouTube tutorial by Audrey Tavares
  <https://www.youtube.com/watch?v=DSIuLnoKLd8>`_

A good starting point would be Eric Holscher's presentations followed
by the reStructuredText primer.

.. _editing_this_user_guide_style:

================
Style guidelines
================

This user guide is written in semantic markup. This is important so
that the user guide remains maintainable. Before contributing to
this documentation, please review our style guidelines
(below). When editing the source, please refrain from using
elements with the wrong semantic meaning for aesthetic
reasons. Aesthetic issues can be addressed by changes to the theme.

Titles and headers
------------------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Element
     - reST Markup
   * - Section header |br| (aka "Heading 1)
     - Overline by :literal:`#` and underline by :literal:`#`
   * - Sub-section header |br| (aka "Heading 2")
     - Overline by :literal:`=` and underline by :literal:`=`
   * - Sub-sub-section header |br| (aka "Heading 3")
     - Underline by :literal:`-`
   * - Sub-sub-sub-section header |br| (aka "Heading 4")
     - Underline by :literal:`~`
   * - Sub-sub-sub-sub-section header |br| (aka "Heading 5")
     - Underline by :literal:`^`

References and links
--------------------

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - Element
     - reST Markup Example
     - Rendered text
   * - Reference to a named anchor
     - ``:ref:`editing_this_user_guide_quickstart```
     - :ref:`editing_this_user_guide_quickstart`
   * - Renamed reference to a named anchor
     - ``:ref:`Getting Started <editing_this_user_guide_quickstart>``
     - :ref:`Getting Started <editing_this_user_guide_quickstart>`
   * - HTML link
     - ```ReadTheDocs <https://geos-chem.readthedocs.io>`_``
     - `GEOS-Chem Manual <https://geos-chem.readthedocs.io>`_

Other common style elements
---------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - Element
     - reST Markup Example
     - Rendered text
   * - File paths
     - ``:file:`myfile.nc```
     - :file:`myfile.nc`
   * - Directories
     - ``:file:`/usr/bin```
     - :file:`/usr/bin`
   * - Program names
     - ``:program:`cmake```
     - :program:`cmake`
   * - OS-level commands
     - ``:program:`rm -rf```
     - :program:`rm -rf`
   * - Environment variables
     - ``:envvar:`$HOME```
     - :envvar:`$HOME`
   * - Inline code or code variables
     - ``:code:`PRINT*, "HELLO!"```
     - :code:`PRINT*, "HELLO!"`
   * - Inline literal text
     - ``:literal:`$```
     - :literal:`$`

Indented code and text blocks
-----------------------------

Code snippets should use :literal:`.. code-block:: <language>`
directives:

Python
~~~~~~

.. code-block:: none

   .. code-block:: python

      import gcpy
      print("hello world")

Renders as:

.. code-block:: python

    import gcpy
    print("hello world")

Fortran
~~~~~~~

.. code-block:: none

   .. code-block:: Fortran

      DO I = 1, 10
         PRINT*, I
      ENDDO

Renders as:

.. code-block:: Fortran

   DO I = 1, 10
      PRINT*, I
   ENDDO

Bash
~~~~

.. code-block:: none

   .. code-block:: bash

      #!/bin/bash

      for f in *.nc; do
          echo $f
      done

Renders as:

.. code-block:: bash

   #!/bin/bash

   for f in *.nc; do
       echo $f
   done

Command line (aka "console")
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: none

   .. code-block:: console

      $ ls -l $HOME

Renders as:

.. code-block:: console

   $ ls -l $HOME

No formatting
~~~~~~~~~~~~~

.. code-block:: none

   .. code-block:: none

      This text renders without any special formatting.

Renders as:

.. code-block:: none

   This text renders without any special formatting.
