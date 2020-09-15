
Editing these docs
==================

This documentation is generated with Sphinx. This page describes how to contribute to the GCPy documentation.

Quick start
-----------

You need the Sphinx Python to build (and therefore edit) this documentation. Assuming you already have Python installed,
install Sphinx:

.. code-block:: console

   $ pip install sphinx

To build the documentation, navigate to :literal:`gcpy/docs` and make the html target:

.. code-block:: shell-session

   gcuser:~$ cd gcpy/docs
   gcuser:~/gcpy/docs$ make html

This will generate the HTML documentation in :literal:`gcpy/docs/build/html` from the reST files in
:literal:`gcpy/docs/source`. You can view this local HTML documentation by opening
:literal:`index.html` in your web-browser.

.. note::

   You can clean the documentation with :code:`make clean`.

Learning reST
-------------

Writing reST can be a bit tricky at first. Whitespace matters (just like in Python), and some directives
can be easily miswritten. Two important things you should know right away are:

* Indents are 3-spaces
* "Things" are separated by 1 blank line (e.g., a list or code-block following a paragraph should be separated from the paragraph by 1 blank line)

You should keep these in mind when you're first getting started. Dedicating an hour to learning reST
will save you time in the long-run. Below are some good resources for learning reST.

* `reStructuredText primer <https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`_: (single best resource; however, it's better read than skimmed)
* Official `reStructuredText reference <https://docutils.sourceforge.io/docs/user/rst/quickref.html>`_ (there is *a lot* of information here)
* `Presentation by Eric Holscher <https://www.youtube.com/watch?v=eWNiwMwMcr4>`_ (co-founder of Read The Docs) at DjangoCon US 2015 (the entire presentation is good, but reST is described from 9:03 to 21:04)
* `YouTube tutorial by Audrey Tavares's <https://www.youtube.com/watch?v=DSIuLnoKLd8>`_

A good starting point would be Eric Holscher's presentations followed by reading the reStructuredText primer.

Style guidelines
----------------

.. important::  

   This documentation is written in semantic markup. This is important so that the documentation
   remains maintainable by the GEOS-Chem Support Team. Before contributing to this documentation,
   please review our style guidelines. When editing the documentation, please refrain from using
   elements with the wrong semantic meaning for aesthetic reasons. Aesthetic issues should be
   addressed by changes to the theme (not changes to reST files).

For **titles and headers**:

* H1 titles should be underlined by :literal:`#` characters
* H2 headers should be underlined by :literal:`-` characters
* H3 headers should be underlined by :literal:`^` characters
* H4 headers should be avoided, but if necessary, they should be underlined by :literal:`"` characters

**File paths** occuring in the text should use the :literal:`:literal:` role.

**Inline code**, or references to variables in code, occuring in the text should use the :literal:`:code:` role.

**Code snippets** should use :literal:`.. code-block:: <language>` directive like so

.. code-block:: none

   .. code-block:: python

      import gcpy
      print("hello world")

The language can be "none" to omit syntax highlighting. 

For command line instructions, the "console" language should be used. The :literal:`$` should be used
to denote the console's prompt. If the current working directory is relevant to the instructions,
a prompt like :literal:`gcuser:~/path1/path2$` should be used.

**Inline literals** (such as the :literal:`$` above) should use the :literal:`:literal:` role.

