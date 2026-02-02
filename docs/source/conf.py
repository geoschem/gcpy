# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'GCPy'
copyright = '2023, GEOS-Chem Support Team'
author = 'GEOS-Chem Support Team'

# The full version, including alpha/beta/rc tags
release = '1.7.0'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx_rtd_theme",   
    "sphinxcontrib.bibtex",
    "recommonmark",
]
bibtex_default_style = 'gcrefstyle'
bibtex_reference_style = "author_year"

from pybtex.style.formatting.alpha import Style as AlphaStyle
from pybtex.style.names.lastfirst import NameStyle as LastFirst
from pybtex.style.template import join, words, optional, sentence
from pybtex.style.labels import BaseLabelStyle

class GCLabelStyle(BaseLabelStyle):
    # Base label definition.  Here we replace text in citations
    # e.g. "et al" to "et al.".  Add others as needed.
    def format_labels(self, sorted_entries):
        for entry in sorted_entries:
            yield entry.key.replace("_", " ").replace("et al.", "et al.,")

class GCRefStyle(AlphaStyle):
    # Sorts authors alphabetically by last name
    default_name_style = LastFirst
    default_sort_style = None
    default_label_style = GCLabelStyle

    def __init__(self):
       super().__init__()
       self.abbreviate_names = True
      #  self.label_style = KeyLabelStyle()
      #  self.format_labels = self.label_style.format_labels

    def format_web_refs(self, e):
       return sentence[ optional[ self.format_doi(e) ], ]

from pybtex.plugin import register_plugin
register_plugin('pybtex.style.formatting', 'gcrefstyle', GCRefStyle)


bibtex_bibliography_header = ".. rubric:: References"
bibtex_footbibliography_header = bibtex_bibliography_header

bibtex_bibfiles = []

# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
#html_theme = 'alabaster'
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = [
    '_static',
]

# CSS files that will override sphinx-rtd-theme default settings
# (paths are relative to _static, which is specified above)
html_css_files = [
    'css/icon_home.css',
    'css/theme_overrides.css',
]

# Display GEOS-Chem favicon and logo
html_favicon = '_static/images/gc-o-logo-favicon.ico'
html_logo = "_static/images/GEOS-Chem_Logo_Light_Background.png"

# More theme settings
html_theme_options = {
    'logo_only': False,                        # Show logo & top text
    'display_version': False,                  # Don't show version number
    'style_nav_header_background': '#FCFCFC',  # 99% white for top left bkgrnd
}
