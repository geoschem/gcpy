# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

import os
import sys

# Point Sphinx to the root of the GCPy repository so that autodoc
# can import gcpy modules directly.
_REPO_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..")
)
sys.path.insert(0, _REPO_ROOT)

# -- Project information -----------------------------------------------------

project = 'GCPy'
copyright = '2023, GEOS-Chem Support Team'
author = 'GEOS-Chem Support Team'

# The full version, including alpha/beta/rc tags
release = '1.8.0'


# -- General configuration ---------------------------------------------------

extensions = [
    # --- Core Sphinx extensions ---
    "sphinx.ext.autodoc",        # Pulls docstrings from source code
    "sphinx.ext.autosummary",    # Generates summary tables + per-item pages
    "sphinx.ext.napoleon",       # Parses NumPy & Google style docstrings
    "sphinx.ext.viewcode",       # Adds [source] links to each function/class
    "sphinx.ext.intersphinx",    # Cross-links to numpy, xarray, etc.
    "sphinx.ext.todo",           # Support for .. todo:: directives
    "sphinx.ext.coverage",       # Reports docstring coverage

    # --- Third-party extensions ---
    "sphinx_rtd_theme",
    "sphinxcontrib.bibtex",
    "myst_parser",               # Parses Markdown (.md) source files
]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}


# -- BibTeX configuration ----------------------------------------------------

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

    def format_web_refs(self, e):
        return sentence[ optional[ self.format_doi(e) ], ]

from pybtex.plugin import register_plugin
register_plugin('pybtex.style.formatting', 'gcrefstyle', GCRefStyle)

bibtex_bibliography_header = ".. rubric:: References"
bibtex_footbibliography_header = bibtex_bibliography_header
bibtex_bibfiles = []


# -- autodoc configuration ---------------------------------------------------

# Show type hints in the description body, not in the signature line
autodoc_typehints = "description"

# Include both class docstring and __init__ docstring
autoclass_content = "both"

# List members in the order they appear in the source file
autodoc_member_order = "bysource"

# Default options applied to every autodoc directive
autodoc_default_options = {
    "members":           True,   # Document all public members
    "undoc-members":     True,   # Include members without docstrings
    "private-members":   False,  # Exclude _private members
    "special-members":   False,  # Exclude __dunder__ members
    "inherited-members": False,  # Exclude inherited members
    "show-inheritance":  True,   # Show base classes in class pages
}


# -- autosummary configuration -----------------------------------------------

# Automatically generate stub .rst pages for each module/class/function
autosummary_generate = True

# Overwrite existing stub files on each build
autosummary_generate_overwrite = True

# Do not document names that were imported from other modules
autosummary_imported_members = False


# -- napoleon configuration ---------------------------------------------------

# We use NumPy-style docstrings throughout GCPy
napoleon_google_docstring            = False
napoleon_numpy_docstring             = True
napoleon_include_init_with_doc       = True
napoleon_include_private_with_doc    = False
napoleon_include_special_with_doc    = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes    = True
napoleon_use_admonition_for_references = False
napoleon_use_ivar                    = False
napoleon_use_param                   = True
napoleon_use_rtype                   = True
napoleon_preprocess_types            = True
napoleon_attr_annotations            = True


# -- intersphinx configuration -----------------------------------------------

# Cross-link to the API docs of packages that GCPy depends on
intersphinx_mapping = {
    "python":     ("https://docs.python.org/3",             None),
    "numpy":      ("https://numpy.org/doc/stable",          None),
    "xarray":     ("https://docs.xarray.dev/en/stable",     None),
    "pandas":     ("https://pandas.pydata.org/docs",        None),
    "scipy":      ("https://docs.scipy.org/doc/scipy",      None),
    "matplotlib": ("https://matplotlib.org/stable",         None),
}


# -- General build options ---------------------------------------------------

# Patterns to ignore when looking for source files
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'

html_theme_options = {
    'logo_only':                       False,
    'style_nav_header_background':     '#FCFCFC',
    'navigation_depth':                4,
    'collapse_navigation':             False,
    'sticky_navigation':               True,
    'includehidden':                   True,
    'titles_only':                     False,
    'prev_next_buttons_location':      'bottom',
    'style_external_links':            True,
}

# Paths containing custom static files (relative to this directory)
html_static_path = [
    '_static',
]

# CSS files that override sphinx-rtd-theme defaults
# (paths are relative to _static, specified above)
html_css_files = [
    'css/icon_home.css',
    'css/theme_overrides.css',
]

# Display GEOS-Chem favicon and logo
html_favicon = '_static/images/gc-o-logo-favicon.ico'
html_logo    = '_static/images/GEOS-Chem_Logo_Light_Background.png'
