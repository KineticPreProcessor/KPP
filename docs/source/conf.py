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

master_doc = "index"

# -- Project information -----------------------------------------------------

project = "KPP: The Kinetic PreProcessor"
copyright = "2023, The KPP Development Team"
author = "A. Sandu, R. Sander, M. Long, H. Lin, R. Yantosca, L. Estrada, L. Shen, and D. J. Jacob"

# The full version, including alpha/beta/rc tags
# (version numbers must be synchronized in CHANGELOG.md, src/gdata.h,
# and docs/source/conf.py)
release = "3.0.1"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx_rtd_theme",   
    "sphinxcontrib.bibtex",
    "recommonmark",
]
bibtex_default_style = "refstyle"
bibtex_reference_style = "author_year"

from pybtex.style.formatting.alpha import Style as AlphaStyle
from pybtex.style.names.lastfirst import NameStyle as LastFirst
from pybtex.style.template import join, words, optional, sentence
from pybtex.style.labels import BaseLabelStyle

class LabelStyle(BaseLabelStyle):
    def format_labels(self, sorted_entries):
        for entry in sorted_entries:
            yield entry.key.replace("_", " ").replace("et al.", "et al.,")

class RefStyle(AlphaStyle):
    # We have to use the Pybtex class for alphabetic sorting,
    # which is imported and renamed to AlphaStyle.  This will
    # sort citations by author last name, and then by year (since
    # we have used the "author_year" reference style above.
    #   -- Bob Yantosca (28 Jun 2022)
    default_name_style = LastFirst
    default_sort_style = None
    default_label_style = LabelStyle
    
    def __init__(self):
       super().__init__()
       self.abbreviate_names = True
      #  self.label_style = KeyLabelStyle()
      #  self.format_labels = self.label_style.format_labels

    def format_web_refs(self, e):
       return sentence[ optional[ self.format_doi(e) ], ]

from pybtex.plugin import register_plugin
register_plugin("pybtex.style.formatting", "refstyle", RefStyle)

bibtex_bibliography_header = ""
bibtex_footbibliography_header = bibtex_bibliography_header
bibtex_bibfiles = ["citations/kpp.bib"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"
#html_theme = "alabaster"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

# These paths are either relative to html_static_path
# or fully qualified paths (e.g. https://...)
html_css_files = ["custom.css"]

# Display KPP logo
html_favicon = "_static/kpp-favicon.ico"
html_logo = "_static/kpp-logo.png"

# RTD theme settings
html_theme_options = {
    'logo_only': True,
    'display_version': False,
    'style_nav_header_background': '#FCFCFC',
}

# -- Options for PDF output via LaTeX ----------------------------------------

# https://www.sphinx-doc.org/en/master/latex.html
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-latex-output

latex_theme = "howto"
latex_logo = "_static/kpp-logo.png"
#latex_show_pagerefs = True
latex_show_urls = "footnote"

latex_elements = {
    "papersize" : "a4paper",
    "pointsize": "12pt"
}
