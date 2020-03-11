# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
import imp
sys.path.insert(0, os.path.abspath('..'))
# sys.path.insert(0, os.path.abspath(os.path.join('..', 'scripts')))
import hybkit
#import hybkit_scripts
hyb_check = imp.load_source('hyb_check', 
                            os.path.abspath(os.path.join('..', 'scripts', 'hyb_check')))
hyb_filter = imp.load_source('hyb_filter',
                             os.path.abspath(os.path.join('..', 'scripts', 'hyb_filter')))
hyb_analysis = imp.load_source('hyb_analysis',
                                os.path.abspath(os.path.join('..', 'scripts', 'hyb_analysis')))
hyb_check = imp.load_source('hyb_check',
                            os.path.abspath(os.path.join('..', 'scripts', 'hyb_check')))


# -- Project information -----------------------------------------------------

project = hybkit.__about__.project_name
copyright = '2020, ' + hybkit.__about__.__author__
author = hybkit.__about__.__author__

# The full version, including alpha/beta/rc tags
version = '.'.join(hybkit.__about__.__version__.split('.'))
release = hybkit.__about__.__version__


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc', 
    'sphinx.ext.napoleon',
    'sphinx.ext.todo',
    'sphinx.ext.autosectionlabel',
    'sphinx.ext.intersphinx',
    'sphinxarg.ext',
    ]

# add_module_names
autodoc_member_order = 'bysource'
autosectionlabel_maxdepth = 1
#autosectionlabel_prefix_document = True
master_doc = 'index'

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', 'setup.py']


intersphinx_mapping = {
                       'python': ('https://docs.python.org/3', None),
                       #'matplotlib': ('https://readthedocs.org/projects/matplotlib/latest/', None)
                      }
napoleon_use_ivar = True

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# html_theme = 'alabaster'
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Fix to RTD table wrapping: https://rackerlabs.github.io/docs-rackspace/tools/rtd-tables.html
html_context = {
    'css_files': [
        '_static/theme_overrides.css',  # override wide tables in RTD theme
        ],
     }

# Define custom variables
rst_epilog = (
"""
.. |3p| replace:: :abbr:`3p (3-Prime)`
.. |5p| replace:: :abbr:`5p (5-Prime)`
.. |spec_version| replace:: %s
""" % (
       hybkit.__about__.spec_version,
      )
)

