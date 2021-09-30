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
import copy

sys.path.insert(0, os.path.abspath('..'))
# sys.path.insert(0, os.path.abspath(os.path.join('..', 'scripts')))
import hybkit

# remove TypeFinder and RegionFinder links from HybRecord class for documentation
hybkit.HybRecord.TypeFinder = None
hybkit.HybRecord.RegionFinder = None

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
copyright = '2021, ' + hybkit.__about__.__author__
author = hybkit.__about__.__author__
# The full version, including alpha/beta/rc tags
version = '.'.join(hybkit.__about__.__version__.split('.'))
release = hybkit.__about__.__version__


# -- Exec Directive -----------------------------------------------------------
# Source: https://stackoverflow.com/a/18143318
from io import StringIO
from docutils.parsers.rst import Directive    
from docutils import nodes, statemachine

class ExecDirective(Directive):
    """Execute the specified python code and insert the output into the document"""
    has_content = True

    def run(self):
        oldStdout, sys.stdout = sys.stdout, StringIO()

        tab_width = self.options.get('tab-width', self.state.document.settings.tab_width)
        source = self.state_machine.input_lines.source(self.lineno - self.state_machine.input_offset - 1)

        try:
            exec('\n'.join(self.content))
            text = sys.stdout.getvalue()
            lines = statemachine.string2lines(text, tab_width, convert_whitespace=True)
            self.state_machine.insert_input(lines, source)
            return []
        except Exception:
            return [nodes.error(None, nodes.paragraph(text = "Unable to execute python code at %s:%d:" % (basename(source), self.lineno)), nodes.paragraph(text = str(sys.exc_info()[1])))]
        finally:
            sys.stdout = oldStdout

def setup(app):
    app.add_directive('expy', ExecDirective)

# -- PPrint Functions -----------------------------------------------------------
import pprint
def pprint_code_block(in_item, prefix_indent=8, obj_indent=1, width=92, item_name=''):
    print('\n.. code-block::\n')
    ptext = pprint.pformat(in_item, indent=obj_indent, compact=False, 
                           width=width)
    if item_name and ptext.startswith('{'):
        ptext = ('%s = {\n' % item_name ) + ptext[1:]
    if ptext.endswith('}'):
       ptext = ptext[:-1] + '\n}'
    for line in ptext.split('\n'):
        print((' '*prefix_indent)+line)

def pprint_settings_info_block(settings_info, prefix_indent=8, obj_indent=1, 
                               width=92, item_name=''):
    use_info = copy.deepcopy(settings_info)
    new_settings_info = {}
    for key in settings_info:
        use_info[key][1] = ' '.join(use_info[key][1].split())
        new_settings_info[key] = {'Def-Val': use_info[key][0],
                                  'Desc.': use_info[key][1], 
                                  'Argp-Type': use_info[key][2], 
                                  'Argp-Flag': use_info[key][3], 
                                  'Argp-Opts': use_info[key][4]
                                 } 
        #new_settings_info[key] = use_info[key]
    pprint_code_block(new_settings_info, prefix_indent=prefix_indent, obj_indent=obj_indent, 
                      width=width, item_name=item_name)


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
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
#html_static_path = ['_static']

# Fix to RTD table wrapping: https://rackerlabs.github.io/docs-rackspace/tools/rtd-tables.html
#html_context = {
#    'css_files': [
#        '_static/theme_overrides.css',  # override wide tables in RTD theme
#        ],
#     }

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

