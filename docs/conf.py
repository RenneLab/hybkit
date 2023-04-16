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

"""Sphinx configuration file for hybkit documentation."""

import os
import sys
import copy
import imp
import pprint
import textwrap
from docutils.parsers.rst import Directive
from docutils import nodes, statemachine

sys.path.insert(0, os.path.abspath('..'))
# sys.path.insert(0, os.path.abspath(os.path.join('..', 'scripts')))
import hybkit
conf_dir = os.path.abspath('.')
main_dir = os.path.abspath('..')
scripts_dir = os.path.join(main_dir, 'scripts')

# remove TypeFinder links from HybRecord class for documentation
hybkit.HybRecord.TypeFinder = None

# import hybkit_scripts
hyb_check = imp.load_source(
    'hyb_check',
    os.path.join(scripts_dir, 'hyb_check'),
)
hyb_filter = imp.load_source(
    'hyb_filter',
    os.path.join(scripts_dir, 'hyb_filter'),
)
hyb_eval = imp.load_source(
    'hyb_eval',
    os.path.join(scripts_dir, 'hyb_eval'),
)
hyb_analyze = imp.load_source(
    'hyb_analyze',
    os.path.join(scripts_dir, 'hyb_analyze'),
)

# -- Project information -----------------------------------------------------

project = hybkit.__about__.project_name
copyright = hybkit.__about__.__copyright__
author = hybkit.__about__.__author__
# The full version, including alpha/beta/rc tags
version = '.'.join(hybkit.__about__.__version__.split('.'))
release = hybkit.__about__.__version__

# -- Create Docs Index File --------------------------------------------------
with open(os.path.join(main_dir, 'README.rst'), 'r') as readme_file, \
     open(os.path.join(conf_dir, 'index.rst'), 'w') as index_file:

    for line in readme_file:
        if line.startswith('.. Github Only'):
            break
        index_file.write(line)
    index_file.write('.. include:: ./index_suffix.rst\n')


# -- PPrint Functions -----------------------------------------------------------
def return_pprint_code_block(in_item, prefix_indent=8, obj_indent=1, width=85, item_name=''):
    """Return a string of a pretty-printed code block for the given item."""
    ret_str = '\n.. code-block:: python\n\n'
    ptext = pprint.pformat(in_item, indent=obj_indent, compact=False,
                           width=width)
    if item_name and ptext.startswith('{'):
        ptext = ('%s = {\n' % item_name) + ptext[1:]
    if ptext.endswith('}'):
        ptext = ptext[:-1] + '\n}'
    for line in ptext.split('\n'):
        ret_str += (' ' * prefix_indent) + line + '\n'
    return ret_str


def return_settings_info_block(settings_name, prefix_indent=8, obj_indent=1, width=85):
    """Return a specially-formatted block for settings information."""
    item_name = settings_name
    settings_obj = getattr(hybkit.settings, settings_name.split('.')[-1])
    use_info = copy.deepcopy(settings_obj)
    new_settings_info = {}
    for key in use_info:
        use_info[key][1] = ' '.join(use_info[key][1].split())
        new_settings_info[key] = {
            'Def-Val': use_info[key][0],
            'Desc.': use_info[key][1],
            'Argp-Type': use_info[key][2],
            'Argp-Flag': use_info[key][3],
            'Argp-Opts': use_info[key][4],
        }
        # new_settings_info[key] = use_info[key]
    return return_pprint_code_block(new_settings_info, prefix_indent=prefix_indent,
                                    obj_indent=obj_indent, width=width, item_name=item_name)


# -- PrettyPrint Directive -----------------------------------------------------------
# Inspriaton: https://stackoverflow.com/a/18143318
class PPrintDictDirective(Directive):
    """Execute the specified python code and insert the output into the document."""

    has_content = True

    def run(self):
        """Run the directive."""
        tab_width = self.options.get('tab-width', self.state.document.settings.tab_width)
        source = self.state_machine.input_lines.source(
            self.lineno - self.state_machine.input_offset - 1)
        try:
            text = return_settings_info_block((''.join(self.content).strip()))
            lines = statemachine.string2lines(text, tab_width, convert_whitespace=True)
            self.state_machine.insert_input(lines, source)
            return []
        except Exception:
            return [
                nodes.error(None, nodes.paragraph(
                    text="Unable to prettyprint dict at %s:%d:" % (source, self.lineno)),
                    nodes.paragraph(text=str(sys.exc_info()[1]))
                )
            ]


def setup(app):
    """Set / setup the directive."""
    app.add_directive('ppdict', PPrintDictDirective)
    app.add_config_value('on_github', False, 'env')

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.todo',
    'sphinx.ext.autosectionlabel',
    'sphinx.ext.intersphinx',
    'sphinx.ext.ifconfig',
    'sphinxarg.ext',
    'sphinx_rtd_theme',
]

# Latex/PDF output options:
latex_toplevel_sectioning = 'chapter'
latex_show_pagerefs = True
latex_show_urls = 'footnote'
latex_elements = {
    'extraclassoptions': 'openany,oneside'
}

# Formatting Options
smartquotes = False

# add_module_names
autodoc_member_order = 'bysource'
autosectionlabel_maxdepth = 1
# autosectionlabel_prefix_document = True
master_doc = 'index'

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', 'setup.py']


intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    # 'matplotlib': ('https://readthedocs.org/projects/matplotlib/latest/', None)
}
napoleon_use_ivar = True

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ['_static']

# Fix to RTD table wrapping: https://rackerlabs.github.io/docs-rackspace/tools/rtd-tables.html
# html_context = {
#    'css_files': [
#        '_static/theme_overrides.css',  # override wide tables in RTD theme
#        ],
#     }

# Define custom variables
rst_epilog = textwrap.dedent(
    """
    .. |3p| replace:: :abbr:`3p (3-Prime)`
    .. |5p| replace:: :abbr:`5p (5-Prime)`
    .. |spec_version| replace:: %s
    .. _UNAFold: http://www.unafold.org/
    .. _ViennaRNA: https://www.tbi.univie.ac.at/RNA/

    """ % (hybkit.__about__.spec_version,)
)
