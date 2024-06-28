#!/usr/bin/env python3
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

"""Package details for the hybkit project."""

# To update:
#   Change Version
#   Change Date


import logging
import os
import sys
from importlib import resources

# Check Python version and throw error if not 3.8+
if sys.version_info < (3, 8):  #noqa: UP036
    message = 'Python 3.8+ is required for hybkit.'
    message += ' Current version is ' + str(sys.version_info.major) + '.'
    message += str(sys.version_info.minor) + '.' + str(sys.version_info.micro)
    message += '.' + str(sys.version_info.releaselevel)
    raise RuntimeError(message)

# Set hybkit module directory
with resources.path('hybkit', '__init__.py') as path_obj:
    module_dir = os.path.dirname(path_obj)

# Hybkit information
project_name = 'hybkit'
version = 'v0.3.5'
python_requires = '>=3.8'
description = 'Toolkit for analysis of chimeric (hybrid) RNA sequence data.'
project_url = 'https://github.com/RenneLab/hybkit'
keywords = 'genetics genomics ribonomics bioinformatics hyb CLASH qCLASH miRNA '
keywords += 'RNA DNA vienna UNAfold'
name_and_version = project_name + '-' + version

prefix_data_dir = os.path.join(sys.prefix, name_and_version)

# Identify paths for data files used in hybkit
try:
    local_data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
except NameError:
    local_data_dir = 'using_with_exec'
local_prefix_data_dir = os.path.join(local_data_dir, name_and_version)
if os.path.isdir(os.path.join(prefix_data_dir, 'ref_data')):
    hybkit_data_dir = prefix_data_dir
elif os.path.isdir(os.path.join(local_data_dir, 'ref_data')):
    hybkit_data_dir = local_data_dir
elif os.path.isdir(os.path.join(local_prefix_data_dir, 'ref_data')):
    hybkit_data_dir = local_prefix_data_dir
else:
    message = 'hybkit_data_dir variable cannot be set, ignore during setup.py.'
    logging.warning(message)
    hybkit_data_dir = ''

ref_data_dir = os.path.join(hybkit_data_dir, 'reference_data')
docs_dir = os.path.join(hybkit_data_dir, 'docs')

# Python package classifiers for PyPI
# For a list of valid classifiers, see https://pypi.org/classifiers/
classifiers = [
    'Development Status :: 5 - Production/Stable',
    'Natural Language :: English',
    'Intended Audience :: Science/Research',
    'Intended Audience :: Developers',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Topic :: Software Development :: Libraries :: Python Modules',
    'Operating System :: OS Independent',
    'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
    'Programming Language :: Python :: 3.8',
    'Programming Language :: Python :: 3.9',
    'Programming Language :: Python :: 3.10',
    'Programming Language :: Python :: 3.11',
]

info_urls = {
    'Download as TAR': ('https://github.com/RenneLab/hybkit/tarball/' + version),
    'Renne Lab Github': 'https://github.com/RenneLab',
    'Renne Lab Mainpage': 'https://www.rennelab.com/',
    'Travis Hyb Format Specification':
    'https://www.sciencedirect.com/science/article/pii/S1046202313004180',
}

keywords = 'genetics, genomics, microRNAs, Ribonomics, Hyb, Hybrids, CLASH, qCLASH, '
keywords += 'CLEAR-CLIP, Chimeric e-CLIP'

__author__ = 'Daniel Stribling'
__contact__ = 'ds@ufl.edu'
__credits__ = ['Daniel Stribling', 'Rolf Renne']
__copyright__ = '2023, ' + __author__
__date__ = '2023/11/20'
__deprecated__ = False
__email__ = 'ds@ufl.edu'
__license__ = 'GPLv3+'
__maintainer__ = 'Renne Group, University of Florida'
__status__ = 'Development'
__version__ = version
# Hyb file specification version.
spec_version = __version__
