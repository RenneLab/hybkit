#!/usr/bin/env python3
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

'''
Package details for the hybkit project.
'''

import os
import sys
if not (sys.version_info.major >= 3 and sys.version_info.minor >= 7):
    message = 'Python 3.7+ is required for hybkit.'
    print(message)
    raise Exception(message)

from importlib import resources
with resources.path('hybkit', '__init__.py') as path_obj:
    module_dir = os.path.dirname(path_obj)

project_name = 'hybkit'
version = "0.3.0a"
description = 'Toolkit for analysis of hybrid genomic sequence data.'
project_url = 'https://github.com/RenneLab/hybkit'
keywords = 'genetics genomics ribonomics bioinformatics hyb CLASH qCLASH miRNA '
keywords += 'RNA DNA vienna viennad unafold'
name_and_version = project_name + '-' + version

prefix_data_dir = os.path.join(sys.prefix, name_and_version)

#Putting in try block to allow use with exec()
try:
    local_data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
except NameError:
    local_data_dir = 'using_with_exec'

local_prefix_data_dir = os.path.join(local_data_dir, name_and_version) 

if os.path.isdir(os.path.join(prefix_data_dir, 'ref_data')):
    hybkit_data_dir = prefix_data_dir
elif os.path.isdir(os.path.join(local_data_dir, 'ref_data')):
    hybkit_data_dir = local_data_dir
else:
    print('WARNING: hybkit_data_dir variable cannot be set, ignore during setup.py.')
    print('prefix_data_dir', prefix_data_dir)
    if os.path.isdir(prefix_data_dir):
        print(os.listdir(prefix_data_dir))
    print('local_data_dir', local_data_dir)
    if os.path.isdir(local_data_dir):
        print(os.listdir(local_data_dir))
    print('local_prefix_data_dir', local_prefix_data_dir)
    if os.path.isdir(local_prefix_data_dir):
        print(os.listdir(local_prefix_data_dir))
    hybkit_data_dir = ''

ref_data_dir = os.path.join(hybkit_data_dir, 'reference_data')
docs_dir = os.path.join(hybkit_data_dir, 'docs')

# For a list of valid classifiers, see https://pypi.org/classifiers/
classifiers = [
    'Development Status :: 4 - Beta',
    'Natural Language :: English',
    'Intended Audience :: Science/Research',
    'Intended Audience :: Developers',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Topic :: Software Development :: Libraries :: Python Modules',
    'Operating System :: OS Independent',
    'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
    'Programming Language :: Python :: 3.9',
    ]

info_urls = {
    'Download as TAR': ('https://github.com/RenneLab/hybkit/tarball/' + version),
    'Renne Lab Github': 'https://github.com/RenneLab',
    'Renne Lab Mainpage': 'https://www.rennelab.com/',
    'Hyb Format Specification':
    'https://www.sciencedirect.com/science/article/pii/S1046202313004180',
    }

keywords = 'genetics genomics ribonomics bioinformatics CLASH qCLASH miRNA'

__author__ = "Daniel Stribling"
__contact__ = "ds@ufl.edu"
__credits__ = ["Daniel Stribling", "Rolf Renne"]
__date__ = "2021/09/30"
__deprecated__ = False
__email__ = "ds@ufl.edu"
__license__ = "GPLv3+"
__maintainer__ = "Renne Group, University of Florida"
__status__ = "Development"
__version__ = version
spec_version = __version__  # Optionally define separate specification version.
