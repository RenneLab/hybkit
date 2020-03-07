#!/usr/bin/env python3
# Daniel Stribling
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

'''
Package details for the hybkit project.
'''

import os

# Set package directory, code directory, and default file locations
package_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
code_dir = os.path.join(package_dir, 'hybkit')
data_dir = os.path.join(package_dir, 'databases')
default_string_match_params = os.path.join(code_dir, 'string_match_params.csv')
# default_id_map_params = os.path.join()

# Preparation for initial PIP release, not yet completed.
project_name = 'hybkit'
description = 'Toolkit for analysis of .hyb format genomic sequence data.'
project_url = 'https://github.com/RenneLab/hybkit'
keywords = 'genetics genomics ribonomics bioinformatics hyb CLASH qCLASH miRNA '
keywords += 'RNA DNA vienna viennad unafold'

# For a list of valid classifiers, see https://pypi.org/classifiers/
classifiers = [
    'Development Status :: 2 - Pre-Alpha',
    'Natural Language :: English',
    'Intended Audience :: Science/Research',
    'Intended Audience :: Developers',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Topic :: Software Development :: Libraries :: Python Modules',
    'Operating System :: OS Independent',
    'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
    ]

info_urls = {
    'Renne Lab Github': 'https://github.com/RenneLab',
    'Renne Lab Mainpage': 'https://www.rennelab.com/',
    'Hyb Format Specification':
    'https://www.sciencedirect.com/science/article/pii/S1046202313004180',
    }

__author__ = "Daniel Stribling"
__contact__ = "ds@ufl.edu"
__credits__ = ["Daniel B. Stribling", "Rolf Renne"]
__date__ = "YYYY/MM/DD"
__deprecated__ = False
__email__ = "ds@ufl.edu"
__license__ = "GPLv3"
__maintainer__ = "Renne Group, University of Florida"
__status__ = "Development"
__version__ = "0.1.0"
spec_version = __version__  # Define separate specification version.
