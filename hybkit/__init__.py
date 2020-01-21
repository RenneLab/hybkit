#!/usr/bin/env python3
# Daniel B. Stribling
# Renne Lab, University of Florida
# Hybkit Project : http://www.github.com/RenneLab/hybkit

'''
Package-level initialization of hybkit.
Public classes and methods of hybkit.hybkit_code are imported so they are accessible
as hybkit.HybRecord() ... etc.
'''
import os

# Set package directory and code directory
package_dir = os.path.abspath(os.path.dirname(__file__))
code_dir = os.path.join(package_dir, 'code_dir')

# Import module-level dunder-names:
from hybkit.__about__ import __author__, __contact__, __credits__, __date__, __depreciated__, \
                             __email__, __license__, __maintainer__, __status__, __version__

# Import public classes and methods of hybkit_code
from hybkit.hybkit_code import HybRecord, HybFile, \
                               FoldRecord, \
                               ViennaFile, HybViennaIter, HybViennaCmbIter, \
                               ViennadFile, HybViennadIter, HybViennadCmbIter, \
                               CtFile, HybCtIter, HybCtCmbIter
