#!/usr/bin/env python3
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

"""
Automatic testing of hybkit ViennaFile Class.
"""

# ruff: noqa: ANN001 ANN201

import os

# import sys
# import copy
# from contextlib import nullcontext as does_not_raise
# import argparse
import pytest

import hybkit

# from auto_tests.test_helper_data import ()
# from auto_tests.test_helper_functions import ()

hybkit.util.set_setting('error_mode', 'raise')
hybkit.util.set_setting('iter_error_mode', 'raise')


# ----- Begin FoldFile Tests -----
# ----- Test FoldFile Base Class Misc -----
def test_foldfile_misc(tmp_path):
    """Test misc methods of FoldFile base class."""
    fold_autotest_file_name = os.path.join(tmp_path, 'vienna_autotest_file.vienna')
    with pytest.raises(NotImplementedError):
        _fold_file = hybkit.FoldFile.open(fold_autotest_file_name, 'w')
