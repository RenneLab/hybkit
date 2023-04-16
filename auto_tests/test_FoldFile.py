#!/usr/bin/env python3
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

"""
Automatic testing of hybkit ViennaFile Class.
"""

import os
import sys
import copy
from contextlib import nullcontext as does_not_raise
import argparse
import pytest
import hybkit


# ----- Import Testing Helper Data -----
from auto_tests.test_helper_data import *
# Includes the following variables:
# TEST_HYBID_STR, TEST_SEQ_STR, TEST_FOLD_STR, TEST_ENERGY_STR
# ART_HYB_PROPS_1, ART_HYB_PROPS_ALL, ART_BAD_HYB_STRS
# ID_ALLOWED_TYPES, SEQ_ALLOWED_TYPES, FOLD_ALLOWED_TYPES, ENERGY_ALLOWED_TYPES
# test_out_dir, vienna_autotest_file_name, hyb_file_name


# ----- Import Testing Helper Functions -----
from auto_tests.test_helper_functions import *
# Includes the following functions:
# get_expected_result_string(is_allowed=False)
# get_expected_result_context(expect_str, error_types = (TypeError, RuntimeError))

hybkit.util.set_setting('error_mode', 'raise')
hybkit.util.set_setting('iter_error_mode', 'raise')


# ----- Begin FoldFile Tests -----
# ----- Test FoldFile Base Class Misc -----
def test_foldfile_misc(tmp_path):
    fold_autotest_file_name = os.path.join(tmp_path, 'vienna_autotest_file.vienna')
    with pytest.raises(NotImplementedError):
        fold_file = hybkit.FoldFile.open(fold_autotest_file_name, 'w')
