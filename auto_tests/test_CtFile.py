#!/usr/bin/env python3
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

"""
Automatic testing of hybkit CtFile Class.
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
# test_out_dir, Ct_autotest_file_name, hyb_file_name


# ----- Import Testing Helper Functions -----
from auto_tests.test_helper_functions import *
# Includes the following functions:
# get_expected_result_string(is_allowed=False)
# get_expected_result_context(expect_str, error_types = (TypeError, RuntimeError))

# ----- Start CtFile Test IO of Ct Strings / FoldRecords -----
#TODO: Implement CT File Tests
# test_parameters = []
# for prop_set in [ART_HYB_Ct_PROPS_1, ART_HYB_Ct_PROPS_2]:
#     if prop_set['overlapping']:
#         test_name = 'Overlapping'
#     else:
#         test_name = 'Static'
#     test_parameters.append(
#         (test_name, 'Pass', prop_set)
#     )

# @pytest.mark.parametrize("test_name,expectation,test_props",[*test_parameters])
# def test_Ctfile_io(test_name, expectation, test_props):
#     if not os.path.isdir(test_out_dir):
#         os.mkdir(test_out_dir)
#     expect_context = get_expected_result_context(expectation)
#     Ct_str = test_props['Ct_str']
#     with open(Ct_autotest_file_name, 'w') as Ct_autotest_file:
#         Ct_autotest_file.write(Ct_str)

#     assert hybkit.util.Ct_exists(Ct_autotest_file_name)

#     with expect_context:
#         with hybkit.CtFile.open(Ct_autotest_file_name, 'r') as Ct_autotest_file:
#             for Ct_record in Ct_autotest_file:
#                 gen_Ct_str = Ct_record.to_Ct_string()
#                 assert gen_Ct_str == Ct_str

#     if expectation.lower() == 'pass':
#         with hybkit.CtFile.open(Ct_autotest_file_name, 'r') as Ct_autotest_file:
#             first_record = Ct_autotest_file.read_record()
#         with hybkit.CtFile.open(Ct_autotest_file_name, 'r') as Ct_autotest_file:
#             all_records = Ct_autotest_file.read_records()
#         assert first_record == all_records[0]
#         with hybkit.CtFile(Ct_autotest_file_name, 'w') as Ct_autotest_file:
#             Ct_autotest_file.write_records([Ct_record, Ct_record])
#             Ct_autotest_file.write_record(Ct_record)
#             Ct_autotest_file.write_fh(Ct_str)
