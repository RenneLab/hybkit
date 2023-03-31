#!/usr/bin/env python3
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

"""
Automatic testing of hybkit code.
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
# test_out_dir, hyb_autotest_file_name, hyb_file_name


# ----- Import Testing Helper Functions -----
from auto_tests.test_helper_functions import *
# Includes the following functions:
# get_expected_result_string(is_allowed=False)
# get_expected_result_context(expect_str, error_types = (TypeError, RuntimeError))

# ----- Start Test HybFoldIter -----
test_param_sets = []
for prop_set in [ART_HYB_VIENNA_PROPS_1, ART_HYB_VIENNA_PROPS_2]:
    if prop_set['overlapping']:
        test_name = 'Good-Overlap'
    else:
        test_name = 'Good-Static'
    test_param_sets.append(
        (test_name, 'Pass', prop_set)
    )
test_param_sets.append(('Mismatch_Seq_Static', 'Raise', ART_BAD_HYB_VIENNA_PROPS_1))
test_param_sets.append(('Disallowed-Overlap', 'Raise', ART_BAD_HYB_VIENNA_PROPS_2))
test_param_sets.append(('Mismatch_Seq_Dynamic', 'Raise', ART_BAD_HYB_VIENNA_PROPS_3))

test_parameters = []
for test_param_set in test_param_sets:
    for combine in ['Combine', 'Separate']:
        test_parameters.append((*test_param_set, combine))

@pytest.mark.parametrize("test_name,expectation,test_props,combine_str",[*test_parameters])
def test_hybfolditer_io(test_name, expectation, test_props, combine_str):
    if not os.path.isdir(test_out_dir):
        os.mkdir(test_out_dir)
    expect_context = get_expected_result_context(expectation)
    combine = (combine_str == 'Combine')
    hyb_str = test_props['hyb_str']
    vienna_str = test_props['vienna_str']
    if test_props['overlapping']:
        seq_type = 'dynamic'
    else:
        seq_type = 'static'
        
    with open(hyb_autotest_file_name, 'w') as hyb_autotest_file:
        hyb_autotest_file.write(hyb_str)
    with open(vienna_autotest_file_name, 'w') as vienna_autotest_file:
        vienna_autotest_file.write(vienna_str)

    assert hybkit.util.hyb_exists(hyb_autotest_file_name)
    assert hybkit.util.vienna_exists(vienna_autotest_file_name)

    with pytest.raises(RuntimeError):
        use_iter= hybkit.HybFoldIter(
            hyb_autotest_file_name,
            vienna_autotest_file_name,
            combine=combine
        )

    hyb_file = hybkit.HybFile(hyb_autotest_file_name, 'r')
    fold_file = hybkit.ViennaFile(vienna_autotest_file_name, 'r', seq_type=seq_type)
    use_iter= hybkit.HybFoldIter(
        hyb_file,
        fold_file,
        combine=combine
    )
    for ret_items in use_iter:
        with expect_context:
            if combine:
                hyb_record = ret_items
                fold_record = hyb_record.fold_record
            else:
                hyb_record, fold_record = ret_items
                hyb_record.set_fold_record(fold_record)

            assert hyb_record.to_line() == hyb_str
            assert fold_record.to_vienna_string() == vienna_str

            assert fold_record.matches_hyb_record(hyb_record)
            assert fold_record.count_hyb_record_mismatches(hyb_record) == 0
            fold_record.ensure_matches_hyb_record(hyb_record)

            print(use_iter.report())
            use_iter.print_report()

        if expectation == 'Raise' and not combine:
            assert not fold_record.matches_hyb_record(hyb_record)
            assert fold_record.count_hyb_record_mismatches(hyb_record) >= test_props['mismatches']
            with pytest.raises(RuntimeError):
                fold_record.ensure_matches_hyb_record(hyb_record)


# Test iter_error_mode values:
raise_error_modes = ['raise']
allow_error_modes = ['warn_return', 'warn_skip', 'skip', 'return']
all_iter_error_modes = [*raise_error_modes, *allow_error_modes]
test_names_allowed = {
    'StaticAllowed': all_iter_error_modes,
    'StaticDisallowed': allow_error_modes,
    'DynamicAllowed': all_iter_error_modes,
    'MismatchSeqStatic': allow_error_modes,
    'MismatchSeqDynamic': allow_error_modes,
    'MismatchEnergy': allow_error_modes,
}
test_param_sets = {
    'StaticAllowed': ART_HYB_VIENNA_PROPS_1,
    'StaticDisallowed': ART_HYB_VIENNA_PROPS_2,
    'DynamicAllowed': ART_HYB_VIENNA_PROPS_2,
    'MismatchSeqStatic': ART_BAD_HYB_VIENNA_PROPS_1,
    'MismatchSeqDynamic': ART_BAD_HYB_VIENNA_PROPS_3,
    'MismatchEnergy': ART_BAD_HYB_VIENNA_PROPS_4,
}
test_seq_types = {
    'StaticAllowed': 'static',
    'StaticDisallowed': 'static',
    'DynamicAllowed': 'dynamic',
    'MismatchSeqStatic': 'static',
    'MismatchSeqDynamic': 'dynamic',
    'MismatchEnergy': 'dynamic',
}
test_parameters = []
for test_name, test_props in test_param_sets.items():
    allowed_errors = test_names_allowed[test_name]
    seq_type = test_seq_types[test_name]
    for iter_error_mode in all_iter_error_modes:
        if iter_error_mode in allowed_errors:
            expectation = 'Pass'
        else:
            expectation = 'Raise'
        test_parameters.append((test_name, iter_error_mode, expectation, seq_type, test_props))

@pytest.mark.parametrize("test_name,iter_error_mode,expectation,seq_type,test_props",[*test_parameters])
def test_hybfolditer_iter_error_mode(test_name, iter_error_mode, expectation, seq_type, test_props):
    expect_context = get_expected_result_context(expectation)
    hyb_str = test_props['hyb_str']
    vienna_str = test_props['vienna_str']
    with open(hyb_autotest_file_name, 'w') as hyb_autotest_file:
        hyb_autotest_file.write(hyb_str)
    with open(vienna_autotest_file_name, 'w') as vienna_autotest_file:
        vienna_autotest_file.write(vienna_str)

    hyb_file = hybkit.HybFile(hyb_autotest_file_name, 'r')
    fold_file = hybkit.ViennaFile(vienna_autotest_file_name, 'r', seq_type=seq_type)
    use_iter = hybkit.HybFoldIter(
        hyb_file,
        fold_file,
        iter_error_mode=iter_error_mode,
    )
    with expect_context:
        for ret_items in use_iter:
            pass


# if __name__ == '__main__':
#     test_type_finder()
