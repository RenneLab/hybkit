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



# ----- FoldRecord Constructor Tests - Minimal -----
def test_foldrecord_constructor_minimal():
    """Test construction of FoldRecord class with minimal information."""
    # Test HybRecord Minimal Constructor:
    test_record = hybkit.FoldRecord(id=TEST_HYBID_STR, seq=TEST_SEQ_STR, fold=TEST_FOLD_STR)
    # Test "seq_id" attribute 
    assert test_record.id == TEST_HYBID_STR
    # Test "seq" attribute
    assert test_record.seq == TEST_SEQ_STR  
    # Test "fold" attribute
    assert test_record.fold == TEST_FOLD_STR
    # Test "fold" attribute
    assert test_record.fold == TEST_FOLD_STR
    # Test "energy" attribute
    assert test_record.energy == None
    
    # Test with Energy
    test_record_2 = hybkit.FoldRecord(id=TEST_HYBID_STR, seq=TEST_SEQ_STR, 
                                      fold=TEST_FOLD_STR, energy=TEST_ENERGY_STR)
    assert test_record_2.energy == TEST_ENERGY_STR

    assert test_record == test_record
    test_record_2.id = 'test2'
    assert test_record != test_record_2
    assert bool(test_record)
    print(str(test_record))
    hash(test_record)

# ----- FoldRecord Type Tests -----
test_seg_props = copy.deepcopy(EMPTY_SEG_PROPS)
default_constructor_args = {
    'id': TEST_HYBID_STR,
    'seq': TEST_SEQ_STR,
    'fold': TEST_FOLD_STR,
    'energy': TEST_ENERGY_STR,
}
field_allowed_types = {
    'id': ID_ALLOWED_TYPES,
    'seq': SEQ_ALLOWED_TYPES,
    'fold': FOLD_ALLOWED_TYPES,
    'energy': ENERGY_ALLOWED_TYPES,
}
test_parameters = []
# Setup test constructor types for each attribute in default_constructor_args:
for constructor_field in default_constructor_args.keys():
    # Setup testing of each possible data type for each field
    for test_name, test_object in TEST_OBJECTS.items():
        # Get types allowed for this field
        allowed_types = field_allowed_types[constructor_field]
        # Setup constructor arguments for this test
        constructor_args = copy.deepcopy(default_constructor_args)
        # Set test data for this field
        constructor_args[constructor_field] = test_object
        # Determine Error vs. Null Context for this test
        expect_result = get_expected_result_string(test_name in allowed_types)
        test_param_set = (
            constructor_field, 
            test_name, 
            expect_result,
            constructor_args, 
        )
        test_parameters.append(test_param_set)

@pytest.mark.parametrize("test_field,test_name,expect_str,test_input",[*test_parameters])
def test_foldrecord_obj_types(test_field, test_name, expect_str, test_input):
    expect_context = get_expected_result_context(expect_str)
    with expect_context:
        print(test_input)
        assert hybkit.FoldRecord(**test_input) is not None


# ----- FoldRecord test reading/writing of vienna-format records. -----
test_parameters = [
    ('NonOverlapping', 'Pass', ART_HYB_VIENNA_PROPS_1),
    ('Overlapping', 'Pass', ART_HYB_VIENNA_PROPS_2),
]
@pytest.mark.parametrize("test_name,expectation,test_props",[*test_parameters])
def test_foldrecord_vienna_io(test_name, expectation, test_props):
    hyb_str = test_props['hyb_str']
    vienna_str = test_props['vienna_str']
    vienna_id = vienna_str.splitlines()[0].split()[0].lstrip('>')
    vienna_lines = vienna_str.splitlines()
    vienna_lines_mod_last = copy.deepcopy(vienna_lines)
    vienna_lines_mod_last[-1] = vienna_lines_mod_last[-1] + '\n'
    full_seq = test_props['seg1_seq'] + test_props['seg2_seq']
    full_fold = test_props['seg1_fold'] + test_props['seg2_fold']
    energy = test_props['vienna_str'].splitlines()[-1].split()[-1].strip('()')

    fold_records = [
         hybkit.FoldRecord.from_vienna_string(vienna_str),
         hybkit.FoldRecord.from_vienna_lines(vienna_str.splitlines())
    ]
    for fold_record in fold_records:
        assert fold_record.id == vienna_id
        assert fold_record.seq == full_seq
        assert fold_record.fold == full_fold
        assert fold_record.energy == energy
    
    assert fold_record.to_vienna_lines() == vienna_lines
    assert ''.join(fold_record.to_vienna_lines(True)) == '\n'.join(vienna_lines_mod_last)
    assert fold_record.to_vienna_string() == vienna_str
    assert fold_record.to_vienna_string(True) == (vienna_str + '\n')


# ----- FoldRecord test reading vienna-format records with errors. -----
test_parameters = [
    ('empty3', {'record_lines': ['', '', '']}),
    ('empty2', {'record_lines': ['', '']}),
    ('empty1', {'record_lines': ['abc', '123', 'singleval']}),
    ('empty0', {'record_lines': ['abc', '123', '(99']}),
]
@pytest.mark.parametrize("test_name,test_kws",[*test_parameters])
def test_foldrecord_vienna_io_errors(test_name, test_kws):
    with pytest.raises(RuntimeError):
        hybkit.FoldRecord.from_vienna_lines(**test_kws)

# ----- FoldRecord test reading vienna-format records with allowed errors. -----
test_parameters = [
    ('warn-99', {'record_lines': ['abc', '123', '(99'], 'error_mode': 'warn_return'}),
]
@pytest.mark.parametrize("test_name,test_kws",[*test_parameters])
def test_foldrecord_vienna_io_warnings(test_name, test_kws):
    hybkit.FoldRecord.from_vienna_lines(**test_kws)
