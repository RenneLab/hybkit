#!/usr/bin/env python3
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

"""
Automatic testing of the hybkit HybRecord class.
"""

import os
import sys
import copy
from contextlib import nullcontext as does_not_raise
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


# Start HybRecord Tests
# ----- HybRecord Constructor Tests - Minimal -----
def test_hybrecord_constructor_minimal():
    """Test construction of HybRecord class with minimal information."""
    # Test HybRecord Minimal Constructor:
    test_record = hybkit.HybRecord(id=TEST_HYBID_STR, seq=TEST_SEQ_STR)
    # Test "seq_id" attribute
    assert test_record.id == TEST_HYBID_STR
    # Test "seq" attribute
    assert test_record.seq == TEST_SEQ_STR
    # Test "seg1_props" attribute
    assert test_record.seg1_props == EMPTY_SEG_PROPS
    # Test "seg2_props" attribute
    assert test_record.seg2_props == EMPTY_SEG_PROPS
    # Test "flags" attribute
    assert test_record.flags == {}
    # Test "fold_record" attribute
    assert test_record.fold_record is None
    # Test "get_seg1_type" method
    assert test_record.get_seg1_type() is None
    # Test "get_seg1_type" method empty error
    with pytest.raises(RuntimeError):
        test_record.get_seg1_type(require=True)
    # Test "get_seg2_type" method
    assert test_record.get_seg2_type() is None
    # Test "get_seg2_type" method empty error
    with pytest.raises(RuntimeError):
        test_record.get_seg2_type(require=True)
    # Test "get_seg_types" method
    assert test_record.get_seg_types() == (None, None)
    # Test "get_seg_types" method empty error
    with pytest.raises(RuntimeError):
        test_record.get_seg_types(require=True)
    # Test "get_read_count" method
    assert test_record.get_read_count() is None
    # Test "get_read_count" method empty error
    with pytest.raises(RuntimeError):
        test_record.get_read_count(require=True)
    # Test "get_record_count" method
    assert test_record.get_record_count() == 1
    # Test "get_record_count" method empty error
    with pytest.raises(RuntimeError):
        test_record.get_record_count(require=True)
    # Test HybRecord Minimal Constructor with missing id
    with pytest.raises(RuntimeError):
        hybkit.HybRecord(id=None, seq=TEST_SEQ_STR)
    # Test HybRecord Minimal Constructor with missing seq
    with pytest.raises(RuntimeError):
        hybkit.HybRecord(id=TEST_HYBID_STR, seq=None)


# ----- HybRecord Constructor Failure Tests -----
default_params = [
    TEST_HYBID_STR, TEST_SEQ_STR, TEST_ENERGY_STR, TEST_SEG_PROPS_STR,
    TEST_SEG_PROPS_STR, TEST_FLAGS_STR, TEST_READ_COUNT_STR
]
test_parameters = []
# TEST_FLAGS_STR is 7th parameter
# Test undefined flag in constructor.
error_params = [''] + copy.deepcopy(default_params)
error_params[7] = {'badflag': True}
test_parameters.append(tuple(error_params))
# Test mismatched read_count and read_count_flag
error_params = [''] + copy.deepcopy(default_params)
error_params[7] = {'read_count': TEST_READ_COUNT + 1}
test_parameters.append(tuple(error_params))
# Test disallowed seg1_type flag in constructor.
error_params = [''] + copy.deepcopy(default_params)
error_params[7] = {'seg1_type': 'badtype'}
test_parameters.append(tuple(error_params))
arg_string = "test_name,test_id,test_seq,test_energy,test_seg1_props,"
arg_string += "test_seg2_props,test_flags,test_read_count"


@pytest.mark.parametrize(arg_string, [*test_parameters])
def test_hybrecord_constructor_errors(test_name, test_id, test_seq, test_energy, test_seg1_props,
                                      test_seg2_props, test_flags, test_read_count):
    """Test construction of HybRecord class with full complement of information."""
    # Test HybRecord Minimal Constructor:
    with pytest.raises((RuntimeError, TypeError)):
        test_record = hybkit.HybRecord(
            id=test_id,
            seq=test_seq,
            energy=test_energy,
            seg1_props=test_seg1_props,
            seg2_props=test_seg2_props,
            flags=test_flags,
            read_count=test_read_count,
        )


# ----- HybRecord Type Tests - Main Attributes -----
test_seg_props = copy.deepcopy(EMPTY_SEG_PROPS)
default_constructor_args = {
    'id': TEST_HYBID_STR,
    'seq': TEST_SEQ_STR,
    'energy': TEST_ENERGY_STR,
    'seg1_props': test_seg_props,
    'seg2_props': test_seg_props,
}
field_allowed_types = {
    'id': ID_ALLOWED_TYPES,
    'seq': SEQ_ALLOWED_TYPES,
    'energy': ENERGY_ALLOWED_TYPES,
    'seg1_props': SEG_PROPS_ALLOWED_TYPES,
    'seg2_props': SEG_PROPS_ALLOWED_TYPES,
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


@pytest.mark.parametrize("test_field,test_name,expect_str,test_input", [*test_parameters])
def test_hybrecord_obj_types(test_field, test_name, expect_str, test_input):
    expect_context = get_expected_result_context(expect_str)
    with expect_context:
        print(test_input)
        assert hybkit.HybRecord(**test_input) is not None


# ----- HybRecord Type Tests - Segment Property Attributes -----
# Setup test constructor types for seg_props attributes:
default_seg_props = {
    'ref_name': 'test_ref',
    'read_start': 1,
    'read_end': 21,
    'ref_start': 1,
    'ref_end': 21,
    'score': 1.0
}
default_constructor_args = {
    'id': TEST_HYBID_STR,
    'seq': TEST_SEQ_STR,
    'energy': TEST_ENERGY_STR,
    'seg1_props': default_seg_props,
    'seg2_props': default_seg_props,
}
props_allowed_types = {
    'ref_name': REF_NAME_ALLOWED_TYPES,
    'read_start': READ_START_ALLOWED_TYPES,
    'read_end': READ_END_ALLOWED_TYPES,
    'ref_start': REF_START_ALLOWED_TYPES,
    'ref_end': REF_END_ALLOWED_TYPES,
    'score': SCORE_ALLOWED_TYPES,
}
test_parameters = []
# Setup test constructor types for each segN_props dict
for prop_set in ['seg1_props', 'seg2_props']:
    for prop_field in props_allowed_types.keys():
        # Setup testing of each possible data type for each field
        for test_name, test_object in TEST_OBJECTS.items():
            # Get types allowed for this field
            allowed_types = props_allowed_types[prop_field]
            # Setup constructor arguments for this test
            constructor_args = copy.deepcopy(default_constructor_args)
            seg_args = copy.deepcopy(default_seg_props)
            seg_args[prop_field] = test_object
            constructor_args[prop_set] = seg_args
            # Determine Error vs. Null Context for this test
            expect_result = get_expected_result_string(test_name in allowed_types)
            test_param_set = (
                prop_set,
                prop_field,
                test_name,
                expect_result,
                constructor_args,
            )
            test_parameters.append(test_param_set)


@pytest.mark.parametrize("prop_set,prop_field,test_name,expect_str,test_input", [*test_parameters])
def test_hybrecord_obj_types_seg_props(prop_set, prop_field, test_name, expect_str, test_input):
    expect_context = get_expected_result_context(expect_str)
    with expect_context:
        print(test_input)
        assert hybkit.HybRecord(**test_input) is not None


# ----- HybRecord eval_type(), eval_mirna(), prop(), and mirna_detail() tests -----
test_parameters = [
    ('miRNA-miRNA', ART_HYB_PROPS_1),
    ('miRNA-coding', ART_HYB_PROPS_2),
    ('coding-miRNA', ART_HYB_PROPS_3),
    ('coding-coding', ART_HYB_PROPS_4),
]


@pytest.mark.parametrize("test_name,test_params", [*test_parameters])
def test_hybrecord_type_mirna(test_name, test_params):
    """Test Hybrecord type_eval(), eval_mirna(), mirna-associated prop(), and mirna_detail()"""
    test_record = hybkit.HybRecord.from_line(
        line=test_params['hyb_str'],
    )
    assert test_record.to_line() == test_params['hyb_str']
    test_record = hybkit.HybRecord.from_line(
        line=test_params['hyb_str'],
        hybformat_id=True,
        hybformat_ref=True,
    )
    test_record = hybkit.HybRecord.from_line(
        line=test_params['hyb_str'],
        hybformat_id=True,
    )
    with pytest.raises((ValueError, RuntimeError)):
        test_record.mirna_detail(detail='all', allow_mirna_dimers=True)

    test_record.eval_types()
    test_record.eval_mirna()
    # Check type setting
    assert test_record.get_seg1_type() == test_params['seg1_type']
    assert test_record.get_seg2_type() == test_params['seg2_type']

    # Check fasta-generation properties - hybrid
    all_fasta_record = test_record.to_fasta_record(mode='hybrid', annotate=False)
    assert str(all_fasta_record.seq) == test_params['seg1_seq'] + test_params['seg2_seq']
    assert str(all_fasta_record.id) == test_params['hybrid_fasta_id']
    all_fasta_record = test_record.to_fasta_record(mode='hybrid', annotate=True)
    assert str(all_fasta_record.id) == test_params['hybrid_fasta_id_annotate']
    all_fasta_record_string = test_record.to_fasta_str(mode='hybrid', annotate=True)
    assert all_fasta_record_string == all_fasta_record.format('fasta')

    # Check fasta-generation properties - seg1
    seg1_fasta_record = test_record.to_fasta_record(mode='seg1', annotate=False)
    assert str(seg1_fasta_record.seq) == test_params['seg1_seq']
    assert str(seg1_fasta_record.id) == test_params['seg1_fasta_id']
    seg1_fasta_record = test_record.to_fasta_record(mode='seg1', annotate=True)
    assert str(seg1_fasta_record.id) == test_params['seg1_fasta_id_annotate']

    # Check fasta-generation properties - seg2
    seg2_fasta_record = test_record.to_fasta_record(mode='seg2', annotate=False)
    assert str(seg2_fasta_record.seq) == test_params['seg2_seq']
    assert str(seg2_fasta_record.id) == test_params['seg2_fasta_id']
    seg2_fasta_record = test_record.to_fasta_record(mode='seg2', annotate=True)
    assert str(seg2_fasta_record.id) == test_params['seg2_fasta_id_annotate']

    # Check properties of hybrids
    for prop in test_params['true_prop_argsets']:
        assert test_record.prop(*prop)
    for prop in test_params['false_prop_argsets']:
        assert not test_record.prop(*prop)
    for prop in test_params['true_is_set_argsets']:
        assert test_record.is_set(*prop)
    for prop in test_params['false_is_set_argsets']:
        assert test_record.not_set(*prop)

    # Check mirna-based properties
    assert test_record.flags['miRNA_seg'] == test_params['miRNA_seg']
    if test_params['one_mirna_error']:
        with pytest.raises((ValueError, RuntimeError)):
            test_record.mirna_detail()
        with pytest.raises((ValueError, RuntimeError)):
            test_record.to_fasta_record(mode='miRNA')
        with pytest.raises((ValueError, RuntimeError)):
            test_record.to_fasta_record(mode='target')
    else:
        # Check fasta-generation properties - mirna
        mirna_fasta = test_record.to_fasta_record('miRNA', annotate=False)
        assert str(mirna_fasta.seq) == test_params['mirna_seq']
        assert str(mirna_fasta.id) == test_params['mirna_fasta_id']
        mirna_fasta = test_record.to_fasta_record('miRNA', annotate=True)
        assert str(mirna_fasta.id) == test_params['mirna_fasta_id_annotate']
        # Check fasta-generation properties - target
        target_fasta = test_record.to_fasta_record('target', annotate=False)
        assert str(target_fasta.seq) == test_params['target_seq']
        assert str(target_fasta.id) == test_params['target_fasta_id']
        target_fasta = test_record.to_fasta_record('target', annotate=True)
        assert str(target_fasta.id) == test_params['target_fasta_id_annotate']

    if test_params['two_mirna_error']:
        with pytest.raises((ValueError, RuntimeError)):
            test_record.mirna_detail(detail='all', allow_mirna_dimers=True)
        with pytest.raises((ValueError, RuntimeError)):
            test_record.to_fasta_record('miRNA')
        with pytest.raises((ValueError, RuntimeError)):
            test_record.to_fasta_record('target')
    else:
        mirna_detail_dict = test_record.mirna_detail(detail='all', allow_mirna_dimers=True)
        assert mirna_detail_dict['mirna_ref'] == test_params['mirna_ref']
        assert mirna_detail_dict['target_ref'] == test_params['target_ref']
        assert mirna_detail_dict['mirna_seg_type'] == test_params['mirna_seg_type']
        assert mirna_detail_dict['target_seg_type'] == test_params['target_seg_type']
        assert mirna_detail_dict['mirna_seq'] == test_params['mirna_seq']
        assert mirna_detail_dict['target_seq'] == test_params['target_seq']


# ----- HybRecord properties tests -----
test_parameters = [
    ('miRNA-coding-props', ART_HYB_STR_PROPS),
]


@pytest.mark.parametrize("test_name,test_params", [*test_parameters])
def test_hybrecord_props(test_name, test_params):
    """Test Hybrecord type_eval(), eval_mirna(), mirna-associated prop(), and mirna_detail()"""
    test_record = hybkit.HybRecord.from_line(
        line=test_params['hyb_str'],
        hybformat_id=True,
    )
    test_record.eval_types()
    test_record.eval_mirna()
    # Check properties of hybrids
    for prop in test_params['true_prop_argsets']:
        assert test_record.prop(*prop)
    for prop in test_params['false_prop_argsets']:
        assert not test_record.prop(*prop)
    for prop in test_params['true_is_set_argsets']:
        assert test_record.is_set(*prop)
    for prop in test_params['false_is_set_argsets']:
        assert test_record.not_set(*prop)


# ----- HybRecord Magic Methods tests -----
def test_hybrecord_magic_methods():
    test_record_1 = hybkit.HybRecord.from_line(
        ART_HYB_PROPS_1['hyb_str'],
        hybformat_id=True,
        hybformat_ref=True
    )
    test_record_2 = hybkit.HybRecord.from_line(
        ART_HYB_PROPS_2['hyb_str'],
        hybformat_id=True,
        hybformat_ref=True
    )
    test_record_2.id = 'NewID'
    print(str(test_record_1))
    assert test_record_1 == test_record_1
    assert not (test_record_1 != test_record_1)
    assert test_record_1 != test_record_2
    assert bool(test_record_1)
    assert str(test_record_1)
    assert hash(test_record_1)
    assert len(test_record_1) == 40


# ----- HybRecord misc disallowed option tests -----
test_parameters = [
    ('to_fasta_record', ('notallowed',)),
    ('is_set', ('badprop',)),
    ('prop', ('badprop',)),
    ('prop', ('any_seg_type_contains', None)),
    ('set_fold_record', (None,)),
    ('set_fold_record', ('not_fold_record',)),
    ('mirna_detail', ('disallowed_detail',)),
    ('_get_flag', ('fake_flag', True)),
    ('_make_flags_dict', ('not_dict',)),
    ('_make_flags_dict', ({'bad_flag': True},)),
    ('_parse_hybformat_id', ('bad_id_name_continues_on',)),
    ('_parse_hybformat_ref', ('bad_ref_name_continues_on',)),
    ('_read_flags', ('bad_flag=B;bad_flag2=C;',)),
]
# TODO: assert not test_record.prop('has_indels')


@pytest.mark.parametrize("method,badval", [*test_parameters])
def test_hybrecord_misc_disallowed(method, badval):
    test_record = hybkit.HybRecord.from_line(
        ART_HYB_PROPS_1['hyb_str'],
        hybformat_id=True,
        hybformat_ref=True,
    )
    with pytest.raises(RuntimeError):
        getattr(test_record, method)(*badval)


test_parameters = [
    ('_ensure_props_read_start_end', tuple()),
    ('to_fasta_record', ('seg1',)),
]


@pytest.mark.parametrize("method,badval", [*test_parameters])
def test_hybrecord_bad_seg_props(method, badval):
    test_record = hybkit.HybRecord.from_line(
        ART_HYB_PROPS_1['hyb_str'],
        hybformat_id=True,
        hybformat_ref=True
    )
    test_record.seg1_props['read_start'] = None
    with pytest.raises(RuntimeError):
        getattr(test_record, method)(*badval)


# ----- HybRecord misc private_function_tests -----
def test_hybrecord_misc_private():
    test_record = hybkit.HybRecord.from_line(
        ART_HYB_PROPS_1['hyb_str'],
        hybformat_id=True,
        hybformat_ref=True
    )
    read_flags = test_record._read_flags('bad_flag=B;', allow_undefined_flags=True)
    assert read_flags['bad_flag'] == 'B'
