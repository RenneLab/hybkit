#!/usr/bin/env python3
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

"""
Automatic testing of the hybkit TypeFinder class.
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

# ----- TypeFinder Misc - Minimal -----
def test_typefinder_misc():
    """Test construction of HybRecord class with minimal information."""
    with pytest.raises((RuntimeError, TypeError)):
        hybkit.type_finder.TypeFinder()

test_parameters = [
    ('make_string_match_params', 'Pass', 
     STRING_MATCH_PARAMS_1['params_str'], STRING_MATCH_PARAMS_1['params_dict']),
    ('make_string_match_params', 'Raise', 
     BAD_STRING_MATCH_PARAMS_1['params_str'], BAD_STRING_MATCH_PARAMS_1['params_dict']),
    ('make_string_match_params', 'Raise', 
     BAD_STRING_MATCH_PARAMS_2['params_str'], BAD_STRING_MATCH_PARAMS_2['params_dict']),
    ('make_id_map_params', 'Pass', 
     ID_MAP_PARAMS_1['params_str'], ID_MAP_PARAMS_1['params_dict']),
    ('make_id_map_params', 'Raise', 
     BAD_ID_MAP_PARAMS_1['params_str'], BAD_ID_MAP_PARAMS_1['params_dict']),
]
arg_string = "method,expect_str,params_str,params_dict"
@pytest.mark.parametrize(arg_string, [*test_parameters])
def test_typefinder_make_params(method, expect_str, params_str, params_dict):
    expect_context = get_expected_result_context(expect_str)
    if not os.path.isdir(test_out_dir):
        os.mkdir(test_out_dir)
    with open(make_params_autotest_file_name, 'w') as make_params_autotest_file:
        make_params_autotest_file.write(params_str)
    with expect_context:
        gen_params = getattr(hybkit.type_finder.TypeFinder, method)(make_params_autotest_file_name)
        assert gen_params == params_dict


test_record_strings = [
    (HYB_STR_1, 'microRNA', 'mRNA'),
]
for hyb_props in ART_HYB_PROPS_ALL + [ART_HYB_NOTYPE_PROPS]:
    test_record_strings.append((hyb_props['hyb_str'], hyb_props['seg1_type'], hyb_props['seg2_type']))
test_records = [
    (hybkit.HybRecord.from_line(hyb_strs[0]), *hyb_strs) for hyb_strs in test_record_strings
]

test_parameter_sets = [
    ('badmethod', 'Raise', None, {}),
    ('hybformat', 'Pass', None, {}),
    ('hybformat', 'Raise', {'badparam': True}, {}),
    ('string_match', 'Pass', STRING_MATCH_PARAMS_1['params_dict'], {}),
    #('string_match', 'Raise', BAD_STRING_MATCH_PARAMS_1['params_dict'], {}),
    #('string_match', 'Raise', BAD_STRING_MATCH_PARAMS_2['params_dict'], {}),
    ('id_map', 'Pass', ID_MAP_PARAMS_1['params_dict'], {}),
    #('id_map', 'Raise', BAD_ID_MAP_PARAMS_1['params_dict'], {}),
]
test_parameters = []
for test_parameter_set in test_parameter_sets:
    for test_record_set in test_records:
        test_parameters.append((*test_parameter_set, test_record_set))

arg_string = "method,expect_str,method_params,test_params,hyb_record_params"
@pytest.mark.parametrize(arg_string, [*test_parameters])
def test_typefinder_methods(method, expect_str, method_params, test_params, hyb_record_params):
    expect_context = get_expected_result_context(expect_str)
    hybkit.type_finder.TypeFinder._reset()
    with expect_context:
        hybkit.type_finder.TypeFinder.set_method(method, method_params)

        test_hyb_record, test_hyb_string, test_seg1_type, test_seg2_type = hyb_record_params
        use_test_hyb_record = copy.deepcopy(test_hyb_record)

        gen_seg_1_type = hybkit.type_finder.TypeFinder.find(use_test_hyb_record.seg1_props)
        assert gen_seg_1_type == test_seg1_type
        gen_seg_2_type = hybkit.type_finder.TypeFinder.find(use_test_hyb_record.seg2_props)
        assert gen_seg_2_type == test_seg2_type
        
        use_test_hyb_record.eval_types(allow_unknown=(test_seg1_type is None))
        if test_seg1_type is None:
            assert use_test_hyb_record.get_seg1_type()  == 'unknown'
        else:
            assert use_test_hyb_record.get_seg1_type() == test_seg1_type
        if test_seg2_type is None:
            assert use_test_hyb_record.get_seg2_type() == 'unknown'
        else:
            assert use_test_hyb_record.get_seg2_type() == test_seg2_type

        if test_seg1_type is None:
            use_test_hyb_record_2 = copy.deepcopy(test_hyb_record)
            with pytest.raises(RuntimeError):
                use_test_hyb_record_2.eval_types()

# ----- Begin Old Tests -----

def old_test_type_finder():
    # Generic Tests
    def do_nothing(*args, **kwargs):
        pass
    with pytest.raises(RuntimeError):
        hybkit.type_finder.TypeFinder()
    with pytest.raises(RuntimeError):
        hybkit.type_finder.TypeFinder.set_method('bad_method')
    hybkit.type_finder.TypeFinder.set_custom_method(do_nothing)

    # Defualt Hybformat
    hyb_record = hybkit.HybRecord.from_line(HYB_STR_1)
    hybkit.type_finder.TypeFinder.set_method('hybformat')
    hyb_record.eval_types()
    with pytest.raises(RuntimeError):
        hybkit.type_finder.TypeFinder()
    with pytest.raises(RuntimeError):
        hybkit.type_finder.TypeFinder.set_method('bad_method')
    # Non-Defualt String-Match
    with pytest.raises(RuntimeError):
        hybkit.type_finder.TypeFinder.make_string_match_params('badfile')
    with pytest.raises(RuntimeError):
        hybkit.type_finder.TypeFinder.make_string_match_params(bad1_match_legend_file_name)
    with pytest.raises(RuntimeError):
        hybkit.type_finder.TypeFinder.make_string_match_params(bad2_match_legend_file_name)
    match_params = hybkit.type_finder.TypeFinder.make_string_match_params(
        bad3_match_legend_file_name)
    hybkit.type_finder.TypeFinder.set_method('string_match', match_params)
    hyb_record = hybkit.HybRecord.from_line(HYB_STR_1)
    hyb_record.settings['check_complete_seg_types'] = True
    with pytest.raises(RuntimeError):
        hyb_record.eval_types()
    # match_params = {'startswith':[('MIMAT', 'MIMAT')], 'endswith':[('microRNA', 'miRNA')]}
    # hybkit.type_finder.TypeFinder.set_method('string_match', match_params)
    # with pytest.raises(RuntimeError):
    #     hyb_record.eval_types()
    match_params = hybkit.type_finder.TypeFinder.make_string_match_params(match_legend_file_name)
    hybkit.type_finder.TypeFinder.set_method('string_match', match_params)
    hyb_record = hybkit.HybRecord.from_line(HYB_STR_1)
    hyb_record.eval_types()
    # Non-Defualt ID-Map
    with pytest.raises(RuntimeError):
        hybkit.type_finder.TypeFinder.make_id_map_params()
    with pytest.raises(RuntimeError):
        hybkit.type_finder.TypeFinder.make_id_map_params('wrong_id_map_type')
    with pytest.raises(RuntimeError):
        hybkit.type_finder.TypeFinder.make_id_map_params(type_file_pairs='wrong_id_map_type')
    with pytest.raises(RuntimeError):
        hybkit.type_finder.TypeFinder.make_id_map_params([bad1_id_map_legend_file_name])
    with pytest.raises(RuntimeError):
        id_map_params = hybkit.type_finder.TypeFinder.make_id_map_params(
            type_file_pairs=[('seq1type', id_map_legend_file_name),
                             ('seq2type', id_map_legend_file_name)])
    id_map_params = hybkit.type_finder.TypeFinder.make_id_map_params(
        type_file_pairs=[('seqtype', id_map_legend_file_name)])
    id_map_params = hybkit.type_finder.TypeFinder.make_id_map_params([id_map_legend_file_name])
    hybkit.type_finder.TypeFinder.set_method('id_map', id_map_params)
    hyb_record = hybkit.HybRecord.from_line(HYB_STR_1)
    hyb_record.eval_types()
    hyb_record = hybkit.HybRecord.from_line(HYB_STR_1)
    hyb_record.seg1_props['ref_name'] = 'not_real_name'
    with pytest.raises(RuntimeError):
        hyb_record.eval_types(allow_unknown=False)