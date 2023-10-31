#!/usr/bin/env python3
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

"""
Automatic testing of the hybkit TypeFinder class.
"""

# import sys
import copy
import os

# from contextlib import nullcontext as does_not_raise
import pytest

import hybkit

# ----- Import Testing Helper Data -----
from auto_tests.test_helper_data import (
    ART_HYB_MATCHTYPE_PROPS,
    ART_HYB_NOTYPE_PROPS,
    ART_HYB_PROPS_ALL,
    BAD_ID_MAP_PARAMS_1,
    BAD_ID_MAP_PARAMS_2,
    BAD_STRING_MATCH_PARAMS_1,
    BAD_STRING_MATCH_PARAMS_2,
    HYB_STR_1,
    ID_MAP_PARAMS_1,
    STRING_MATCH_PARAMS_1,
)

# ----- Import Testing Helper Functions -----
from auto_tests.test_helper_functions import get_expected_result_context

# ----- Linting Directives -----
# ruff: noqa: SLF001 ARG001

# ----- TypeFinder Misc - Minimal -----
def test_typefinder_misc():
    """Test construction of HybRecord class with minimal information."""
    type_finder = hybkit.type_finder.TypeFinder
    with pytest.raises((RuntimeError, TypeError)):
        type_finder()
    type_finder.set_custom_method(print)
    type_finder._reset()
    with pytest.raises((RuntimeError, TypeError)):
        type_finder.find({})


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
    ('make_id_map_params', 'Raise',
     BAD_ID_MAP_PARAMS_2['params_str'], BAD_ID_MAP_PARAMS_2['params_dict']),
]
arg_string = 'method,expect_str,params_str,params_dict'


@pytest.mark.parametrize(arg_string, [*test_parameters])
def test_typefinder_make_params(method, expect_str, params_str, params_dict, tmp_path):
    """Test TypeFinder.make_params methods."""
    make_params_autotest_file_name = os.path.join(tmp_path, 'make_params_autotest_file.txt')
    expect_context = get_expected_result_context(expect_str)
    with pytest.raises(TypeError):
        gen_params = getattr(hybkit.type_finder.TypeFinder, method)(23)
    with pytest.raises(FileNotFoundError):
        gen_params = getattr(hybkit.type_finder.TypeFinder, method)(make_params_autotest_file_name)
    with open(make_params_autotest_file_name, 'w') as make_params_autotest_file:
        make_params_autotest_file.write(params_str)
    with expect_context:
        gen_params = getattr(hybkit.type_finder.TypeFinder, method)(make_params_autotest_file_name)
        assert gen_params == params_dict


test_parameter_sets = [
    ('badmethod', 'Raise', None, {}),
    ('hybformat', 'Pass', None, {}),
    ('hybformat', 'Raise', {'badparam': True}, {}),
    ('string_match', 'Pass', STRING_MATCH_PARAMS_1['params_dict'], {}),
    ('string_match', 'Raise', {}, {}),
    # ('string_match', 'Raise', BAD_STRING_MATCH_PARAMS_1['params_dict'], {}),
    # ('string_match', 'Raise', BAD_STRING_MATCH_PARAMS_2['params_dict'], {}),
    ('id_map', 'Pass', ID_MAP_PARAMS_1['params_dict'], {}),
    # ('id_map', 'Raise', BAD_ID_MAP_PARAMS_1['params_dict'], {}),
]
test_record_strings = [
    (HYB_STR_1, 'microRNA', 'mRNA'),
]
for hyb_props in ART_HYB_PROPS_ALL + [ART_HYB_NOTYPE_PROPS, ART_HYB_MATCHTYPE_PROPS]:
    test_record_strings.append(
        (hyb_props['hyb_str'], hyb_props['seg1_type'], hyb_props['seg2_type'])
    )

test_records = [
    (hybkit.HybRecord.from_line(hyb_strs[0]), *hyb_strs) for hyb_strs in test_record_strings
]
test_parameters = []
for test_parameter_set in test_parameter_sets:
    for test_record_set in test_records:
        test_parameters.append((*test_parameter_set, test_record_set))
arg_string = 'method,expect_str,method_params,test_params,hyb_record_params'


@pytest.mark.parametrize(arg_string, [*test_parameters])
def test_typefinder_methods(method, expect_str, method_params, test_params, hyb_record_params):
    """Test TypeFinder methods."""
    expect_context = get_expected_result_context(expect_str)
    hybkit.type_finder.TypeFinder._reset()
    with expect_context:
        hybkit.type_finder.TypeFinder.set_method(method, method_params)

        test_hyb_record, test_hyb_string, test_seg1_type, test_seg2_type = hyb_record_params
        use_test_hyb_record = copy.deepcopy(test_hyb_record)

        gen_seg_1_type = hybkit.type_finder.TypeFinder.find(use_test_hyb_record.seg1_props)
        gen_seg_2_type = hybkit.type_finder.TypeFinder.find(use_test_hyb_record.seg2_props)

        if not (method == 'hybformat' and test_seg1_type == 'MatchType'):
            assert gen_seg_1_type == test_seg1_type
            assert gen_seg_2_type == test_seg2_type

            use_test_hyb_record.eval_types(allow_unknown=(test_seg1_type is None))
            if test_seg1_type is None:
                assert use_test_hyb_record.get_seg1_type() == 'unknown'
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
