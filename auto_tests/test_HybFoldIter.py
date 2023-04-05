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


# ----- Begin HybFoldIter Tests -----
# ----- Constructor Tests -----
def test_hybfolditer_constructor_misc(tmp_path):
    hyb_autotest_file_name = os.path.join(tmp_path, 'hyb_autotest_file.hyb')
    vienna_autotest_file_name = os.path.join(tmp_path, 'vienna_autotest_file.vienna')
    hyb_file = hybkit.HybFile(hyb_autotest_file_name, 'w')
    fold_file = hybkit.ViennaFile(vienna_autotest_file_name, 'w')
    with pytest.raises(RuntimeError):
        use_iter = hybkit.HybFoldIter(hyb_file, None)
    with pytest.raises(RuntimeError):
        use_iter = hybkit.HybFoldIter(None, fold_file)
    with pytest.raises(RuntimeError):
        use_iter = hybkit.HybFoldIter(hyb_file, fold_file, iter_error_mode='BadMode')




# ----- Direct Match Tests -----
def test_hybfolditer_direct_match_misc():
    """Test misc HybRecord / FoldRecord Match Properties"""
    use_props = ART_BAD_HYB_VIENNA_PROPS_4
    hyb_record = hybkit.HybRecord.from_line(use_props['hyb_str'])
    fold_record = hybkit.FoldRecord.from_vienna_string(use_props['vienna_str'])
    copy.deepcopy(hyb_record).set_fold_record(
        copy.deepcopy(fold_record),
        allow_energy_mismatch=True
    )
    with pytest.raises(RuntimeError):
        copy.deepcopy(hyb_record).set_fold_record(
            copy.deepcopy(fold_record),
            allow_energy_mismatch=False
        )


# ----- Start Test HybFoldIter -----
test_param_sets = []
for prop_set in [ART_HYB_VIENNA_PROPS_1, ART_HYB_VIENNA_PROPS_2]:
    if prop_set['overlapping']:
        test_name = 'Good-Overlap'
    else:
        test_name = 'Good-Static'
    test_param_sets.append(
        (test_name, 'Pass', prop_set, False)
    )
test_param_sets.append(('Mismatch_Seq_Static', 'Raise', ART_BAD_HYB_VIENNA_PROPS_1, False))
test_param_sets.append(('Disallowed-Overlap', 'Raise', ART_BAD_HYB_VIENNA_PROPS_2, False))
test_param_sets.append(('Mismatch_Seq_Dynamic', 'Raise', ART_BAD_HYB_VIENNA_PROPS_3, False))
test_param_sets.append(('Mismatch_Energy', 'Raise', ART_BAD_HYB_VIENNA_PROPS_4, True))
test_param_sets.append(('Missing_Vienna', 'Raise', ART_BAD_HYB_VIENNA_PROPS_5, False))
test_param_sets.append(('Bad_Vienna', 'Raise', ART_BAD_HYB_VIENNA_PROPS_6, True))
test_param_sets.append(('Error99', 'Raise', ART_BAD_HYB_VIENNA_PROPS_7, True))
test_param_sets.append(('Insertion', 'Raise', ART_BAD_HYB_VIENNA_PROPS_8, True))
test_param_sets.append(('Deletion', 'Raise', ART_BAD_HYB_VIENNA_PROPS_9, True))


test_parameters = []
for test_param_set in test_param_sets:
    for combine in ['Combine', 'Separate']:
        test_parameters.append((*test_param_set, combine))
arg_string = 'test_name,expectation,test_props,combine_str,skip_match_check'


@pytest.mark.parametrize(arg_string, [*test_parameters])
def test_hybfolditer_io(test_name, expectation, test_props, combine_str,
                        skip_match_check, tmp_path):
    hyb_autotest_file_name = os.path.join(tmp_path, 'hyb_autotest_file.hyb')
    vienna_autotest_file_name = os.path.join(tmp_path, 'vienna_autotest_file.vienna')
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
        vienna_autotest_file.write(vienna_str + '\n')

    assert hybkit.util.hyb_exists(hyb_autotest_file_name)
    assert hybkit.util.vienna_exists(vienna_autotest_file_name)

    with pytest.raises(RuntimeError):
        use_iter = hybkit.HybFoldIter(
            hyb_autotest_file_name,
            vienna_autotest_file_name,
            combine=combine
        )

    hyb_file = hybkit.HybFile(hyb_autotest_file_name, 'r')
    fold_file = hybkit.ViennaFile(vienna_autotest_file_name, 'r', seq_type=seq_type)
    use_iter = hybkit.HybFoldIter(
        hyb_file,
        fold_file,
        combine=combine,
        iter_error_mode='raise',
    )
    with expect_context:
        for ret_items in use_iter:
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

        # Print iteration report
        print(use_iter.report())
        use_iter.print_report()

        # Ensure iteration has occurred
        assert hyb_record
        assert fold_record

    if expectation == 'Raise':
        hyb_record = hybkit.HybFile(hyb_autotest_file_name, 'r').read_record()
        fold_file = hybkit.ViennaFile(
            vienna_autotest_file_name,
            'r',
            seq_type=seq_type,
        )
        fold_record = fold_file.read_record(override_error_mode='warn_return')
        if not skip_match_check:
            assert not "made it here"
            assert not fold_record.matches_hyb_record(hyb_record)
            assert fold_record.count_hyb_record_mismatches(hyb_record) >= test_props['mismatches']
            with pytest.raises(RuntimeError):
                fold_record.ensure_matches_hyb_record(hyb_record)

    else:
        hyb_file = hybkit.HybFile(hyb_autotest_file_name, 'r')
        fold_file = hybkit.ViennaFile(vienna_autotest_file_name, 'r', seq_type=seq_type)
        use_iter = hybkit.HybFoldIter(
            hyb_file,
            fold_file,
            combine=True,
            iter_error_mode='raise',
        )
        for hyb_record in use_iter:
            fold_record = hyb_record.fold_record
            hyb_record.eval_types()
            hyb_record.eval_mirna()
            mirna_detail_dict = hyb_record.mirna_detail(detail='all')
            assert mirna_detail_dict['mirna_seq'] == test_props['seg1_seq']
            assert mirna_detail_dict['target_seq'] == test_props['seg2_seq']
            assert mirna_detail_dict['mirna_fold'] == test_props['seg1_fold']
            assert mirna_detail_dict['target_fold'] == test_props['seg2_fold']


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
    'Missing_Vienna': allow_error_modes,
    'Bad_Vienna': allow_error_modes,
    'Error99': allow_error_modes,
    'Insertion': allow_error_modes,
    'Deletion': allow_error_modes,
}
test_param_sets = {
    'StaticAllowed': ART_HYB_VIENNA_PROPS_1,
    'StaticDisallowed': ART_HYB_VIENNA_PROPS_2,
    'DynamicAllowed': ART_HYB_VIENNA_PROPS_2,
    'MismatchSeqStatic': ART_BAD_HYB_VIENNA_PROPS_1,
    'MismatchSeqDynamic': ART_BAD_HYB_VIENNA_PROPS_3,
    'MismatchEnergy': ART_BAD_HYB_VIENNA_PROPS_4,
    'Missing_Vienna': ART_BAD_HYB_VIENNA_PROPS_5,
    'Bad_Vienna': ART_BAD_HYB_VIENNA_PROPS_6,
    'Error99': ART_BAD_HYB_VIENNA_PROPS_7,
    'Insertion': ART_BAD_HYB_VIENNA_PROPS_8,
    'Deletion': ART_BAD_HYB_VIENNA_PROPS_9,
}
test_seq_types = {
    'StaticAllowed': 'static',
    'StaticDisallowed': 'static',
    'DynamicAllowed': 'dynamic',
    'MismatchSeqStatic': 'static',
    'MismatchSeqDynamic': 'dynamic',
    'MismatchEnergy': 'dynamic',
    'Missing_Vienna': 'dynamic',
    'Bad_Vienna': 'dynamic',
    'Error99': 'dynamic',
    'Insertion': 'dynamic',
    'Deletion': 'dynamic',
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
arg_string = 'test_name,iter_error_mode,expectation,seq_type,test_props'


@pytest.mark.parametrize(arg_string, [*test_parameters])
def test_hybfolditer_iter_error_mode(test_name, iter_error_mode,
                                     expectation, seq_type, test_props, tmp_path):
    hyb_autotest_file_name = os.path.join(tmp_path, 'hyb_autotest_file.hyb')
    vienna_autotest_file_name = os.path.join(tmp_path, 'vienna_autotest_file.vienna')
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


default_skips = hybkit.settings.HybFoldIter_settings_info['max_sequential_skips'][0]
test_parameters = [
    ('Allowed_Skips', 'warn_skip', 'Pass', ART_BAD_HYB_VIENNA_PROPS_1, (default_skips - 1)),
    ('Disallowed_Skips', 'warn_skip', 'Raise', ART_BAD_HYB_VIENNA_PROPS_1, (default_skips + 2))
]
arg_string = 'test_name,iter_error_mode,expectation,test_props,num_skips'


@pytest.mark.parametrize(arg_string, [*test_parameters])
def test_hybfolditer_iter_num_skips(test_name, iter_error_mode,
                                    expectation, test_props, num_skips, tmp_path):
    hyb_autotest_file_name = os.path.join(tmp_path, 'hyb_autotest_file.hyb')
    vienna_autotest_file_name = os.path.join(tmp_path, 'vienna_autotest_file.vienna')
    expect_context = get_expected_result_context(expectation)
    hyb_str = test_props['hyb_str']
    vienna_str = test_props['vienna_str']
    with open(hyb_autotest_file_name, 'w') as hyb_autotest_file:
        for i in range(num_skips):
            hyb_autotest_file.write(hyb_str)
    with open(vienna_autotest_file_name, 'w') as vienna_autotest_file:
        for i in range(num_skips):
            vienna_autotest_file.write(vienna_str)

    hyb_file = hybkit.HybFile(hyb_autotest_file_name, 'r')
    fold_file = hybkit.ViennaFile(vienna_autotest_file_name, 'r')
    use_iter = hybkit.HybFoldIter(
        hyb_file,
        fold_file,
        iter_error_mode=iter_error_mode,
    )
    with expect_context:
        ret_items = None
        for ret_items in use_iter:
            pass
        if 'skip' not in iter_error_mode:
            assert ret_items

# # ----- Start CT-Format HybFoldIter Tests -----
# test_param_sets = []
# for prop_set in [ART_HYB_CT_PROPS_1, ART_HYB_CT_PROPS_2]:
#     if prop_set['overlapping']:
#         test_name = 'Good-Overlap'
#     else:
#         test_name = 'Good-Static'
#     test_param_sets.append(
#         (test_name, 'Pass', prop_set, False)
#     )
# test_param_sets.append(('Mismatch_Seq_Static', 'Raise', ART_BAD_HYB_CT_PROPS_1, False))
# test_param_sets.append(('Disallowed-Overlap', 'Raise', ART_BAD_HYB_CT_PROPS_2, False))
# test_param_sets.append(('Mismatch_Seq_Dynamic', 'Raise', ART_BAD_HYB_CT_PROPS_3, False))
# test_param_sets.append(('Mismatch_Energy', 'Raise', ART_BAD_HYB_CT_PROPS_4, True))
# test_param_sets.append(('Missing_CT', 'Raise', ART_BAD_HYB_CT_PROPS_5, False))
# test_param_sets.append(('Bad_CT', 'Raise', ART_BAD_HYB_CT_PROPS_6, True))
# test_param_sets.append(('Error99', 'Raise', ART_BAD_HYB_CT_PROPS_7, True))
# test_param_sets.append(('InDel', 'Raise', ART_BAD_HYB_CT_PROPS_8, True))


# test_parameters = []
# for test_param_set in test_param_sets:
#     for combine in ['Combine', 'Separate']:
#         test_parameters.append((*test_param_set, combine))

# @pytest.mark.parametrize("test_name,expectation,test_props,
# combine_str,skip_match_check", [*test_parameters])
# def test_hybfolditer_io(test_name, expectation, test_props, combine_str, skip_match_check):
#     if not os.path.isdir(test_out_dir):
#         os.mkdir(test_out_dir)
#     expect_context = get_expected_result_context(expectation)
#     combine = (combine_str == 'Combine')
#     hyb_str = test_props['hyb_str']
#     ct_str = test_props['ct_str']
#     if test_props['overlapping']:
#         seq_type = 'dynamic'
#     else:
#         seq_type = 'static'

#     with open(hyb_autotest_file_name, 'w') as hyb_autotest_file:
#         hyb_autotest_file.write(hyb_str)
#     with open(ct_autotest_file_name, 'w') as ct_autotest_file:
#         ct_autotest_file.write(ct_str + '\n')

#     assert hybkit.util.hyb_exists(hyb_autotest_file_name)
#     assert hybkit.util.ct_exists(ct_autotest_file_name)

#     with pytest.raises(RuntimeError):
#         use_iter= hybkit.HybFoldIter(
#             hyb_autotest_file_name,
#             ct_autotest_file_name,
#             combine=combine
#         )

#     hyb_file = hybkit.HybFile(hyb_autotest_file_name, 'r')
#     fold_file = hybkit.CTFile(ct_autotest_file_name, 'r', seq_type=seq_type)
#     use_iter= hybkit.HybFoldIter(
#         hyb_file,
#         fold_file,
#         combine=combine,
#         iter_error_mode='raise',
#     )
#     with expect_context:
#         for ret_items in use_iter:
#             if combine:
#                 hyb_record = ret_items
#                 fold_record = hyb_record.fold_record
#             else:
#                 hyb_record, fold_record = ret_items
#                 hyb_record.set_fold_record(fold_record)

#             assert hyb_record.to_line() == hyb_str
#             assert fold_record.to_ct_string() == ct_str

#             assert fold_record.matches_hyb_record(hyb_record)
#             assert fold_record.count_hyb_record_mismatches(hyb_record) == 0
#             fold_record.ensure_matches_hyb_record(hyb_record)

#         #Print iteration report
#         print(use_iter.report())
#         use_iter.print_report()

#         # Ensure iteration has occurred
#         assert hyb_record
#         assert fold_record

#     if expectation == 'Raise':
#         hyb_record = hybkit.HybFile(hyb_autotest_file_name, 'r').read_record()
#         fold_record = hybkit.CTFile(
#             ct_autotest_file_name,
#             'r',
#             seq_type=seq_type,
#         ).read_record(override_error_mode='return')
#         if not isinstance(fold_record, tuple) and not skip_match_check:
#             assert not fold_record.matches_hyb_record(hyb_record)
#             assert (fold_record.count_hyb_record_mismatches(hyb_record)
#                       >= test_props['mismatches'])
#             with pytest.raises(RuntimeError):
#                 fold_record.ensure_matches_hyb_record(hyb_record)


# # Test iter_error_mode values:
# raise_error_modes = ['raise']
# allow_error_modes = ['warn_return', 'warn_skip', 'skip', 'return']
# all_iter_error_modes = [*raise_error_modes, *allow_error_modes]
# test_names_allowed = {
#     'StaticAllowed': all_iter_error_modes,
#     'StaticDisallowed': allow_error_modes,
#     'DynamicAllowed': all_iter_error_modes,
#     'MismatchSeqStatic': allow_error_modes,
#     'MismatchSeqDynamic': allow_error_modes,
#     'MismatchEnergy': allow_error_modes,
# }
# test_param_sets = {
#     'StaticAllowed': ART_HYB_CT_PROPS_1,
#     'StaticDisallowed': ART_HYB_CT_PROPS_2,
#     'DynamicAllowed': ART_HYB_CT_PROPS_2,
#     'MismatchSeqStatic': ART_BAD_HYB_CT_PROPS_1,
#     'MismatchSeqDynamic': ART_BAD_HYB_CT_PROPS_3,
#     'MismatchEnergy': ART_BAD_HYB_CT_PROPS_4,
# }
# test_seq_types = {
#     'StaticAllowed': 'static',
#     'StaticDisallowed': 'static',
#     'DynamicAllowed': 'dynamic',
#     'MismatchSeqStatic': 'static',
#     'MismatchSeqDynamic': 'dynamic',
#     'MismatchEnergy': 'dynamic',
# }
# test_parameters = []
# for test_name, test_props in test_param_sets.items():
#     allowed_errors = test_names_allowed[test_name]
#     seq_type = test_seq_types[test_name]
#     for iter_error_mode in all_iter_error_modes:
#         if iter_error_mode in allowed_errors:
#             expectation = 'Pass'
#         else:
#             expectation = 'Raise'
#         test_parameters.append((test_name, iter_error_mode, expectation, seq_type, test_props))

# @pytest.mark.parametrize("test_name,iter_error_mode,expectation,
# seq_type,test_props", [*test_parameters])
# def test_hybfolditer_iter_error_mode(test_name, iter_error_mode,
# expectation, seq_type, test_props):
#     expect_context = get_expected_result_context(expectation)
#     hyb_str = test_props['hyb_str']
#     ct_str = test_props['ct_str']
#     with open(hyb_autotest_file_name, 'w') as hyb_autotest_file:
#         hyb_autotest_file.write(hyb_str)
#     with open(ct_autotest_file_name, 'w') as ct_autotest_file:
#         ct_autotest_file.write(ct_str)

#     hyb_file = hybkit.HybFile(hyb_autotest_file_name, 'r')
#     fold_file = hybkit.CTFile(ct_autotest_file_name, 'r', seq_type=seq_type)
#     use_iter = hybkit.HybFoldIter(
#         hyb_file,
#         fold_file,
#         iter_error_mode=iter_error_mode,
#     )
#     with expect_context:
#         for ret_items in use_iter:
#             pass


# default_skips = hybkit.settings.HybFoldIter_settings_info['max_sequential_skips'][0]
# test_parameters = [
#     ('Allowed_Skips', 'warn_skip', 'Pass',
# ART_BAD_HYB_CT_PROPS_1, (default_skips - 1)),
#     ('Disallowed_Skips', 'warn_skip', 'Raise',
# ART_BAD_HYB_CT_PROPS_1, (default_skips + 2))
# ]
# @pytest.mark.parametrize("test_name,iter_error_mode,
# expectation,test_props,num_skips", [*test_parameters])
# def test_hybfolditer_iter_num_skips(test_name, iter_error_mode,
# expectation, test_props, num_skips):
#     expect_context = get_expected_result_context(expectation)
#     hyb_str = test_props['hyb_str']
#     ct_str = test_props['ct_str']
#     with open(hyb_autotest_file_name, 'w') as hyb_autotest_file:
#         for i in range(num_skips):
#             hyb_autotest_file.write(hyb_str + '\n')
#     with open(ct_autotest_file_name, 'w') as ct_autotest_file:
#         for i in range(num_skips):
#             ct_autotest_file.write(ct_str + '\n')

#     hyb_file = hybkit.HybFile(hyb_autotest_file_name, 'r')
#     fold_file = hybkit.CTFile(ct_autotest_file_name, 'r')
#     use_iter = hybkit.HybFoldIter(
#         hyb_file,
#         fold_file,
#         iter_error_mode=iter_error_mode,
#     )
#     with expect_context:
#         for ret_items in use_iter:
#             pass


# if __name__ == '__main__':
#     test_type_finder()
