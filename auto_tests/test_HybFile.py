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

# ----- HybFile test reading/writing of hyb records. -----
test_parameters = [
    ('One-Record', 'Pass', [ART_HYB_PROPS_1['hyb_str']], ),
    ('All-Records', 'Pass', [props['hyb_str'] for props in ART_HYB_PROPS_ALL]),
    ('All-Records-Multi', 'Pass', [props['hyb_str'] for props in ART_HYB_PROPS_ALL]),
]
for i, bad_hyb_str in enumerate(ART_BAD_HYB_STRS, start=1):
    test_parameters.append(
        ('Bad-Record-' + str(i), 'Raise', [bad_hyb_str])
    )

@pytest.mark.parametrize("test_name,expectation,hyb_strs",[*test_parameters])
def test_hybfile_io(test_name, expectation, hyb_strs):
    if not os.path.isdir(test_out_dir):
        os.mkdir(test_out_dir)
    expect_context = get_expected_result_context(expectation)
    all_hyb_strs = '\n'.join(hyb_strs)
    with open(hyb_autotest_file_name, 'w') as hyb_autotest_file:
        hyb_autotest_file.write(all_hyb_strs + '\n')

    assert hybkit.util.hyb_exists(hyb_autotest_file_name)

    with expect_context:
        with hybkit.HybFile.open(hyb_autotest_file_name, 'r') as hyb_autotest_file:
            for hyb_record in hyb_autotest_file:
                hyb_record_str = hyb_record.to_line()
                assert hyb_record_str in hyb_strs

    if expectation.lower() == 'pass':
        with hybkit.HybFile.open(hyb_autotest_file_name, 'r') as hyb_autotest_file:
            first_record = hyb_autotest_file.read_record()
        with hybkit.HybFile.open(hyb_autotest_file_name, 'r') as hyb_autotest_file:
            all_records = hyb_autotest_file.read_records()
        assert first_record == all_records[0]
        with hybkit.HybFile(hyb_autotest_file_name, 'w') as hyb_autotest_file:
            hyb_autotest_file.write_records([hyb_record, hyb_record])
            hyb_autotest_file.write_record(hyb_record)
            hyb_autotest_file.write_fh(all_hyb_strs)

