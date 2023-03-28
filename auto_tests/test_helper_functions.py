#!/usr/bin/env python3
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

"""
Helper functions for automatic testing of hybkit code.
"""

import copy
import os
import hybkit
import pytest
from contextlib import nullcontext as does_not_raise

import test_helper_data

# ----- Test Assistance Functions -----
# Generate HybRecord and FoldRecord objects for testing.

# def default_hyb_records():
#     """Generate two HybRecord objects for tests"""
#     hyb_record_1 = hybkit.HybRecord.from_line(
#         HYB_STR_1,
#         hybformat_id=True,
#         hybformat_ref=True,
#     )
#     fold_record = None
#     #fold_record = hybkit.DynamicFoldRecord.from_vienna_string(VIENNA_STR_1)

#     hyb_record_2 = hybkit.HybRecord(
#         id=hyb_record_1.id,
#         seq=hyb_record_1.seq,
#         seg1_props=copy.deepcopy(hyb_record_1.seg1_props),
#         seg2_props=copy.deepcopy(hyb_record_1.seg2_props),
#         fold_record=fold_record,
#     )
#     return hyb_record_1, hyb_record_2, fold_record

# Get expected result string for exception testing.
def get_expected_result_string(is_allowed=False):
    """Return string identifying expected pass/error result"""
    if is_allowed:
        return 'Pass'
    else:
        return 'Raise'
    
# Get expected result context for exception testing.
def get_expected_result_context(expect_str, error_types = (TypeError, RuntimeError)):
    """Return context for testing allowed types."""
    if expect_str.lower() == 'pass':
        return does_not_raise()
    elif expect_str.lower() == 'raise':
        if isinstance(error_types, list):
            error_types = tuple(error_types)
        elif not isinstance(error_types, tuple):
            error_types = (error_types,)
        return pytest.raises(error_types)
