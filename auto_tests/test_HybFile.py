#!/usr/bin/env python3
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

"""
Automatic testing of the hybkit HybRecord class.
"""

# import copy
import os

# import sys
# from contextlib import nullcontext as does_not_raise
import pytest

import hybkit
from auto_tests.test_helper_data import ART_BAD_HYB_STRS, ART_HYB_PROPS_1, ART_HYB_PROPS_ALL
from auto_tests.test_helper_functions import get_expected_result_context

# ----- Linting Directives:
# ruff: noqa: SLF001 ARG001

# ----- HybFile test misc disallowed commands -----
def test_hybfile_misc_disallowed(tmp_path):
    """Test misc methods of HybFile class that should raise errors."""
    hyb_autotest_file_name = os.path.join(tmp_path, 'hyb_autotest_file.hyb')
    with hybkit.HybFile(hyb_autotest_file_name, 'w') as hyb_autotest_file:
        with pytest.raises(RuntimeError):
            hyb_autotest_file.write('test_str')
        with pytest.raises(RuntimeError):
            hyb_autotest_file._ensure_hybrecord(None)


# ----- HybFile test reading/writing of hyb records. -----
test_parameters = [
    ('One-Record', 'Pass', [ART_HYB_PROPS_1['hyb_str']]),
    ('All-Records', 'Pass', [props['hyb_str'] for props in ART_HYB_PROPS_ALL]),
    ('All-Records-Multi', 'Pass', [props['hyb_str'] for props in ART_HYB_PROPS_ALL]),
]
for i, bad_hyb_str in enumerate(ART_BAD_HYB_STRS, start=1):
    test_parameters.append(
        ('Bad-Record-' + str(i), 'Raise', [bad_hyb_str])
    )


@pytest.mark.parametrize(('test_name', 'expectation', 'hyb_strs'), [*test_parameters])
def test_hybfile_io(test_name, expectation, hyb_strs, tmp_path):
    """Test reading/writing of hyb records."""
    hyb_autotest_file_name = os.path.join(tmp_path, 'hyb_autotest_file.hyb')
    expect_context = get_expected_result_context(expect_str=expectation)
    all_hyb_strs = ''.join(hyb_strs)
    with open(hyb_autotest_file_name, mode='w') as hyb_autotest_file:
        hyb_autotest_file.write(all_hyb_strs)

    assert hybkit.util.hyb_exists(hyb_autotest_file_name)

    with expect_context:
        with hybkit.HybFile.open(hyb_autotest_file_name, 'r') as hyb_autotest_file:
            for hyb_record in hyb_autotest_file:
                hyb_record_str = hyb_record.to_line(newline=True)
                assert hyb_record_str in hyb_strs

    if expectation.lower() == 'pass':
        with hybkit.HybFile.open(hyb_autotest_file_name, 'r') as hyb_autotest_file:
            first_record = hyb_autotest_file.read_record()
        with hybkit.HybFile.open(hyb_autotest_file_name, 'r') as hyb_autotest_file:
            all_records = hyb_autotest_file.read_records()
        assert first_record == all_records[0]
        with hybkit.HybFile(hyb_autotest_file_name, 'w') as hyb_autotest_file:
            hyb_autotest_file.write_records(write_records=[hyb_record, hyb_record])
            hyb_autotest_file.write_record(write_record=hyb_record)
            hyb_autotest_file.write_fh(all_hyb_strs)
